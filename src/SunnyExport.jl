# Export a fitted SCE model to a Sunny.jl linear-spin-wave-theory script.
#
# The SCE energy is `E = j0 + Σ_ν jphi[ν] Φ_ν({ê_i})`, with Φ_ν built from
# products of real tesseral harmonics Z_{l,m} of the unit spin directions.
# This file decomposes the supported low-order SALCs into conventional
# spin-model interactions and emits a runnable Sunny.jl script.
#
# Supported terms:
#   - 2-body, l1 = l2 = 1  → 3×3 bilinear exchange tensor (Heisenberg + DM + Γ).
#   - 1-body, l = 2        → single-ion anisotropy (symmetric traceless 3×3).
# Spin-independent (all-l = 0) SALCs fold into a constant offset that is
# dropped (Sunny carries no scalar energy term, exactly as j0 is dropped).
# Every other SALC (higher-l pairs, 3-body+) is recorded as skipped.
#
# The decomposition is convention-pinned by the tesseral harmonics:
#   Z_{1,m}(ê) = √(3/4π) e_μ,  μ(-1)=y, μ(0)=z, μ(+1)=x.
#   Z_{2,m}(ê) = ê' Q^(m) ê    (quadratic forms below), with the same
#   normalization. Both are cross-checked against `TesseralHarmonics.Zₗₘ`
#   in the component tests.

# ── intermediate representation ─────────────────────────────────────────────

# Decomposed interactions in supercell-atom-index space. `pairs` is keyed by
# an ordered atom pair `(a, b)` with `a < b`; its value `M` is accumulated so
# that the pair energy is `ê_a' M ê_b`. `onsites[a]` gives `A` with single-ion
# energy `ê_a' A ê_a`. `const_offset` collects spin-independent contributions
# (dropped on emission). `skipped` describes SALCs that cannot be represented.
struct _SunnyTerms
	pairs::Dict{Tuple{Int, Int}, SMatrix{3, 3, Float64, 9}}
	onsites::Dict{Int, SMatrix{3, 3, Float64, 9}}
	const_offset::Float64
	skipped::Vector{String}
end

# ── tesseral → Cartesian linear maps ────────────────────────────────────────

# Cartesian axis (x=1, y=2, z=3) → folded-tensor m-index for l = 1.
# folded m-index 1,2,3 ↔ m = -1,0,+1 ↔ axis y,z,x, so the inverse map is
# x→3, y→1, z→2.
const _SUNNY_L1_AXIS_TO_MIDX = (3, 1, 2)

# Build the 3×3 bilinear matrix `M` (Cartesian, xyz) for one l1 = l2 = 1
# coupled basis from its folded tensor (already Mf-collapsed), such that
# `Σ_{m1,m2} folded[m1,m2] Z_{1,m1}(ê_a) Z_{1,m2}(ê_b) = ê_a' M ê_b`.
# The Z prefactor (3/4π) = (√(3/4π))² is factored out here.
function _sunny_l1_pair_matrix(folded::AbstractMatrix{Float64})::SMatrix{3, 3, Float64, 9}
	pref = 3 / (4π)
	p = _SUNNY_L1_AXIS_TO_MIDX
	return SMatrix{3, 3, Float64}(
		ntuple(9) do k
			row = (k - 1) % 3 + 1
			col = (k - 1) ÷ 3 + 1
			pref * folded[p[row], p[col]]
		end...,
	)
end

# Build the 3×3 single-ion matrix `A` (Cartesian, xyz, symmetric traceless)
# for one l = 2 coupled basis from its folded tensor (length 5, m = -2..+2),
# such that `Σ_m folded[m] Z_{2,m}(ê) = ê' A ê`. Constants verified against
# `TesseralHarmonics.Zₗₘ`: a = √(15/16π), b = √(5/16π).
function _sunny_l2_onsite_matrix(folded::AbstractVector{Float64})::SMatrix{3, 3, Float64, 9}
	a = sqrt(15 / (16π))
	b = sqrt(5 / (16π))
	# folded index 1..5 ↔ m = -2,-1,0,+1,+2.
	Bxy = a * folded[1]   # m = -2
	Byz = a * folded[2]   # m = -1
	d0 = folded[3]        # m =  0
	Bxz = a * folded[4]   # m = +1
	d2 = folded[5]        # m = +2
	Bxx = -b * d0 + a * d2
	Byy = -b * d0 - a * d2
	Bzz = 2b * d0
	return SMatrix{3, 3, Float64}(
		Bxx, Bxy, Bxz,
		Bxy, Byy, Byz,
		Bxz, Byz, Bzz,
	)
end

# ── decomposition ───────────────────────────────────────────────────────────

# Walk the fitted SALCs and accumulate the supported interactions. The outer
# index ν aligns with `model.jphi[ν]`. Each key group holds one or more
# `CoupledBasis_with_coefficient`; each carries its own `multiplicity` (a
# scalar weight) and `clusters` (translation images, atom indices in site
# order). This mirrors the energy oracle `Fitting.design_matrix_energy_element`.
_sunny_build_terms(model::SCEModel)::_SunnyTerms =
	_sunny_decompose(model.basis.salcbasis.salc_list, model.jphi)

# Core SALC walk. Split from `_sunny_build_terms` so it can be exercised with
# synthetic key groups (skip-logic tests) without building a full `SCEModel`.
function _sunny_decompose(
	salc_list::AbstractVector{<:AbstractVector{<:CoupledBases.CoupledBasis_with_coefficient}},
	jphi::AbstractVector{<:Real},
)::_SunnyTerms
	pairs = Dict{Tuple{Int, Int}, SMatrix{3, 3, Float64, 9}}()
	onsites = Dict{Int, SMatrix{3, 3, Float64, 9}}()
	skipped = String[]
	const_offset = 0.0
	Z3 = zero(SMatrix{3, 3, Float64})

	for (ν, group) in enumerate(salc_list)
		isempty(group) && continue
		ls = group[1].ls
		Lf = group[1].Lf
		n_C = length(group[1].atoms)
		scaling = Fitting._cluster_scaling(n_C)
		coeff = Float64(jphi[ν])

		if all(==(0), ls)
			# Spin-independent SALC. The per-site Z_{0,0} = (4π)^(-1/2) factors
			# exactly cancel the (4π)^(n_C/2) scaling, leaving a pure constant.
			for cbc in group
				fval = first(cbc.folded_tensor)
				const_offset += coeff * cbc.multiplicity * length(cbc.clusters) * fval
			end
			continue
		end

		if n_C == 2 && ls == [1, 1]
			for cbc in group
				M0 = _sunny_l1_pair_matrix(cbc.folded_tensor)
				w = coeff * scaling * cbc.multiplicity
				for cluster in cbc.clusters
					a, b = cluster[1], cluster[2]
					a == b && error(
						"unexpected self-pair in a 2-body cluster (atoms = $cluster)")
					if a < b
						key = (a, b)
						pairs[key] = get(pairs, key, Z3) + w * M0
					else
						key = (b, a)
						pairs[key] = get(pairs, key, Z3) + w * transpose(M0)
					end
				end
			end
		elseif n_C == 1 && ls == [2]
			for cbc in group
				A0 = _sunny_l2_onsite_matrix(cbc.folded_tensor)
				w = coeff * scaling * cbc.multiplicity
				for cluster in cbc.clusters
					a = cluster[1]
					onsites[a] = get(onsites, a, Z3) + w * A0
				end
			end
		else
			push!(skipped,
				"SALC #$ν: n_body=$n_C ls=$ls Lf=$Lf " *
				"(unsupported; only l=l=1 pairs and l=2 single-ion are exported)")
		end
	end

	return _SunnyTerms(pairs, onsites, const_offset, skipped)
end

# Reconstruct the spin-dependent energy from the decomposed terms. Equals
# `predict_energy(model, spin_directions) - model.j0` when the model contains
# no skipped (unsupported) SALCs. Used by the round-trip component tests.
function _sunny_reconstruct_energy(
	terms::_SunnyTerms,
	spin_directions::AbstractMatrix{<:Real},
)::Float64
	e = terms.const_offset
	for ((a, b), M) in terms.pairs
		ea = SVector{3, Float64}(spin_directions[1, a], spin_directions[2, a], spin_directions[3, a])
		eb = SVector{3, Float64}(spin_directions[1, b], spin_directions[2, b], spin_directions[3, b])
		e += dot(ea, M * eb)
	end
	for (a, A) in terms.onsites
		ea = SVector{3, Float64}(spin_directions[1, a], spin_directions[2, a], spin_directions[3, a])
		e += dot(ea, A * ea)
	end
	return e
end
