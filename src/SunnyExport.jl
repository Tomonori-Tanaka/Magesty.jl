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

# ── physical spin (effective spin length) ────────────────────────────────────
#
# The SCE couplings are fit with unit spin directions, so they absorb the spin
# magnitude: the per-bond matrix M reconstructs `ê_a' M ê_b` regardless of the
# moment size, and the classical energy is independent of the spin length chosen
# in Sunny. The *magnon dispersion*, however, scales as ħω ∝ 1/S for a fixed
# energy landscape (precession of an angular momentum ħS), so the emitted script
# must encode the physical effective spin S_eff = m/(g μ_B) to be correct.
#
# Two emission routes (`scaling`) carry S_eff differently:
#
#   :moment   — Sunny's `Moment` spin IS the physical S_eff. The static energy is
#               preserved (`energy(sys) == predict_energy - j0`) AND the dispersion
#               is physical. Requires S_eff to be a half-integer (Sunny rejects
#               other `Moment` spins).
#     bilinear : J = M / (s_i s_j)          (s factors cancel in the energy)
#     single-ion : mode factor, with s = S_eff
#
#   :coupling — Sunny's `Moment` spin is a fixed half-integer placeholder s₀ = 1;
#               the couplings carry S_eff so only the dispersion is physical. This
#               works for ANY positive real S_eff (itinerant / non-half-integer
#               moments). The dispersion is invariant under an overall spin scale
#               c (s_i → c S_i, J → J/c), so picking c = s₀/S keeps the Moment
#               half-integer while the dispersion stays physical. The price: the
#               whole represented landscape is scaled by c, so `energy(sys)` is no
#               longer the SCE energy (LSWT dispersion is unaffected). Exact for a
#               uniform S_eff; off-diagonal-exact / Larmor-approximate otherwise.
#     bilinear : J = M / (√(s_i s_j) √(S_i S_j))   → M/(s₀ S) for uniform S
#     single-ion : classical factor 1/(s₀ S_i) (scales the gap as A/S, matching
#                  the bilinear's 1/S); the quantum :dipole renorm needs the real
#                  S and cannot be represented by the placeholder, so :coupling +
#                  :dipole + single-ion is rejected.

# Prefactor converting the classical quadratic form `ê'Aê` (single-ion) into the
# Sunny spin operator whose mode expectation reproduces the physical magnon. The
# Moment spin `smom` may differ from the physical `sphys` (the :coupling route),
# so both are passed; for the :moment route they coincide.
#   :dipole_uncorrected — classical limit. The operator's classical value scales
#       as smom², so factor = 1/(smom·sphys) makes the landscape scale by
#       smom/sphys (consistent with the bilinear rescale). Equals 1/S² when
#       smom = sphys = S.
#   :dipole — Sunny renormalizes a rank-2 operator by smom(2smom-1)/2 (the quantum
#       coherent-state quadrupole; it vanishes at smom = 1/2). The inverse
#       2/(smom(2smom-1)) reproduces the classical form, but this only yields the
#       physical single-ion dispersion when smom = sphys (the :moment route);
#       callers forbid the placeholder (:coupling) case.
function _sunny_onsite_factor(smom::Real, sphys::Real, mode::Symbol)::Float64
    sm = Float64(smom)
    S = Float64(sphys)
    if mode === :dipole_uncorrected
        return 1 / (sm * S)
    elseif mode === :dipole
        # The quantum rank-2 renormalization depends on the *actual* Moment spin, so
        # it only reproduces the physical single-ion dispersion when the Moment spin
        # is the physical spin (the :moment route). Reject a placeholder Moment here
        # so the invariant is enforced at the point of use, not only by `_sunny_check_plan`.
        sm ≈ S || throw(ArgumentError(
            "single-ion anisotropy in :dipole mode needs the Moment spin to be the " *
            "physical spin (s = $smom, S_eff = $sphys); a placeholder Moment cannot " *
            "carry the quantum quadrupole. Use mode = :dipole_uncorrected or scaling = :moment."))
        sm > 0.5 || throw(ArgumentError(
            "single-ion anisotropy with spin s = $smom in :dipole mode: a spin s ≤ 1/2 " *
            "carries no rank-2 quadrupole (the s(2s-1)/2 renormalization vanishes). " *
            "Pass mode = :dipole_uncorrected for this model."))
        return 2 / (sm * (2sm - 1))
    else
        throw(ArgumentError("mode must be :dipole or :dipole_uncorrected; got :$mode"))
    end
end

# Convenience for the :moment route, where the Moment spin equals the physical
# spin. Keeps `1/s²` (uncorrected) and `2/(s(2s-1))` (:dipole) at s = S_eff.
_sunny_onsite_factor(s::Real, mode::Symbol)::Float64 = _sunny_onsite_factor(s, s, mode)

# Placeholder Moment spin used by the :coupling route. The dispersion is invariant
# under the overall spin scale, so any half-integer works; s₀ = 1 is the simplest.
const _SUNNY_COUPLING_S0 = 1

# Fully resolved per-species spin plan for one route. `sphys` is the physical S_eff
# (drives the coupling rescale and the mode auto-selection); `smom` is what is
# written to Sunny's `Moment` (= sphys for :moment, = s₀ for :coupling). Defined
# here, before the emission functions that take it in their signatures.
struct _SunnySpinPlan
    sphys::Dict{String, Real}
    smom::Dict{String, Real}
    gmap::Dict{String, Real}
    mode::Symbol     # resolved :dipole | :dipole_uncorrected
    route::Symbol    # :moment | :coupling
end

# A spin length is "half-integer" when 2s is an integer (positivity is checked
# separately by `_sunny_resolve_spin_maps`).
_sunny_is_half_integer(s::Real)::Bool =
    isapprox(Float64(2s), round(Float64(2s)); atol = 1e-9)

# Literal for a spin / g value: keep an exact rational as `p//q`, otherwise %.12g.
_sunny_num_literal(x::Rational)::String = string(numerator(x), "//", denominator(x))
_sunny_num_literal(x::Integer)::String = string(x)
_sunny_num_literal(x::Real)::String = _sunny_fmt(x)

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

# Classify a SALC by its per-site angular momenta into the interaction kind the
# Sunny export supports. `ls` is the per-site `l` list, so `length(ls)` equals
# the cluster's site count. Both export routes (supercell and primitive) share
# this mapping.
#   :spin_independent — all l = 0 (pure scalar; a constant)
#   :pair_l1          — two sites with l = l = 1 (off-diagonal exchange matrix)
#   :onsite_l2        — one site with l = 2 (single-ion anisotropy)
#   :unsupported      — anything else (not representable in the Sunny model)
function _classify_salc(ls::AbstractVector{<:Integer})::Symbol
    if all(==(0), ls)
        return :spin_independent
    elseif ls == [1, 1]
        return :pair_l1
    elseif ls == [2]
        return :onsite_l2
    else
        return :unsupported
    end
end

# Shared message for a SALC whose interaction kind the Sunny export cannot
# represent. `n_body` is the cluster's site count.
_unsupported_salc_string(ν::Integer, n_body::Integer, ls, Lf)::String =
    "SALC #$ν: n_body=$n_body ls=$ls Lf=$Lf " *
    "(unsupported; only l=l=1 pairs and l=2 single-ion are exported)"

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
        kind = _classify_salc(ls)

        if kind === :spin_independent
            # Spin-independent SALC. The per-site Z_{0,0} = (4π)^(-1/2) factors
            # exactly cancel the (4π)^(n_C/2) scaling, leaving a pure constant.
            for cbc in group
                fval = first(cbc.folded_tensor)
                const_offset += coeff * cbc.multiplicity * length(cbc.clusters) * fval
            end
            continue
        end

        if kind === :pair_l1
            for cbc in group
                M0 = _sunny_l1_pair_matrix(cbc.folded_tensor)
                w = coeff * scaling * cbc.multiplicity
                for cluster in cbc.clusters
                    a, b = cluster[1], cluster[2]
                    a == b && error(
                        "internal error: self-pair in a 2-body cluster (atoms = $cluster); " *
                        "this indicates corrupt SCEModel data")
                    if a < b
                        key = (a, b)
                        pairs[key] = get(pairs, key, Z3) + w * M0
                    else
                        key = (b, a)
                        pairs[key] = get(pairs, key, Z3) + w * transpose(M0)
                    end
                end
            end
        elseif kind === :onsite_l2
            for cbc in group
                A0 = _sunny_l2_onsite_matrix(cbc.folded_tensor)
                w = coeff * scaling * cbc.multiplicity
                for cluster in cbc.clusters
                    a = cluster[1]
                    onsites[a] = get(onsites, a, Z3) + w * A0
                end
            end
        else
            push!(skipped, _unsupported_salc_string(ν, n_C, ls, Lf))
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

# ── primitive-cell route ─────────────────────────────────────────────────────
#
# The supercell route above bakes the orbit multiplicity into each coupling and
# places one bond per supercell atom pair; it reproduces the energy exactly but
# the magnetic cell equals the training supercell, so the spin-wave dispersion
# is folded. To obtain an unfolded dispersion we map the interactions onto the
# chemical primitive cell.
#
# This is exact only when the fitted model respects the modeling assumption that
# the interaction range is below half the supercell, equivalently when every
# pair satisfies `multiplicity * n_clusters == n_translations`. Then each primitive
# bond is set once (WITHOUT the multiplicity; the periodic replication of the
# cell reproduces it) and the offset is the minimum image. When the assumption
# is violated (cutoff >= L/2 with non-collinear bonds folding onto one atom pair)
# the primitive model cannot reproduce the supercell energy; `clean` is false and
# callers fall back to the supercell route.

# Primitive-cell geometry derived from the model's symmetry (spglib) data.
struct _SunnyPrimitive
    latvecs::SMatrix{3, 3, Float64, 9}      # columns = primitive lattice vectors (Å)
    positions::Vector{SVector{3, Float64}}  # sublattice fractional positions
    types::Vector{String}                   # sublattice species
    reshape_matrix::SMatrix{3, 3, Int, 9}   # supercell = primitive * reshape_matrix
end

# Decomposed primitive-cell interactions.
struct _SunnyPrimitiveModel
    prim::_SunnyPrimitive
    # bonds keyed by (sublattice_i, sublattice_j, n1, n2, n3); coupling M with
    # pair energy ê_i' M ê_j (sublattice i in any cell, j in cell + n).
    bonds::Dict{NTuple{5, Int}, SMatrix{3, 3, Float64, 9}}
    onsites::Dict{Int, SMatrix{3, 3, Float64, 9}}  # per sublattice, ê' A ê
    clean::Bool                             # primitive unfolding is exact
    skipped::Vector{String}
end

# Reconstruct the primitive lattice vectors and sublattice data from the
# supercell + symmetry. The primitive lattice vectors are the three shortest
# independent pure-translation vectors (falling back to the supercell vectors
# along non-periodic directions); this matches the translations used to build
# the SALC orbits, so the bond mapping stays consistent.
function _sunny_primitive(model::SCEModel)::_SunnyPrimitive
    struc = model.basis.structure
    sym = model.basis.symmetry
    L = SMatrix{3, 3, Float64}(struc.supercell.lattice_vectors)
    xf = struc.supercell.x_frac
    rep0 = sym.map_p2s[1, 1]

    cands = SVector{3, Float64}[]
    for k = 1:sym.ntran
        df = xf[:, sym.map_p2s[1, k]] .- xf[:, rep0]
        df = df .- round.(df)
        push!(cands, L * SVector{3, Float64}(df))
    end
    for j = 1:3
        push!(cands, SVector{3, Float64}(L[1, j], L[2, j], L[3, j]))
    end
    sort!(cands; by = norm)

    basis = SVector{3, Float64}[]
    for v in cands
        norm(v) < 1e-8 && continue
        if length(basis) == 0
            push!(basis, v)
        elseif length(basis) == 1
            norm(cross(basis[1], v)) > 1e-6 && push!(basis, v)
        elseif length(basis) == 2
            abs(dot(cross(basis[1], basis[2]), v)) > 1e-6 && push!(basis, v)
        end
        length(basis) == 3 && break
    end
    length(basis) == 3 || error(
        "could not determine three primitive lattice vectors from the symmetry data")

    Lp = hcat(basis[1], basis[2], basis[3])
    det(Lp) < 0 && (Lp = hcat(basis[1], basis[2], -basis[3]))
    Lp = SMatrix{3, 3, Float64}(Lp)
    Lpi = inv(Lp)

    natp = sym.nat_prim
    prim_frac(i) = SVector{3, Float64}(
        mod.(Lpi * (L * SVector{3, Float64}(xf[:, sym.atoms_in_prim[i]])), 1.0))
    positions = [prim_frac(i) for i = 1:natp]
    types = [struc.kd_name[struc.supercell.kd_int_list[sym.atoms_in_prim[i]]] for i = 1:natp]
    reshape_matrix = SMatrix{3, 3, Int}(round.(Int, Lpi * L))
    return _SunnyPrimitive(Lp, positions, types, reshape_matrix)
end

# Distinct primitive-cell offsets `n` such that sublattice `subl(b)` in cell `n`
# sits at the minimum image distance from sublattice `subl(a)` in the home cell.
# Replicates the equal-minimum-distance selection of `Clusters._set_mindist_pairs`
# from the stored 27-image arrays, so every degenerate (equal-distance) lattice
# vector connecting the pair is recovered — not only one minimum image. The
# second return value is `true` when all offsets are integer to tolerance (the
# pair maps cleanly onto the primitive cell).
function _sunny_equal_distance_offsets(
    struc, sym, prim::_SunnyPrimitive, Lpi, a::Int, b::Int,
)::Tuple{Vector{NTuple{3, Int}}, Bool}
    XC = struc.x_image_cart
    exist = struc.exist_image
    tol = sym.tol
    ncell = size(XC, 3)
    a0 = SVector{3, Float64}(XC[1, a, 1], XC[2, a, 1], XC[3, a, 1])
    dfp = prim.positions[sym.map_s2p[b].atom] .- prim.positions[sym.map_s2p[a].atom]

    mind = Inf
    @inbounds for c = 1:ncell
        exist[c] || continue
        bc = SVector{3, Float64}(XC[1, b, c], XC[2, b, c], XC[3, b, c])
        d = norm(bc - a0)
        d < mind && (mind = d)
    end

    offs = NTuple{3, Int}[]
    ok = true
    @inbounds for c = 1:ncell
        exist[c] || continue
        bc = SVector{3, Float64}(XC[1, b, c], XC[2, b, c], XC[3, b, c])
        disp = bc - a0
        isapprox(norm(disp), mind; atol = tol) || continue
        nf = Lpi * disp .- dfp
        n = round.(Int, nf)
        maximum(abs.(nf .- n)) < 1e-6 || (ok = false)
        push!(offs, (n[1], n[2], n[3]))
    end
    return offs, ok
end

# Map the supercell interactions onto primitive-cell bonds and single-ion terms.
function _sunny_build_primitive(model::SCEModel)::_SunnyPrimitiveModel
    prim = _sunny_primitive(model)
    struc = model.basis.structure
    sym = model.basis.symmetry
    Lpi = inv(prim.latvecs)

    subl(a) = sym.map_s2p[a].atom

    bonds = Dict{NTuple{5, Int}, SMatrix{3, 3, Float64, 9}}()
    onsites = Dict{Int, SMatrix{3, 3, Float64, 9}}()
    skipped = String[]
    clean = true
    Z3 = zero(SMatrix{3, 3, Float64})

    for (ν, group) in enumerate(model.basis.salcbasis.salc_list)
        isempty(group) && continue
        ls = group[1].ls
        kind = _classify_salc(ls)
        if kind === :pair_l1
            for cbc in group
                M0 = _sunny_l1_pair_matrix(cbc.folded_tensor)
                # Per-bond coupling WITHOUT multiplicity. Each equal-distance lattice
                # vector of the pair is placed as its own primitive bond; together
                # (after periodic replication) they reproduce the orbit's
                # multiplicity. ± / translation duplicates collapse to one canonical
                # bond and are placed once (dedup within the cbc).
                W = model.jphi[ν] * Fitting._cluster_scaling(2) * M0
                a, b = cbc.clusters[1][1], cbc.clusters[1][2]
                i, j = subl(a), subl(b)
                offs, ok = _sunny_equal_distance_offsets(struc, sym, prim, Lpi, a, b)
                ok || (clean = false)
                seen = Set{NTuple{5, Int}}()
                for n in offs
                    # Canonical orientation: lower sublattice first, ties broken by
                    # the lexicographically non-negative offset (Julia tuple `>=`).
                    if i < j || (i == j && n >= (0, 0, 0))
                        key = (i, j, n...)
                        val = W
                    else
                        key = (j, i, -n[1], -n[2], -n[3])
                        val = SMatrix{3, 3, Float64}(transpose(W))
                    end
                    key in seen && continue   # ± / duplicate of this same cbc
                    push!(seen, key)
                    bonds[key] = get(bonds, key, Z3) + val
                end
            end
        elseif kind === :onsite_l2
            for cbc in group
                A0 = _sunny_l2_onsite_matrix(cbc.folded_tensor)
                W = model.jphi[ν] * Fitting._cluster_scaling(1) * A0
                i = sym.map_s2p[cbc.clusters[1][1]].atom
                onsites[i] = get(onsites, i, Z3) + W
            end
        elseif kind === :spin_independent
            # spin-independent constant; dropped (Sunny carries no scalar term)
        else
            push!(skipped, _unsupported_salc_string(ν, length(ls), ls, group[1].Lf))
        end
    end

    return _SunnyPrimitiveModel(prim, bonds, onsites, clean, skipped)
end

# ── script emission ──────────────────────────────────────────────────────────

_sunny_fmt(x::Real)::String = @sprintf("%.12g", x)

# 3×3 matrix as a Julia literal "[a b c; d e f; g h i]".
function _sunny_mat_literal(M::AbstractMatrix)::String
    rows = [join((_sunny_fmt(M[r, c]) for c = 1:3), " ") for r = 1:3]
    return "[" * join(rows, "; ") * "]"
end

_sunny_vec_literal(v::AbstractVector{<:Real})::String =
    "[" * join((_sunny_fmt(x) for x in v), ", ") * "]"

# Single-ion operator line shared by both routes. `setter` is the Sunny call name,
# `target` its trailing argument (sublattice index or site tuple), and `factor` the
# mode-dependent prefactor (from `_sunny_onsite_factor`) that converts the classical
# quadratic form ê'Aê into the matching Sunny spin operator.
function _sunny_onsite_line(setter::String, A::AbstractMatrix, target::String,
                            factor::Float64)::String
    return "    $setter(sys, S -> $(_sunny_fmt(factor))*(" *
        join(("$(_sunny_fmt(A[p, q]))*S[$p]*S[$q]" for p = 1:3, q = 1:3), " + ") *
        "), $target)"
end

# `cell_route` is the spatial placement (`:primitive` / `:explicit`); `plan.route`
# is the spin-scaling route (`:moment` / `:coupling`). They are independent.
function _sunny_header(model::SCEModel, skipped::Vector{String}, cell_route::Symbol,
                       plan::_SunnySpinPlan)::String
    io = IOBuffer()
    println(io, "# Auto-generated by Magesty `sce_to_sunny` (do not edit by hand;")
    println(io, "# regenerate from the SCEModel instead).")
    println(io, "#")
    println(io, "# The SCE couplings are fit with unit spin directions, so they absorb")
    println(io, "# the spin magnitude (J_SCE = J_phys·S²). The magnon dispersion needs")
    println(io, "# the physical effective spin S_eff = m/(g μ_B); s = 1 inflates it by ~S.")
    if plan.route === :moment
        println(io, "# Spin scaling: :moment — Sunny's Moment spin IS S_eff (see the Moment")
        println(io, "# lines below); mode :$(plan.mode). Bilinear bonds are rescaled by")
        println(io, "# 1/(s_i s_j) and single-ion terms by the :$(plan.mode) factor, so the")
        println(io, "# static energy is unchanged (energy(sys) == predict_energy - j0) and")
        println(io, "# the dispersion is physical.")
    elseif plan.route === :coupling
        println(io, "# Spin scaling: :coupling — Sunny's Moment spin is the placeholder")
        println(io, "# s₀ = $(_SUNNY_COUPLING_S0); the physical S_eff is carried by the couplings")
        println(io, "# (J = M/(s₀·√(S_i S_j)), single-ion 1/(s₀ S_i)). Mode :$(plan.mode).")
        println(io, "# Only the magnon DISPERSION is physical; the static energy is rescaled")
        println(io, "# and is NOT the SCE energy. Exact for a uniform S_eff.")
    else
        error("unreachable scaling route: $(plan.route)")
    end
    println(io, "# Energies are in the unit of the fit (typically eV). The reference")
    println(io, "# energy j0 = $(_sunny_fmt(model.j0)) and any spin-independent terms are")
    println(io, "# dropped (Sunny has no constant energy term); the dispersion is unaffected.")
    println(io, "# Cell route: $(cell_route) " *
        (cell_route === :primitive ? "(chemical primitive cell; unfolded dispersion)." :
         "(training supercell; the dispersion is folded into the supercell BZ)."))
    if !isempty(skipped)
        println(io, "#")
        println(io, "# WARNING: the following SALCs are not representable in Sunny and were")
        println(io, "# dropped from this script:")
        for s in skipped
            println(io, "#   - $s")
        end
    end
    return String(take!(io))
end

# Emit the LSWT / magnetic-structure tail shared by both routes.
function _sunny_emit_tail(io::IO; folded::Bool)::Nothing
    println(io)
    println(io, "# ---- Magnetic structure (EDIT THIS BLOCK) ----")
    if folded
        println(io, "# This system is the training supercell, so the dispersion is folded.")
        println(io, "# `to_inhomogeneous` systems cannot be reshaped to a smaller cell.")
    else
        println(io, "# Set up the magnetic unit cell for your ground state. For a ferromagnet")
        println(io, "# (k = 0) nothing is needed. For a larger magnetic cell, reshape first, e.g.")
        println(io, "#   sys = reshape_supercell(sys, [2 0 0; 0 2 0; 0 0 2])")
    end
    println(io, "randomize_spins!(sys)")
    println(io, "minimize_energy!(sys)")
    println(io)
    println(io, "# ---- Linear spin-wave dispersion ----")
    println(io, "# EDIT: high-symmetry path in reciprocal lattice units of `cryst`.")
    println(io, "qs = [[0, 0, 0], [1/2, 0, 0], [1/2, 1/2, 0], [0, 0, 0]]")
    println(io, "path = q_space_path(cryst, qs, 400)")
    println(io, "swt = SpinWaveTheory(sys; measure = ssf_perp(sys))")
    println(io, "disp = dispersion(swt, path)")
    println(io)
    println(io, "# ---- Plot (requires a Makie backend) ----")
    println(io, "# `disp` is a (bands × q) matrix; plot each band as its own line.")
    println(io, "# `path.xticks` supplies the high-symmetry q-point positions and labels")
    println(io, "# (e.g. \"[0, 0, 0]\"); vlines mark the interior path corners.")
    println(io, "# Energy uses the fit's unit (typically eV); EDIT the label if not eV.")
    println(io, "# using GLMakie")
    println(io, "# fig = Figure()")
    println(io, "# ax = Axis(fig[1, 1]; xlabel = \"Wavevector\", ylabel = \"Energy (eV)\",")
    println(io, "#           xticks = path.xticks)")
    println(io, "# vlines!(ax, path.xticks[1][2:end-1]; color = :gray, linewidth = 0.75)")
    println(io, "# for b in axes(disp, 1); lines!(ax, disp[b, :]); end")
    println(io, "# fig")
    return nothing
end

# Emit the primitive-cell (unfolded) script.
function _sunny_emit_primitive(model::SCEModel, pm::_SunnyPrimitiveModel,
                               plan::_SunnySpinPlan)::String
    io = IOBuffer()
    print(io, _sunny_header(model, pm.skipped, :primitive, plan))
    println(io)
    println(io, "using Sunny")
    println(io)
    natp = length(pm.prim.positions)
    types = pm.prim.types
    println(io, "# ---- Crystal (chemical primitive cell) ----")
    println(io, "# P1 is forced so the fitted couplings are placed exactly as given.")
    println(io, "latvecs = $(_sunny_mat_literal(pm.prim.latvecs))")
    println(io, "positions = [")
    for p in pm.prim.positions
        println(io, "    $(_sunny_vec_literal(p)),")
    end
    println(io, "]")
    println(io, "types = [", join(("\"$t\"" for t in types), ", "), "]")
    println(io, "cryst = Crystal(latvecs, positions, 1; types = types)")
    println(io)
    moments = join(("$i => Moment(s = $(_sunny_num_literal(plan.smom[types[i]])), " *
                    "g = $(_sunny_num_literal(plan.gmap[types[i]])))" for i = 1:natp), ", ")
    println(io, "sys = System(cryst, [$moments], :$(plan.mode))")
    println(io)
    println(io, "# ---- Bilinear exchange (per primitive bond; ê_i' J ê_j) ----")
    for key in sort!(collect(keys(pm.bonds)))
        i, j, n1, n2, n3 = key
        M = pm.bonds[key] ./ _sunny_bond_denom(plan, types[i], types[j])
        println(io, "set_exchange!(sys, $(_sunny_mat_literal(M)), Bond($i, $j, [$n1, $n2, $n3]))")
    end
    if !isempty(pm.onsites)
        println(io)
        println(io, "# ---- Single-ion anisotropy (ê' A ê) ----")
        for i in sort!(collect(keys(pm.onsites)))
            factor = _sunny_plan_onsite_factor(plan, types[i])
            println(io, "let")
            println(io, _sunny_onsite_line("set_onsite_coupling!", pm.onsites[i], string(i), factor))
            println(io, "end")
        end
    end
    _sunny_emit_tail(io; folded = false)
    return String(take!(io))
end

# Minimum-image cell offset (supercell lattice units) for the bond a→b, used by
# the explicit (folded) route where each supercell atom is its own sublattice.
function _sunny_supercell_offset(model::SCEModel, a::Int, b::Int)::NTuple{3, Int}
    xf = model.basis.structure.supercell.x_frac
    isper = model.basis.structure.is_periodic
    o = MVector{3, Int}(0, 0, 0)
    for d = 1:3
        isper[d] || continue
        df = xf[d, b] - xf[d, a]
        # Minimum image; `floor(df + 0.5)` rounds the +0.5 boundary deterministically
        # (unlike `round`, which uses banker's rounding there). The two atoms carry
        # the same dipole under periodicity, so the choice does not affect energy.
        o[d] = -floor(Int, df + 0.5)
    end
    return (o[1], o[2], o[3])
end

# Emit the explicit supercell (folded) script; exact for any model.
function _sunny_emit_explicit(model::SCEModel, terms::_SunnyTerms,
                              plan::_SunnySpinPlan)::String
    io = IOBuffer()
    print(io, _sunny_header(model, terms.skipped, :explicit, plan))
    println(io)
    println(io, "using Sunny")
    println(io)
    struc = model.basis.structure
    L = SMatrix{3, 3, Float64}(struc.supercell.lattice_vectors)
    xf = struc.supercell.x_frac
    nat = size(xf, 2)
    species(a) = _sunny_atom_species(struc, a)
    println(io, "# ---- Crystal (training supercell; each atom is its own sublattice) ----")
    println(io, "latvecs = $(_sunny_mat_literal(L))")
    println(io, "positions = [")
    for a = 1:nat
        println(io, "    $(_sunny_vec_literal(SVector{3, Float64}(xf[1, a], xf[2, a], xf[3, a]))),")
    end
    println(io, "]")
    println(io, "types = [", join(("\"S$a\"" for a = 1:nat), ", "), "]")
    println(io, "cryst = Crystal(latvecs, positions, 1; types = types)")
    println(io)
    moments = join(("$a => Moment(s = $(_sunny_num_literal(plan.smom[species(a)])), " *
                    "g = $(_sunny_num_literal(plan.gmap[species(a)])))" for a = 1:nat), ", ")
    println(io, "sys = System(cryst, [$moments], :$(plan.mode))")
    println(io, "sys = to_inhomogeneous(sys)")
    println(io)
    println(io, "# ---- Bilinear exchange (ê_i' J ê_j) ----")
    for key in sort!(collect(keys(terms.pairs)))
        a, b = key
        M = terms.pairs[key] ./ _sunny_bond_denom(plan, species(a), species(b))
        o = _sunny_supercell_offset(model, a, b)
        println(io, "set_exchange_at!(sys, $(_sunny_mat_literal(M)), " *
            "(1, 1, 1, $a), (1, 1, 1, $b); offset = ($(o[1]), $(o[2]), $(o[3])))")
    end
    if !isempty(terms.onsites)
        println(io)
        println(io, "# ---- Single-ion anisotropy ----")
        for a in sort!(collect(keys(terms.onsites)))
            factor = _sunny_plan_onsite_factor(plan, species(a))
            line = _sunny_onsite_line("set_onsite_coupling_at!", terms.onsites[a], "(1, 1, 1, $a)", factor)
            println(io, "let")
            println(io, line)
            println(io, "end")
        end
    end
    _sunny_emit_tail(io; folded = true)
    return String(take!(io))
end

# Emit a `@warn` to the caller when SALCs were dropped (the generated script also
# lists them in its header, but that is easy to miss).
function _sunny_warn_skipped(skipped::Vector{String})::Nothing
    isempty(skipped) && return nothing
    @warn "sce_to_sunny: $(length(skipped)) SALC(s) cannot be represented in Sunny " *
        "and were dropped from the script (see its header for the list)."
    return nothing
end

# Look up `val` (a scalar or a `species => value` Dict) for one species. Magnetic
# species (those carrying couplings) must be present when a Dict is given;
# non-magnetic species fall back to `placeholder` (inert in LSWT).
function _sunny_lookup_species(val::Union{Real, AbstractDict}, sp::AbstractString,
                               magnetic::AbstractSet, placeholder::Real,
                               name::AbstractString)::Real
    if val isa AbstractDict
        if haskey(val, sp)
            return val[sp]
        elseif sp in magnetic
            throw(ArgumentError(
                "`$name` is missing magnetic species \"$sp\"; provide `$name` for " *
                "every magnetic species, or pass a scalar to cover all of them."))
        else
            return placeholder
        end
    else
        return sp in magnetic ? val : placeholder
    end
end

# Resolve the per-species effective spin length and g-factor. `spin` / `g` are a
# scalar (applied to every magnetic species) or a `species => value` Dict.
# Non-magnetic species get a neutral placeholder (s = 1, g = 2). `spin === nothing`
# (omitted) is a hard error explaining the physics.
function _sunny_resolve_spin_maps(spin::Union{Real, AbstractDict, Nothing},
                                  g::Union{Real, AbstractDict},
                                  all_species::AbstractVector{<:AbstractString},
                                  magnetic::AbstractSet)::Tuple{Dict{String, Real}, Dict{String, Real}}
    spin === nothing && throw(ArgumentError(
        "sce_to_sunny requires the keyword `spin`: the effective spin length " *
        "S_eff = m/(g μ_B) of each magnetic species (its local-moment magnitude). " *
        "The SCE couplings absorb the spin magnitude (J_SCE = J_phys·S²), so the " *
        "magnon dispersion needs the physical S_eff — without it the frequencies " *
        "are off by a factor ~S. Pass a scalar for all magnetic sites or a Dict, " *
        "e.g. `spin = 5//2` or `spin = Dict(\"Mn\" => 5//2)`. It need not be a " *
        "half-integer; use m/(g μ_B) for itinerant moments."))
    smap = Dict{String, Real}()
    gmap = Dict{String, Real}()
    for sp in all_species
        smap[sp] = _sunny_lookup_species(spin, sp, magnetic, 1, "spin")
        gmap[sp] = _sunny_lookup_species(g, sp, magnetic, 2, "g")
    end
    # The effective spin length is a magnitude, so it must be positive. The
    # half-integer restriction is route-specific (Sunny's `Moment` accepts only
    # multiples of 1/2 directly, the :moment route) and is checked in `_sunny_plan_spins`.
    for sp in magnetic
        smap[sp] > 0 || throw(ArgumentError(
            "spin for species \"$sp\" is $(smap[sp]); the effective spin length " *
            "S_eff = m/(g μ_B) must be positive."))
    end
    return smap, gmap
end

# Choose the Sunny system mode. `:auto` uses `:dipole` when every magnetic spin is
# a half-integer (a genuine quantum spin), otherwise `:dipole_uncorrected` (the
# classical limit, appropriate for non-half-integer / itinerant moments).
function _sunny_select_mode(mode::Symbol, magnetic_spins::AbstractVector{<:Real})::Symbol
    if mode === :auto
        return all(_sunny_is_half_integer, magnetic_spins) ? :dipole : :dipole_uncorrected
    elseif mode in (:dipole, :dipole_uncorrected)
        return mode
    else
        throw(ArgumentError(
            "mode must be :auto, :dipole, or :dipole_uncorrected; got :$mode"))
    end
end

# Choose the spin-scaling route. `:auto` uses `:moment` (Sunny's `Moment` spin = the
# physical S_eff; static energy and dispersion both exact) when every magnetic
# S_eff is a half-integer, otherwise `:coupling` (placeholder Moment, couplings
# carry S_eff; dispersion-only) so itinerant / non-half-integer moments still work.
function _sunny_select_scaling(scaling::Symbol,
                               magnetic_spins::AbstractVector{<:Real})::Symbol
    if scaling === :auto
        return all(_sunny_is_half_integer, magnetic_spins) ? :moment : :coupling
    elseif scaling === :moment
        all(_sunny_is_half_integer, magnetic_spins) || throw(ArgumentError(
            "scaling = :moment requires every magnetic spin to be a half-integer " *
            "(Sunny's `Moment` rejects other values), but got $(collect(magnetic_spins)). " *
            "Use scaling = :coupling (or :auto) to keep the physical S_eff in the " *
            "couplings with a placeholder Moment spin."))
        return :moment
    elseif scaling === :coupling
        return :coupling
    else
        throw(ArgumentError("scaling must be :auto, :moment, or :coupling; got :$scaling"))
    end
end

# Resolve spin/g, pick the scaling route and the system mode, and build the Moment
# spin map. `mode` auto-selection keys off the *physical* S_eff (so a :coupling
# placeholder Moment of s₀ = 1 does not masquerade as a quantum spin). Warns when
# the :coupling route is used with a non-uniform magnetic S_eff (off-diagonal
# exchange stays exact; the on-site/Larmor term becomes approximate).
function _sunny_plan_spins(spin::Union{Real, AbstractDict, Nothing},
                           g::Union{Real, AbstractDict},
                           all_species::AbstractVector{<:AbstractString},
                           magnetic::AbstractSet, mode::Symbol, scaling::Symbol)::_SunnySpinPlan
    sphys, gmap = _sunny_resolve_spin_maps(spin, g, all_species, magnetic)
    mag_spins = [sphys[sp] for sp in magnetic]
    route = _sunny_select_scaling(scaling, mag_spins)
    sel = _sunny_select_mode(mode, mag_spins)
    if route === :moment
        smom = copy(sphys)
    else
        smom = Dict{String, Real}(sp => _SUNNY_COUPLING_S0 for sp in all_species)
        if length(unique(round.(Float64.(mag_spins); digits = 9))) > 1
            per_species = join(("$sp = $(sphys[sp])" for sp in sort!(collect(magnetic))), ", ")
            @warn(
                "sce_to_sunny: scaling = :coupling with a non-uniform magnetic S_eff " *
                "($per_species); the dispersion is exact only for a uniform S_eff. The " *
                "off-diagonal exchange stays exact, but the on-site (Larmor) term is " *
                "approximate. For an exact multi-sublattice dispersion use scaling = :moment " *
                "with half-integer spins.")
        end
    end
    return _SunnySpinPlan(sphys, smom, gmap, sel, route)
end

# Per-bond denominator D so the emitted coupling is M/D (pair energy ê_i' (M/D) ê_j
# evaluated with Moment spins). :moment → s_i s_j (energy-preserving). :coupling →
# √(smom_i smom_j) √(S_i S_j), which gives the physical off-diagonal exchange
# √(S_i S_j)·J_phys for any spins and reduces to M/(s₀ S) for a uniform S.
function _sunny_bond_denom(plan::_SunnySpinPlan, spi::AbstractString,
                           spj::AbstractString)::Float64
    if plan.route === :moment
        return Float64(plan.sphys[spi]) * Float64(plan.sphys[spj])
    elseif plan.route === :coupling
        return sqrt(Float64(plan.smom[spi]) * Float64(plan.smom[spj])) *
               sqrt(Float64(plan.sphys[spi]) * Float64(plan.sphys[spj]))
    else
        error("unreachable scaling route: $(plan.route)")
    end
end

# Single-ion factor for one species under a plan (route- and mode-aware).
_sunny_plan_onsite_factor(plan::_SunnySpinPlan, sp::AbstractString)::Float64 =
    _sunny_onsite_factor(plan.smom[sp], plan.sphys[sp], plan.mode)

# Guard the one single-ion case the :coupling route cannot represent: a quantum
# (:dipole) rank-2 operator needs the real spin, not the s₀ placeholder.
function _sunny_check_plan(plan::_SunnySpinPlan, has_onsite::Bool)::Nothing
    if plan.route === :coupling && plan.mode === :dipole && has_onsite
        throw(ArgumentError(
            "scaling = :coupling with mode = :dipole and single-ion anisotropy is " *
            "unsupported: the quantum rank-2 renormalization s(2s-1)/2 depends on the " *
            "physical spin, which the placeholder Moment (s = $(_SUNNY_COUPLING_S0)) does " *
            "not carry, so the single-ion dispersion would be wrong. Use " *
            "mode = :dipole_uncorrected (classical single-ion, consistent with the " *
            "coupling route) or scaling = :moment with a half-integer spin."))
    end
    return nothing
end

# Real element of supercell atom `a` (the explicit route labels each atom with a
# fake "S$a" P1 type, so the physical spin keys off the actual species instead).
_sunny_atom_species(struc, a::Integer)::String =
    struc.kd_name[struc.supercell.kd_int_list[a]]

# Per-route species sets: all sublattice/atom species and the subset that carries
# couplings (magnetic). Drives spin resolution and mode auto-selection.
function _sunny_primitive_species(pm::_SunnyPrimitiveModel)::Tuple{Vector{String}, Set{String}}
    types = pm.prim.types
    mag = Set{Int}()
    for k in keys(pm.bonds)
        push!(mag, k[1]); push!(mag, k[2])
    end
    for i in keys(pm.onsites)
        push!(mag, i)
    end
    return unique(types), Set(types[i] for i in mag)
end

function _sunny_explicit_species(model::SCEModel, terms::_SunnyTerms)::Tuple{Vector{String}, Set{String}}
    struc = model.basis.structure
    nat = size(struc.supercell.x_frac, 2)
    mag = Set{Int}()
    for (a, b) in keys(terms.pairs)
        push!(mag, a); push!(mag, b)
    end
    for a in keys(terms.onsites)
        push!(mag, a)
    end
    return unique(_sunny_atom_species(struc, a) for a = 1:nat),
           Set(_sunny_atom_species(struc, a) for a in mag)
end

# Resolve per-species spin/g and mode for a route, then emit. Split into two named
# helpers (one per route) so the resolve+select+emit sequence is written once and
# the tuple unpacking is in a typed top-level function (a closure tripped JET).
function _sunny_primitive_script(model::SCEModel, pm::_SunnyPrimitiveModel,
                                 spin::Union{Real, AbstractDict, Nothing},
                                 g::Union{Real, AbstractDict}, mode::Symbol,
                                 scaling::Symbol)::String
    all_sp, mag = _sunny_primitive_species(pm)
    plan = _sunny_plan_spins(spin, g, all_sp, mag, mode, scaling)
    _sunny_check_plan(plan, !isempty(pm.onsites))
    return _sunny_emit_primitive(model, pm, plan)
end

function _sunny_explicit_script(model::SCEModel, terms::_SunnyTerms,
                                spin::Union{Real, AbstractDict, Nothing},
                                g::Union{Real, AbstractDict}, mode::Symbol,
                                scaling::Symbol)::String
    all_sp, mag = _sunny_explicit_species(model, terms)
    plan = _sunny_plan_spins(spin, g, all_sp, mag, mode, scaling)
    _sunny_check_plan(plan, !isempty(terms.onsites))
    return _sunny_emit_explicit(model, terms, plan)
end

"""
    sce_to_sunny(model::SCEModel; spin, g=2, mode=:auto, scaling=:auto, output=nothing, placement=:auto) -> String

Export a fitted spin-cluster-expansion model to a runnable [Sunny.jl](https://github.com/SunnySuite/Sunny.jl)
script that computes a linear spin-wave-theory (LSWT) magnon dispersion, and
return the script text.

The lowest-order SALCs are converted to a conventional spin Hamiltonian:
two-site `l₁ = l₂ = 1` terms become 3×3 bilinear exchange matrices (Heisenberg,
Dzyaloshinskii–Moriya, and anisotropic symmetric parts together), and single-site
`l = 2` terms become single-ion anisotropy. Higher-order SALCs (higher-`l` pairs,
three-body and beyond) cannot be represented in Sunny and are skipped, with a
warning listing what was dropped. The reference energy `j0` and spin-independent
terms are dropped (Sunny carries no constant energy term); the dispersion is
unaffected. Energies are in the unit of the fit (typically eV).

# Physical spin

The SCE couplings are fit with unit spin directions, so they absorb the spin
magnitude (`J_SCE = J_phys·S²`). The classical energy is therefore independent of
the spin length, but the magnon dispersion scales as `ħω ∝ 1/S` for a fixed energy
landscape. You must pass the physical effective spin `S_eff = m/(g μ_B)` (the
local-moment magnitude); using `s = 1` would inflate the dispersion by a factor
`~S` (for MnTe, `S = 5/2` ⇒ ~2.5× too high).

`scaling` selects how `S_eff` is encoded, since Sunny's `Moment` only accepts spin
lengths that are exact multiples of `1/2`:

- `:moment` — put `S_eff` directly into `Moment`. Each bilinear bond is rescaled by
  `1/(s_i s_j)` and each single-ion term by a mode-dependent factor, so `energy(sys)`
  still reproduces `predict_energy(model, …) - j0` *and* the dispersion is physical.
  Requires `S_eff` to be a half-integer.
- `:coupling` — keep `Moment` at a fixed half-integer placeholder (`s₀ = 1`) and let
  the couplings carry `S_eff` (`J = M/(s₀·√(S_i S_j))`, single-ion `1/(s₀ S_i)`).
  Accepts **any positive real `S_eff`** (itinerant / non-half-integer moments). The
  dispersion is invariant under an overall spin scale, so this still yields the
  physical magnon dispersion, but the represented energy landscape is rescaled —
  `energy(sys)` is then *not* the SCE energy. Exact for a uniform `S_eff`; for a
  non-uniform `S_eff` the off-diagonal exchange stays exact while the on-site
  (Larmor) term is approximate (a warning is emitted).

`:auto` (default) uses `:moment` when every magnetic `S_eff` is a half-integer and
`:coupling` otherwise, so half-integer moments keep the exact (energy-preserving)
behavior while itinerant moments still produce a physical dispersion.

# Arguments
- `model::SCEModel`: a fitted model (e.g. from `SCEModel(fit)` or
  `Magesty.load(SCEModel, path)`).

# Keyword arguments
- `spin::Union{Real, AbstractDict}` (required): the effective spin length
  `S_eff = m/(g μ_B)`. A scalar applies to every magnetic species; a
  `Dict(species => S_eff)` sets it per species (every magnetic species must be
  present). Must be positive. With `scaling = :moment` it must also be a
  half-integer; `scaling = :coupling` (or `:auto`) accepts any positive real.
  Omitting it is an error.
- `g::Union{Real, AbstractDict} = 2`: the `g`-factor passed to `Moment`; scalar or
  per-species `Dict`. It does not affect the bare dispersion (only an external
  field or neutron intensities would use it).
- `mode::Symbol = :auto`: Sunny system mode. `:auto` selects `:dipole` when every
  magnetic spin is a half-integer, otherwise `:dipole_uncorrected` (the classical
  limit, appropriate for non-half-integer / itinerant moments). `:dipole` applies
  the quantum single-ion renormalization `s(2s-1)/2` (undefined for `s ≤ 1/2`);
  `:dipole_uncorrected` uses the classical `s²`. Note: `mode = :dipole` combined
  with `scaling = :coupling` and single-ion anisotropy is rejected — the quantum
  renormalization cannot ride the placeholder Moment (use `:dipole_uncorrected`).
- `scaling::Symbol = :auto`: how the physical `S_eff` is encoded (see *Physical
  spin*). `:moment` puts `S_eff` in `Moment` (energy-preserving, half-integer only);
  `:coupling` keeps `Moment` at a placeholder `s₀ = 1` and rescales the couplings
  (dispersion-only, any positive real `S_eff`); `:auto` picks `:moment` for
  all-half-integer magnetic spins and `:coupling` otherwise.
- `output::Union{AbstractString, Nothing} = nothing`: when given, the script is
  also written to this file (`.jl` is appended if absent).
- `placement::Symbol = :auto`: `:primitive` maps interactions onto the chemical
  primitive cell for an unfolded dispersion; `:explicit` keeps the training
  supercell (the dispersion is folded into the supercell Brillouin zone) but is
  exact for any model. `:auto` chooses `:primitive` when the model is cleanly
  unfoldable (interaction range below half the supercell) and `:explicit`
  otherwise.

# Returns
- `String`: the full Sunny.jl script.

# Throws
- `ArgumentError`: if `spin` is omitted; if a `spin` / `g` `Dict` is missing a
  magnetic species; if any magnetic `spin` is not positive; if `scaling = :moment`
  (or `:auto` resolving to it) is used with a non-half-integer spin; if `mode`,
  `scaling`, or `placement` is invalid; if a model with single-ion anisotropy uses
  `mode = :dipole` with `s ≤ 1/2`; or if `scaling = :coupling` is combined with
  `mode = :dipole` and single-ion anisotropy.

# Examples
```julia
model = Magesty.load(SCEModel, "model.xml")
# MnTe: Mn²⁺ has S = 5/2 (half-integer ⇒ :moment route, energy-preserving).
script = sce_to_sunny(model; spin = 5//2, output = "lswt.jl")
# Itinerant moment, e.g. Fe 2.2 μB ⇒ S_eff = 1.1 (:auto picks the :coupling route).
script = sce_to_sunny(model; spin = 1.1)
# Per-species:
script = sce_to_sunny(model; spin = Dict("Mn" => 5//2, "Fe" => 1.1))
```
"""
function sce_to_sunny(
    model::SCEModel;
    spin::Union{Real, AbstractDict, Nothing} = nothing,
    g::Union{Real, AbstractDict} = 2,
    mode::Symbol = :auto,
    scaling::Symbol = :auto,
    output::Union{AbstractString, Nothing} = nothing,
    placement::Symbol = :auto,
)::String
    placement in (:auto, :primitive, :explicit) || throw(ArgumentError(
        "placement must be :auto, :primitive, or :explicit; got :$placement"))

    # The decomposition (`_sunny_build_*`) returns the unscaled M^SCE; the per-route
    # `_sunny_*_script` helpers resolve per-species spin/g, pick the scaling route and
    # system mode, and apply the spin rescale at emission.
    #
    # Build the primitive route only when it may be used; a forced :explicit request
    # skips it and avoids a redundant second SALC walk.
    if placement === :explicit
        terms = _sunny_build_terms(model)
        text = _sunny_explicit_script(model, terms, spin, g, mode, scaling)
        _sunny_warn_skipped(terms.skipped)
    else
        pm = _sunny_build_primitive(model)
        if pm.clean
            text = _sunny_primitive_script(model, pm, spin, g, mode, scaling)
            _sunny_warn_skipped(pm.skipped)
        else
            placement === :primitive && @warn(
                "Model is not cleanly unfoldable to the primitive cell " *
                "(interaction range reaches >= half the supercell). Falling back to " *
                ":explicit; the dispersion will be folded into the supercell BZ.")
            terms = _sunny_build_terms(model)
            text = _sunny_explicit_script(model, terms, spin, g, mode, scaling)
            _sunny_warn_skipped(terms.skipped)
        end
    end

    if output !== nothing
        outfile = endswith(output, ".jl") ? output : output * ".jl"
        write(outfile, text)
    end
    return text
end
