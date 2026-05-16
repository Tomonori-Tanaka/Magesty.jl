"""
	Fitting.jl

This module contains functions for optimizing the SCE coefficients.
"""
module Fitting

using Base.Threads
using LinearAlgebra
using Printf
using MultivariateStats
using Statistics
using StaticArrays
using ..TesseralHarmonics
using ..Structures
using ..Symmetries
using ..Clusters
using ..CoupledBases
using ..SALCBases
using ..SpinConfigs

export AbstractEstimator, OLS, Ridge

"""
	_cluster_scaling(n_sites::Integer) -> Float64

Return the SCE basis normalization factor `(4π)^(n_sites/2)` for an
`n_sites`-body cluster contribution. The factor compensates the
`1/√(4π)` carried by each per-site tesseral harmonic so that fitted
coefficients `Jφ` are in the input energy unit (typically eV).

See the Magesty.jl technical notes for the derivation. This helper is
internal; callers in this module invoke it directly as
`_cluster_scaling(n_C)` to keep the scaling convention in one place.
"""
@inline _cluster_scaling(n_sites::Integer)::Float64 = (4π)^(n_sites / 2)

"""
	AbstractEstimator

Abstract type for SCE coefficient estimation methods.
"""
abstract type AbstractEstimator end

"""
	OLS <: AbstractEstimator

Ordinary Least Squares estimator (no regularization).
"""
struct OLS <: AbstractEstimator end

"""
	Ridge <: AbstractEstimator

L2-regularized least-squares (ridge) estimator. The bias column is
excluded from the penalty by the dispatch boundary (see
`solve_coefficients`).

# Fields
- `lambda::Float64`: Regularization strength. `λ = 0` reduces to OLS.

# Examples
```julia
estimator = Ridge(lambda = 0.1)
```
"""
struct Ridge <: AbstractEstimator
	lambda::Float64
end

Ridge(; lambda::Real = 0.0) = Ridge(Float64(lambda))

"""
	EnergyWorkspace

Internal scratch state reused across calls to `design_matrix_energy_element`.
Pools the heap-allocated buffers that would otherwise be re-created on every
call. Each thread in `build_design_matrix_energy` owns one workspace; the
workspace itself never leaves the threaded loop.

# Fields

- `searched_pairs::Set{UInt}` — de-duplication set keyed by
  `_atoms_hash_key(atoms_sorted)`. Inside the
  `for itrans in symmetry.symnum_translation` loop, each translation
  produces a sorted-atom tuple; translations that yield the same sorted
  configuration would contribute identically, so we keep only one
  representative per orbit. The Set is `empty!`'d at the start of every
  element call but keeps its hash-table capacity, so the second-and-later
  calls in the same workspace incur no rehashing.

- `sh_values::Vector{Float64}` — flattened tesseral spherical-harmonic
  values for all sites in a cluster. The slice for site `i` is
  `sh_values[sh_offsets[i] + 1 : sh_offsets[i + 1]]` and has length
  `2 * ls[i] + 1`; entry `sh_values[sh_offsets[i] + m_idx]` holds
  `Zₗₘ_unsafe(ls[i], m, spin_direction)` for `m = m_idx - ls[i] - 1`.
  A single contiguous `Vector{Float64}` removes the outer-vector pointer
  indirection and keeps all sites' values cache-adjacent. Sized via
  `_ensure_sh_buffer!`; the buffer only grows so capacity is preserved
  across calls.

- `sh_offsets::Vector{Int}` — cumulative offsets for `sh_values`. Length
  `N + 1` after `_ensure_sh_buffer!`, with `sh_offsets[1] = 0` and
  `sh_offsets[i + 1] = sh_offsets[i] + 2 * ls[i] + 1`. Built once per
  element call from `cbc.ls`.

- `legendre_buf::Vector{Float64}` — Legendre-polynomial cache passed to
  the buffered overload `Zₗₘ_unsafe(l, m, uvec, buf)`. Sized to
  `max(ls) + 1` (the upper bound required by `_required_buf_size`).
  Reusing this buffer instead of letting `Zₗₘ_unsafe` allocate one per
  `(l, m)` call eliminates the dominant remaining allocation source.

Not exported. Construct with `EnergyWorkspace()`.
"""
mutable struct EnergyWorkspace
	searched_pairs::Set{UInt}
	sh_values::Vector{Float64}
	sh_offsets::Vector{Int}
	legendre_buf::Vector{Float64}

	EnergyWorkspace() = new(
		Set{UInt}(),
		Vector{Float64}(),
		Vector{Int}(),
		Vector{Float64}(),
	)
end

"""
	GradWorkspace

Internal scratch state reused across calls to `calc_∇ₑu!`. Each thread in
`build_design_matrix_torque` owns one workspace.

# Fields

- `searched_pairs::Set{UInt}` — same role as in
  [`EnergyWorkspace`](@ref): de-duplicates symmetry translations that
  yield identical sorted-atom configurations.

- `sh_values::Vector{Float64}` — same flattened layout as in
  [`EnergyWorkspace`](@ref): the slice for site `i` is
  `sh_values[sh_offsets[i] + 1 : sh_offsets[i + 1]]` with length
  `2 * ls[i] + 1`. Holds tesseral spherical-harmonic values for every
  site (including the differentiated one — the gradient site uses
  these for the product over non-differentiated sites).

- `sh_offsets::Vector{Int}` — cumulative offsets for `sh_values`. See
  [`EnergyWorkspace`](@ref) for the construction rule.

- `atom_grad_values::Vector{SVector{3, Float64}}` — gradient values for
  the differentiated site only. Length equals `2*l_atom + 1` where
  `l_atom = cbc.ls[atom_site_idx]`. Each entry is the 3-component
  vector `∂ᵢZlm_unsafe(l_atom, m, spin_direction)`. Stored as
  `SVector{3, Float64}` (not `Vector{Float64}`) to match the producer's
  return type and avoid a per-entry `convert`.

- `legendre_buf::Vector{Float64}` — same role as in
  [`EnergyWorkspace`](@ref): shared cache for the buffered
  `Zₗₘ_unsafe(l, m, uvec, buf)` and `∂ᵢZlm_unsafe(l, m, uvec, buf)`
  overloads. The two SH calls per site share this buffer because their
  buffer-size requirements `(l - |m| + 1)` are bounded by the same
  `max(ls) + 1`.

Not exported. Construct with `GradWorkspace()`.
"""
mutable struct GradWorkspace
	searched_pairs::Set{UInt}
	sh_values::Vector{Float64}
	sh_offsets::Vector{Int}
	atom_grad_values::Vector{SVector{3, Float64}}
	legendre_buf::Vector{Float64}

	GradWorkspace() = new(
		Set{UInt}(),
		Vector{Float64}(),
		Vector{Int}(),
		Vector{SVector{3, Float64}}(),
		Vector{Float64}(),
	)
end

# Ensure `legendre_buf` has at least `max(ls) + 1` entries — the upper bound
# required by `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe`'s buffered overloads
# (`length(buf) >= l - |m| + 1`, bounded above by `max(ls) + 1`).
@inline function _ensure_legendre_buf!(buf::Vector{Float64}, ls::Vector{Int})
	needed = maximum(ls) + 1
	if length(buf) < needed
		resize!(buf, needed)
	end
	return buf
end

# Size a flattened SH-value buffer to match a coupled-basis `ls` vector.
# After the call, `offsets` has length `N + 1` with `offsets[1] = 0` and
# `offsets[i + 1] = offsets[i] + 2*ls[i] + 1`, so site `i`'s slice is
# `buf[offsets[i] + 1 : offsets[i + 1]]`. `buf` is grown to at least the
# total size `offsets[N + 1]` and never shrunk, so capacity is reused
# across calls within the same workspace.
function _ensure_sh_buffer!(
	buf::Vector{Float64},
	offsets::Vector{Int},
	ls::Vector{Int},
)
	N = length(ls)
	resize!(offsets, N + 1)
	@inbounds offsets[1] = 0
	@inbounds for i in 1:N
		offsets[i + 1] = offsets[i] + 2 * ls[i] + 1
	end
	total = @inbounds offsets[N + 1]
	if length(buf) < total
		resize!(buf, total)
	end
	return buf
end

# `alpha` is a deprecated remnant of the old ElasticNet(alpha, lambda) API.
@inline function _atoms_hash_key(atoms_sorted::AbstractVector{Int})::UInt
	h = UInt(0x9e3779b97f4a7c15)
	@inbounds for a in atoms_sorted
		x = reinterpret(UInt, Int64(a))
		h ⊻= x + UInt(0x9e3779b97f4a7c15) + (h << 6) + (h >> 2)
	end
	return h ⊻ UInt(length(atoms_sorted))
end

function build_design_matrix_energy(
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
	spinconfig_list::AbstractVector{SpinConfig},
	symmetry::Symmetry,
)::Matrix{Float64}
	num_salcs = length(salc_list)  # Number of key groups
	num_spinconfigs = length(spinconfig_list)

	# construct design matrix A in Ax = b
	design_matrix = zeros(Float64, num_spinconfigs, num_salcs + 1)

	# set first column to 1 (reference_energy term)
	design_matrix[:, 1] .= 1.0

	# Parallel over key-group columns; each thread writes to a disjoint column.
	# Rank-erasing annotations on `key_group` / `cbc` are intentionally absent so
	# that Julia specializes `design_matrix_energy_element` on each element's
	# concrete `CoupledBasis_with_coefficient{R}` type via call-site dispatch.
	@threads for i = 1:num_salcs
		key_group = salc_list[i]
		n_C = length(key_group[1].atoms)  # Number of sites in the cluster
		scaling_factor = _cluster_scaling(n_C)
		# One workspace per @threads iteration. The Set/Vector buffers grow on
		# first use and are reused across all (j, cbc) calls in this column.
		ws = EnergyWorkspace()
		@inbounds for j in 1:num_spinconfigs
			# Sum contributions from all CoupledBasis_with_coefficient in this key group
			group_value = 0.0
			for cbc in key_group
				group_value += design_matrix_energy_element(
					cbc,
					spinconfig_list[j].spin_directions,
					symmetry,
					ws,
				)
			end
			design_matrix[j, i+1] = group_value * scaling_factor
		end
	end

	return design_matrix
end


"""
	design_matrix_energy_element(cbc, spin_directions, symmetry, ws) -> Float64

Compute one energy-design feature for a given CoupledBasis_with_coefficient and spin directions.

# Description
- Contracts coupled angular momentum basis tensor with spherical harmonics over atoms following symmetry translations.
- Equivalent to one column entry (excluding bias) in the energy design matrix.

# Arguments
- `cbc::CoupledBases.CoupledBasis_with_coefficient`: CoupledBasis_with_coefficient object
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)
- `symmetry::Symmetry`: Symmetry information of the structure
- `ws::EnergyWorkspace`: Reusable scratch state (cleared on entry). Each
  thread should own one; see `build_design_matrix_energy` for the
  recommended allocation pattern.

# Returns
- `Float64`: Feature value for the CoupledBasis_with_coefficient
"""
function design_matrix_energy_element(
	cbc::CoupledBases.CoupledBasis_with_coefficient{R},
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
	ws::EnergyWorkspace,
)::Float64 where {R}
	# N = number of sites; R = tensor rank = N + 1 (compile-time constant).
	result::Float64 = 0.0
	dims_t = ntuple(i -> 2 * cbc.ls[i] + 1, Val(R - 1))
	Mf_size = size(cbc.coeff_tensor, R)

	# cbc.ls is fixed in this function, so de-dup by sorted translated atoms only.
	empty!(ws.searched_pairs)
	searched_pairs = ws.searched_pairs
	_ensure_sh_buffer!(ws.sh_values, ws.sh_offsets, cbc.ls)
	sh_values = ws.sh_values
	sh_offsets = ws.sh_offsets
	_ensure_legendre_buf!(ws.legendre_buf, cbc.ls)
	legendre_buf = ws.legendre_buf
	translated_atoms = MVector{R - 1, Int}(undef)
	atoms_sorted_buf = MVector{R - 1, Int}(undef)

	# These quantities depend only on cbc.ls (constant within the call), so hoist
	# them out of the translation loop instead of rebuilding per `itrans`.
	other_dims_t = ntuple(i -> dims_t[i], Val(R - 2))
	other_site_indices = CartesianIndices(other_dims_t)
	idx_buf = MVector{R - 1, Int}(undef)
	base_last = @inbounds sh_offsets[R - 1]

	@inbounds for itrans in symmetry.symnum_translation
		# Translate atoms
		for site_idx in 1:(R - 1)
			translated_atoms[site_idx] = symmetry.map_sym[cbc.atoms[site_idx], itrans]
		end
		# Sort atoms for comparison, but keep ls in original order
		copyto!(atoms_sorted_buf, translated_atoms)
		sort!(atoms_sorted_buf)
		atoms_key = _atoms_hash_key(atoms_sorted_buf)
		if atoms_key in searched_pairs
			continue
		end
		push!(searched_pairs, atoms_key)


		# Compute spherical harmonics for each site into the workspace buffer.
		# `Zₗₘ_unsafe(l, m, uvec, legendre_buf)` reuses the Legendre cache
		# rather than reallocating it on every (l, m, atom).
		for site_idx in 1:(R - 1)
			atom = translated_atoms[site_idx]
			l = cbc.ls[site_idx]
			base = sh_offsets[site_idx]
			for m_idx in 1:(2*l+1)
				# Convert tesseral index to m value: m = m_idx - l - 1
				m = m_idx - l - 1
				sh_values[base + m_idx] =
					@views Zₗₘ_unsafe(l, m, spin_directions[:, atom], legendre_buf)
			end
		end

		# Contract coeff_tensor with spherical harmonics
		# coeff_tensor has shape (d1, d2, ..., dN, Mf_size) where di = 2*li + 1.
		# Reuse product over first N-1 sites and iterate the last-site index
		# separately; `idx_buf::MVector{R-1, Int}` makes the splat indexing
		# `coeff_tensor[idx_buf..., mf_idx]` statically resolvable.
		tensor_result = 0.0

		# Iterate over all Mf values
		for mf_idx in 1:Mf_size
			mf_contribution = 0.0

			for other_tuple in other_site_indices
				product_other = 1.0
				for site_idx in 1:(R - 2)
					m_idx_other = other_tuple.I[site_idx]
					idx_buf[site_idx] = m_idx_other
					product_other *= sh_values[sh_offsets[site_idx] + m_idx_other]
				end
				for m_idx_last in 1:dims_t[R - 1]
					idx_buf[R - 1] = m_idx_last
					mf_contribution +=
						cbc.coeff_tensor[idx_buf..., mf_idx] *
						(product_other * sh_values[base_last + m_idx_last])
				end
			end

			tensor_result += cbc.coefficient[mf_idx] * mf_contribution
		end

		result += tensor_result * cbc.multiplicity
	end

	return result
end

"""
	build_design_matrix_torque(salc_list, spinconfig_list, num_atoms, symmetry) -> Matrix{Float64}

Build the torque design matrix used for regression.

# Description
- For each spin configuration, constructs a block of size (3·num_atoms × num_salcs)
  whose rows are per-atom XYZ components of `cross(spin_dir, ∇ₑu)`.
- Blocks are vertically concatenated across configurations.

# Arguments
- `salc_list`: List of CoupledBasis_with_coefficient
- `spinconfig_list`: Vector of spin configurations
- `num_atoms`: Number of atoms in the structure
- `symmetry`: Symmetry information

# Returns
- `Matrix{Float64}`: Torque design matrix
"""
function build_design_matrix_torque(
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
	spinconfig_list::AbstractVector{SpinConfig},
	num_atoms::Integer,
	symmetry::Symmetry,
)::Matrix{Float64}
	num_salcs = length(salc_list)  # Number of key groups
	num_spinconfigs = length(spinconfig_list)
	scaling_factors = Vector{Float64}(undef, num_salcs)
	for (salc_idx, key_group) in enumerate(salc_list)
		n_C = length(key_group[1].atoms)  # Number of sites in the cluster
		scaling_factors[salc_idx] = _cluster_scaling(n_C)
	end

	# Preallocate the full design matrix and let each thread write into its
	# disjoint row block (sc_idx → rows [block_size*(sc_idx-1)+1 : block_size*sc_idx]).
	# Avoids per-thread block allocation and the final vcat copy.
	block_size = 3 * num_atoms
	design_matrix = Matrix{Float64}(undef, num_spinconfigs * block_size, num_salcs)

	@threads for sc_idx in 1:num_spinconfigs
		spinconfig = spinconfig_list[sc_idx]
		row_offset = block_size * (sc_idx - 1)
		grad_u_buf = MVector{3, Float64}(0.0, 0.0, 0.0)
		# One workspace per @threads iteration; reused across all (iatom,
		# salc_idx, cbc) calls in this row block.
		ws = GradWorkspace()
		@inbounds for iatom in 1:num_atoms
			@views dir_iatom = spinconfig.spin_directions[:, iatom]
			dir_iatom_svec = SVector{3, Float64}(dir_iatom)
			@inbounds for (salc_idx, key_group) in enumerate(salc_list)
				# Sum contributions from all CoupledBasis_with_coefficient in this key group
				group_grad = MVector{3, Float64}(0.0, 0.0, 0.0)
				for cbc in key_group
					calc_∇ₑu!(
						grad_u_buf,
						cbc,
						iatom,
						spinconfig.spin_directions,
						symmetry,
						ws,
					)
					group_grad .+= grad_u_buf
				end
				scaling_factor = scaling_factors[salc_idx]
				torque_vec = cross(dir_iatom_svec, SVector{3, Float64}(group_grad)) * scaling_factor
				row_base = row_offset + 3 * (iatom - 1)
				design_matrix[row_base + 1, salc_idx] = torque_vec[1]
				design_matrix[row_base + 2, salc_idx] = torque_vec[2]
				design_matrix[row_base + 3, salc_idx] = torque_vec[3]
			end
		end
	end

	return design_matrix
end


"""
	calc_∇ₑu(cbc, atom, spin_directions, symmetry) -> Vector{Float64}

Compute the gradient of the coupled angular momentum basis for a CoupledBasis_with_coefficient
with respect to the spin direction of a specific atom.

# Description
- Returns a 3-vector corresponding to (∂/∂x, ∂/∂y, ∂/∂z) components.
- Applies symmetry translations before accumulation.

# Arguments
- `cbc::CoupledBases.CoupledBasis_with_coefficient`: CoupledBasis_with_coefficient object
- `atom::Integer`: Target atom index (1-based)
- `spin_directions::AbstractMatrix{<:Real}`: 3×N spin directions
- `symmetry::Symmetry`: Symmetry information

# Returns
- `Vector{Float64}`: Gradient vector (length 3)
"""
function calc_∇ₑu!(
	result::MVector{3, Float64},
	cbc::CoupledBases.CoupledBasis_with_coefficient{R},
	atom::Integer,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
	ws::GradWorkspace,
) where {R}
	# N = number of sites; R = tensor rank = N + 1 (compile-time constant).
	result[1] = 0.0
	result[2] = 0.0
	result[3] = 0.0
	dims_t = ntuple(i -> 2 * cbc.ls[i] + 1, Val(R - 1))
	Mf_size = size(cbc.coeff_tensor, R)
	translated_atoms = MVector{R - 1, Int}(undef)
	atoms_sorted_buf = MVector{R - 1, Int}(undef)

	# Hoist scratch buffers out of the translation loop. `other_sites_buf` and
	# `other_dims_buf` get refilled per iteration based on `atom_site_idx`, but
	# the storage is allocated once up front.
	idx_buf = MVector{R - 1, Int}(undef)
	other_sites_buf = MVector{R - 2, Int}(undef)
	other_dims_buf = MVector{R - 2, Int}(undef)

	# cbc.ls is fixed in this function, so de-dup by sorted translated atoms only.
	empty!(ws.searched_pairs)
	searched_pairs = ws.searched_pairs
	_ensure_sh_buffer!(ws.sh_values, ws.sh_offsets, cbc.ls)
	sh_values = ws.sh_values
	sh_offsets = ws.sh_offsets
	_ensure_legendre_buf!(ws.legendre_buf, cbc.ls)
	legendre_buf = ws.legendre_buf

	@inbounds for itrans in symmetry.symnum_translation
		# Translate atoms and identify the differentiated site index.
		atom_site_idx = 0
		for site_idx in 1:(R - 1)
			ta = symmetry.map_sym[cbc.atoms[site_idx], itrans]
			translated_atoms[site_idx] = ta
			if ta == atom
				atom_site_idx = site_idx
			end
		end
		# Sort buffer for de-dup, preserving translated_atoms order for computation.
		copyto!(atoms_sorted_buf, translated_atoms)
		sort!(atoms_sorted_buf)
		atoms_key = _atoms_hash_key(atoms_sorted_buf)
		if atoms_key in searched_pairs
			continue
		end
		push!(searched_pairs, atoms_key)

		if atom_site_idx == 0
			continue
		end

		# Compute spherical harmonics and derivatives for the target site.
		# `sh_values` is shared scratch; `atom_grad_values` is sized to the
		# differentiated site's `2*l_atom+1` m-values.
		l_atom = cbc.ls[atom_site_idx]
		resize!(ws.atom_grad_values, 2 * l_atom + 1)
		atom_grad_values = ws.atom_grad_values

		# Both buffered overloads use `legendre_buf` as a write-only scratch:
		# the Legendre cache they need is fully recomputed on every call from
		# the `(l, m, z)` arguments, so the prior contents are not consulted.
		# Calling `Zₗₘ_unsafe` and then `∂ᵢZlm_unsafe` against the same buffer
		# (the `site_idx == atom_site_idx` case below) is therefore safe; the
		# second call simply overwrites the first call's cache state.
		for site_idx in 1:(R - 1)
			translated_atom = translated_atoms[site_idx]
			l = cbc.ls[site_idx]
			base = sh_offsets[site_idx]
			for m_idx in 1:(2*l+1)
				# Convert tesseral index to m value: m = m_idx - l - 1
				m = m_idx - l - 1

				sh_values[base + m_idx] =
					@views Zₗₘ_unsafe(l, m, spin_directions[:, translated_atom], legendre_buf)
				if site_idx == atom_site_idx
					# Keep gradients only for the differentiated site.
					atom_grad_values[m_idx] =
						@views ∂ᵢZlm_unsafe(l, m, spin_directions[:, translated_atom], legendre_buf)
				end
			end
		end

		# Contract coeff_tensor with spherical harmonics and gradients.
		# coeff_tensor has shape (d1, d2, ..., dN, Mf_size) where di = 2*li + 1.
		grad_result = MVector{3, Float64}(0.0, 0.0, 0.0)
		# Refill the other_sites/other_dims buffers for this `atom_site_idx`.
		ki = 0
		for s in 1:(R - 1)
			if s != atom_site_idx
				ki += 1
				other_sites_buf[ki] = s
				other_dims_buf[ki] = dims_t[s]
			end
		end
		# `Tuple(::MVector{R-2,Int})` is `NTuple{R-2,Int}`, so the
		# resulting CartesianIndices has statically known rank.
		other_site_indices = CartesianIndices(Tuple(other_dims_buf))

		# Iterate over all Mf values
		for mf_idx in 1:Mf_size
			mf_grad_contribution = MVector{3, Float64}(0.0, 0.0, 0.0)

			# Reuse product over non-differentiated sites for each m on target site.
			for other_tuple in other_site_indices
				product_other = 1.0
				for k in 1:(R - 2)
					site_idx = other_sites_buf[k]
					m_idx_other = other_tuple.I[k]
					idx_buf[site_idx] = m_idx_other
					product_other *= sh_values[sh_offsets[site_idx] + m_idx_other]
				end

				for m_idx_atom in 1:(2*l_atom+1)
					idx_buf[atom_site_idx] = m_idx_atom
					grad_atom = atom_grad_values[m_idx_atom]
					coeff_val = cbc.coeff_tensor[idx_buf..., mf_idx]
					mf_grad_contribution .+= coeff_val * product_other .* grad_atom
				end
			end

			grad_result .+= cbc.coefficient[mf_idx] .* mf_grad_contribution
		end

		result .+= grad_result .* cbc.multiplicity
	end

	return result
end

function calc_∇ₑu(
	cbc::CoupledBases.CoupledBasis_with_coefficient,
	atom::Integer,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)::Vector{Float64}
	result = MVector{3, Float64}(0.0, 0.0, 0.0)
	ws = GradWorkspace()
	calc_∇ₑu!(result, cbc, atom, spin_directions, symmetry, ws)
	return Vector{Float64}(result)
end


"""
	_predict_energy(j0, jphi, salc_list, symmetry, spin_directions) -> Float64

Internal SCE energy evaluation from raw pieces. The `SCEModel`-typed
`predict_energy` method (in `Magesty`) unpacks a model and delegates
here. Returns `j0 + Σ jphi · ϕ(spin)` in eV.
"""
function _predict_energy(
	j0::Float64,
	jphi::AbstractVector{<:Real},
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
	symmetry::Symmetry,
	spin_directions::AbstractMatrix{<:Real},
)::Float64
	num_salcs = length(salc_list)
	design_vector = Vector{Float64}(undef, num_salcs)
	# Single-threaded inference path; allocate one workspace and reuse.
	ws = EnergyWorkspace()

	for i = 1:num_salcs
		key_group = salc_list[i]
		n_C = length(key_group[1].atoms)
		scaling_factor = _cluster_scaling(n_C)
		group_value = 0.0
		for cbc in key_group
			group_value += design_matrix_energy_element(cbc, spin_directions, symmetry, ws)
		end
		design_vector[i] = group_value * scaling_factor
	end

	return dot(design_vector, jphi) + j0
end

"""
	_predict_torque(jphi, salc_list, symmetry, spin_directions) -> Matrix{Float64}

Internal SCE per-atom torque evaluation from raw pieces. The
`SCEModel`-typed `predict_torque` method (in `Magesty`) unpacks a model
and delegates here. Returns a `3×N` torque matrix (eV); column `i` is
the torque on atom `i`. The sign convention matches
`build_design_matrix_torque` (`cross(spin, ∇ₑu)`).

Throws `ArgumentError` if `spin_directions` is not a `3×N` matrix.
"""
function _predict_torque(
	jphi::AbstractVector{<:Real},
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
	symmetry::Symmetry,
	spin_directions::AbstractMatrix{<:Real},
)::Matrix{Float64}
	if size(spin_directions, 1) != 3
		throw(ArgumentError("spin_directions must be a 3xN matrix"))
	end
	num_atoms = size(spin_directions, 2)
	torque = zeros(Float64, 3, num_atoms)
	grad_u_buf = MVector{3, Float64}(0.0, 0.0, 0.0)
	# Single-threaded inference path; allocate one workspace and reuse.
	ws = GradWorkspace()

	@inbounds for iatom = 1:num_atoms
		dir_iatom = SVector{3, Float64}(@view spin_directions[:, iatom])
		for (salc_idx, key_group) in enumerate(salc_list)
			group_grad = MVector{3, Float64}(0.0, 0.0, 0.0)
			for cbc in key_group
				calc_∇ₑu!(grad_u_buf, cbc, iatom, spin_directions, symmetry, ws)
				group_grad .+= grad_u_buf
			end
			n_C = length(key_group[1].atoms)
			scaling_factor = _cluster_scaling(n_C)
			torque[:, iatom] .+=
				cross(dir_iatom, SVector{3, Float64}(group_grad)) .*
				scaling_factor .* jphi[salc_idx]
		end
	end

	return torque
end


"""
	assemble_weighted_problem(design_matrix_energy, design_matrix_torque,
	                          observed_energy_list, observed_torque, weight)
	    -> (X::Matrix{Float64}, y::Vector{Float64}, bias_col::Int)

Build the augmented `(X, y)` linear system that all estimators solve.

`observed_torque` may be passed either as a vector of `3×n_atoms`
matrices (one per config) or as the already-flattened
`3 * n_atoms * n_config` vector; the matrix form is flattened internally
and both produce identical output.

Energy rows are scaled by `√((1 - weight) / n_E)` and torque rows
(after flattening) by `√(weight / n_T)`, where `n_E` is the number of
energy samples (configs) and `n_T = 3 * n_atoms * n_E` the number of
torque components. This per-sample normalization makes the augmented
least-squares objective equal to `(1 - weight) * MSE_energy + weight *
MSE_torque`, so `weight` selects a convex combination of the two mean
squared errors independent of how many rows each block contributes
(mirroring the energy/force loss convention used by ML interatomic
potentials). The bias column of the energy block is reset to `1.0`
after scaling so that the bias coefficient is not absorbed into the
weight factor. A zero column is prepended to the torque block so the
two blocks share `bias_col = 1`.

# Arguments
- `design_matrix_energy`: Energy design matrix (bias column at column 1).
- `design_matrix_torque`: Torque design matrix (no bias column).
- `observed_energy_list`: Observed energies.
- `observed_torque`: Observed torques, either as `3×n_atoms` matrices per
  config or as the flattened `3 * n_atoms * n_config` vector.
- `weight`: Trade-off; energy weight = `1 - weight`, torque weight = `weight`.

# Returns
- `X`: Augmented design matrix (energy rows stacked above flattened torque rows).
- `y`: Augmented observation vector.
- `bias_col`: Column index of the bias term (always `1`).
"""
function assemble_weighted_problem(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_torque::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_torque_list::AbstractVector{<:AbstractMatrix{<:Real}},
	weight::Real,
)
	return assemble_weighted_problem(
		design_matrix_energy,
		design_matrix_torque,
		observed_energy_list,
		convert(Vector{Float64}, vcat(vec.(observed_torque_list)...)),
		weight,
	)
end

function assemble_weighted_problem(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_torque::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_torque_flattened::AbstractVector{<:Real},
	weight::Real,
)
	# weight parameters
	w_e = 1 - weight
	w_m = weight

	# Per-sample MSE normalization. The objective we want to minimize is
	#   L = (1 - weight) * MSE_energy + weight * MSE_torque
	#     = (1 - weight) / n_E * Σ ε²  +  weight / n_T * Σ τ²
	# but a linear least-squares solver minimizes only the plain residual
	# sum of squares ‖X β - y‖². Scaling a block's rows by a factor c
	# turns its residual sum of squares into c² * Σ(residual²), so to give
	# a block the weight w / n we scale its rows (and the matching entries
	# of y) by √(w / n). This is the standard weighted-least-squares
	# whitening transform W^(1/2) applied to (X, y); the √ appears because
	# the solver squares the residuals.
	# n_E = number of energy samples (configs),
	# n_T = 3 * n_atoms * n_E = number of torque components.
	n_E = length(observed_energy_list)
	n_T = length(observed_torque_flattened)
	scale_e = sqrt(w_e / n_E)
	scale_m = sqrt(w_m / n_T)

	# Normalize the design matrices
	normalized_design_matrix_energy =
		design_matrix_energy .* scale_e
	normalized_design_matrix_energy[:, 1] .= 1.0
	normalized_design_matrix_torque =
		design_matrix_torque .* scale_m

	# Also normalize the observed vectors
	normalized_observed_energy_list =
		observed_energy_list .* scale_e
	normalized_observed_torque_flattened =
		observed_torque_flattened .* scale_m

	# Add 0 bias column to the torque design matrix
	# to align with the energy design matrix
	with_bias_design_matrix_torque =
		hcat(zeros(size(normalized_design_matrix_torque, 1)), normalized_design_matrix_torque)

	# Construct the augmented design matrix
	X = vcat(
		normalized_design_matrix_energy,
		with_bias_design_matrix_torque,
	)
	y = vcat(
		normalized_observed_energy_list,
		normalized_observed_torque_flattened,
	)

	return X, y, 1
end

"""
	extract_j0_jphi(j_values, design_matrix_energy, observed_energy_list)
	    -> (j0::Float64, jphi::Vector{Float64})

Split the augmented coefficient vector into the bias `j0` and the SCE
coefficients `jphi`. The bias is re-estimated from the unscaled energy
residual `mean(ye .- Xe[:, 2:end] * jphi)` to undo the √weight scaling
that affected the augmented bias column.

# Arguments
- `j_values`: Augmented coefficient vector returned by `solve_coefficients`.
- `design_matrix_energy`: Unscaled energy design matrix (bias column at 1).
- `observed_energy_list`: Unscaled observed energies.

# Returns
- `(j0, jphi)`: Bias and SCE coefficients.
"""
function extract_j0_jphi(
	j_values::AbstractVector{<:Real},
	design_matrix_energy::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
)
	jphi = j_values[2:end]
	j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)
	return j0, jphi
end

"""
	solve_coefficients(estimator::AbstractEstimator, X, y; bias_col=1)
	    -> Vector{Float64}

Solve the augmented linear system `(X, y)` for the coefficient vector.
The element at index `bias_col` is the bias (j0) and MUST be excluded
from any regularization penalty by the estimator's method.

Each `AbstractEstimator` subtype defines exactly one method of this
function. Unknown estimator types raise `MethodError`.
"""
function solve_coefficients end

function solve_coefficients(
	::OLS,
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real};
	bias_col::Int = 1,
)
	return X \ y
end

"""
	solve_coefficients(estimator::Ridge, X, y; bias_col=1) -> Vector{Float64}

L2-regularized solve. Uses `MultivariateStats.ridge` with a per-column
penalty vector that zeros out the bias column at `bias_col`. When
`estimator.lambda ≈ 0` falls back to `X \\ y` to avoid building the
penalty vector for an effectively unregularized problem.
"""
function solve_coefficients(
	e::Ridge,
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real};
	bias_col::Int = 1,
)
	if e.lambda ≈ 0.0
		return X \ y
	end
	lambda_vec = fill(e.lambda, size(X, 2))
	lambda_vec[bias_col] = 0.0  # exclude bias term from regularization
	return ridge(X, y, lambda_vec; bias = false)
end

function calc_rmse(list1::AbstractVector{<:Real}, list2::AbstractVector{<:Real})::Float64
	# Calculate the Root Mean Square Error (RMSE) between two lists
	if length(list1) != length(list2)
		throw(ArgumentError("The lengths of the two lists must be equal."))
	end
	return sqrt(mean((list1 .- list2) .^ 2))
end

function calc_r2score(
	observed_list::AbstractVector{<:Real},
	predicted_list::AbstractVector{<:Real},
)::Float64
	# R² = 1 - (SS_res / SS_tot)
	# SS_res = Σ(y_observed - y_predicted)²
	# SS_tot = Σ(y_observed - y_mean)²
	ss_res = sum((observed_list .- predicted_list) .^ 2)
	ss_tot = sum((observed_list .- mean(observed_list)) .^ 2)
	return 1 - ss_res / ss_tot
end

end

