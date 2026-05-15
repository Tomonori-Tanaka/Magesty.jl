"""
	Optimize.jl

This module contains functions for optimizing the SCE coefficients.
"""
module Optimize

using Base.Threads
using LinearAlgebra
using Printf
using MultivariateStats
using Statistics
using StaticArrays
using ..MySphericalHarmonics
using ..ConfigParser
using ..Structures
using ..Symmetries
using ..Clusters
using ..Basis
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
	salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
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
	@threads for i = 1:num_salcs
		key_group::Vector{Basis.CoupledBasis_with_coefficient} = salc_list[i]
		n_C = length(key_group[1].atoms)  # Number of sites in the cluster
		scaling_factor = _cluster_scaling(n_C)
		@inbounds for j in 1:num_spinconfigs
			# Sum contributions from all CoupledBasis_with_coefficient in this key group
			group_value = 0.0
			for cbc::Basis.CoupledBasis_with_coefficient in key_group
				group_value += design_matrix_energy_element(
					cbc,
					spinconfig_list[j].spin_directions,
					symmetry,
				)
			end
			design_matrix[j, i+1] = group_value * scaling_factor
		end
	end

	return design_matrix
end


"""
	design_matrix_energy_element(cbc, spin_directions, symmetry) -> Float64

Compute one energy-design feature for a given CoupledBasis_with_coefficient and spin directions.

# Description
- Contracts coupled angular momentum basis tensor with spherical harmonics over atoms following symmetry translations.
- Equivalent to one column entry (excluding bias) in the energy design matrix.

# Arguments
- `cbc::Basis.CoupledBasis_with_coefficient`: CoupledBasis_with_coefficient object
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)
- `symmetry::Symmetry`: Symmetry information of the structure

# Returns
- `Float64`: Feature value for the CoupledBasis_with_coefficient
"""
function design_matrix_energy_element(
	cbc::Basis.CoupledBasis_with_coefficient,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)::Float64
	result::Float64 = 0.0
	N = length(cbc.atoms)
	dims = [2*l + 1 for l in cbc.ls]
	site_indices = CartesianIndices(Tuple(dims))
	Mf_size = size(cbc.coeff_tensor, N+1)

	# cbc.ls is fixed in this function, so de-dup by sorted translated atoms only.
	searched_pairs = Set{UInt}()
	translated_atoms = Vector{Int}(undef, N)
	atoms_sorted_buf = Vector{Int}(undef, N)

	# These quantities depend only on cbc.ls (constant within the call), so hoist
	# them out of the translation loop instead of rebuilding per `itrans`.
	last_site = N
	other_sites = 1:(N-1)
	other_dims = dims[other_sites]
	other_site_indices = CartesianIndices(Tuple(other_dims))
	idx_buf = Vector{Int}(undef, N)

	@inbounds for itrans in symmetry.symnum_translation
		# Translate atoms
		for site_idx in 1:N
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


		# Compute spherical harmonics for each site
		sh_values = Vector{Vector{Float64}}(undef, N)
		for (site_idx, atom) in enumerate(translated_atoms)
			l = cbc.ls[site_idx]
			sh_values[site_idx] = Vector{Float64}(undef, 2*l+1)
			for m_idx in 1:(2*l+1)
				# Convert tesseral index to m value: m = m_idx - l - 1
				m = m_idx - l - 1
				sh_values[site_idx][m_idx] = @views Zₗₘ_unsafe(l, m, spin_directions[:, atom])
			end
		end

		# Contract coeff_tensor with spherical harmonics
		# coeff_tensor has shape (d1, d2, ..., dN, Mf_size)
		# where di = 2*li + 1
		tensor_result = 0.0

		# Iterate over all Mf values
		for mf_idx in 1:Mf_size
			mf_contribution = 0.0

			# Reuse product over first N-1 sites and iterate last-site index separately.
			for other_tuple in other_site_indices
				product_other = 1.0
				for site_idx in other_sites
					m_idx_other = other_tuple.I[site_idx]
					idx_buf[site_idx] = m_idx_other
					product_other *= sh_values[site_idx][m_idx_other]
				end
				for m_idx_last in 1:dims[last_site]
					idx_buf[last_site] = m_idx_last
					mf_contribution +=
						cbc.coeff_tensor[idx_buf..., mf_idx] *
						(product_other * sh_values[last_site][m_idx_last])
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
	salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
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
- `cbc::Basis.CoupledBasis_with_coefficient`: CoupledBasis_with_coefficient object
- `atom::Integer`: Target atom index (1-based)
- `spin_directions::AbstractMatrix{<:Real}`: 3×N spin directions
- `symmetry::Symmetry`: Symmetry information

# Returns
- `Vector{Float64}`: Gradient vector (length 3)
"""
function calc_∇ₑu!(
	result::MVector{3, Float64},
	cbc::Basis.CoupledBasis_with_coefficient,
	atom::Integer,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)
	result[1] = 0.0
	result[2] = 0.0
	result[3] = 0.0
	N = length(cbc.atoms)
	dims = [2 * l + 1 for l in cbc.ls]
	Mf_size = size(cbc.coeff_tensor, N + 1)
	translated_atoms = Vector{Int}(undef, N)
	atoms_sorted_buf = Vector{Int}(undef, N)

	# Hoist scratch buffers out of the translation loop. `other_sites_buf` and
	# `other_dims_buf` get refilled per iteration based on `atom_site_idx`, but
	# the storage is allocated once up front.
	idx_buf = Vector{Int}(undef, N)
	other_sites_buf = Vector{Int}(undef, N - 1)
	other_dims_buf = Vector{Int}(undef, N - 1)

	# cbc.ls is fixed in this function, so de-dup by sorted translated atoms only.
	searched_pairs = Set{UInt}()

	@inbounds for itrans in symmetry.symnum_translation
		# Translate atoms and identify the differentiated site index.
		atom_site_idx = 0
		for site_idx in 1:N
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
		sh_values = Vector{Vector{Float64}}(undef, N)
		l_atom = cbc.ls[atom_site_idx]
		atom_grad_values = Vector{Vector{Float64}}(undef, 2 * l_atom + 1)

		for (site_idx, translated_atom) in enumerate(translated_atoms)
			l = cbc.ls[site_idx]
			sh_values[site_idx] = Vector{Float64}(undef, 2*l+1)
			for m_idx in 1:(2*l+1)
				# Convert tesseral index to m value: m = m_idx - l - 1
				m = m_idx - l - 1

				# Use regular spherical harmonic for all sites
				sh_values[site_idx][m_idx] =
					@views Zₗₘ_unsafe(l, m, spin_directions[:, translated_atom])
				if site_idx == atom_site_idx
					# Keep gradients only for the differentiated site.
					atom_grad_values[m_idx] =
						@views ∂ᵢZlm_unsafe(l, m, spin_directions[:, translated_atom])
				end
			end
		end

		# Contract coeff_tensor with spherical harmonics and gradients
		# coeff_tensor has shape (d1, d2, ..., dN, Mf_size)
		# where di = 2*li + 1
		grad_result = MVector{3, Float64}(0.0, 0.0, 0.0)
		# Refill the other_sites/other_dims buffers for this `atom_site_idx`.
		ki = 0
		for s in 1:N
			if s != atom_site_idx
				ki += 1
				other_sites_buf[ki] = s
				other_dims_buf[ki] = dims[s]
			end
		end
		other_site_indices = CartesianIndices(Tuple(other_dims_buf))

		# Iterate over all Mf values
		for mf_idx in 1:Mf_size
			mf_grad_contribution = MVector{3, Float64}(0.0, 0.0, 0.0)

			# Reuse product over non-differentiated sites for each m on target site.
			for other_tuple in other_site_indices
				product_other = 1.0
				for k in 1:(N-1)
					site_idx = other_sites_buf[k]
					m_idx_other = other_tuple.I[k]
					idx_buf[site_idx] = m_idx_other
					product_other *= sh_values[site_idx][m_idx_other]
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
	cbc::Basis.CoupledBasis_with_coefficient,
	atom::Integer,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)::Vector{Float64}
	result = MVector{3, Float64}(0.0, 0.0, 0.0)
	calc_∇ₑu!(result, cbc, atom, spin_directions, symmetry)
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
	salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
	symmetry::Symmetry,
	spin_directions::AbstractMatrix{<:Real},
)::Float64
	num_salcs = length(salc_list)
	design_vector = Vector{Float64}(undef, num_salcs)

	for i = 1:num_salcs
		key_group = salc_list[i]
		n_C = length(key_group[1].atoms)
		scaling_factor = _cluster_scaling(n_C)
		group_value = 0.0
		for cbc in key_group
			group_value += design_matrix_energy_element(cbc, spin_directions, symmetry)
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
	salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
	symmetry::Symmetry,
	spin_directions::AbstractMatrix{<:Real},
)::Matrix{Float64}
	if size(spin_directions, 1) != 3
		throw(ArgumentError("spin_directions must be a 3xN matrix"))
	end
	num_atoms = size(spin_directions, 2)
	torque = zeros(Float64, 3, num_atoms)
	grad_u_buf = MVector{3, Float64}(0.0, 0.0, 0.0)

	@inbounds for iatom = 1:num_atoms
		dir_iatom = SVector{3, Float64}(@view spin_directions[:, iatom])
		for (salc_idx, key_group) in enumerate(salc_list)
			group_grad = MVector{3, Float64}(0.0, 0.0, 0.0)
			for cbc in key_group
				calc_∇ₑu!(grad_u_buf, cbc, iatom, spin_directions, symmetry)
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

