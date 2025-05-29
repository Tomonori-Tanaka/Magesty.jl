"""
	Optimize.jl

This module contains functions for optimizing the SCE coefficients.
"""
module Optimize

using Base.Threads
using LinearAlgebra
using GLMNet
using Optim
using Printf
using StatsBase
using Statistics
using StaticArrays
using ..MySphericalHarmonics
using ..SALCs
using ..AtomicIndices
using ..ConfigParser
using ..Structures
using ..Symmetries
using ..Clusters
using ..BasisSets
using ..SpinConfigs

export Optimizer

struct Optimizer
	spinconfig_list::Vector{SpinConfig}
	reference_energy::Float64
	SCE::Vector{Float64}
	predicted_energy_list::Vector{Float64}
	observed_energy_list::Vector{Float64}
	predicted_magfield_vertical_flattened_list::Vector{Float64}
	observed_magfield_vertical_flattened_list::Vector{Float64}
	elapsed_time::Float64  # Time taken to create the optimizer in seconds
end

function Optimizer(
	structure::Structure,
	symmetry::Symmetry,
	basisset::BasisSet,
	alpha::Real,
	weight::Real,
	spinconfig_list::AbstractVector{SpinConfig},
)
	# Start timing
	start_time = time_ns()

	# construct design matrix for energy and magfield_vertical
	design_matrix_energy = construct_design_matrix_energy(
		basisset.salc_list,
		spinconfig_list,
		symmetry,
	)
	design_matrix_magfield_vertical = construct_design_matrix_magfield_vertical(
		basisset.salc_list,
		spinconfig_list,
		structure.supercell.num_atoms,
		symmetry,
	)

	# construct observed_energy_list and observed_magfield_vertical_list
	observed_energy_list = [spinconfig.energy for spinconfig in spinconfig_list]
	observed_magfield_vertical_list =
		Vector{Vector{Float64}}(undef, length(spinconfig_list))
	@threads for i in eachindex(spinconfig_list)
		observed_magfield_vertical_list[i] =
			calc_magfield_vertical_list_of_spinconfig(
				spinconfig_list[i],
				structure.supercell.num_atoms,
			)
	end

	j0, jphi = ridge_regression(
		design_matrix_energy,
		design_matrix_magfield_vertical,
		observed_energy_list,
		observed_magfield_vertical_list,
		alpha,
		weight,
	)


	predicted_energy_list = design_matrix_energy[:, 2:end] * jphi .+ j0
	predicted_magfield_vertical_flattened_list = design_matrix_magfield_vertical * jphi
	observed_magfield_vertical_flattened_list =
		-1 * vcat(observed_magfield_vertical_list...)
	# End timing
	elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds

	return Optimizer(
		spinconfig_list,
		j0,
		jphi,
		predicted_energy_list,
		observed_energy_list,
		predicted_magfield_vertical_flattened_list,
		observed_magfield_vertical_flattened_list,
		elapsed_time,
	)
end

function Optimizer(
	structure::Structure,
	symmetry::Symmetry,
	basisset::BasisSet,
	config::Config4Optimize,
)
	return Optimizer(structure, symmetry, basisset, config.alpha, config.weight, config.datafile)
end

function Optimizer(	structure::Structure,
	symmetry::Symmetry,
	basisset::BasisSet,
	alpha::Real,
	weight::Real,
	datafile::AbstractString,
)
	# read datafile
	spinconfig_list = read_embset(datafile, structure.supercell.num_atoms)

	return Optimizer(structure, symmetry, basisset, alpha, weight, spinconfig_list)
end

"""
	construct_design_matrix_energy(salc_list, spinconfig_list, symmetry) -> Matrix{Float64}

Construct the design matrix for energy prediction.

# Arguments
- `salc_list::AbstractVector{SALC}`: List of SALC objects
- `spinconfig_list:AbstractVector{SpinConfig}`: Vector of SpinConfig objects
- `symmetry::Symmetry`: Symmetry information

# Returns
- `Matrix{Float64}`: Design matrix for energy prediction

"""
function construct_design_matrix_energy(
	salc_list::AbstractVector{SALC},
	spinconfig_list::AbstractVector{SpinConfig},
	symmetry::Symmetry,
)::Matrix{Float64}
	num_salcs = length(salc_list)
	num_spinconfigs = length(spinconfig_list)

	# construct design matrix A in Ax = b
	design_matrix = zeros(Float64, num_spinconfigs, num_salcs + 1)
	initialize_check = falses(num_spinconfigs, num_salcs + 1)

	# set first column to 1 (reference_energy term)
	design_matrix[:, 1] .= 1.0
	initialize_check[:, 1] .= true

	@threads for i in 1:num_salcs
		for j in 1:num_spinconfigs
			design_matrix[j, i+1] =
				calc_X_element_energy(
					salc_list[i],
					spinconfig_list[j].spin_directions,
					symmetry,
				)
			initialize_check[j, i+1] = true
		end
	end

	if false in initialize_check
		error("Failed to initialize the design matrix.")
	end

	return design_matrix
end

"""
	calc_X_element_energy(salc, spin_directions, symmetry) -> Float64

calculate an element of the design matrix X in the case of using the energy information.

# Arguments
- `salc::SALC`: Symmetry-Adapted Linear Combination object
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)
- `symmetry::Symmetry`: Symmetry information of the structure

# Returns
- `Float64`: Element of the design matrix X
"""
function calc_X_element_energy(
	salc::SALC,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
)::Float64

	result::Float64 = 0.0
	for itrans in symmetry.symnum_translation
		for (basis_idx, basis::IndicesUniqueList) in enumerate(salc.basisset)
			product_tmp::Float64 = 1.0
			for ibasis::Indices in basis
				atom::Int = symmetry.map_sym[ibasis.atom, itrans]
				l::Int = ibasis.l
				m::Int = ibasis.m
				product_tmp *= Sₗₘ(l, m, spin_directions[:, atom])
			end
			result += salc.coeffs[basis_idx] * salc.multiplicity[basis_idx] * product_tmp
		end
	end

	return result
end

function construct_design_matrix_magfield_vertical(
	salc_list,
	spinconfig_list,
	num_atoms,
	symmetry,
)::Matrix{Float64}
	# dimensions
	num_salcs = length(salc_list)
	num_spinconfigs = length(spinconfig_list)

	# construct design matrix A in Ax = b
	# [num_spindconif][3*num_atoms, num_salcs]
	design_matrix_list = Vector{Matrix{Float64}}(undef, num_spinconfigs)

	@threads for i in eachindex(spinconfig_list)
		design_matrix_list[i] = zeros(Float64, 3 * num_atoms, num_salcs)
		for row_idx in 1:(3*num_atoms)
			for isalc in eachindex(salc_list)
				design_matrix_list[i][row_idx, isalc] =
					calc_X_element_magfield_vertical(
						salc_list[isalc],
						spinconfig_list[i].spin_directions,
						symmetry,
						row_idx,
					)
			end
		end
	end

	design_matrix = vcat(design_matrix_list...)

	return design_matrix
end

"""
	calc_X_element_magfield_vertical(salc, spin_directions, symmetry, row_idx) -> Float64

Calculate an element of the design matrix X in the case of using the derivative of SALCs.

# Arguments
- `salc::SALC`: Symmetry-Adapted Linear Combination object
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)
- `symmetry::Symmetry`: Symmetry information of the structure
- `row_idx::Integer`: Row index corresponding to atom and direction (3*(atom-1) + dir)

"""
function calc_X_element_magfield_vertical(
	salc::SALC,
	spin_directions::AbstractMatrix{<:Real},
	symmetry::Symmetry,
	row_idx::Int,
)::Float64
	atom_idx = (row_idx - 1) ÷ 3 + 1
	direction = mod1(row_idx, 3)# 1, 2, 3 -> x, y, z

	result::Float64 = 0.0
	for itrans in symmetry.symnum_translation
		translated_salc = translate_atom_idx_of_salc(salc, symmetry.map_sym, itrans)
		result +=
			calc_derivative_of_salc(translated_salc, atom_idx, direction, spin_directions)
	end

	return result
end

"""
	ols_energy(design_matrix_energy, design_matrix_magfield_vertical, observed_energy_list, observed_magfield_vertical_list)
	-> Tuple{Vector{Float64}, Float64, Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

Optimize SCE coefficients by using ordinary least squares method.

# Arguments
- `design_matrix_energy::AbstractMatrix{<:Real}`: Design matrix for energy prediction
- `design_matrix_magfield_vertical::AbstractMatrix{<:Real}`: Design matrix for magfield_vertical prediction
- `observed_energy_list::Vector{Float64}`: Observed energy list
- `observed_magfield_vertical_list::Vector{Vector{Float64}}`: Observed magfield_vertical list

# Returns
- `Vector{Float64}`: Optimized SCE coefficients
- `Float64`: Bias term
- `Float64`: Relative error of magfield_vertical
- `Float64`: Relative error of energy
- `Vector{Float64}`: Predicted energy list
- `Vector{Float64}`: Observed energy list
- `Vector{Float64}`: Predicted magfield_vertical flattened list
- `Vector{Float64}`: Observed magfield_vertical flattened list
"""
function ols_energy(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_magfield_vertical::AbstractMatrix{<:Real},
	observed_energy_list::Vector{Float64},
	observed_magfield_vertical_list::Vector{Vector{Float64}},
)
	ols_coeffs = design_matrix_energy \ observed_energy_list

	# predict energy using SCE coefficients from energy information
	predicted_energy_list::Vector{Float64} = design_matrix_energy * ols_coeffs

	observed_magfield_vertical_flattened::Vector{Float64} =
		vcat(observed_magfield_vertical_list...)
	observed_magfield_vertical_flattened = -1 * observed_magfield_vertical_flattened

	reference_energy::Float64 = ols_coeffs[1]
	ols_coeffs_wo_ref_energy::Vector{Float64} = ols_coeffs[2:end]
	predicted_magfield_vertical_list::Vector{Float64} =
		design_matrix_magfield_vertical * ols_coeffs_wo_ref_energy
	relative_error_energy::Float64 =
		√(
		sum((observed_energy_list - predicted_energy_list) .^ 2) /
		sum(observed_energy_list .^ 2),
	)
	relative_error_magfield_vertical::Float64 =
		√(
		sum(
			(observed_magfield_vertical_flattened - predicted_magfield_vertical_list) .^ 2,
		) /
		sum(observed_magfield_vertical_flattened .^ 2),
	)

	return ols_coeffs_wo_ref_energy,
	reference_energy,
	relative_error_magfield_vertical,
	relative_error_energy,
	predicted_energy_list,
	observed_energy_list,
	predicted_magfield_vertical_list,
	observed_magfield_vertical_flattened
end

"""
	ols_magfield_vertical(design_matrix_energy, design_matrix_magfield_vertical, observed_energy_list, observed_magfield_vertical_list)
	-> Tuple{Vector{Float64}, Float64, Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

Optimize SCE coefficients by using ordinary least squares method.

# Arguments
- `design_matrix_energy::AbstractMatrix{<:Real}`: Design matrix for energy prediction
- `design_matrix_magfield_vertical::AbstractMatrix{<:Real}`: Design matrix for magfield_vertical prediction
- `observed_energy_list::Vector{Float64}`: Observed energy list
- `observed_magfield_vertical_list::Vector{Vector{Float64}}`: Observed magfield_vertical list

# Returns
- `Vector{Float64}`: Optimized SCE coefficients
- `Float64`: Bias term
- `Float64`: Relative error of magfield_vertical
- `Float64`: Relative error of energy
- `Vector{Float64}`: Predicted energy list
- `Vector{Float64}`: Observed energy list
- `Vector{Float64}`: Predicted magfield_vertical flattened list
- `Vector{Float64}`: Observed magfield_vertical flattened list
"""
function ols_magfield_vertical(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_magfield_vertical::AbstractMatrix{<:Real},
	observed_energy_list::Vector{Float64},
	observed_magfield_vertical_list::Vector{Vector{Float64}},
)
	observed_magfield_vertical_flattened = vcat(observed_magfield_vertical_list...)
	observed_magfield_vertical_flattened = -1 * observed_magfield_vertical_flattened

	# calculate SCE coefficients from magfield_vertical information
	ols_coeffs = design_matrix_magfield_vertical \ observed_magfield_vertical_flattened
	predicted_magfield_vertical_flattened = design_matrix_magfield_vertical * ols_coeffs
	relative_error_magfield_vertical =
		√(
		sum(
			(observed_magfield_vertical_flattened - predicted_magfield_vertical_flattened) .^
			2,
		) /
		sum(observed_magfield_vertical_flattened .^ 2),
	)

	# calculate reference energy term

	reference_energy =
		mean(observed_energy_list .- design_matrix_energy[:, 2:end] * ols_coeffs)
	relative_error_energy =
		√(
		sum(
			(
				observed_energy_list .-
				(design_matrix_energy[:, 2:end] * ols_coeffs .+ reference_energy)
			) .^
			2,
		) /
		sum(observed_energy_list .^ 2),
	)

	predicted_energy_list = design_matrix_energy[:, 2:end] * ols_coeffs .+ reference_energy

	return ols_coeffs,
	reference_energy,
	relative_error_magfield_vertical,
	relative_error_energy,
	predicted_energy_list,
	observed_energy_list,
	predicted_magfield_vertical_flattened,
	observed_magfield_vertical_flattened
end

"""
	calc_magfield_vertical_list_of_spinconfig(spinconfig, num_atoms) -> Vector{Float64}

Calculate the magfield_vertical vectors for each atom in the spin configuration.

# Arguments
- `spinconfig::SpinConfig`: Spin configuration containing magnetic moments and fields
- `num_atoms::Integer`: Number of atoms in the structure

# Returns
A flattened vector of length 3*num_atoms containing the magfield_vertical components:
[B₁ₓ, B₁ᵧ, B₁ᵣ, B₂ₓ, B₂ᵧ, B₂ᵣ, ...]

where Bᵢ = μᵢ × Hᵢ (magfield_vertical = magnetic moment × local magnetic field)
"""
function calc_magfield_vertical_list_of_spinconfig(
	spinconfig::SpinConfig,
	num_atoms::Int,
)::Vector{Float64}
	# Preallocate the result vector
	magfield_vertical_list = zeros(3 * num_atoms)

	for iatom in 1:num_atoms
		# Get local magnetic field using SVector for better performance
		magfield_vertical =
			SVector{3, Float64}(spinconfig.local_magfield_vertical[:, iatom])

		# Calculate magfield_vertical and store in preallocated vector
		idx = (iatom - 1) * 3 + 1
		# Use SVector for the multiplication
		result = spinconfig.magmom_size[iatom] * magfield_vertical
		magfield_vertical_list[idx:(idx+2)] = result
	end

	return magfield_vertical_list
end


"""
	calc_derivative_of_salc(basislist, coeffs, atom, direction, spin_directions)
	calc_derivative_of_salc(salc, atom, direction, spin_directions)

Calculate the derivative of a Symmetry-Adapted Linear Combination (SALC) with respect 
to the spin direction of a specific atom.

# Arguments
- `basislist::AbstractVector{IndicesUniqueList}`: List of basis functions
- `coeffs::AbstractVector{<:Real}`: Coefficients for each basis
- `atom::Integer`: Target atom index
- `direction::Integer`: Direction of derivative (1=x, 2=y, 3=z)
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)

"""
function calc_derivative_of_salc(
	basislist::AbstractVector{IndicesUniqueList},
	coeffs::AbstractVector{<:Real},
	multiplicity_list::AbstractVector{<:Real},
	atom::Integer,
	direction::Integer,
	spin_directions::AbstractMatrix{<:Real},
)::Float64
	# Input validation
	if !(1 ≤ direction ≤ 3)
		throw(ArgumentError("direction must be 1, 2, or 3"))
	end

	derivative_function = d_Slm[direction]

	result::Float64 = 0.0

	# Iterate through each basis and coefficient
	for (basis, coeff, multiplicity) in zip(basislist, coeffs, multiplicity_list)
		# Skip if atom is not in the basis
		atom ∉ [indices.atom for indices in basis] && continue

		# Calculate product of spherical harmonics and their derivatives
		product = 1.0
		for indices in basis
			# Use SVector for better performance with vector operations
			spin_dir = SVector{3, Float64}(spin_directions[:, indices.atom])
			if indices.atom == atom
				product *= derivative_function(indices.l, indices.m, spin_dir)
			else
				product *= Sₗₘ(indices.l, indices.m, spin_dir)
			end
		end

		result += coeff * product * multiplicity
	end

	return result
end

function calc_derivative_of_salc(
	salc::SALC,
	atom::Integer,
	direction::Integer,
	spin_directions::AbstractMatrix{<:Real},
)::Float64
	return calc_derivative_of_salc(
		salc.basisset,
		salc.coeffs,
		salc.multiplicity,
		atom,
		direction,
		spin_directions,
	)
end


function ridge_regression(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_magfield::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_magfield_list::AbstractVector{<:AbstractVector{<:Real}},
	alpha::Real,
	weight::Real,
)

	w1 = weight
	w2 = 1 - weight
	# Flatten observed magfield_vertical
	observed_magfield_flattened = vcat(observed_magfield_list...)
	observed_magfield_flattened = -1 * observed_magfield_flattened

	# Normalize the design matrices by using factor of 1/2N_data and √weight
	num_data = size(design_matrix_energy, 1)
	normalized_design_matrix_energy =
		design_matrix_energy ./ (2 * num_data) .* sqrt(w1)
	normalized_design_matrix_magfield =
		design_matrix_magfield ./ (2 * num_data) .* sqrt(w2)

	# Also normalise the observed vectors
	normalized_observed_energy_list =
		observed_energy_list ./ (2 * num_data) .* sqrt(w1)
	normalized_observed_magfield_flattened =
		observed_magfield_flattened ./ (2 * num_data) .* sqrt(w2)

	# Add 0 bias term to the design matrix for magfield
	# to align with the energy design matrix
	with_bias_design_matrix_magfield =
		hcat(zeros(size(normalized_design_matrix_magfield, 1)), normalized_design_matrix_magfield)

	# Construct the augmented design matrix
	X = vcat(
		normalized_design_matrix_energy,
		with_bias_design_matrix_magfield,
	)

	y = vcat(
		normalized_observed_energy_list,
		normalized_observed_magfield_flattened
	)

	# Ridge regression solution
	fit = glmnet(X, y; alpha=0, lambda=[alpha], standardize=false)
	# Extract coefficients
	j0 = fit.betas[1, 1]  # Extract coefficients for the first lambda
	jphi = fit.betas[2:end, 1]  # Extract coefficients for the second lambda
	# j = (X' * X + alpha * I) \ (X' * y)
	# jphi = j[2:end]  # SCE coefficients without reference energy
	# j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)
	j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)
	return j0, jphi
end

"""
	calc_relative_error(
		sce_coeffs_with_ref_energy::AbstractVector{<:Real},
		design_matrix_energy::AbstractMatrix{<:Real},
		design_matrix_magfield_vertical::AbstractMatrix{<:Real},
		observed_energy_list::AbstractVector{<:Real},
		observed_magfield_vertical_list::AbstractVector{<:Real},
	) -> Tuple{Float64, Float64}

Compute the relative errors of energy and magfield_vertical.

# Arguments
- `sce_coeffs_with_ref_energy::AbstractVector{<:Real}`: SCE coefficients with reference energy
- `design_matrix_energy::AbstractMatrix{<:Real}`: Design matrix for energy
- `design_matrix_magfield_vertical::AbstractMatrix{<:Real}`: Design matrix for magfield_vertical
- `observed_energy_list::AbstractVector{<:Real}`: Observed energy list
- `observed_magfield_vertical_list::AbstractVector{<:Real}`: Observed magfield_vertical list
"""
function calc_relative_error(
	sce_coeffs_with_ref_energy::AbstractVector{<:Real},
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_magfield_vertical::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_magfield_vertical_list::AbstractVector{<:Real},
)::Tuple{Float64, Float64}
	predicted_energy_list = design_matrix_energy * sce_coeffs_with_ref_energy
	predicted_magfield_vertical_list =
		design_matrix_magfield_vertical * sce_coeffs_with_ref_energy[2:end]
	relative_error_energy::Float64 =
		√(
		sum((observed_energy_list .- predicted_energy_list) .^ 2) /
		sum(observed_energy_list .^ 2),
	)
	relative_error_magfield_vertical::Float64 =
		√(
		sum((observed_magfield_vertical_list .- predicted_magfield_vertical_list) .^ 2) /
		sum(observed_magfield_vertical_list .^ 2),
	)
	return relative_error_energy, relative_error_magfield_vertical
end


function translate_atom_idx_of_salc(
	salc::SALC,
	map_sym::AbstractMatrix{<:Integer},
	itrans::Integer,
)::SALC
	translated_basisset = Vector{IndicesUniqueList}()
	for basis::IndicesUniqueList in salc.basisset
		translated_basis = IndicesUniqueList()
		for indices::Indices in basis
			translated_atom_idx = map_sym[indices.atom, itrans]
			push!(
				translated_basis,
				Indices(translated_atom_idx, indices.l, indices.m, indices.cell),
			)
		end
		push!(translated_basisset, translated_basis)
	end

	return SALC(translated_basisset, salc.coeffs, salc.multiplicity)
end

function print_info(optimizer::Optimizer)
	println(
		"""
		============
		OPTIMIZATION
		============
		""",
	)
	println(@sprintf("E_ref: %.10f", optimizer.reference_energy))
	for (i, sce) in enumerate(optimizer.SCE)
		println(@sprintf("%9d: %15.10f", i, sce))
	end
	# println(
	# 	@sprintf(
	# 		"fitting error of force: %.4f %%",
	# 		optimizer.relative_error_magfield_vertical * 100
	# 	)
	# )
	# println(@sprintf("fitting error of energy: %.4e %%", optimizer.relative_error_energy * 100))
	println("")
	println(@sprintf("Elapsed time: %.6f seconds", optimizer.elapsed_time))

	println("-------------------------------------------------------------------")
end

end

