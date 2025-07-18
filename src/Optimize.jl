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
	predicted_magfield_list::Vector{Float64}
	observed_magfield_list::Vector{Float64}
	rmse_energy::Float64  # Root Mean Square Error for energy

	function Optimizer(
		structure::Structure,
		symmetry::Symmetry,
		basisset::BasisSet,
		alpha::Real,
		lambda::Real,
		weight::Real,
		spinconfig_list::AbstractVector{SpinConfig},
		;
		verbosity::Bool = true,
	)
		# Start timing
		start_time = time_ns()

		if verbosity
			println("""

			OPTIMIZATION
			============

			""")
		end

		# construct design matrix for energy and magfield
		if verbosity
			println("Constructing design matrix for energy...")
		end
		design_matrix_energy = construct_design_matrix_energy(
			basisset.salc_list,
			spinconfig_list,
			symmetry,
		)
		if verbosity
			println("Constructing design matrix for magfield...")
		end
		design_matrix_magfield = construct_design_matrix_magfield(
			basisset.salc_list,
			spinconfig_list,
			structure.supercell.num_atoms,
			symmetry,
		)

		# construct observed_energy_list and observed_magfield_list
		if verbosity
			println("Constructing observed data for energy and magfield...")
		end
		observed_energy_list = [spinconfig.energy for spinconfig in spinconfig_list]
		observed_magfield_list =
			Vector{Vector{Float64}}(undef, length(spinconfig_list))
		@threads for i in eachindex(spinconfig_list)
			observed_magfield_list[i] =
				calc_magfield_list_of_spinconfig(
					spinconfig_list[i],
					structure.supercell.num_atoms,
				)
		end

		if verbosity
			println("Fitting SCE coefficients...")
		end
		j0, jphi = elastic_net_regression(
			design_matrix_energy,
			design_matrix_magfield,
			observed_energy_list,
			observed_magfield_list,
			alpha,
			lambda,
			weight,
		)

		predicted_energy_list = design_matrix_energy[:, 2:end] * jphi .+ j0
		predicted_magfield_list = design_matrix_magfield * jphi
		observed_magfield_list =
			-1 * vcat(observed_magfield_list...)

		# Calculater RMSE for energy
		rmse_energy = sqrt(
			mean((observed_energy_list .- predicted_energy_list) .^ 2),
		)

		if verbosity
			print_optimize_stdout(
				j0,
				jphi,
				observed_energy_list,
				predicted_energy_list,
				observed_magfield_list,
				predicted_magfield_list,
			)
			elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds
			println(@sprintf(
				"""

				Time Elapsed: %.6f sec.
				""",
				elapsed_time,
			))
		end


		return new(
			spinconfig_list,
			j0,
			jphi,
			predicted_energy_list,
			observed_energy_list,
			predicted_magfield_list,
			observed_magfield_list,
			rmse_energy,
		)
	end
end

function Optimizer(
	structure::Structure,
	symmetry::Symmetry,
	basisset::BasisSet,
	alpha::Real,
	lambda::Real,
	weight::Real,
	datafile::AbstractString,
	;
	verbosity::Bool = true,
)
	# read datafile
	spinconfig_list = read_embset(datafile, structure.supercell.num_atoms)

	return Optimizer(
		structure,
		symmetry,
		basisset,
		alpha,
		lambda,
		weight,
		spinconfig_list,
		verbosity = verbosity,
	)
end

function Optimizer(
	structure::Structure,
	symmetry::Symmetry,
	basisset::BasisSet,
	config::Config4Optimize,
	;
	verbosity::Bool = true,
)
	return Optimizer(
		structure,
		symmetry,
		basisset,
		config.alpha,
		config.lambda,
		config.weight,
		config.datafile,
		verbosity = verbosity,
	)
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

	for i in 1:num_salcs
		for j in 1:num_spinconfigs
			design_matrix[j, i+1] = calc_X_element_energy(
				salc_list[i],
				spinconfig_list[j].spin_directions,
				symmetry,
			)
			initialize_check[j, i+1] = true
		end
	end

	if false in initialize_check
		false_indices = findall(x -> x == false, initialize_check)
		error("""
			Failed to initialize the design matrix.
			False values found at indices: $false_indices
			Full initialize_check array: $initialize_check
			""")
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

function construct_design_matrix_magfield(
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
	init_flags = fill(false, num_spinconfigs)

	@threads for i in eachindex(spinconfig_list)
		M = zeros(Float64, 3 * num_atoms, num_salcs)
		for row_idx in 1:(3*num_atoms)
			for isalc in eachindex(salc_list)
				M[row_idx, isalc] =
					calc_X_element_magfield(
						salc_list[isalc],
						spinconfig_list[i].spin_directions,
						symmetry,
						row_idx,
					)
			end
		end
		design_matrix_list[i] = M
		init_flags[i] = true
	end

	if false in init_flags
		false_indices = findall(x -> x == false, init_flags)
		error("""
			Failed to initialize the design matrix for magfield.
			False values found at indices: $false_indices
			Full init_flags array: $init_flags
			""")
	end

	design_matrix = vcat(design_matrix_list...)

	return design_matrix
end

"""
	calc_X_element_magfield(salc, spin_directions, symmetry, row_idx) -> Float64

Calculate an element of the design matrix X in the case of using the derivative of SALCs.

# Arguments
- `salc::SALC`: Symmetry-Adapted Linear Combination object
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)
- `symmetry::Symmetry`: Symmetry information of the structure
- `row_idx::Integer`: Row index corresponding to atom and direction (3*(atom-1) + dir)

"""
function calc_X_element_magfield(
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

# """
# 	ols_energy(design_matrix_energy, design_matrix_magfield,
# 	-> Tuple{Vector{Float64}, Float64, Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

# Optimize SCE coefficients by using ordinary least squares method.

# # Arguments
# - `design_matrix_energy::AbstractMatrix{<:Real}`: Design matrix for energy prediction
# - `design_matrix_magfield::AbstractMatrix{<:Real}`: Design matrix for magfield prediction
# - `observed_energy_list::Vector{Float64}`: Observed energy list
# - `observed_magfield_list::Vector{Vector{Float64}}`: Observed magfield list

# # Returns
# - `Vector{Float64}`: Optimized SCE coefficients
# - `Float64`: Bias term
# - `Float64`: Relative error of magfield
# - `Float64`: Relative error of energy
# - `Vector{Float64}`: Predicted energy list
# - `Vector{Float64}`: Observed energy list
# - `Vector{Float64}`: Predicted magfield flattened list
# - `Vector{Float64}`: Observed magfield flattened list
# """
# function ols_energy(
# 	design_matrix_energy::AbstractMatrix{<:Real},
# 	observed_energy_list::Vector{Float64},
# )
# 	# fit = glmnet(
# 	# 	design_matrix_energy,
# 	# 	observed_energy_list;
# 	# 	alpha = 0.0,
# 	# 	lambda = [0.0],
# 	# 	standardize = true,
# 	# )
# 	# Extract coefficients
# 	# jphi = fit.betas[2:end, 1]
# 	# j0 = fit.a0[1]

# 	# predict energy using SCE coefficients from energy information
# 	ols_coeffs = design_matrix_energy \ observed_energy_list
# 	j0::Float64 = ols_coeffs[1]
# 	jphi::Vector{Float64} = ols_coeffs[2:end]

# 	return j0, jphi
# end

# """
# 	ols_magfield(design_matrix_energy, design_matrix_magfield, observed_energy_list, observed_magfield_list)
# 	-> Tuple{Vector{Float64}, Float64, Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

# Optimize SCE coefficients by using ordinary least squares method.

# # Arguments
# - `design_matrix_energy::AbstractMatrix{<:Real}`: Design matrix for energy prediction
# - `design_matrix_magfield::AbstractMatrix{<:Real}`: Design matrix for magfield prediction
# - `observed_energy_list::Vector{Float64}`: Observed energy list
# - `observed_magfield_list::Vector{Vector{Float64}}`: Observed magfield list

# # Returns
# - `Vector{Float64}`: Optimized SCE coefficients
# - `Float64`: Bias term
# - `Float64`: Relative error of magfield
# - `Float64`: Relative error of energy
# - `Vector{Float64}`: Predicted energy list
# - `Vector{Float64}`: Observed energy list
# - `Vector{Float64}`: Predicted magfield flattened list
# - `Vector{Float64}`: Observed magfield flattened list
# """
# function ols_magfield(
# 	design_matrix_energy::AbstractMatrix{<:Real},
# 	design_matrix_magfield::AbstractMatrix{<:Real},
# 	observed_energy_list::Vector{Float64},
# 	observed_magfield_list::Vector{Vector{Float64}},
# )
# 	observed_magfield_flattened = vcat(observed_magfield_list...)
# 	observed_magfield_flattened = -1 * observed_magfield_flattened

# 	# calculate SCE coefficients from magfield information
# 	ols_coeffs = design_matrix_magfield \ observed_magfield_flattened

# 	# calculate reference energy term
# 	reference_energy =
# 		mean(observed_energy_list .- design_matrix_energy[:, 2:end] * ols_coeffs)

# 	return reference_energy, ols_coeffs
# end

"""
	calc_magfield_vertical_list_of_spinconfig(spinconfig, num_atoms) -> Vector{Float64}

Calculate the vertical component of magfield vectors for each atom in the spin configuration.


# Arguments
- `spinconfig::SpinConfig`: Spin configuration containing magnetic moments and fields
- `num_atoms::Integer`: Number of atoms in the structure

# Returns
A flattened vector of length 3*num_atoms containing the vertical components:
[B₁ₓ, B₁ᵧ, B₁ᵣ, B₂ₓ, B₂ᵧ, B₂ᵣ, ...]

where Bᵢ = μᵢ × Hᵢ (magfield_vertical = magnetic moment × local magnetic field)
"""
function calc_magfield_list_of_spinconfig(
	spinconfig::SpinConfig,
	num_atoms::Int,
)::Vector{Float64}
	# Preallocate the result vector
	magfield_list = zeros(3 * num_atoms)

	for iatom in 1:num_atoms

		# Calculate magfield_vertical and store in preallocated vector
		idx = (iatom - 1) * 3 + 1
		magfield_list[idx:(idx+2)] =
			spinconfig.magmom_size[iatom] * spinconfig.local_magfield_vertical[:, iatom]
	end

	return magfield_list
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

	spin_dirs = [SVector{3, Float64}(spin_directions[:, i]) for i in axes(spin_directions, 2)]

	derivative_function = d_Slm[direction]

	result::Float64 = 0.0

	# Iterate through each basis and coefficient
	for (basis, coeff, multiplicity) in zip(basislist, coeffs, multiplicity_list)
		# Skip if atom is not in the basis
		atom ∉ [indices.atom for indices in basis] && continue

		# Calculate product of spherical harmonics and their derivatives
		@fastmath begin
			product = 1.0
			for indices in basis
				if indices.atom == atom
					product *=
						derivative_function(indices.l, indices.m, @inbounds spin_dirs[indices.atom])
				else
					product *= Sₗₘ(indices.l, indices.m, @inbounds spin_dirs[indices.atom])
				end
			end
			result += coeff * product * multiplicity
		end
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


function elastic_net_regression(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_magfield::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_magfield_list::AbstractVector{<:AbstractVector{<:Real}},
	alpha::Real,
	lambda::Real,
	weight::Real,
)
	# weight parameters
	w_e = 1 - weight
	w_m = weight

	# Flatten observed magfield
	observed_magfield_flattened = vcat(observed_magfield_list...)
	observed_magfield_flattened = -1 * observed_magfield_flattened

	# Normalize the design matrices by using factor of 1/2N_data and √weight
	num_data = size(design_matrix_energy, 1)
	normalized_design_matrix_energy =
		design_matrix_energy .* sqrt(w_e)
	normalized_design_matrix_energy[:, 1] .= 1.0
	normalized_design_matrix_magfield =
		design_matrix_magfield .* sqrt(w_m)

	# Also normalise the observed vectors
	normalized_observed_energy_list =
		observed_energy_list .* sqrt(w_e)
	normalized_observed_magfield_flattened =
		observed_magfield_flattened .* sqrt(w_m)

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
		normalized_observed_magfield_flattened,
	)

	# betas = X \ y
	# jphi = betas[2:end]
	# j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)
	# println("j0: $j0")
	# println("jphi: $jphi")

	# Elastic net regression solution
	lambda_vec = fill(lambda, size(X, 2))
	lambda_vec[1] = 0.0  # exclude bias term from regularization
	j_values = begin
		if lambda ≈ 0.0
			X \ y
		else
			ridge(X, y, lambda_vec; bias = false)
		end
	end
	jphi = j_values[2:end]
	j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)


	# fit = glmnet(X, y; alpha = alpha, lambda = [lambda], standardize = true)
	# fit = glmnet(normalized_design_matrix_energy, normalized_observed_energy_list; alpha = alpha, lambda = [lambda], standardize = true)
	# Extract coefficients
	# jphi = fit.betas[2:end, 1]
	# j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)
	# j0 = fit.a0[1] * sqrt(2 * num_data) * num_atoms / sqrt(w_e)


	return j0, jphi
end

"""
	calc_relative_error(
		sce_coeffs_with_ref_energy::AbstractVector{<:Real},
		design_matrix_energy::AbstractMatrix{<:Real},
		design_matrix_magfield::AbstractMatrix{<:Real},
		observed_energy_list::AbstractVector{<:Real},
		observed_magfield_list::AbstractVector{<:Real},
	) -> Tuple{Float64, Float64}

Compute the relative errors of energy and magfield

# Arguments
- `sce_coeffs_with_ref_energy::AbstractVector{<:Real}`: SCE coefficients with reference energy
- `design_matrix_energy::AbstractMatrix{<:Real}`: Design matrix for energy
- `design_matrix_magfield::AbstractMatrix{<:Real}`: Design matrix for magfield
- `observed_energy_list::AbstractVector{<:Real}`: Observed energy list
- `observed_magfield_list::AbstractVector{<:Real}`: Observed magfield list
"""
function calc_relative_error(
	sce_coeffs_with_ref_energy::AbstractVector{<:Real},
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_magfield::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_magfield_list::AbstractVector{<:Real},
)::Tuple{Float64, Float64}
	predicted_energy_list = design_matrix_energy * sce_coeffs_with_ref_energy
	predicted_magfield_list =
		design_matrix_magfield * sce_coeffs_with_ref_energy[2:end]
	relative_error_energy::Float64 =
		√(
		sum((observed_energy_list .- predicted_energy_list) .^ 2) /
		sum(observed_energy_list .^ 2),
	)
	relative_error_magfield::Float64 =
		√(
		sum((observed_magfield_list .- predicted_magfield_list) .^ 2) /
		sum(observed_magfield_list .^ 2),
	)
	return relative_error_energy, relative_error_magfield
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

function print_optimize_stdout(
	reference_energy::Float64,
	sce_list::AbstractVector{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	predicted_energy_list::AbstractVector{<:Real},
	observed_magfield_list::AbstractVector{<:Real},
	predicted_magfield_list::AbstractVector{<:Real},
)

	rmse_energy = calc_rmse(observed_energy_list, predicted_energy_list)
	rmse_magfield = calc_rmse(observed_magfield_list, predicted_magfield_list)
	relative_error_energy = calc_relative_error(observed_energy_list, predicted_energy_list)
	relative_error_magfield = calc_relative_error(observed_magfield_list, predicted_magfield_list)

	delta_energy = (1/2)*mean((observed_energy_list - predicted_energy_list) .^ 2)
	delta_magfield = (1/2)*mean((observed_magfield_list - predicted_magfield_list) .^ 2)

	println(@sprintf("   E_ref: %.10f", reference_energy))
	for (i, sce) in enumerate(sce_list)
		println(@sprintf("%8d: %15.10f", i, sce))
	end

	println(@sprintf(
		"""

		Loss function part
		------------------
		L = (1-w)*Delta_E + w*Delta_B
		Delta_E: %.10f (meV)^2
		Delta_B: %.10f (meV)^2
		""",
		delta_energy * 1e6,
		delta_magfield * 1e6,
	))

	println(@sprintf(
		"""
		Root Mean Square Error (RMSE)
		-----------------------------
		RMSE for energy: %.4f meV
		RMSE for magnetic field: %.4f meV
		""",
		rmse_energy * 1000,
		rmse_magfield * 1000,
	))
	println(
		@sprintf(
			"""
			Relative Errors
			---------------
			Relative error for energy: %.4f %%
			Relative error for magnetic field: %.4f %%
			""",
			relative_error_energy * 100,
			relative_error_magfield * 100,
		)
	)
end

function calc_rmse(list1::AbstractVector{<:Real}, list2::AbstractVector{<:Real})::Float64
	# Calculate the Root Mean Square Error (RMSE) between two lists
	if length(list1) != length(list2)
		throw(ArgumentError("The lengths of the two lists must be equal."))
	end
	return sqrt(mean((list1 .- list2) .^ 2))
end

function calc_relative_error(
	observed_list::AbstractVector{<:Real},
	predicted_list::AbstractVector{<:Real},
)::Float64
	# Calculate the relative error between observed and predicted lists
	if length(observed_list) != length(predicted_list)
		throw(ArgumentError("The lengths of the two lists must be equal."))
	end
	return sqrt(sum((observed_list .- predicted_list) .^ 2) / sum(observed_list .^ 2))
end

end

