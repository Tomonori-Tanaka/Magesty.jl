"""
	Optimize.jl

This module contains functions for optimizing the SCE coefficients.
"""
module Optimize

using LinearAlgebra
using ..MySphericalHarmonics
using ..SALCs
using ..AtomicIndices
using ..Systems
using ..Symmetries
using ..Clusters
using ..BasisSets
using ..SpinConfigs

export SCEOptimizer

struct SCEOptimizer
	SCE::Vector{Float64}
	energy_list::Vector{Float64}
end

function SCEOptimizer(
	system::System,
	symmetry::Symmetry,
	basisset::BasisSet,
	j_zero_thr::Real,
	weight::Real,
	spinconfig_list::AbstractVector{SpinConfig},
)
	SCE, energy_list = ols_energy(basisset.salc_list, spinconfig_list, symmetry)
	SCE_torque = ols_torque(basisset.salc_list, spinconfig_list, system.supercell.num_atoms, symmetry)

	display(SCE)
	display(SCE_torque)

	return SCEOptimizer(SCE, energy_list)
end

function SCEOptimizer(
	system::System,
	symmetry::Symmetry,
	basisset::BasisSet,
	j_zero_thr::Real,
	weight::Real,
	datafile::AbstractString,
)
	# read datafile
	spinconfig_list::Vector{SpinConfig} = read_embset(datafile, system.supercell.num_atoms)

	return SCEOptimizer(system, symmetry, basisset, j_zero_thr, weight, spinconfig_list)
end


function ols_energy(
	salc_list::AbstractVector{SALC},
	spinconfig_list::AbstractVector{SpinConfig},
	symmetry::Symmetry,
)::Tuple{Vector{Float64}, Vector{Float64}}
	# dimensions
	num_salcs = length(salc_list)
	num_spinconfigs = length(spinconfig_list)

	# construct design matrix A in Ax = b
	design_matrix = zeros(Float64, num_spinconfigs, num_salcs + 1)
	initialize_check = falses(num_spinconfigs, num_salcs + 1)

	# set first column to 1 (bias term)
	design_matrix[:, 1] .= 1.0
	initialize_check[:, 1] .= true

	for i in 1:num_salcs
		for j in 1:num_spinconfigs
			design_matrix[j, i+1] =
				calc_design_matrix_element(
					salc_list[i],
					spinconfig_list[j].spin_directions,
					symmetry,
				)
			# if i == 1
			# 	design_matrix[j, i+1] = 0.0
			# end
			initialize_check[j, i+1] = true
		end
	end

	if false in initialize_check
		error("Failed to initialize the design matrix.")
	end


	energy_list::Vector{Float64} = [spinconfig.energy for spinconfig in spinconfig_list]
	# solve Ax = b
	ols_coeffs = design_matrix \ energy_list

	predicted_energy_list::Vector{Float64} = design_matrix * ols_coeffs

	# end
	open("energy_comparison.txt", "w") do file
		for i in 1:num_spinconfigs
			println(file, "$(energy_list[i])\t$(predicted_energy_list[i])")
		end
	end

	return ols_coeffs, energy_list
end

function calc_design_matrix_element(
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

function ols_torque(
	salc_list::AbstractVector{SALC},
	spinconfig_list::AbstractVector{SpinConfig},
	num_atoms::Integer,
	symmetry::Symmetry,
)
	# dimensions
	num_salcs = length(salc_list)
	num_spinconfigs = length(spinconfig_list)

	# observed torque
	observed_torque_list = Vector{Vector{Float64}}(undef, num_spinconfigs)
	for i in 1:num_spinconfigs
		observed_torque_list[i] =
			calc_torque_list_of_spinconfig(spinconfig_list[i], num_atoms)
	end
	observed_torque_flattened = vcat(observed_torque_list...)


	# construct design matrix A in Ax = b
	# [num_spindconif][3*num_atoms, num_salcs]
	design_matrix_list = Vector{Matrix{Float64}}(undef, num_spinconfigs)

	for i in 1:num_spinconfigs
		design_matrix_list[i] = zeros(Float64, 3 * num_atoms, num_salcs)
		for row_idx in 1:3*num_atoms
			for isalc in 1:num_salcs
				design_matrix_list[i][row_idx, isalc] =
					calc_X_matrix_element(
						salc_list[isalc],
						spinconfig_list[i].spin_directions,
						symmetry,
						row_idx,
					)
			end
		end
	end

	# display(design_matrix_list[end])
	design_matrix = vcat(design_matrix_list...)

	ols_coeffs = design_matrix \ observed_torque_flattened

	return ols_coeffs

end

"""
	calc_torque_list_of_spinconfig(spinconfig, num_atoms) -> Vector{Float64}

Calculate the torque vectors for each atom in the spin configuration.

# Arguments
- `spinconfig::SpinConfig`: Spin configuration containing magnetic moments and fields
- `num_atoms::Integer`: Number of atoms in the system

# Returns
A flattened vector of length 3*num_atoms containing the torque components:
[τ₁ₓ, τ₁ᵧ, τ₁ᵣ, τ₂ₓ, τ₂ᵧ, τ₂ᵣ, ...]

where τᵢ = μᵢ × Bᵢ (torque = magnetic moment × local magnetic field)
"""
function calc_torque_list_of_spinconfig(
	spinconfig::SpinConfig,
	num_atoms::Int,
)::Vector{Float64}
	# Preallocate the result vector
	torque_list = zeros(3 * num_atoms)

	for iatom in 1:num_atoms
		# Calculate magnetic moment vector
		magmom = @view(spinconfig.spin_directions[:, iatom]) * spinconfig.magmom_size[iatom]

		# Get local magnetic field
		magfield = @view spinconfig.local_magfield[:, iatom]

		# Calculate torque and store in preallocated vector
		idx = (iatom - 1) * 3 + 1
		torque_list[idx:idx+2] .= cross(magmom, magfield)
	end

	return torque_list
end



"""
	calc_X_matrix_element(salc, spin_directions, symmetry, row_idx) -> Float64

Calculate an element of the design matrix X in the case of using the derivative of SALCs.

# Arguments
- `salc::SALC`: Symmetry-Adapted Linear Combination object
- `spin_directions::AbstractMatrix{<:Real}`: Matrix of spin directions (3×N)
- `symmetry::Symmetry`: Symmetry information of the system
- `row_idx::Integer`: Row index corresponding to atom and direction (3*(atom-1) + dir)

"""
function calc_X_matrix_element(
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
		atom ∉ get_atomlist(basis) && continue

		# Calculate product of spherical harmonics and their derivatives
		product = mapreduce(*, basis) do indices
			spin_dir = @view spin_directions[:, indices.atom]

			if indices.atom == atom
				derivative_function(indices.l, indices.m, spin_dir)
			else
				Sₗₘ(indices.l, indices.m, spin_dir)
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

end

