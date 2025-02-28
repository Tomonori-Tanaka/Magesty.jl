"""
	Optimize.jl

This module contains functions for optimizing the SCE coefficients.
"""
module Optimize

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

	display(design_matrix)

	energy_list::Vector{Float64} = [spinconfig.energy for spinconfig in spinconfig_list]
	# solve Ax = b
	ols_coeffs = design_matrix \ energy_list

	predicted_energy_list::Vector{Float64} = design_matrix * ols_coeffs

	# for i in 1:num_spinconfigs
	# 	println(
	# 		"predicted: $(predicted_energy_list[i]) actual: $(energy_list[i])",
	# 	)
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
	num_salc_basis::Int = length(salc.basisset)

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

	# result = num_salc_basis * result
	return result
end

end
