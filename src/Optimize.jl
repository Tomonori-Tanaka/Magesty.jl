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
end

function SCEOptimizer(
	system::System,
	symmetry::Symmetry,
	cluster::Cluster,
	basisset::BasisSet,
	j_zero_thr::Real,
	weight::Real,
	datafile::AbstractString,
)
	# read datafile
	spinconfig_list::Vector{SpinConfig} = read_embset(datafile, system.supercell.num_atoms)

	SCE = ols_energy(basisset.salc_list, spinconfig_list, symmetry)
	@show SCE
	return SCEOptimizer(SCE)
end

function ols_energy(
	salc_list::AbstractVector{SALC},
	spinconfig_list::AbstractVector{SpinConfig},
	symmetry::Symmetry,
)
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
				calc_design_matrix_element(salc_list[i], spinconfig_list[j], symmetry)
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

	# for i in 1:num_spinconfigs
	# 	println(
	# 		"predicted: $(predicted_energy_list[i]) actual: $(energy_list[i])",
	# 		" diff: $(predicted_energy_list[i] - energy_list[i])",
	# 	)
	# end

	return ols_coeffs * 1000 # convert to meV
end

function calc_design_matrix_element(
	salc::SALC,
	spinconfig::SpinConfig,
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
				product_tmp *= Sₗₘ(l, m, spinconfig.spin_directions[:, atom])
			end
			result += salc.coeffs[basis_idx] * salc.multiplicity[basis_idx] * product_tmp
		end
	end

	result = num_salc_basis * result
	return result
end

end
