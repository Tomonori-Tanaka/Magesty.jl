"""
	Module SALCs

This module provides a data structure for symmetry-adapted linear combinations (SALCs).
"""
module SALCs

using LinearAlgebra
using ..SortedContainer
using ..AtomicIndices

import Base: show

export SALC

struct SALC
	basisset::Vector{IndicesUniqueList}
	coeffs::Vector{Float64}
	multiplicity::Vector{Int}
end

function SALC(
	basislist::SortedCountingUniqueVector{IndicesUniqueList},
	coeffs::Vector{<:Real},
)
	if length(basislist) != length(coeffs)
		throw(ArgumentError("The length of basislist and coeffs must be the same."))
	end

	result_basisset = Vector{IndicesUniqueList}()
	result_coeffs = Vector{Float64}()
	result_multiplicity = Vector{Int}()

	for (idx, basis) in enumerate(basislist)
		count::Int = basislist.counts[basis]
		coeff::Float64 = coeffs[idx] * count
		if !isapprox(coeff, 0.0, atol = 1e-8)
			push!(result_basisset, basis)
			push!(result_coeffs, coeff)
			push!(result_multiplicity, count)
		end
	end

	# normalize coefficient vector
	norm_coeffs = norm(result_coeffs)
	if isapprox(norm_coeffs, 0.0, atol = 1e-8)
		throw(ArgumentError("The norm of the coefficient vector is zero."))
	end
	result_coeffs ./= norm_coeffs

	return SALC(result_basisset, result_coeffs, result_multiplicity)
end

function show(io::IO, salc::SALC)
	println(io, "number of SALCs: ", length(salc.basisset))
	for (basis, coeff, multiplicity) in zip(salc.basisset, salc.coeffs, salc.multiplicity)
		println(io, multiplicity, "\t", coeff, "\t", basis)
	end
end

end
