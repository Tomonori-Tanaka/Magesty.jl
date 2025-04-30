"""
	Module SALCs

This module provides a data structure for symmetry-adapted linear combinations (SALCs).

# Types
- `SALC`: A structure representing a symmetry-adapted linear combination

# Examples
```julia
julia> basis = SortedCountingUniqueVector([IndicesUniqueList([1, 2]), IndicesUniqueList([1, 2])])
julia> coeffs = [1.0, -1.0]
julia> salc = SALC(basis, coeffs)
number of terms: 1
2  -1.4142135624  IndicesUniqueList([1, 2])
```
"""
module SALCs

using Printf
using LinearAlgebra
using ..SortedContainer
using ..AtomicIndices

import Base: show

export SALC

"""
	SALC

A structure representing a symmetry-adapted linear combination.

# Fields
- `basisset::Vector{IndicesUniqueList}`: The basis set of atomic indices
- `coeffs::Vector{Float64}`: The coefficients of the linear combination
- `multiplicity::Vector{Int}`: The multiplicity of each basis element

# Constructors
- `SALC(basisset::Vector{IndicesUniqueList}, coeffs::Vector{Float64}, multiplicity::Vector{Int})`:
  Create a SALC from a basis set, coefficients, and multiplicity
"""
struct SALC
	basisset::Vector{IndicesUniqueList}
	coeffs::Vector{Float64}
	multiplicity::Vector{Int}

	function SALC(
		basisset::Vector{IndicesUniqueList},
		coeffs::Vector{Float64},
		multiplicity::Vector{Int},
	)
		length(basisset) == length(coeffs) == length(multiplicity) ||
			throw(ArgumentError("All vectors must have the same length"))
		all(x -> x > 0, multiplicity) ||
			throw(ArgumentError("Multiplicity must be positive"))
		new(basisset, coeffs, multiplicity)
	end
end

"""
	SALC(
		basislist::SortedCountingUniqueVector{IndicesUniqueList},
		coeffs::Vector{<:Real},
	) -> SALC

Create a SALC from a basis list and coefficients.

# Arguments
- `basislist::SortedCountingUniqueVector{IndicesUniqueList}`: The basis list with counts
- `coeffs::Vector{<:Real}`: The coefficients for each basis element

# Returns
- `SALC`: A new symmetry-adapted linear combination

# Throws
- `ArgumentError` if the lengths of `basislist` and `coeffs` differ
- `ArgumentError` if the resulting coefficient vector has zero norm
"""
function SALC(
	basislist::SortedCountingUniqueVector{IndicesUniqueList},
	coeffs::Vector{<:Real},
)
	length(basislist) == length(coeffs) ||
		throw(ArgumentError("The length of basislist and coeffs must be the same"))

	result_basisset = Vector{IndicesUniqueList}()
	result_coeffs = Vector{Float64}()
	result_multiplicity = Vector{Int}()

	for (idx, basis) in enumerate(basislist)
		count = basislist.counts[basis]
		coeff = coeffs[idx] * count
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

"""
	show(io::IO, salc::SALC)

Display a SALC in a human-readable format.

# Arguments
- `io::IO`: The output stream
- `salc::SALC`: The SALC to display

# Output Format
```
number of terms: N
M  COEFFICIENT  BASIS
```
where:
- N is the number of terms
- M is the multiplicity
- COEFFICIENT is the normalized coefficient
- BASIS is the basis element
"""
function show(io::IO, salc::SALC)
	println(io, "number of terms: ", length(salc.basisset))
	for (basis, coeff, multiplicity) in zip(salc.basisset, salc.coeffs, salc.multiplicity)
		println(io, @sprintf("%2d  % 15.10f  %s", multiplicity, coeff, basis))
	end
end

end
