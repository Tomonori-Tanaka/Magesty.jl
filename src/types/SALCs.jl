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

using EzXML
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

function SALC(
	basislist::SortedCountingUniqueVector{SHProduct},
	coeffs::Vector{<:Real},
)
	length(basislist) == length(coeffs) ||
		throw(ArgumentError("The length of basislist and coeffs must be the same"))

	result_basisset = Vector{IndicesUniqueList}()
	result_coeffs = Vector{Float64}()
	result_multiplicity = Vector{Int}()
	basis_index_map = Dict{IndicesUniqueList, Int}()

	for (idx, basis::SHProduct) in enumerate(basislist)
		iul::IndicesUniqueList = convert2indices(basis)
		count = basislist.counts[basis]
		coeff = coeffs[idx] * count
		if haskey(basis_index_map, iul)
			acc_idx = basis_index_map[iul]
			result_coeffs[acc_idx] += coeff
			result_multiplicity[acc_idx] += count
		else
			push!(result_basisset, iul)
			push!(result_coeffs, coeff)
			push!(result_multiplicity, count)
			basis_index_map[iul] = length(result_basisset)
		end
	end

	# remove near-zero coefficients
	keep_mask = .!isapprox.(result_coeffs, 0.0, atol = 1e-8)
	result_basisset = result_basisset[keep_mask]
	result_coeffs = result_coeffs[keep_mask]
	result_multiplicity = result_multiplicity[keep_mask]

	# normalize coefficient vector
	norm_coeffs = norm(result_coeffs)
	if isapprox(norm_coeffs, 0.0, atol = 1e-8)
		throw(ArgumentError("The norm of the coefficient vector is zero."))
	end
	result_coeffs ./= norm_coeffs

	return SALC(result_basisset, result_coeffs, result_multiplicity)
end


"""
	SALC(xml_node::EzXML.Node)

Create a SALC from an XML (SALC) node.

# Arguments
- `xml_node::EzXML.Node`: An XML node representing a SALC

# Returns
- `SALC`: A new symmetry-adapted linear combination

# Throws
- `ArgumentError` if the XML node is invalid
- `ArgumentError` if multiplicity attribute is missing or invalid
- `ArgumentError` if coefficient is missing or invalid
- `ArgumentError` if indices are not continuous or missing
- `ArgumentError` if the resulting coefficient vector has zero norm
"""
function SALC(xml_salc_node::EzXML.Node)
	multiplicity_list = Int[]
	basis_list = IndicesUniqueList[]
	coeff_list = Float64[]

	for basis_node in eachelement(xml_salc_node)
		# Check and parse multiplicity
		if !haskey(basis_node, "multiplicity")
			throw(ArgumentError("Multiplicity attribute is missing"))
		end
		multiplicity = try
			parse(Int, basis_node["multiplicity"])
		catch e
			throw(ArgumentError("Invalid multiplicity value: $(basis_node["multiplicity"])"))
		end
		if multiplicity <= 0
			throw(ArgumentError("Multiplicity must be positive"))
		end
		push!(multiplicity_list, multiplicity)

		# Parse indices
		basis_indices = Indices[]

		# Count number of index attributes
		num_indices = count(attr -> startswith(attr.name, "index-"), attributes(basis_node))
		if num_indices == 0
			throw(ArgumentError("No indices found in basis node"))
		end

		# Check index continuity and parse values
		for i in 1:num_indices
			name = "index-$i"
			# "1 1 -1 1" → ["1","1","-1","1"] → [1,1,-1,1]
			idx = parse.(Int, split(basis_node[name]))
			push!(basis_indices, Indices(idx...))
		end

		basis = IndicesUniqueList(basis_indices)
		push!(basis_list, basis)

		# Parse coefficient
		try
			coeff = parse(Float64, nodecontent(basis_node))
			push!(coeff_list, coeff)
		catch e
			throw(ArgumentError("Invalid coefficient value: $(nodecontent(basis_node))"))
		end
	end

	# Normalize coefficients
	norm_coeffs = norm(coeff_list)
	if isapprox(norm_coeffs, 0.0, atol = 1e-8)
		throw(ArgumentError("The norm of the coefficient vector is zero"))
	end
	coeff_list ./= norm_coeffs

	return SALC(basis_list, coeff_list, multiplicity_list)
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
