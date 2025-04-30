module CountingContainers

import Base: append!, copy, getindex, in, isless, iterate, length, push!, show, size, ==

export CountingUniqueVector, getcounts

"""
	CountingUniqueVector{T} <: AbstractVector{T}

A mutable vector of unique elements with counts of each element's occurrences.

# Fields
- `data::Vector{T}`: The underlying array of unique elements.
- `counts::Dict{T, Int}`: A dictionary mapping each element to its insertion count.

# Constructors
- `CountingUniqueVector()`: Create an empty `CountingUniqueVector{T}`.
- `CountingUniqueVector(data::AbstractVector{T})`: Create a counting vector from `data`,
  counting each element's occurrences.

# Examples
```julia
julia> cuv = CountingUniqueVector([1, 2, 2, 3, 3, 3])
counts: 1 data: 1
counts: 2 data: 2
counts: 3 data: 3

julia> getcounts(cuv, 2)
2
```
"""
mutable struct CountingUniqueVector{T} <: AbstractVector{T}
	data::Vector{T}
	counts::Dict{T, Int}

	function CountingUniqueVector{T}() where T
		new{T}(Vector{T}(), Dict{T, Int}())
	end

	function CountingUniqueVector{T}(data::Vector{T}, counts::Dict{T, Int}) where T
		length(data) == length(counts) || throw(ArgumentError("Data and counts must have the same number of unique elements"))
		new{T}(data, counts)
	end
end

function CountingUniqueVector(data::AbstractVector{T}) where T
	vec = unique(data)
	cnts = Dict{T, Int}()
	for value in data
		cnts[value] = get(cnts, value, 0) + 1
	end
	CountingUniqueVector{T}(vec, cnts)
end

# Base interface implementations
Base.getindex(cuv::CountingUniqueVector, idx::Int) = cuv.data[idx]
Base.length(cuv::CountingUniqueVector) = length(cuv.data)
Base.size(cuv::CountingUniqueVector) = size(cuv.data)
Base.iterate(cuv::CountingUniqueVector) = iterate(cuv.data)     # Support iteration
Base.iterate(cuv::CountingUniqueVector, state) = iterate(cuv.data, state)
function Base.:(==)(cuv1::CountingUniqueVector, cuv2::CountingUniqueVector)
	return cuv1.data == cuv2.data && cuv1.counts == cuv2.counts
end
Base.isless(cuv1::CountingUniqueVector, cuv2::CountingUniqueVector) = cuv1.data < cuv2.data

function Base.push!(cuv::CountingUniqueVector{T}, val::T) where T
	if haskey(cuv.counts, val)
		cuv.counts[val] += 1
	else
		push!(cuv.data, val)
		cuv.counts[val] = 1
	end
	return cuv
end

function Base.append!(
	cuv::CountingUniqueVector{T},
	new_vals::AbstractVector{T},
) where T
	for val in new_vals
		push!(cuv, val)
	end
	return cuv
end

function Base.in(val::T, cuv::CountingUniqueVector{T}) where T
	return haskey(cuv.counts, val)
end

function Base.copy(cuv::CountingUniqueVector{T}) where T
	CountingUniqueVector{T}(copy(cuv.data), copy(cuv.counts))
end

function Base.show(io::IO, cuv::CountingUniqueVector{T}) where T
	for val in cuv.data
		print(io, "counts: ", cuv.counts[val], " ")
		println(io, "data: ", val)
	end
end

"""
	getcounts(cuv::CountingUniqueVector{T}, val::T) -> Int

Get the count of occurrences of `val` in the counting vector.

# Arguments
- `cuv::CountingUniqueVector{T}`: The counting vector to query
- `val::T`: The value to count

# Returns
- `Int`: The number of times `val` appears in the vector

# Examples
```julia
julia> cuv = CountingUniqueVector([1, 2, 2, 3, 3, 3])
julia> getcounts(cuv, 2)
2
```
"""
function getcounts(cuv::CountingUniqueVector{T}, val::T)::Int where T
	return get(cuv.counts, val, 0)
end
end