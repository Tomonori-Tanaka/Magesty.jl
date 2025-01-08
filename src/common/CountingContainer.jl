module CountingContainers

import Base: append!, copy, getindex, in, isless, iterate, length, push!, show, size, ==

export CountingUniqueVector, getcounts
# ─────────────────────────────────────────────────────────────────────────────
# CountingUniqueVector
# ─────────────────────────────────────────────────────────────────────────────
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
"""
mutable struct CountingUniqueVector{T} <: AbstractVector{T}
	data::Vector{T}
	counts::Dict{T, Int}
end

function CountingUniqueVector{T}() where T
	CountingUniqueVector{T}(Vector{T}(), Dict{T, Int}())
end

function CountingUniqueVector(data::AbstractVector{T}) where T
	vec::Vector{T} = unique(data)
	cnts = Dict{T, Int}()
	for value in data
		if haskey(cnts, value)
			cnts[value] += 1
		else
			cnts[value] = 1
		end
	end
	CountingUniqueVector{T}(vec, cnts)
end

getindex(cuv::CountingUniqueVector, idx::Int) = cuv.data[idx]
length(cuv::CountingUniqueVector) = length(cuv.data)
size(cuv::CountingUniqueVector) = size(cuv.data)
iterate(cuv::CountingUniqueVector) = iterate(cuv.data)     # Support iteration
iterate(cuv::CountingUniqueVector, state) = iterate(cuv.data, state)
function ==(cuv1::CountingUniqueVector, cuv2::CountingUniqueVector)
	return cuv1.data == cuv2.data && cuv1.counts == cuv2.counts
end
isless(cuv1::CountingUniqueVector, cuv2::CountingUniqueVector) = cuv1.data < cuv2.data

function push!(cuv::CountingUniqueVector{T}, val::T) where T
	if haskey(cuv.counts, val)
		cuv.counts[val] += 1
	else
		push!(cuv.data, val)
		cuv.counts[val] = 1
	end
	return cuv
end

function append!(
	cuv::CountingUniqueVector{T},
	new_vals::AbstractVector{T},
) where T
	for val in new_vals
		push!(cuv, val)
	end
	return cuv
end

function in(val::T, cuv::CountingUniqueVector{T}) where T
	return haskey(cuv.counts, val)
end

function copy(cuv::CountingUniqueVector{T}) where T
	data_copy = copy(cuv.data)
	counts_copy = copy(cuv.counts)
	return CountingUniqueVector(data_copy, counts_copy)
end

function show(io::IO, cuv::CountingUniqueVector{T}) where T
	for val in cuv.data
		print(io, "counts: ", cuv.counts[val], " ")
		println(io, "data: ", val)
	end
end

function getcounts(cuv::CountingUniqueVector{T}, val::T)::Int where T
	return get(cuv.counts, val, 0)
end
end