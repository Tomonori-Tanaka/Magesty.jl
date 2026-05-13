"""
	Module SortedContainer

This module provides various sorted container data structures, including:
- `SortedVector`: A mutable sorted vector allowing duplicate elements.
- `SortedUniqueVector`: A mutable sorted vector that avoids duplicates.
- `SortedCountingUniqueVector`: A mutable sorted vector that also keeps track of element counts.

All data structures implement `AbstractVector` and support operations such as indexing,
iteration, insertion, deletion, and clearing in a sorted manner.
"""
module SortedContainer

import Base: append!, copy, findall, findfirst, delete!, deleteat!, getindex, in, isless,
	isempty, iterate, length, push!, show, size, ==

export SortedVector
export SortedUniqueVector
export SortedCountingUniqueVector, getcount, addcount!

"""
	AbstractSortedVector{T} <: AbstractVector{T}

An abstract type representing a sorted vector of element type `T`.
Concrete subtypes maintain a sorted order of their elements and must expose
a `data` field that itself supports indexing, iteration, length, and sorted
search via `searchsortedfirst` (i.e. either a `Vector{T}` or another
`AbstractSortedVector{T}` for nested containers).
"""
abstract type AbstractSortedVector{T} <: AbstractVector{T} end

# ─────────────────────────────────────────────────────────────────────────────
# Generic fallbacks (apply to all AbstractSortedVector subtypes that expose
# a sorted `data` field). Type-specific behavior (push! semantics, == across
# different concrete types, counted-container overrides) is defined per type.
# ─────────────────────────────────────────────────────────────────────────────

getindex(asv::AbstractSortedVector, i::Int) = asv.data[i]
length(asv::AbstractSortedVector) = length(asv.data)
size(asv::AbstractSortedVector) = size(asv.data)
iterate(asv::AbstractSortedVector) = iterate(asv.data)
iterate(asv::AbstractSortedVector, state) = iterate(asv.data, state)
isempty(asv::AbstractSortedVector) = isempty(asv.data)

function append!(asv::AbstractSortedVector{T}, new_vec::AbstractVector{T}) where T
	for value in new_vec
		push!(asv, value)
	end
	return asv
end

function findfirst(asv::AbstractSortedVector{T}, value::T) where T
	idx = searchsortedfirst(asv.data, value)
	return idx <= length(asv.data) && asv.data[idx] == value ? idx : nothing
end

function in(value::T, asv::AbstractSortedVector{T}) where T
	idx = searchsortedfirst(asv.data, value)
	return idx <= length(asv.data) && asv.data[idx] == value
end



# ─────────────────────────────────────────────────────────────────────────────
# SortedVector
# ─────────────────────────────────────────────────────────────────────────────

"""
	SortedVector{T} <: AbstractSortedVector{T}

A mutable sorted vector of elements of type `T`.  
This structure maintains its internal `data::Vector{T}` in sorted order.

# Fields
- `data::Vector{T}`: The underlying sorted array storing all elements.

# Constructors
- `SortedVector(data)`: Create a `SortedVector` by sorting and storing a copy of `data`.
- `SortedVector{T}()`: Create an empty `SortedVector{T}`.
- `SortedVector()`: Create an empty `SortedVector{Any}`.
"""
mutable struct SortedVector{T} <: AbstractSortedVector{T}
	data::Vector{T}  # Underlying sorted vector

	# Constructor from AbstractVector
	function SortedVector(data)
		sorted_data = sort(collect(data))
		new{eltype(data)}(sorted_data)
	end

	function SortedVector(data::AbstractVector{T}) where T
		sorted_data = sort(collect(data))
		new{T}(sorted_data)
	end

	# Empty constructors
	SortedVector{T}() where T = new{T}(Vector{T}())
	SortedVector{T}(data::AbstractVector{T}) where T = new{T}(sort(collect(data)))
	SortedVector() = new{Any}(Vector())
end

# Type-specific operations (semantics differ from generic fallbacks)
==(sv1::SortedVector, sv2::SortedVector) = sv1.data == sv2.data
isless(sv1::SortedVector, sv2::SortedVector) = sv1.data < sv2.data

"""
	push!(sv::SortedVector{T}, value::T) -> SortedVector{T}

Insert `value` into the `SortedVector` while maintaining sorted order.
Returns the updated `SortedVector`.
"""
function push!(sv::SortedVector{T}, value::T) where T
	idx = searchsortedfirst(sv.data, value)  # Find the insertion position
	insert!(sv.data, idx, value)            # Insert the value at the calculated position
	return sv
end

"""
	findall(sv::SortedVector{T}, value::T) -> Vector{Int}

Return all indices where `value` occurs in `sv`. If `value` is not present,
an empty vector is returned.
"""
function findall(sv::SortedVector{T}, value::T) where T
	indices = []
	start_idx = searchsortedfirst(sv.data, value)
	while start_idx <= length(sv.data) && sv.data[start_idx] == value
		push!(indices, start_idx)
		start_idx += 1
	end
	return indices
end

function delete!(sv::SortedVector{T}, value::T) where T
	idx = Base.findfirst(sv.data .== value)
	if !isnothing(idx)
		deleteat!(sv.data, idx)
	end
	return sv
end

function deleteat!(sv::SortedVector{T}, idx::Int) where T
	if idx < 1 || idx > length(sv.data)
		throw(BoundsError(sv, idx))  # Throw an error if index is out of bounds
	end
	deleteat!(sv.data, idx)      # Remove the element at the specified index
	return sv
end

function deleteall!(sv::SortedVector{T}, value::T) where T
	indices = Base.findall(sv.data .== value)
	for idx in reverse(indices)
		deleteat!(sv.data, idx)
	end
	return sv
end

function clear!(sv::SortedVector)
	empty!(sv.data)
	return sv
end

function copy(sv::SortedVector{T}) where T
	return SortedVector{T}(copy(sv.data))
end


# ─────────────────────────────────────────────────────────────────────────────
# SortedUniqueVector
# ─────────────────────────────────────────────────────────────────────────────

"""
	SortedUniqueVector{T} <: AbstractSortedVector{T}

A mutable sorted vector that maintains unique elements of type `T`.

# Fields
- `data::Vector{T}`: The underlying array, always sorted and without duplicates.

# Constructors
- `SortedUniqueVector(data)`: Create a `SortedUniqueVector` by sorting, removing duplicates, and storing a copy of `data`.
- `SortedUniqueVector{T}()`: Create an empty `SortedUniqueVector{T}`.
"""
mutable struct SortedUniqueVector{T} <: AbstractSortedVector{T}
	data::Vector{T}  # Embed SortedVector inside

	# Constructors
	function SortedUniqueVector(data::AbstractVector)
		return new{eltype(data)}(unique(sort(data)))
	end
	SortedUniqueVector{T}() where T = new{T}(Vector{T}())
	SortedUniqueVector() = new{Any}(Vector{Any}())
end

# Type-specific operations (semantics differ from generic fallbacks)
==(suv1::SortedUniqueVector, suv2::SortedUniqueVector) = suv1.data == suv2.data
==(suv::SortedUniqueVector, vec::AbstractVector) = suv.data == vec
==(vec::AbstractVector, suv::SortedUniqueVector) = vec == suv.data
isless(suv1::SortedUniqueVector, suv2::SortedUniqueVector) = suv1.data < suv2.data

# Overload push! to avoid duplicates
function push!(suv::SortedUniqueVector{T}, value::T) where T
	idx = searchsortedfirst(suv.data, value)
	if idx > length(suv.data) || suv.data[idx] != value
		insert!(suv.data, idx, value)
	end
	return suv
end

# Remove a specific element if it exists
function delete!(suv::SortedUniqueVector{T}, value::T) where T
	idx = Base.findfirst(x -> x == value, suv.data)
	if !isnothing(idx)
		deleteat!(suv.data, idx)
	end
	return suv
end

function deleteat!(suv::SortedUniqueVector{T}, idx::Int) where T
	if idx < 1 || idx > length(suv.data)
		throw(BoundsError(suv, idx))  # Throw an error if index is out of bounds
	end
	deleteat!(suv.data, idx)      # Remove the element at the specified index
	return suv
end

function copy(suv::SortedUniqueVector{T}) where T
	data = copy(suv.data)
	return SortedUniqueVector(data)
end

function clear!(suv::SortedUniqueVector)
	empty!(suv.data)
	return suv
end


# ─────────────────────────────────────────────────────────────────────────────
# SortedCountingUniqueVector
# ─────────────────────────────────────────────────────────────────────────────
"""
	SortedCountingUniqueVector{T} <: AbstractSortedVector{T}

A mutable sorted vector of unique elements (like `SortedUniqueVector`) that also
tracks how many times each element was inserted.

# Fields
- `data::SortedUniqueVector{T}`: The underlying sorted collection of unique elements.
- `counts::Dict{T,Int}`: A dictionary mapping each element to its insertion count.

# Constructors
- `SortedCountingUniqueVector{T}()`: Create an empty counting vector.
- `SortedCountingUniqueVector(data::AbstractVector{T})`: Create a counting vector from `data`,
  counting each element's occurrences.
"""
mutable struct SortedCountingUniqueVector{T} <: AbstractSortedVector{T}
	data::SortedUniqueVector{T}
	counts::Dict{T, Int}

	function SortedCountingUniqueVector{T}() where T
		new{T}(SortedUniqueVector{T}(), Dict{T, Int}())
	end

	function SortedCountingUniqueVector(data::AbstractVector{T}) where T
		suv = SortedUniqueVector(data)
		cnts = Dict{T, Int}()
		for value in data
			if haskey(cnts, value)
				cnts[value] += 1
			else
				cnts[value] = 1
			end
		end
		new{T}(suv, cnts)
	end

	function SortedCountingUniqueVector(
		data::SortedUniqueVector{T},
		counts::Dict{T, Int},
	) where {T}
		new{T}(data, counts)
	end
end

# Type-specific operations (semantics differ from generic fallbacks)
function ==(suv1::SortedCountingUniqueVector, suv2::SortedCountingUniqueVector)
	return suv1.data == suv2.data && suv1.counts == suv2.counts
end
==(suv::SortedCountingUniqueVector, vec::AbstractVector) = suv.data == vec
==(vec::AbstractVector, suv::SortedCountingUniqueVector) = vec == suv.data
isless(suv1::SortedCountingUniqueVector, suv2::SortedCountingUniqueVector) =
	suv1.data < suv2.data

function push!(scv::SortedCountingUniqueVector{T}, val::T) where T
	if haskey(scv.counts, val)
		scv.counts[val] += 1
	else
		push!(scv.data, val)
		scv.counts[val] = 1
	end
	return scv
end

function push!(scv::SortedCountingUniqueVector{T}, val::T, count::Integer) where T
	if haskey(scv.counts, val)
		scv.counts[val] += count
	else
		push!(scv.data, val)
		scv.counts[val] = count
	end
	return scv
end

function delete!(scv::SortedCountingUniqueVector{T}, val::T) where T
	if haskey(scv.counts, val)
		delete!(scv.counts, val)
		delete!(scv.data, val)
	end
	return scv
end

# Override generic `in` fallback: counts-dict lookup is O(1) vs searchsortedfirst's
# O(log N), and stays accurate even if `data` and `counts` ever drift.
function in(val::T, scv::SortedCountingUniqueVector{T}) where T
	return haskey(scv.counts, val)
end

"""
	copy(scv::SortedCountingUniqueVector{T}) -> SortedCountingUniqueVector{T}

Return a new `SortedCountingUniqueVector` that is a deep copy of `scv`.
It copies both the underlying `SortedUniqueVector` and the `counts` dictionary.
"""
function copy(scv::SortedCountingUniqueVector{T}) where T
	data_copy = copy(scv.data)
	counts_copy = copy(scv.counts)
	return SortedCountingUniqueVector(data_copy, counts_copy)
end

function show(io::IO, scv::SortedCountingUniqueVector{T}) where T
	for val in scv.data
		print(io, "counts: ", scv.counts[val], " ")
		println(io, "data: ", val)
	end
end

function getcount(scv::SortedCountingUniqueVector{T}, val::T)::Int where T
	return get(scv.counts, val, 0)
end

function addcount!(scv::SortedCountingUniqueVector{T}, val::T, count::Integer) where T
	scv.counts[val] += count
end

end
