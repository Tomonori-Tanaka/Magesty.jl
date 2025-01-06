"""
	Module SortedContainer

This module provides various sorted container data structures, including:
- `SortedVector`: A mutable sorted vector allowing duplicate elements.
- `SortedUniqueVector`: A mutable sorted vector that avoids duplicates.
- `SortedCountingVector`: A mutable sorted vector that also keeps track of element counts.

All data structures implement `AbstractVector` and support operations such as indexing,
iteration, insertion, deletion, and clearing in a sorted manner.
"""
module SortedContainer

export SortedVector
export SortedUniqueVector
export SortedCountingVector, getcount

"""
	AbstractSortedVector{T} <: AbstractVector{T}

An abstract type representing a sorted vector of element type `T`.
This type indicates that any concrete subtype maintains a sorted order of its elements.
"""
abstract type AbstractSortedVector{T} <: AbstractVector{T} end



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

	# Empty constructors
	SortedVector{T}() where T = new{T}(Vector{T}())
	SortedVector() = new{Any}(Vector())
end

# Support array-like operations
Base.getindex(sv::SortedVector, i::Int) = sv.data[i]  # Allow index-based access
Base.length(sv::SortedVector) = length(sv.data)       # Return the length of the vector
Base.iterate(sv::SortedVector) = iterate(sv.data)     # Support iteration
Base.iterate(sv::SortedVector, state) = iterate(sv.data, state)
Base.isempty(sv::SortedVector) = isempty(sv.data)
Base.size(sv::SortedVector) = size(sv.data)
Base.:(==)(sv1::SortedVector, sv2::SortedVector) = sv1.data == sv2.data
Base.isless(sv1::SortedVector, sv2::SortedVector) = sv1.data < sv2.data

"""
	push!(sv::SortedVector{T}, value::T) -> SortedVector{T}

Insert `value` into the `SortedVector` while maintaining sorted order.
Returns the updated `SortedVector`.
"""
function Base.push!(sv::SortedVector{T}, value::T) where T
	idx = searchsortedfirst(sv.data, value)  # Find the insertion position
	insert!(sv.data, idx, value)            # Insert the value at the calculated position
	return sv
end

function Base.append!(sv::SortedVector{T}, new_vec::AbstractVector{T}) where T
	for value in new_vec
		push!(sv, value)
	end
	return sv
end

"""
	findfirst(sv::SortedVector{T}, value::T) -> Union{Int, Nothing}

Return the index of the first occurrence of `value` in `sv`, or `nothing`
if `value` is not found.
"""
function Base.findfirst(sv::SortedVector{T}, value::T) where T
	idx = searchsortedfirst(sv.data, value)
	return idx <= length(sv.data) && sv.data[idx] == value ? idx : nothing
end

"""
	findall(sv::SortedVector{T}, value::T) -> Vector{Int}

Return all indices where `value` occurs in `sv`. If `value` is not present,
an empty vector is returned.
"""
function Base.findall(sv::SortedVector{T}, value::T) where T
	indices = []
	start_idx = searchsortedfirst(sv.data, value)
	while start_idx <= length(sv.data) && sv.data[start_idx] == value
		push!(indices, start_idx)
		start_idx += 1
	end
	return indices
end

function Base.delete!(sv::SortedVector{T}, value::T) where T
	idx = findfirst(sv.data .== value)
	if !isnothing(idx)
		deleteat!(sv.data, idx)
	end
	return sv
end

function Base.deleteat!(sv::SortedVector{T}, idx::Int) where T
	if idx < 1 || idx > length(sv.data)
		throw(BoundsError(sv, idx))  # Throw an error if index is out of bounds
	end
	deleteat!(sv.data, idx)      # Remove the element at the specified index
	return sv
end

function deleteall!(sv::SortedVector{T}, value::T) where T
	indices = findall(sv.data .== value)
	for idx in reverse(indices)
		deleteat!(sv.data, idx)
	end
	return sv
end

function Base.in(value::T, sv::SortedVector{T}) where T
	idx = searchsortedfirst(sv.data, value)
	return idx <= length(sv.data) && sv.data[idx] == value
end

function clear!(sv::SortedVector)
	empty!(sv.data)
	return sv
end

function Base.copy(sv::SortedVector{T}) where T
	return SortedVector{T}(data = copy(sv.data))
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

Base.getindex(suv::SortedUniqueVector, i::Int) = suv.data[i]  # Allow index-based access
Base.length(suv::SortedUniqueVector) = length(suv.data)       # Return the length of the vector
Base.iterate(suv::SortedUniqueVector) = iterate(suv.data)     # Support iteration
Base.iterate(suv::SortedUniqueVector, state) = iterate(suv.data, state)
Base.isempty(suv::SortedUniqueVector) = isempty(suv.data)
Base.size(suv::SortedUniqueVector) = size(suv.data)
Base.:(==)(suv1::SortedUniqueVector, suv2::SortedUniqueVector) = suv1.data == suv2.data
Base.:(==)(suv::SortedUniqueVector, vec::AbstractVector) = suv.data == vec
Base.:(==)(vec::AbstractVector, suv::SortedUniqueVector) = vec == suv.data
Base.isless(suv1::SortedUniqueVector, suv2::SortedUniqueVector) = suv1.data < suv2.data

# Overload push! to avoid duplicates
function Base.push!(suv::SortedUniqueVector{T}, value::T) where T
	idx = searchsortedfirst(suv.data, value)
	if idx > length(suv.data) || suv.data[idx] != value
		insert!(suv.data, idx, value)
	end
	return suv
end

function Base.append!(suv::SortedUniqueVector{T}, new_vec::AbstractVector{T}) where T
	for value in new_vec
		push!(suv, value)
	end
	return suv
end

function Base.findfirst(suv::SortedUniqueVector{T}, value::T) where T
	idx = searchsortedfirst(suv.data, value)
	return idx <= length(suv.data) && suv.data[idx] == value ? idx : nothing
end

# Remove a specific element if it exists
function Base.delete!(suv::SortedUniqueVector{T}, value::T) where T
	idx = findfirst(suv.data .== value)
	if !isnothing(idx)
		deleteat!(suv.data, idx)
	end
	return suv
end

function Base.deleteat!(suv::SortedUniqueVector{T}, idx::Int) where T
	if idx < 1 || idx > length(suv.data)
		throw(BoundsError(suv, idx))  # Throw an error if index is out of bounds
	end
	deleteat!(suv.data, idx)      # Remove the element at the specified index
	return suv
end

function Base.in(value::T, suv::SortedUniqueVector{T}) where T
	idx = searchsortedfirst(suv.data, value)
	return idx <= length(suv.data) && suv.data[idx] == value
end

function Base.copy(suv::SortedUniqueVector{T}) where T
	data = copy(suv.data)
	return SortedUniqueVector(data)
end

function clear!(suv::SortedUniqueVector)
	empty!(suv.data)
	return suv
end

# ─────────────────────────────────────────────────────────────────────────────
# SortedCountingVector
# ─────────────────────────────────────────────────────────────────────────────

"""
	SortedCountingVector{T} <: AbstractSortedVector{T}

A mutable sorted vector of unique elements (like `SortedUniqueVector`) that also
tracks how many times each element was inserted.

# Fields
- `data::SortedUniqueVector{T}`: The underlying sorted collection of unique elements.
- `counts::Dict{T,Int}`: A dictionary mapping each element to its insertion count.

# Constructors
- `SortedCountingVector{T}()`: Create an empty counting vector.
- `SortedCountingVector(data::AbstractVector{T})`: Create a counting vector from `data`,
  counting each element's occurrences.
"""
mutable struct SortedCountingVector{T} <: AbstractSortedVector{T}
	data::SortedUniqueVector{T}
	counts::Dict{T, Int}

	function SortedCountingVector{T}() where T
		new{T}(SortedUniqueVector{T}(), Dict{T, Int}())
	end

	function SortedCountingVector(data::AbstractVector{T}) where T
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

	function SortedCountingVector(
		data::SortedUniqueVector{T},
		counts::Dict{T, Int},
	) where {T}
		new{T}(data, counts)
	end
end

Base.getindex(suv::SortedCountingVector, i::Int) = suv.data[i]  # Allow index-based access
Base.length(suv::SortedCountingVector) = length(suv.data)       # Return the length of the vector
Base.size(suv::SortedCountingVector) = size(suv.data)
Base.iterate(suv::SortedCountingVector) = iterate(suv.data)     # Support iteration
Base.iterate(suv::SortedCountingVector, state) = iterate(suv.data, state)
function Base.:(==)(suv1::SortedCountingVector, suv2::SortedCountingVector)
	return suv1.data == suv2.data && suv1.counts == suv2.counts
end
Base.:(==)(suv::SortedCountingVector, vec::AbstractVector) = suv.data == vec
Base.:(==)(vec::AbstractVector, suv::SortedCountingVector) = vec == suv.data
Base.isless(suv1::SortedCountingVector, suv2::SortedCountingVector) = suv1.data < suv2.data

function Base.push!(scv::SortedCountingVector{T}, val::T) where T
	if haskey(scv.counts, val)
		scv.counts[val] += 1
	else
		push!(scv.data, val)
		scv.counts[val] = 1
	end
	return scv
end


function Base.append!(scv::SortedCountingVector{T}, new_vals::AbstractVector{T}) where T
	for val in new_vals
		push!(scv, val)
	end
	return scv
end

function Base.delete!(scv::SortedCountingVector{T}, val::T) where T
	if haskey(scv.counts, val)
		delete!(scv.counts, val)
		delete!(scv.data, val)
	end
	return scv
end

function Base.in(val, scv::SortedCountingVector{T}) where T
	return haskey(scv.counts, val)
end

"""
	copy(scv::SortedCountingVector{T}) -> SortedCountingVector{T}

Return a new `SortedCountingVector` that is a deep copy of `scv`.
It copies both the underlying `SortedUniqueVector` and the `counts` dictionary.
"""
function Base.copy(scv::SortedCountingVector{T}) where T
	data_copy = copy(scv.data)
	counts_copy = copy(scv.counts)
	return SortedCountingVector(data_copy, counts_copy)
end

function Base.show(io::IO, scv::SortedCountingVector{T}) where T
	for val in scv.data
		print(io, "counts: ", scv.counts[val], " ")
		println(io, "data: ", val)
	end
end

function getcount(scv::SortedCountingVector{T}, val::T)::Int where T
	return get(scv.counts, val, 0)
end

end
