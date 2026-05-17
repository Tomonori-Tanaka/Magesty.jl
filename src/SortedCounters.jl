"""
	module SortedCounters

Internal helper that pairs an unordered `Dict{T, Int}` count store with a
lazily-sorted view of its keys. Used by `Clusters.jl`, `SALCBases.jl`, and
`XMLIO.jl` for the irreducible cluster / coupled basis containers that
must iterate in sorted key order but also support O(1) count lookup.

Build phase: `push!(sc, k)` and `push!(sc, k, n)` mutate `counts` directly
(O(1) amortized). Read phase: `iterate`, `getindex`, and `length` see the
sorted-unique-key view; the sort runs once after the next push and stays
cached until the following mutation.

The container is **not** an `AbstractVector`; only the small subset of the
`AbstractVector` API actually used by callers (`length`, `isempty`,
`getindex(::Int)`, `iterate`) is implemented.
"""
module SortedCounters

import Base: ==, copy, eltype, firstindex, getindex, isempty, isless, iterate,
	lastindex, length, push!, show

export SortedCounter

"""
	SortedCounter{T}

Counter mapping unique keys of type `T` to insertion multiplicities, with
ordered iteration over keys.

# Constructors
- `SortedCounter{T}()`: create an empty counter.

# Returns
- `SortedCounter{T}`: empty counter whose `counts::Dict{T, Int}` field
  accumulates per-key multiplicities and whose iteration yields the unique
  keys in `isless` order.

# Examples
```julia
sc = SortedCounter{Int}()
push!(sc, 3)
push!(sc, 1, 4)
push!(sc, 3)
collect(sc) == [1, 3]
sc.counts[3] == 2
sc.counts[1] == 4
```
"""
mutable struct SortedCounter{T}
	counts::Dict{T, Int}
	_sorted_keys::Vector{T}
	_dirty::Bool

	SortedCounter{T}() where {T} = new{T}(Dict{T, Int}(), T[], false)

	function SortedCounter{T}(counts::Dict{T, Int}) where {T}
		return new{T}(counts, T[], true)
	end
end

"""
	push!(sc::SortedCounter{T}, key::T) -> SortedCounter{T}

Increment the count of `key` by 1, inserting it if absent.
"""
function push!(sc::SortedCounter{T}, key::T) where {T}
	if haskey(sc.counts, key)
		sc.counts[key] += 1
	else
		sc.counts[key] = 1
		sc._dirty = true
	end
	return sc
end

"""
	push!(sc::SortedCounter{T}, key::T, n::Integer) -> SortedCounter{T}

Increment the count of `key` by `n`, inserting it with count `n` if absent.
"""
function push!(sc::SortedCounter{T}, key::T, n::Integer) where {T}
	if haskey(sc.counts, key)
		sc.counts[key] += n
	else
		sc.counts[key] = Int(n)
		sc._dirty = true
	end
	return sc
end

"""
	_refresh_sorted!(sc::SortedCounter{T}) -> Vector{T}

Refill `_sorted_keys` with the current unique keys in `isless` order if the
cache is dirty, then return it. Internal helper.
"""
function _refresh_sorted!(sc::SortedCounter{T})::Vector{T} where {T}
	if sc._dirty
		resize!(sc._sorted_keys, length(sc.counts))
		i = 0
		for k in keys(sc.counts)
			i += 1
			sc._sorted_keys[i] = k
		end
		sort!(sc._sorted_keys)
		sc._dirty = false
	end
	return sc._sorted_keys
end

length(sc::SortedCounter) = length(sc.counts)
isempty(sc::SortedCounter) = isempty(sc.counts)
eltype(::Type{SortedCounter{T}}) where {T} = T
firstindex(sc::SortedCounter) = 1
lastindex(sc::SortedCounter) = length(sc)
getindex(sc::SortedCounter, i::Int) = _refresh_sorted!(sc)[i]
iterate(sc::SortedCounter) = iterate(_refresh_sorted!(sc))
iterate(sc::SortedCounter, state) = iterate(_refresh_sorted!(sc), state)

==(a::SortedCounter, b::SortedCounter) = a.counts == b.counts

function isless(a::SortedCounter{T}, b::SortedCounter{T}) where {T}
	return _refresh_sorted!(a) < _refresh_sorted!(b)
end

function copy(sc::SortedCounter{T}) where {T}
	out = SortedCounter{T}(copy(sc.counts))
	if !sc._dirty
		out._sorted_keys = copy(sc._sorted_keys)
		out._dirty = false
	end
	return out
end

function show(io::IO, sc::SortedCounter{T}) where {T}
	for k in _refresh_sorted!(sc)
		print(io, "counts: ", sc.counts[k], " ")
		println(io, "data: ", k)
	end
end

end
