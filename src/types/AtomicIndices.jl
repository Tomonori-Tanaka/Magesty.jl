module AtomicIndices

using ..SortedContainer

export Indices,
	IndicesUniqueList, getatoms, getls, gettotall

# ─────────────────────────────────────────────────────────────────────────────
# Indices
# ─────────────────────────────────────────────────────────────────────────────
"""
	Indices(atom::Int, cell::Int, l::Int, m::Int)

Represents the indices of an atom along with its associated labels (l and m) in a real spherical harmonics basis.

### Fields
- `atom::Int`: The index of the atom.
- `cell::Int` : The index of the imaginary cell
- `l::Int`: The quantum number ( l ).
- `m::Int`: The quantum number ( m ).

### Constructor
- `Indices(atom::Int, cell::Int, l::Int, m::Int)`: Creates a new `Indices` instance after validating that ( l ) is positive and ( |m| \\leq l ).
"""
struct Indices
	atom::Int
	cell::Int
	l::Int
	m::Int

	function Indices(atom::Integer, cell::Integer, l::Integer, m::Integer)
		if l < 1
			error("negative or 0 l value is detected.")
		elseif abs(m) > l
			error("|m| exceeds l.")
		end
		return new(atom, cell, l, m)
	end

	function Indices(tuple::NTuple{4, Integer})
		return new(tuple...)
	end
end

Base.isless(atom_i::Indices, atom_j::Indices) =
	(atom_i.atom, atom_i.cell, atom_i.l, atom_i.m) <
	(atom_j.atom, atom_j.cell, atom_j.l, atom_j.m)
Base.isless(atom_i::Indices, tuple::NTuple{4, Integer}) =
	(atom_i.atom, atom_i.cell, atom_i.l, atom_i.m) < tuple
Base.isless(tuple::NTuple{4, Integer}, atom_j::Indices) =
	tuple < (atom_j.atom, atom_j.cell, atom_j.l, atom_j.m)
Base.:(==)(atom_i::Indices, atom_j::Indices) =
	(atom_i.atom, atom_i.cell, atom_i.l, atom_i.m) ==
	(atom_j.atom, atom_j.cell, atom_j.l, atom_j.m)
Base.:(==)(atom_i::Indices, tuple::NTuple{4, Integer}) =
	(atom_i.atom, atom_i.cell, atom_i.l, atom_i.m) == tuple
Base.:(==)(tuple::NTuple{4, Integer}, atom_j::Indices) =
	tuple == (atom_j.atom, atom_j.cell, atom_j.l, atom_j.m)
Base.hash(atom_i::Indices, h::UInt) =
	hash((atom_i.atom, atom_i.cell, atom_i.l, atom_i.m), h)
function Base.show(io::IO, indices::Indices)
	print(
		io,
		"(atom: $(indices.atom), cell: $(indices.cell), l: $(indices.l), m: $(indices.m))",
	)
end

# ─────────────────────────────────────────────────────────────────────────────
# IndicesUniqueList
# ─────────────────────────────────────────────────────────────────────────────
mutable struct IndicesUniqueList <: AbstractVector{Indices}
	data::SortedUniqueVector{Indices}
end

function IndicesUniqueList()
	return IndicesUniqueList(SortedUniqueVector{Indices}())
end

function IndicesUniqueList(indices::Indices)
	return IndicesUniqueList(SortedUniqueVector([indices]))
end

function IndicesUniqueList(data::AbstractVector{Indices})
	return IndicesUniqueList(SortedUniqueVector(data))
end

"""
	getatoms(iul::IndicesUniqueList) -> Vector{Int}

Extracts the atom indices from `IndicesUniqueList` and returns them as a tuple.

# Arguments
- `iul::IndicesUniqueList`: The `IndicesUniqueList` instance.

# Returns
- A tuple containing the atom indices.
"""
function getatoms(iul::IndicesUniqueList)::Vector{Int}
	return [indices.atom for indices in iul]
end

function getls(iul::IndicesUniqueList)::Vector{Int}
	return [indices.l for indices in iul]
end

"""
	gettotall(iul::AbstractVector) -> Int

Extracts the total angular momentum index (L) from `IndicesUniqueList`

# Arguments
- `iul::AbstractVector`: The `IndicesUniqueList` instance.

# Returns
- total angular momentum index (L)::Int
"""
function gettotall(iul::IndicesUniqueList)::Int
	return sum([indices.l for indices in iul])
end

Base.:(==)(iul1::IndicesUniqueList, iul2::IndicesUniqueList) = iul1.data == iul2.data

function Base.isless(iul1::IndicesUniqueList, iul2::IndicesUniqueList)
	return length(iul1.data) < length(iul2.data) || iul1.data < iul2.data
end

function Base.push!(iul::IndicesUniqueList, indices::Indices)
	push!(iul.data, indices)
	return iul
end

function Base.append!(iul::IndicesUniqueList, elems::AbstractVector{Indices})
	append!(iul.data, elems)
	return iul
end

function Base.show(io::IO, iul::IndicesUniqueList)
	for indices in iul
		print(
			io,
			"(atom: $(indices.atom), cell: $(indices.cell), l: $(indices.l), m: $(indices.m))",
		)
	end
end

Base.getindex(iul::IndicesUniqueList, i::Int) = iul.data[i]
Base.length(iul::IndicesUniqueList) = length(iul.data)
Base.eltype(::Type{IndicesUniqueList}) = Indices
Base.hash(iul::IndicesUniqueList, h::UInt) = hash(iul.data, h)
Base.iterate(iul::IndicesUniqueList, state::Int = 1) = iterate(iul.data, state)
Base.isempty(iul::IndicesUniqueList) = isempty(iul.data)
Base.size(iul::IndicesUniqueList) = size(iul.data)


"""
	product_indices(
		atoms::AbstractVector{Int},
		cells::AbstractVector{Int},
		maxl::AbstractVector{Int},
	) :: Vector{IndicesUniqueList}

Generates all possible combinations (Cartesian product) of `Indices` objects across multiple
atoms, cells, and maximum angular momentum (l) values. For each triple `(atom, cell, l)` in
the input arrays, a corresponding `IndicesUniqueList` is built via `indices_singleatom`.
The function then takes the Cartesian product of all these lists, forming a new
`IndicesUniqueList` for each unique combination.

# Arguments
- `atoms::AbstractVector{Int}`: A vector of atom indices.
- `cells::AbstractVector{Int}`: A vector of cell indices (same length as `atoms`).
- `maxl::AbstractVector{Int}`: A vector of maximum `l` values (same length as `atoms`).

# Returns
- A `Vector{IndicesUniqueList}` where each element contains a unique combination of
  `Indices` formed from the Cartesian product.

# Throws
- An `error` if `atoms`, `cells`, and `maxl` do not have the same length.

# Examples
```julia-repl
julia> atoms = [1, 2]
julia> cells = [2, 3]
julia> maxl  = [1, 2]

julia> combos = product_indices(atoms, cells, maxl)
# This will generate every combination of IndicesUniqueList from the input arrays.
# [(atom: 1, cell: 2, l: 1, m: -1)(atom: 2, cell: 3, l: 2, m: -2), (atom: 1, cell: 2, l: 1, m: -1)(atom: 2, cell: 3, l: 2, m: -1), ... ]

julia> length(combos)
15  # Number of unique (l,m) combinations across the product
"""
function product_indices(
	atoms::AbstractVector{Int},
	cells::AbstractVector{Int},
	maxl::AbstractVector{Int},
)::Vector{IndicesUniqueList}
	if (length(atoms) != length(cells)) || (length(atoms) != length(maxl))
		error(
			"Different vector lengths detected. atoms, maxl, and cells must have the same length.",
		)
	end

	# 1. Build a temporary vector of IndicesUniqueList,
	#    one for each triple (atom, cell, l) in zip(atoms, cells, maxl).
	vec_tmp = [indices_singleatom(a, c, l) for (a, c, l) in zip(atoms, cells, maxl)]

	# 2. Take the Cartesian product of all IndicesUniqueList objects in vec_tmp.
	#    Each element in the product is a Tuple{Vararg{Indices}}.
	combined_vec = Vector{IndicesUniqueList}()
	prod_iter = Iterators.product(vec_tmp...)
	for comb::Tuple{Vararg{Indices}} in prod_iter
		# 3. Build a new IndicesUniqueList from each tuple of Indices.
		iul_tmp = IndicesUniqueList()
		for ind in comb
			push!(iul_tmp, ind)
		end
		push!(combined_vec, iul_tmp)
	end
	return combined_vec
end

"""
	indices_singleatom(atom::Int, cell::Int, maxl::Int) :: IndicesUniqueList

Constructs an `IndicesUniqueList` for a single atom with index `atom`, using quantum
numbers `l` and `m` in the range `1 <= l <= maxl` and `-l <= m <= l`. Each
`Indices` object includes the specified atom index (`atom`), cell index (`cell`),
and the `(l, m)` pair.

# Arguments
- `atom::Int`: The index of the atom.
- `cell::Int`: The imaginary cell index.
- `maxl::Int`: The angular momentum quantum number `l`.

# Returns
- An `IndicesUniqueList` containing all valid `Indices` for the specified ranges.

 Examples

```julia-repl
julia> # Example 1: atom = 2, cell = 1, lmax = 1
	   iul1 = indices_singleatom(2, 1, 1)
(atom: 2, cell: 1, l: 1, m: -1)(atom: 2, cell: 1, l: 1, m: 0)(atom: 2, cell: 1, l: 1, m: 1)

julia> length(iul1)
3

julia> # Example 2: atom = 3, cell = 1, lmax = 2
	   iul2 = indices_singleatom(3, 1, 2)
(atom: 3, cell: 1, l: 2, m: -2)(atom: 3, cell: 1, l: 2, m: -1)(atom: 3, cell: 1, l: 2, m: 0)(atom: 3, cell: 1, l: 2, m: 1)(atom: 3, cell: 1, l: 2, m: 2)

julia> length(iul2)
5
"""
function indices_singleatom(atom::Int, cell::Int, maxl::Int)::IndicesUniqueList
	iul = IndicesUniqueList()
	for l in 1:maxl
		for m in -l:l
			push!(iul, Indices(atom, cell, l, m))
		end
	end
	return iul
end

end
