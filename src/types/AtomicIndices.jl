module AtomicIndices

import Base:
	append!, eltype, getindex, hash, in, isempty, isless, iterate, length, push!, show,
	size, sort, ==

export Indices, IndicesUniqueList, get_atomlist, get_llist, get_celllist, get_totalL, get_atom_l_list,
	equivalent,
	product_indices,
	indices_singleatom

# ─────────────────────────────────────────────────────────────────────────────
# Indices
# ─────────────────────────────────────────────────────────────────────────────
"""
	Indices(atom::Int, l::Int, m::Int)

Represents the indices of an atom along with its associated labels (l and m) in a real spherical harmonics basis.

### Fields
- `atom::Int`: The index of the atom.
- `l::Int`: The quantum number ( l ).
- `m::Int`: The quantum number ( m ).

### Constructor
- `Indices(atom::Int,  l::Int, m::Int)`: Creates a new `Indices` instance after validating that ( l ) is positive and ( |m| \\leq l ).
"""
struct Indices
	atom::Int
	l::Int
	m::Int
	cell::Int
end

function Indices(atom::Integer, l::Integer, m::Integer, cell::Integer)
	if l < 1
		throw(DomainError("l: $l, l should be positive."))
	elseif abs(m) > l
		throw(DomainError("l: $l, m: $m, |m| must not exceed l."))
	elseif cell <= 0 || cell > 27
		throw(DomainError("cell: $cell, it should be in 1 <= cell <= 27. "))
	end
	return Indices(atom, l, m, cell)
end

function Indices(tuple::NTuple{4, Integer})
	return Indices(tuple...)
end

isless(atom_1::Indices, atom_2::Indices) =
	(atom_1.atom, atom_1.l, atom_1.m, atom_1.cell) <
	(atom_2.atom, atom_2.l, atom_2.m, atom_2.cell)
isless(atom_1::Indices, tuple::NTuple{4, Integer}) =
	(atom_1.atom, atom_1.l, atom_1.m, atom_1.cell) < tuple
isless(tuple::NTuple{4, Integer}, atom_2::Indices) =
	tuple < (atom_2.atom, atom_2.l, atom_2.m, atom_2.cell)
==(atom_1::Indices, atom_2::Indices) =
	(atom_1.atom, atom_1.l, atom_1.m, atom_1.cell) ==
	(atom_2.atom, atom_2.l, atom_2.m, atom_2.cell)
==(atom_1::Indices, tuple::NTuple{4, Integer}) =
	(atom_1.atom, atom_1.l, atom_1.m, atom_1.cell) == tuple
==(tuple::NTuple{4, Integer}, atom_2::Indices) =
	tuple == (atom_2.atom, atom_2.l, atom_2.m, atom_2.cell)
hash(atom_1::Indices, h::UInt) =
	hash((atom_1.atom, atom_1.l, atom_1.m, atom_1.cell), h)
function show(io::IO, indices::Indices)
	print(
		io,
		"(atom: $(indices.atom), l: $(indices.l), m: $(indices.m), cell: $(indices.cell))",
	)
end


# ─────────────────────────────────────────────────────────────────────────────
# IndicesUniqueList
# ─────────────────────────────────────────────────────────────────────────────
mutable struct IndicesUniqueList <: AbstractVector{Indices}
	data::Vector{Indices}

	function IndicesUniqueList()
		new(Vector{Indices}())
	end

	function IndicesUniqueList(vec::AbstractVector{Indices})
		new(unique(vec))
	end

	function IndicesUniqueList(indices::Indices)
		new(Vector{Indices}([indices]))
	end
end

getindex(iul::IndicesUniqueList, idx::Int) = iul.data[idx]
length(iul::IndicesUniqueList) = length(iul.data)
eltype(::Type{IndicesUniqueList}) = Indices
hash(iul::IndicesUniqueList, h::UInt) = hash(iul.data, h)
iterate(iul::IndicesUniqueList, state::Int = 1) = iterate(iul.data, state)
isempty(iul::IndicesUniqueList) = isempty(iul.data)
size(iul::IndicesUniqueList) = size(iul.data)
==(iul1::IndicesUniqueList, iul2::IndicesUniqueList) = iul1.data == iul2.data
in(val::Indices, iul::IndicesUniqueList) = val in iul.data

"""
	isless(iul1::IndicesUniqueList, iul2::IndicesUniqueList)

Compare two `IndicesUniqueList`s based on their lengths first and then lexicographically.
"""
function isless(iul1::IndicesUniqueList, iul2::IndicesUniqueList)
	return length(iul1.data) < length(iul2.data) || iul1.data < iul2.data
end

function sort(iul::IndicesUniqueList)
	IndicesUniqueList(sort(iul.data))
end

"""
	push!(iul::IndicesUniqueList, indices::Indices)

Add `indices` to `iul` only if it does not already exist in `iul`.
"""
function push!(iul::IndicesUniqueList, indices::Indices)
	# guarantees the uniqueness
	if indices in iul
		return iul
	else
		push!(iul.data, indices)
		return iul
	end
end

"""
	append!(iul::IndicesUniqueList, vec::AbstractVector{Indices})

Append all elements of `vec` to `iul`, keeping them unique.
"""
function append!(iul::IndicesUniqueList, vec::AbstractVector{Indices})
	for val in vec
		push!(iul, val)
	end
	return iul
end

# Show method for IndicesUniqueList
function show(io::IO, iul::IndicesUniqueList)
	for indices in iul
		print(
			io,
			"(atom: $(indices.atom), l: $(indices.l), m: $(indices.m), cell: $(indices.cell))",
		)
	end
end

# ─────────────────────────────────────────────────────────────────────────────
# Getter functions
# ─────────────────────────────────────────────────────────────────────────────
"""
	get_atomlist(iul::IndicesUniqueList) -> Vector{Int}

Extracts the atom indices from `IndicesUniqueList` and returns them as a Vector.

# Arguments
- `iul::IndicesUniqueList`: The `IndicesUniqueList` instance.

# Returns
- A Vector containing the atom indices.
"""
function get_atomlist(iul::IndicesUniqueList)::Vector{Int}
	[indices.atom for indices in iul]
end

"""
	get_llist(iul::IndicesUniqueList) -> Vector{Int}

Extracts the l indices from `IndicesUniqueList` and returns them as a Vector.

# Arguments
- `iul::IndicesUniqueList`: The `IndicesUniqueList` instance.

# Returns
- A Vector containing the l indices.
"""
function get_llist(iul::IndicesUniqueList)::Vector{Int}
	[indices.l for indices in iul]
end

"""
	get_celllist(iul::IndicesUniqueList) -> Vector{Int}

Extracts the cell indices from `IndicesUniqueList` and returns them as a Vector.

# Arguments
- `iul::IndicesUniqueList`: The `IndicesUniqueList` instance.

# Returns
- A Vector containing the cell indices.
"""
function get_celllist(iul::IndicesUniqueList)::Vector{Int}
	[indices.cell for indices in iul]
end


"""
	get_totalL(iul::IndicesUniqueList) -> Int

Returns the sum of all l values in `IndicesUniqueList`.
"""
function get_totalL(iul::IndicesUniqueList)::Int
	return sum(get_llist(iul))
end

function get_atom_l_list(iul::IndicesUniqueList)::Vector{Vector{Int}}
	atom_list = get_atomlist(iul)
	l_list = get_llist(iul)
	vec = Vector{Vector{Int}}()
	for (atom, l) in zip(atom_list, l_list)
		push!(vec, Int[atom, l])
	end

	return vec
end

"""
	equivalent(iul1::IndicesUniqueList, iul2::IndicesUniqueList) -> Bool
judge whether given 2 IndicesUniqueList are equivalent or not.
"""
function equivalent(iul1::IndicesUniqueList, iul2::IndicesUniqueList)::Bool
	return sort(iul1) == sort(iul2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Functions for generating Indices and their combinations
# ─────────────────────────────────────────────────────────────────────────────

function product_indices(
	atom_list::AbstractVector{<:Integer},
	lmax_list::AbstractVector{<:Integer},
	cell_list::AbstractVector{<:Integer},
)::Vector{IndicesUniqueList}
	if length(atom_list) != length(lmax_list) != length(cell_list)
		error(
			"Different vector lengths detected. atom_list, l_list, and cell_list must have the same length.",
		)
	end

	list_tmp = Vector{Vector{Indices}}()
	for (atom, lmax, cell) in zip(atom_list, lmax_list, cell_list)
		singleatom_list = Vector{Indices}()
		for l in 1:lmax
			append!(singleatom_list, indices_singleatom(atom, l, cell))
		end
		# Examle for (atom=3, lmax=2),
		# singleatom_list = [Indices(3, 1, -1), Indices(3, 1, 0), Indices(3, 1, 1), Indices(3, 2, -2), Indices(3, 2, -1), ..., Indices(3, 2, 2)]
		push!(list_tmp, singleatom_list)
	end

	# 2. Take the Cartesian product of all IndicesUniqueList objects in vec_tmp.
	#    Each element in the product is a Tuple{Vararg{Indices}}.
	combined_vec = Vector{IndicesUniqueList}()
	prod_iter = Iterators.product(list_tmp...)
	for comb::Tuple{Vararg{Indices}} in prod_iter
		# 3. Build a new IndicesUniqueList from each tuple of Indices.
		iul_tmp = IndicesUniqueList()
		for ind in comb
			push!(iul_tmp, ind)
		end
		push!(combined_vec, iul_tmp)
	end
	return sort(combined_vec)
end

function product_indices_fixed_l(
	atom_list::AbstractVector{<:Integer},
	l_list::AbstractVector{<:Integer},
	cell_list::AbstractVector{<:Integer},
)
	if length(atom_list) != length(l_list) != length(cell_list)
		error(
			"Different vector lengths detected. atom_list, l_list, and cell_list must have the same length.",
		)
	end
	list_tmp = Vector{Vector{Indices}}()
	for (atom, l, cell) in zip(atom_list, l_list, cell_list)
		singleatom_list = Vector{Indices}()
		append!(singleatom_list, indices_singleatom(atom, l, cell))
		push!(list_tmp, singleatom_list)
	end

	combined_vec = Vector{IndicesUniqueList}()
	prod_iter = Iterators.product(list_tmp...)
	for comb::Tuple{Vararg{Indices}} in prod_iter
		iul_tmp = IndicesUniqueList()
		for ind in comb
			push!(iul_tmp, ind)
		end
		push!(combined_vec, iul_tmp)
	end
	return sort(combined_vec)
end

"""
	indices_singleatom(atom::Integer,  l::Integer, cell::Integer) :: Vector{Indices}

"""
function indices_singleatom(atom::Integer, l::Integer, cell::Integer)::Vector{Indices}
	vec = Vector{Indices}()
	for m in -l:l
		push!(vec, Indices(atom, l, m, cell))
	end
	return sort(vec)
end

end
