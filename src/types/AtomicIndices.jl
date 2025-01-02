module AtomicIndices

include("../common/SortedContainer.jl")
using .SortedContainer

export Indices,
	IndicesUniqueList, product_indices, indices_singleatom, getatoms, gettotall

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

function IndicesUniqueList(data::AbstractVector{Indices})
	return IndicesUniqueList(SortedUniqueVector{Indices}(data))
end

"""
	getatoms(iul::IndicesUniqueList) -> Tuple{Vararg{Int}}

Extracts the atom indices from `IndicesUniqueList` and returns them as a tuple.

# Arguments
- `iul::IndicesUniqueList`: The `IndicesUniqueList` instance.

# Returns
- A tuple containing the atom indices.
"""
function getatoms(iul::AbstractVector)::Tuple{Vararg{Int}}
	return Tuple([indices.atom for indices in iul])
end

"""
	gettotall(iul::AbstractVector) -> Int

Extracts the total angular momentum index (L) from `IndicesUniqueList`

# Arguments
- `iul::AbstractVector`: The `IndicesUniqueList` instance.

# Returns
- total angular momentum index (L)::Int
"""
function total_l(iul::AbstractVector)::Int
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

function product_indices(
	iat_vec::AbstractVector{Int},
	cells::AbstractVector{Int},
	maxl::AbstractVector{Int},
)::Vector{IndicesUniqueList}
	if length(iat_vec) != length(maxl) != length(cells)
		error(
			"Different vector lengths detected. iat_vec, maxl, and cells must have the same length.",
		)
	end
	# Create a vector to hold IndicesUniqueList of Indices for each atom
	vec_tmp = Vector{IndicesUniqueList}()
	for (atom, l, cell) in zip(iat_vec, maxl, cells)
		push!(vec_tmp, indices_singleatom(atom, l, cell))
	end


	combined_vec = Vector{IndicesUniqueList}()
	prod_iter = Iterators.product(vec_tmp...)
	for comb::Tuple{Vararg{Indices}} in prod_iter
		isvu_tmp = IndicesUniqueList()
		for ind in comb
			push!(isvu_tmp, ind)
		end
		push!(combined_vec, isvu_tmp)
	end
	return combined_vec
end

function indices_singleatom(atom::Int, maxl::Int, cell::Int)::IndicesUniqueList
	indices_set = IndicesUniqueList()
	for l in 1:maxl
		for m in -l:l
			push!(indices_set, Indices(atom, cell, l, m))
		end
	end
	return indices_set
end

end
