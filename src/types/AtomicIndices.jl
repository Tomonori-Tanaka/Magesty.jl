module AtomicIndices

using LinearAlgebra
import Base:
	append!, eltype, getindex, hash, in, insert!, isempty, isless, iterate, length, push!,
	searchsortedfirst, show, size, sort, ==

export Indices, IndicesUniqueList, get_total_L, get_atom_l_list,
	equivalent,
	product_indices_of_all_comb,
	product_indices,
	product_shsiteindex,
	indices_singleatom,
	shsiteindex_singleatom,
	SHSiteIndex, SHProduct, replace_atom, LinearCombo, inner_product

# ─────────────────────────────────────────────────────────────────────────────
# Indices
# ─────────────────────────────────────────────────────────────────────────────
"""
	Indices(atom::Int, l::Int, m::Int, cell::Int)

Represents the indices of an atom along with its associated quantum numbers (l and m) and cell index
in a real spherical harmonics basis.

### Fields
- `atom::Int`: The index of the atom.
- `l::Int`: The quantum number (l).
- `m::Int`: The quantum number (m).
- `cell::Int`: The cell index (1 ≤ cell ≤ 27).

### Constructor
- `Indices(atom::Int, l::Int, m::Int, cell::Int)`: Creates a new `Indices` instance after validating:
  - l must be positive
  - |m| must not exceed l
  - cell must be between 1 and 27 (inclusive)
"""
struct Indices
	atom::Int
	l::Int
	m::Int
	cell::Int

	function Indices(atom::Integer, l::Integer, m::Integer, cell::Integer)
		if l < 1
			throw(DomainError("l must be positive, got l = $l"))
		elseif abs(m) > l
			throw(DomainError("|m| must not exceed l, got m = $m for l = $l"))
		elseif cell <= 0 || cell > 27
			throw(DomainError("cell must be between 1 and 27, got cell = $cell"))
		end
		new(atom, l, m, cell)
	end
end

function Indices(tuple::NTuple{4, Integer})
	return Indices(tuple...)
end

# Comparison operations
isless(atom_1::Indices, atom_2::Indices) =
	(atom_1.atom, atom_1.l, atom_1.m, atom_1.cell) <
	(atom_2.atom, atom_2.l, atom_2.m, atom_2.cell)

isless(atom_1::Indices, tuple::NTuple{4, Integer}) =
	(atom_1.atom, atom_1.l, atom_1.m, atom_1.cell) < tuple

isless(tuple::NTuple{4, Integer}, atom_2::Indices) =
	tuple < (atom_2.atom, atom_2.l, atom_2.m, atom_2.cell)

# Equality operations
==(atom_1::Indices, atom_2::Indices) =
	(atom_1.atom, atom_1.l, atom_1.m, atom_1.cell) ==
	(atom_2.atom, atom_2.l, atom_2.m, atom_2.cell)

==(atom_1::Indices, tuple::NTuple{4, Integer}) =
	(atom_1.atom, atom_1.l, atom_1.m, atom_1.cell) == tuple

==(tuple::NTuple{4, Integer}, atom_2::Indices) =
	tuple == (atom_2.atom, atom_2.l, atom_2.m, atom_2.cell)

function Base.convert(::Type{NTuple{4, Int}}, indices::Indices)
	try
		return (
			convert(Int, indices.atom),
			convert(Int, indices.l),
			convert(Int, indices.m),
			convert(Int, indices.cell),
		)
	catch e
		throw(InexactError(:convert, NTuple{4, Int}, indices))
	end
end

# Hashing
hash(atom_1::Indices, h::UInt) =
	hash((atom_1.atom, atom_1.l, atom_1.m, atom_1.cell), h)

# Display
function show(io::IO, indices::Indices)
	print(
		io,
		"(atom: $(indices.atom), l: $(indices.l), m: $(indices.m), cell: $(indices.cell))",
	)
end

struct SHSiteIndex
	i::Int
	l::Int
	m::Int

	function SHSiteIndex(i::Integer, l::Integer, m::Integer)
		if l < 1
			throw(DomainError("l must be positive, got l = $l"))
		elseif abs(m) > l
			throw(DomainError("|m| must not exceed l, got m = $m for l = $l"))
		end
		new(convert(Int, i), convert(Int, l), convert(Int, m))
	end
end

function SHSiteIndex(tuple::NTuple{3, Integer})
	return SHSiteIndex(tuple...)
end

isless(shsi1::SHSiteIndex, shsi2::SHSiteIndex) =
	(shsi1.i, shsi1.l, shsi1.m) < (shsi2.i, shsi2.l, shsi2.m)
isless(shsi1::SHSiteIndex, tuple::NTuple{3, Integer}) = (shsi1.i, shsi1.l, shsi1.m) < tuple
isless(tuple::NTuple{3, Integer}, shsi2::SHSiteIndex) = tuple < (shsi2.i, shsi2.l, shsi2.m)

==(shsi1::SHSiteIndex, shsi2::SHSiteIndex) =
	(shsi1.i, shsi1.l, shsi1.m) == (shsi2.i, shsi2.l, shsi2.m)
==(shsi1::SHSiteIndex, tuple::NTuple{3, Integer}) = (shsi1.i, shsi1.l, shsi1.m) == tuple
==(tuple::NTuple{3, Integer}, shsi2::SHSiteIndex) = tuple == (shsi2.i, shsi2.l, shsi2.m)

hash(shsi::SHSiteIndex, h::UInt) = hash((shsi.i, shsi.l, shsi.m), h)

function Base.convert(::Type{NTuple{3, Int}}, shsi::SHSiteIndex)::NTuple{3, Int}
	try
		return (convert(Int, shsi.i), convert(Int, shsi.l), convert(Int, shsi.m))
	catch e
		throw(InexactError(:convert, NTuple{3, Int}, shsi))
	end
end

function show(io::IO, shsi::SHSiteIndex)
	print(io, "(i: $(shsi.i), l: $(shsi.l), m: $(shsi.m))")
end

# ─────────────────────────────────────────────────────────────────────────────
# IndicesUniqueList
# ─────────────────────────────────────────────────────────────────────────────
"""
	IndicesUniqueList

A mutable collection of unique `Indices` objects that implements the `AbstractVector` interface.
Maintains uniqueness of elements and provides specialized operations for handling atomic indices.
"""
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

# AbstractVector interface implementation
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
	isless(iul1::IndicesUniqueList, iul2::IndicesUniqueList) -> Bool

Compare two `IndicesUniqueList`s based on their lengths first and then lexicographically.
"""
function isless(iul1::IndicesUniqueList, iul2::IndicesUniqueList)
	return length(iul1.data) < length(iul2.data) || iul1.data < iul2.data
end

"""
	sort(iul::IndicesUniqueList) -> IndicesUniqueList

Return a new sorted `IndicesUniqueList` containing the same elements as `iul`.
"""
function sort(iul::IndicesUniqueList)
	IndicesUniqueList(sort(iul.data))
end

"""
	push!(iul::IndicesUniqueList, indices::Indices) -> IndicesUniqueList

Add `indices` to `iul` only if it does not already exist in `iul`.
"""
function push!(iul::IndicesUniqueList, indices::Indices)
	if indices in iul
		return iul
	else
		push!(iul.data, indices)
		return iul
	end
end

"""
	append!(iul::IndicesUniqueList, vec::AbstractVector{Indices}) -> IndicesUniqueList

Append all elements of `vec` to `iul`, keeping them unique.
"""
function append!(iul::IndicesUniqueList, vec::AbstractVector{Indices})
	for val in vec
		push!(iul, val)
	end
	return iul
end

# Display
function show(io::IO, iul::IndicesUniqueList)
	for indices in iul
		# (atom, l, m, cell) with fixed width
		print(
			io,
			"($(lpad(indices.atom, 5)),$(lpad(indices.l, 2)),$(lpad(indices.m, 3)), $(lpad(indices.cell, 2)))",
		)
	end
end

# ─────────────────────────────────────────────────────────────────────────────
# Getter functions
# ─────────────────────────────────────────────────────────────────────────────
"""
	get_total_L(iul::IndicesUniqueList) -> Int

Return the sum of all l values in `IndicesUniqueList`.
"""
function get_total_L(iul::IndicesUniqueList)::Int
	return sum([indices.l for indices in iul])
end

function get_atom_l_list(iul::IndicesUniqueList)::Vector{Vector{Int}}
	atom_list = [indices.atom for indices in iul]
	l_list = [indices.l for indices in iul]
	vec = Vector{Vector{Int}}()
	for (atom, l) in zip(atom_list, l_list)
		push!(vec, Int[atom, l])
	end

	return vec
end

function get_atom_l_m_list(iul::IndicesUniqueList)::Vector{Vector{Int}}
	atom_list = [indices.atom for indices in iul]
	l_list = [indices.l for indices in iul]
	m_list = [indices.m for indices in iul]
	vec = Vector{Vector{Int}}()
	for (atom, l, m) in zip(atom_list, l_list, m_list)
		push!(vec, Int[atom, l, m])
	end
	return vec
end

"""
	equivalent(iul1::IndicesUniqueList, iul2::IndicesUniqueList) -> Bool

Determine whether two `IndicesUniqueList`s are equivalent, i.e., contain the same elements
regardless of order.

# Arguments
- `iul1::IndicesUniqueList`: First `IndicesUniqueList` to compare.
- `iul2::IndicesUniqueList`: Second `IndicesUniqueList` to compare.

# Returns
- `true` if the lists contain the same elements, `false` otherwise.
"""
function equivalent(iul1::IndicesUniqueList, iul2::IndicesUniqueList)::Bool
	return sort(iul1) == sort(iul2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Functions for generating Indices and their combinations
# ─────────────────────────────────────────────────────────────────────────────
"""
	product_indices_of_all_comb(atom_list::AbstractVector{<:Integer},
				   lmax_list::AbstractVector{<:Integer},
				   cell_list::AbstractVector{<:Integer}) -> Vector{IndicesUniqueList}

Generate all possible combinations of `Indices` for given atoms, maximum l values, and cells.

# Arguments
- `atom_list::AbstractVector{<:Integer}`: List of atom indices.
- `lmax_list::AbstractVector{<:Integer}`: List of maximum l values for each atom.
- `cell_list::AbstractVector{<:Integer}`: List of cell indices for each atom.

# Returns
- A Vector of `IndicesUniqueList` containing all possible combinations.

# Throws
- `ErrorException` if the input vectors have different lengths.
"""
function product_indices_of_all_comb(
	atom_list::AbstractVector{<:Integer},
	lmax_list::AbstractVector{<:Integer},
	cell_list::AbstractVector{<:Integer},
)::Vector{IndicesUniqueList}
	if length(atom_list) != length(lmax_list) != length(cell_list)
		throw(
			ErrorException(
				"Input vectors must have the same length. Got lengths: " *
				"atom_list=$(length(atom_list)), " *
				"lmax_list=$(length(lmax_list)), " *
				"cell_list=$(length(cell_list))",
			),
		)
	end

	list_tmp = Vector{Vector{Indices}}()
	for (atom, lmax, cell) in zip(atom_list, lmax_list, cell_list)
		singleatom_list = Vector{Indices}()
		for l in 1:lmax
			append!(singleatom_list, indices_singleatom(atom, l, cell))
		end
		# Examle for (atom=3, lmax=2, cell=1),
		# singleatom_list = [Indices(3, 1, -1, 1), Indices(3, 1, 0, 1), Indices(3, 1, 1, 1), Indices(3, 2, -2, 1), Indices(3, 2, -1, 1), ..., Indices(3, 2, 2, 1)]
		# the length of singleatom_list in this case is 8
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

function product_indices(
	atom_list::AbstractVector{<:Integer},
	l_list::AbstractVector{<:Integer},
	cell_list::AbstractVector{<:Integer},
)::Vector{IndicesUniqueList}
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

function product_shsiteindex(
	atom_list::AbstractVector{<:Integer},
	l_list::AbstractVector{<:Integer},
)::Vector{SHProduct}
	if length(atom_list) != length(l_list)
		throw(
			ErrorException(
				"Input vectors must have the same length. Got lengths: " *
				"atom_list=$(length(atom_list)), " *
				"l_list=$(length(l_list))",
			),
		)
	end

	list_tmp = Vector{Vector{SHSiteIndex}}()
	for (atom, l) in zip(atom_list, l_list)
		push!(list_tmp, shsiteindex_singleatom(atom, l))
	end

	combined_vec = Vector{SHProduct}()
	prod_iter = Iterators.product(list_tmp...)
	for comb::Tuple{Vararg{SHSiteIndex}} in prod_iter
		shp_tmp = SHProduct()
		for shsi in comb
			push!(shp_tmp, shsi)
		end
		push!(combined_vec, shp_tmp)
	end
	return combined_vec
end

"""
	indices_singleatom(atom::Integer, l::Integer, cell::Integer) -> Vector{Indices}

Generate all possible `Indices` for a single atom with given l value and cell index.

# Arguments
- `atom::Integer`: The atom index.
- `l::Integer`: The quantum number l.
- `cell::Integer`: The cell index.

# Returns
- A Vector of `Indices` containing all possible combinations for the given parameters.
"""
function indices_singleatom(atom::Integer, l::Integer, cell::Integer)::Vector{Indices}
	vec = Vector{Indices}()
	for m in (-l):l
		push!(vec, Indices(atom, l, m, cell))
	end
	return sort(vec)
end

function shsiteindex_singleatom(atom::Integer, l::Integer)::Vector{SHSiteIndex}
	vec = Vector{SHSiteIndex}()
	for m in (-l):l
		push!(vec, SHSiteIndex(atom, l, m))
	end
	return sort(vec)
end





struct SHProduct
	data::Vector{SHSiteIndex}

	function SHProduct()
		new(Vector{SHSiteIndex}())
	end

	function SHProduct(vec::AbstractVector{SHSiteIndex})
		if length(vec) != length(unique(vec))
			throw(ErrorException("SHSiteIndex vector must be unique."))
		end
		new(vec)
	end

	function SHProduct(shsi::SHSiteIndex)
		new(Vector{SHSiteIndex}([shsi]))
	end
end

getindex(shp::SHProduct, idx::Int) = shp.data[idx]
length(shp::SHProduct) = length(shp.data)
eltype(::Type{SHProduct}) = SHSiteIndex
hash(shp::SHProduct, h::UInt) = hash(shp.data, h)
iterate(shp::SHProduct, state::Int = 1) = iterate(shp.data, state)
isempty(shp::SHProduct) = isempty(shp.data)
size(shp::SHProduct) = size(shp.data)
==(shp1::SHProduct, shp2::SHProduct) = shp1.data == shp2.data
in(val::SHSiteIndex, shp::SHProduct) = val in shp.data

function isless(shp1::SHProduct, shp2::SHProduct)
	return length(shp1.data) < length(shp2.data) || shp1.data < shp2.data
end

function sort(shp::SHProduct)
	return SHProduct(sort(shp.data))
end

function push!(shp::SHProduct, shsi::SHSiteIndex)
	if shsi in shp
		throw(ErrorException("SHSiteIndex must be unique."))
	end
	push!(shp.data, shsi)
	return shp
end

function append!(shp::SHProduct, vec::AbstractVector{SHSiteIndex})
	if length(vec) != length(unique(vec))
		throw(ErrorException("SHSiteIndex vector must be unique."))
	end
	for shsi in vec
		push!(shp, shsi)
	end
	return shp
end

function get_atom_l_list(shp::SHProduct)::Vector{Vector{Int}}
	atom_list = [shsi.i for shsi in shp.data]
	l_list = [shsi.l for shsi in shp.data]
	vec = Vector{Vector{Int}}()
	for (atom, l) in zip(atom_list, l_list)
		push!(vec, Int[atom, l])
	end
	return vec
end

function replace_atom(shp::SHProduct, atom_list::AbstractVector{<:Integer})
	if length(shp) != length(atom_list)
		throw(ErrorException("Lengths of shp and atom_list must be the same."))
	end
	new_shp = SHProduct()
	for (i, shsi) in enumerate(shp)
		push!(new_shp, SHSiteIndex(atom_list[i], shsi.l, shsi.m))
	end
	return new_shp
end


struct LinearCombo
	data::Vector{SHProduct}
	coeffs::Vector{Float64}
	multiplicity::Vector{Int}

	function LinearCombo(data::Vector{SHProduct}, coeffs::Vector{Float64}, multiplicity::Vector{Int})
		if length(data) != length(coeffs) != length(multiplicity)
			throw(ErrorException("Lengths of data and coeffs must be the same."))
		end
		if !isapprox(norm(coeffs), 1.0, atol = 1e-8)
			throw(ErrorException("The norm of the coefficient vector is not 1."))
		end
		new(Vector{SHProduct}(data), Vector{Float64}(coeffs), Vector{Int}(multiplicity))
	end
end

function LinearCombo(data::Vector{SHProduct}, coeffs::Vector{<:Real})
	if length(data) != length(coeffs)
		throw(ErrorException("Lengths of data and coeffs must be the same."))
	end
	multiplicity = [1 for _ in data]
	return LinearCombo(data, coeffs, multiplicity)
end

function inner_product(shp::SHProduct, lco::LinearCombo)::Float64
	for shp_i in lco.data
		if sort(shp) == sort(shp_i)
			return lco.coeffs[findfirst(x -> x == shp_i, lco.data)]
		end
	end
	return 0.0

end

end
