module AtomicIndices

using LinearAlgebra
import Base:
	append!, eltype, getindex, hash, in, isempty, isless, iterate, length, push!,
	show, size, sort, ==

export get_atom_l_list,
	product_shsiteindex,
	shsiteindex_singleatom,
	SHSiteIndex, SHProduct, replace_atom, LinearCombo, inner_product

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
# Functions for generating SHSiteIndex and their combinations
# ─────────────────────────────────────────────────────────────────────────────

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
	# Iterators.product iterates with first factor varying fastest
	# But kron orders elements with last factor varying fastest
	# So we need to reverse the order to match kron
	prod_iter = Iterators.product(reverse(list_tmp)...)
	for comb::Tuple{Vararg{SHSiteIndex}} in prod_iter
		shp_tmp = SHProduct()
		# comb is in reversed order, so reverse it to restore original atom order
		# but this makes the iteration order match kron (last factor varies fastest)
		for shsi in reverse(collect(comb))
			push!(shp_tmp, shsi)
		end
		push!(combined_vec, shp_tmp)
	end
	return combined_vec
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
		new(Vector{SHSiteIndex}(vec))
	end

	function SHProduct(shsi::SHSiteIndex)
		new(Vector{SHSiteIndex}([shsi]))  # single element is already canonical
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
	for (idx, shp_i) in enumerate(lco.data)
		if sort(shp) == sort(shp_i)
			return lco.coeffs[idx]
		end
	end
	return 0.0
end

end
