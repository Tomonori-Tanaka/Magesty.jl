module UnitaryMatrixCl

import Base:
	getindex, setindex!, *, transpose, inv, show, adjoint, conj, size, eltype, ==, isapprox

using LinearAlgebra

export UniMatCl, getindex_m

"""
	UniMatCl(l::Int)

A structure storing unitary matrix Cˡ, where Cˡ is defined as
	Sₗ = CˡYₗ
where Sₗ and Yₗ are vectors consisting of real and complex spherical harmonics, respectively.
The dimension of Cˡ is (2l+1) * (2l+1).

# References
M.A. Blanco et al., Journal of Molecular Structure (Theochem) 419 19-27 (1997).
"""
struct UniMatCl{T <: Complex} <: AbstractMatrix{T}
	umat_cl::Matrix{T}
	l::Int
end

UniMatCl(l::Integer) = UniMatCl{Complex}(Int(l))

function UniMatCl{T}(l::Integer) where T <: Complex
	# Initialize a (2l + 1) × (2l + 1) matrix of zeros with complex entries
	umat_cl = zeros(T, 2l + 1, 2l + 1)
	# Set the element corresponding to m = m' = 0 to 1
	umat_cl[l+1, l+1] = 1

	# Populate the matrix for each j from 1 to l
	for j in 1:l
		i_rev_pos = j + l + 1      # m = +i
		j_rev_pos = j + l + 1      # m' = +j
		i_rev_neg = -j + l + 1     # m = -i
		j_rev_neg = -j + l + 1     # m' = -j

		# Assign values to the matrix elements
		# The specific values are based on the properties of unitary matrices and spherical harmonics
		# `im` represents the imaginary unit in Julia
		umat_cl[i_rev_pos, j_rev_pos] = (-1)^j / √2          # Cˡ(m=i, m'=j)
		umat_cl[i_rev_pos, j_rev_neg] = 1 / √2               # Cˡ(m=i, m'=-j)
		umat_cl[i_rev_neg, j_rev_pos] = -1im * (-1)^j / √2    # Cˡ(m=-i, m'=j)
		umat_cl[i_rev_neg, j_rev_neg] = 1im / √2              # Cˡ(m=-i, m'=-j)

	end
	return UniMatCl{T}(umat_cl, Int(l))
end

function getindex_m(cl::UniMatCl, m1::Integer, m2::Integer)
	# Convert m and m' to matrix indices
	i = m1 + cl.l + 1
	j = m2 + cl.l + 1
	matrix_size = 2 * cl.l + 1

	# Check if the calculated indices are within the valid range
	if i < 1 || i > matrix_size || j < 1 || j > matrix_size
		throw(BoundsError(cl.umat_cl, (i, j)))
	end

	return cl.umat_cl[i, j]
end

getindex(cl::UniMatCl, i::Integer, j::Integer) = cl.umat_cl[i, j]

function getindex(cl::UniMatCl, i::Integer)
	throw(
		ArgumentError(
			"Single indexing is not supported for UniMatCl. Use two indices (i, j) instead.",
		),
	)
end

# Override setindex! to prohibit element modification
function setindex!(cl::UniMatCl, value, args...)
	throw(
		ArgumentError(
			"UniMatCl is immutable. Modification is not allowed.",
		),
	)
end

*(cl::UniMatCl, m::AbstractArray{<:Number}) = cl.umat_cl * m
*(m::AbstractArray{<:Number}, cl::UniMatCl) = m * cl.umat_cl
function *(cl1::UniMatCl, cl2::UniMatCl)
	if cl1.l != cl2.l
		throw(
			ArgumentError(
				"Cannot multiply UniMatCl instances with different l values: cl1.l=$(cl1.l), cl2.l=$(cl2.l)",
			),
		)
	end
	return UniMatCl(cl1.umat_cl * cl2.umat_cl, cl1.l)
end

function transpose(cl::UniMatCl)
	transposed_mat = Matrix(transpose(cl.umat_cl))
	return UniMatCl(transposed_mat, cl.l)
end

# Implement inverse as the adjoint for unitary matrices
function inv(cl::UniMatCl)
	# For unitary matrices, the inverse is equal to the adjoint
	return UniMatCl(Matrix(adjoint(cl.umat_cl)), cl.l)
end

function show(io::IO, cl::UniMatCl)
	println(io, "UniMatCl(l=$(cl.l)):")
	show(io, cl.umat_cl)
end

adjoint(cl::UniMatCl) = UniMatCl(Matrix(adjoint(cl.umat_cl)), cl.l)
conj(cl::UniMatCl) = UniMatCl(conj(cl.umat_cl), cl.l)
size(cl::UniMatCl) = size(cl.umat_cl)
eltype(::Type{UniMatCl{T}}) where T <: Complex = T
==(cl1::UniMatCl, cl2::UniMatCl) = cl1.umat_cl == cl2.umat_cl
isapprox(cl1::UniMatCl, cl2::UniMatCl) = isapprox(cl1.umat_cl, cl2.umat_cl)

end