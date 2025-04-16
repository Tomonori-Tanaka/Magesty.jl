"""
	Module UnitaryMatrixCl

This module provides a data structure for unitary matrices Cˡ that transform between
real and complex spherical harmonics.

# Types
- `UniMatCl`: A structure representing the unitary matrix Cˡ

# References
M.A. Blanco et al., Journal of Molecular Structure (Theochem) 419 19-27 (1997).
"""
module UnitaryMatrixCl

import Base:
	getindex, setindex!, *, transpose, inv, show, adjoint, conj, size, eltype, ==, isapprox

using LinearAlgebra

export UniMatCl, getindex_m

"""
	UniMatCl{T <: Complex} <: AbstractMatrix{T}

A structure storing unitary matrix Cˡ, where Cˡ is defined as
	Sₗ = CˡYₗ
where Sₗ and Yₗ are vectors consisting of real and complex spherical harmonics, respectively.
The dimension of Cˡ is (2l+1) × (2l+1).

# Fields
- `umat_cl::Matrix{T}`: The unitary matrix Cˡ
- `l::Int`: The angular momentum quantum number

# Constructors
- `UniMatCl(l::Integer)`: Create a unitary matrix Cˡ for a given l
- `UniMatCl(mat::AbstractMatrix{<:Complex})`: Create a unitary matrix Cˡ from an existing matrix
"""
struct UniMatCl{T <: Complex} <: AbstractMatrix{T}
	umat_cl::Matrix{T}
	l::Int

	function UniMatCl{T}(l::Integer) where T <: Complex
		l ≥ 0 || throw(ArgumentError("Angular momentum l must be non-negative"))
		
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
		new{T}(umat_cl, Int(l))
	end

	function UniMatCl{T}(mat::AbstractMatrix{<:Complex}) where T <: Complex
		# Check if the matrix is square
		n, m = size(mat)
		n == m || throw(ArgumentError("Matrix must be square"))
		
		# Calculate l from matrix dimension
		l = (n - 1) ÷ 2
		(2l + 1) == n || throw(ArgumentError(
			"Matrix dimension $(n) is not valid for a unitary matrix Cˡ. " *
			"Expected dimension of the form 2l+1 where l is a non-negative integer"
		))
		
		# Create new instance
		new{T}(convert(Matrix{T}, mat), l)
	end
end

# Outer constructors
UniMatCl(l::Integer) = UniMatCl{Complex}(Int(l))
UniMatCl(mat::AbstractMatrix{<:Complex}) = UniMatCl{Complex}(mat)

"""
	getindex_m(cl::UniMatCl, m1::Integer, m2::Integer) -> Complex

Get the matrix element Cˡ(m1, m2) using quantum numbers m1 and m2.

# Arguments
- `cl::UniMatCl`: The unitary matrix
- `m1::Integer`: The first quantum number
- `m2::Integer`: The second quantum number

# Returns
- `Complex`: The matrix element Cˡ(m1, m2)

# Throws
- `BoundsError` if m1 or m2 is outside the valid range [-l, l]
"""
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

# Base interface implementations
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

# Matrix operations
function *(cl1::UniMatCl, cl2::UniMatCl)
	if cl1.l != cl2.l
		throw(
			ArgumentError(
				"Cannot multiply UniMatCl instances with different l values: cl1.l=$(cl1.l), cl2.l=$(cl2.l)",
			),
		)
	end
	return UniMatCl(cl1.umat_cl * cl2.umat_cl)
end

*(cl::UniMatCl, m::AbstractMatrix{<:Number}) = cl.umat_cl * m
*(m::AbstractMatrix{<:Number}, cl::UniMatCl) = m * cl.umat_cl

function transpose(cl::UniMatCl)
	transposed_mat = Matrix(transpose(cl.umat_cl))
	return UniMatCl(transposed_mat)
end

# Implement inverse as the adjoint for unitary matrices
function inv(cl::UniMatCl)
	# For unitary matrices, the inverse is equal to the adjoint
	return UniMatCl(Matrix(adjoint(cl.umat_cl)))
end

function show(io::IO, cl::UniMatCl)
	println(io, "UniMatCl(l=$(cl.l)):")
	show(io, cl.umat_cl)
end

# Additional operations
adjoint(cl::UniMatCl) = UniMatCl(Matrix(adjoint(cl.umat_cl)))
conj(cl::UniMatCl) = UniMatCl(conj(cl.umat_cl))
size(cl::UniMatCl) = size(cl.umat_cl)
eltype(::Type{UniMatCl{T}}) where T <: Complex = T
==(cl1::UniMatCl, cl2::UniMatCl) = cl1.umat_cl == cl2.umat_cl
isapprox(cl1::UniMatCl, cl2::UniMatCl) = isapprox(cl1.umat_cl, cl2.umat_cl)

end