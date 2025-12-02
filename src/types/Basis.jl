module Basis

"""
	LinearCombo{T,N}

Container type for linear combinations of N-body (site) angular-momentum couplings.

- `ls`          : orbital angular momenta for each site (length N)
- `Lf`          : total angular momentum of the final multiplet
- `Lseq`        : intermediate L values along a left-coupling tree (length N-1)
- `atoms`       : atom indices associated with each site (stored as a resizable vector)
- `coeff_list`  : list of scalar coefficients
- `coeff_tensor`: (N+1)-way tensor over m₁,…,m_N and M_f
"""
struct LinearCombo{T <: Number, N}
	ls::NTuple{N, Int}
	Lf::Int
	Lseq::NTuple{N-1, Int}
	atoms::Vector{Int}
	coeff_list::Vector{T}
	coeff_tensor::Array{T, N+1}

	"""
		LinearCombo(ls, Lf, Lseq, atoms, coeff_list, coeff_tensor)

User-facing constructor.

- `ls`           : `AbstractVector{<:Integer}` (length N)
- `Lf`           : `Integer`
- `Lseq`         : `AbstractVector{<:Integer}` (length N-1)
- `atoms`        : `AbstractVector{<:Integer}` (length N), stored internally as `Vector{Int}`
- `coeff_list`   : `AbstractVector{T}` with `T<:Number`
- `coeff_tensor` : `AbstractArray{T,Np}` with `Np == N+1`
	"""
	function LinearCombo(
		ls::AbstractVector{<:Integer},
		Lf::Integer,
		Lseq::AbstractVector{<:Integer},
		atoms::AbstractVector{<:Integer},
		coeff_list::AbstractVector{T},
		coeff_tensor::AbstractArray{T, Np},
	) where {T <: Number, Np}
		N = length(ls)

		length(Lseq) == N-1 ||
			throw(ArgumentError("length(Lseq) must be N-1; got $(length(Lseq)), N=$N"))

		length(atoms) == N ||
			throw(ArgumentError("length(atoms) must equal N; got $(length(atoms)), N=$N"))

		Np == N+1 ||
			throw(ArgumentError("ndims(coeff_tensor) must be N+1; got $Np, N=$N"))

		return new{T, N}(
			NTuple{N, Int}(Int.(ls)),
			Int(Lf),
			NTuple{N-1, Int}(Int.(Lseq)),
			Int.(atoms),
			collect(coeff_list),
			Array{T, N+1}(coeff_tensor),
		)
	end
end

function Base.:*(mat::AbstractMatrix{<:Number}, lc::LinearCombo{T, N}) where {T <: Number, N}
	@assert size(mat, 1) == size(mat, 2) "Matrix must be square"
	@assert (length(lc.coeff_list) == size(mat, 1)) "Dimension mismatch: length(coeff_list) = $(length(lc.coeff_list)), size(mat, 1) = $(size(mat, 1))"
	return LinearCombo(lc.ls, lc.Lf, lc.Lseq, lc.atoms, mat * lc.coeff_list, lc.coeff_tensor)
end

"""
	permute_atoms!(lc::LinearCombo, perm::AbstractVector{<:Integer})

Permute the atoms in the linear combination.
This is used at the step of symmetry operations.
"""
function permute_atoms!(lc::LinearCombo, perm::AbstractVector{<:Integer})
	lc.atoms .= perm[lc.atoms]
	return lc
end

function permute_atoms(lc::LinearCombo, perm::AbstractVector{<:Integer})
	return LinearCombo(
		lc.ls,
		lc.Lf,
		lc.Lseq,
		perm[lc.atoms],
		lc.coeff_list,
		lc.coeff_tensor,
	)
end




end # module Basis
