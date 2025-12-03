module Basis

using ..UnitaryMatrixCl
using ..AngularMomentumCoupling
using LinearAlgebra
using DataStructures

export LinearCombo,
	permute_atoms!,
	permute_atoms,
	*,
	linear_combos_from_complex_bases,
	tesseral_linear_combos_from_tesseral_bases,
	tesseral_linear_combos_from_complex_bases

"""
	LinearCombo{T,N}

Container type for linear combinations of N-body (site) angular-momentum couplings.

- `ls`          : orbital angular momenta for each site (length N)
- `Lf`          : total angular momentum of the final multiplet
- `Lseq`        : intermediate L values along a left-coupling tree (length N-1)
- `atoms`       : atom indices associated with each site (stored as a resizable vector)
- `coeff_list`  : list of scalar coefficients
- `coeff_tensor`: N-way tensor over m₁,…,m_N (one tensor per fixed final M_f)
"""
struct LinearCombo{T <: Number, N}
	ls::NTuple{N, Int}
	Lf::Int
	Lseq::Vector{Int}
	atoms::Vector{Int}
	coeff_list::Vector{T}
	coeff_tensor::Array{T}

	"""
			LinearCombo(ls, Lf, Lseq, atoms, coeff_list, coeff_tensor)

	User-facing constructor.

	- `ls`           : `AbstractVector{<:Integer}` or `NTuple{N, <:Integer}` (length N)
	- `Lf`           : `Integer`
	- `Lseq`         : `AbstractVector{<:Integer}` (length N-1)
	- `atoms`        : `AbstractVector{<:Integer}` (length N), stored internally as `Vector{Int}`
	- `coeff_list`   : `AbstractVector{T}` with `T<:Number`
	- `coeff_tensor` : `AbstractArray{T}` with `ndims(coeff_tensor) == N`
	"""
	# Constructor accepting NTuple for ls
	function LinearCombo(
		ls::NTuple{N, <:Integer},
		Lf::Integer,
		Lseq::AbstractVector{<:Integer},
		atoms::AbstractVector{<:Integer},
		coeff_list::AbstractVector{T},
		coeff_tensor::AbstractArray{T},
	) where {T <: Number, N}
		return LinearCombo(collect(Int.(ls)), Lf, Lseq, atoms, coeff_list, coeff_tensor)
	end

	# Constructor accepting AbstractVector for ls
	function LinearCombo(
		ls::AbstractVector{<:Integer},
		Lf::Integer,
		Lseq::AbstractVector{<:Integer},
		atoms::AbstractVector{<:Integer},
		coeff_list::AbstractVector{T},
		coeff_tensor::AbstractArray{T},
	) where {T <: Number}
		N = length(ls)

		length(Lseq) == N - 2 ||
			throw(ArgumentError("length(Lseq) must be N-2; got $(length(Lseq)), N=$N"))

		length(atoms) == N ||
			throw(ArgumentError("length(atoms) must equal N; got $(length(atoms)), N=$N"))

		ndims(coeff_tensor) == N ||
			throw(ArgumentError("ndims(coeff_tensor) must be N; got $(ndims(coeff_tensor)), N=$N"))

		return new{T, N}(
			NTuple{N, Int}(Int.(ls)),
			Int(Lf),
			collect(Int.(Lseq)),
			collect(Int.(atoms)),
			collect(coeff_list),
			Array{T}(coeff_tensor),
		)
	end
end

function Base.show(io::IO, lc::LinearCombo{T, N}) where {T, N}
	print(io, "LinearCombo{$T, $N}(")
	print(io, "ls=$(lc.ls), ")
	print(io, "Lf=$(lc.Lf), ")
	print(io, "Lseq=$(lc.Lseq), ")
	print(io, "atoms=$(lc.atoms), ")
	print(io, "coeff_list=$(lc.coeff_list), ")
	print(io, "coeff_tensor=$(size(lc.coeff_tensor))")
	print(io, ")")
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

# Helper function to find a permutation that matches (atom, l) pairs
function _find_matching_permutation(
	pairs1::Vector{Tuple{Int, Int}},
	pairs2::Vector{Tuple{Int, Int}},
)::Union{Vector{Int}, Nothing}
	N = length(pairs1)
	N != length(pairs2) && return nothing

	# Quick check: count occurrences of each pair
	count1 = counter(pairs1)
	count2 = counter(pairs2)
	count1 != count2 && return nothing

	# Backtracking to find a permutation
	perm = zeros(Int, N)
	used = falses(N)

	function backtrack(i::Int)::Bool
		if i > N
			return true
		end
		target_pair = pairs1[i]
		for j in 1:N
			if !used[j] && pairs2[j] == target_pair
				perm[i] = j
				used[j] = true
				if backtrack(i + 1)
					return true
				end
				used[j] = false
			end
		end
		return false
	end

	return backtrack(1) ? perm : nothing
end

"""
	dot(lc1::LinearCombo, lc2::LinearCombo)

Compute the inner product of two `LinearCombo` objects.

The inner product is zero if:
- `Lf` values differ
- `Lseq` values differ
- `ls` values differ (even if permuted)
- `(atom, l)` pairs cannot be matched by any permutation

If the `(atom, l)` pairs can be matched by a permutation, returns the dot product
of `coeff_list` vectors.
"""
function LinearAlgebra.dot(lc1::LinearCombo{T1, N1}, lc2::LinearCombo{T2, N2}) where {T1, T2, N1, N2}
	# Different number of sites
	N1 != N2 && return zero(promote_type(T1, T2))

	# Different Lf
	lc1.Lf != lc2.Lf && return zero(promote_type(T1, T2))

	# Different Lseq
	lc1.Lseq != lc2.Lseq && return zero(promote_type(T1, T2))

	# Check if (atom, l) pairs match (with possible permutation)
	# This automatically ensures ls values match after permutation
	pairs1 = collect(zip(lc1.atoms, collect(lc1.ls)))
	pairs2 = collect(zip(lc2.atoms, collect(lc2.ls)))
	perm = _find_matching_permutation(pairs1, pairs2)
	perm === nothing && return zero(promote_type(T1, T2))

	# All checks passed: compute dot product of coeff_list
	return dot(lc1.coeff_list, lc2.coeff_list)
end

"""
	build_linear_combos_from_complex_bases(ls, atoms; normalize=:none, isotropy::Bool=false)

Construct a flat list of `LinearCombo` objects from the output of
`AngularMomentumCoupling.build_all_complex_bases`.

Each coupled basis tensor (for a given coupling path and final `Lf`) is wrapped
into a `LinearCombo` with a single scalar coefficient `1.0`.
"""
function linear_combos_from_complex_bases(
	ls::AbstractVector{<:Integer},
	atoms::AbstractVector{<:Integer};
	normalize::Symbol = :none,
	isotropy::Bool = false,
)
	bases_by_L, paths_by_L =
		build_all_complex_bases(
			collect(Int.(ls));
			normalize = normalize,
			isotropy = isotropy,
		)

	combo_list = LinearCombo[]

	for Lf in sort(collect(keys(bases_by_L)))
		tensors::Vector{Array{Float64}} = bases_by_L[Lf]
		Lseqs::Vector{Vector{Int}} = paths_by_L[Lf]

		for (tensor, Lseq) in zip(tensors, Lseqs)
			nd = ndims(tensor)                  # should be N+1
			@assert nd == length(ls) + 1 "tensor rank must be N+1"
			site_axes = ntuple(_ -> Colon(), nd - 1)

			for Mf_idx in 1:(2*Lf+1)
				coeffs = zeros(Float64, 2*Lf+1)
				coeffs[Mf_idx] = 1.0
				coeff_tensor = Array(tensor[site_axes..., Mf_idx])
				push!(combo_list, LinearCombo(ls, Lf, Lseq, atoms, coeffs, coeff_tensor))
			end
		end
	end

	return combo_list
end


"""
	tesseral_linear_combos_from_complex_bases(ls, atoms; normalize=:none, isotropy::Bool=false)
	Construct a flat list of `LinearCombo` objects from the output of
	`AngularMomentumCoupling.build_all_complex_bases`.
	After constructing the complex bases, the tensors are converted to real (tesseral) bases.

"""
function tesseral_linear_combos_from_complex_bases(
	ls::AbstractVector{<:Integer},
	atoms::AbstractVector{<:Integer};
	normalize::Symbol = :none,
	isotropy::Bool = false,
)
	bases_by_L, paths_by_L =
		build_all_complex_bases(
			collect(Int.(ls));
			normalize = normalize,
			isotropy = isotropy,
		)

	combo_list = LinearCombo[]

	for Lf in sort(collect(keys(bases_by_L)))
		tensors::Vector{Array{Float64}} = bases_by_L[Lf]
		Lseqs::Vector{Vector{Int}} = paths_by_L[Lf]

		Cl_matrix::Matrix{ComplexF64} = UniMatCl(Lf).umat_cl

		for (tensor, Lseq) in zip(tensors, Lseqs)
			nd = ndims(tensor)                  # should be N+1
			@assert nd == length(ls) + 1 "tensor rank must be N+1"

			# Apply transformation to Mf axis (last axis) once
			tensor_Mf_tesseral = nmode_mul(tensor, Cl_matrix, nd)

			site_axes = ntuple(_ -> Colon(), nd - 1)
			for Mf_tesseral_idx in 1:(2*Lf+1)
				coeffs = zeros(Float64, 2*Lf+1)
				coeffs[Mf_tesseral_idx] = 1.0
				coeff_tensor = Array(tensor_Mf_tesseral[site_axes..., Mf_tesseral_idx])
				push!(combo_list, LinearCombo(ls, Lf, Lseq, atoms, coeffs, coeff_tensor))
			end
		end
	end

	return combo_list
end

"""
	tesseral_linear_combos_from_tesseral_bases(ls, atoms; normalize=:none, isotropy::Bool=false)

Construct a flat list of `LinearCombo` objects from the output of
`AngularMomentumCoupling.build_all_real_bases`.

The input tensors are already in real (tesseral) basis for all sites and the final multiplet,
so no additional transformation is needed. Each coupled basis tensor (for a given coupling path
and final `Lf`) is wrapped into a `LinearCombo` with a single scalar coefficient `1.0` for
each Mf tesseral index.
"""
function tesseral_linear_combos_from_tesseral_bases(
	ls::AbstractVector{<:Integer},
	atoms::AbstractVector{<:Integer};
	normalize::Symbol = :none,
	isotropy::Bool = false,
)
	bases_by_L, paths_by_L =
		build_all_real_bases(
			collect(Int.(ls));
			normalize = normalize,
			isotropy = isotropy,
		)

	combo_list = LinearCombo[]

	for Lf in sort(collect(keys(bases_by_L)))
		tensors::Vector{Array{Float64}} = bases_by_L[Lf]
		Lseqs::Vector{Vector{Int}} = paths_by_L[Lf]

		for (tensor, Lseq) in zip(tensors, Lseqs)
			nd = ndims(tensor)                  # should be N+1
			@assert nd == length(ls) + 1 "tensor rank must be N+1"

			site_axes = ntuple(_ -> Colon(), nd - 1)
			for Mf_tesseral_idx in 1:(2*Lf+1)
				coeffs = zeros(Float64, 2*Lf+1)
				coeffs[Mf_tesseral_idx] = 1.0
				coeff_tensor = Array(tensor[site_axes..., Mf_tesseral_idx])
				push!(combo_list, LinearCombo(ls, Lf, Lseq, atoms, coeffs, coeff_tensor))
			end
		end
	end

	return combo_list
end

end # module Basis
