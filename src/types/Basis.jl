module Basis

using ..UnitaryMatrixCl
using ..AngularMomentumCoupling
using LinearAlgebra
using DataStructures

export CoupledBasis,
	CoupledBasis_with_coefficient,
	reorder_atoms,
	tesseral_coupled_bases_from_tesseral_bases

"""
	CoupledBasis{T,N}

Container type for coupled angular momentum bases for N-body (site) angular-momentum couplings.

- `ls`          : orbital angular momenta for each site (length N)
- `Lf`          : total angular momentum of the final multiplet
- `Lseq`        : intermediate L values along a left-coupling tree (length N-1)
- `atoms`       : atom indices associated with each site (stored as a resizable vector)
- `coeff_tensor`: N-way tensor over m₁,…,m_N (one tensor per fixed final M_f)
"""
struct CoupledBasis
	ls::Vector{Int}
	Lf::Int
	Lseq::Vector{Int}
	atoms::Vector{Int}
	coeff_tensor::AbstractArray

	"""
			CoupledBasis(ls, Lf, Lseq, atoms, coeff_tensor)

	User-facing constructor.

- `ls`           : `AbstractVector{<:Integer}` or `NTuple{N, <:Integer}` (length N)
- `Lf`           : `Integer`
- `Lseq`         : `AbstractVector{<:Integer}` (length N-2 for N≥2, empty for N=1)
- `atoms`        : `AbstractVector{<:Integer}` (length N), stored internally as `Vector{Int}`
- `coeff_tensor` : `AbstractArray{T}` with `ndims(coeff_tensor) == N+1` (N dimensions for sites, 1 dimension for final Mf)
	"""
	# Constructor accepting NTuple for ls
	function CoupledBasis(
		ls::AbstractVector{<:Integer},
		Lf::Integer,
		Lseq::AbstractVector{<:Integer},
		atoms::AbstractVector{<:Integer},
		coeff_tensor::AbstractArray{<:Number},
	)
		N = length(ls)
		if N == 1
			return new(
				collect(Int.(ls)),
				Int(Lf),
				Int[],
				Int[atoms[1]],
				coeff_tensor,
			)
		end

		length(Lseq) == N - 2 ||
			throw(ArgumentError("length(Lseq) must be length(ls)-2; got $(length(Lseq)), N=$N"))

		length(atoms) == N ||
			throw(ArgumentError("length(atoms) must equal length(ls); got $(length(atoms)), N=$N"))

		ndims(coeff_tensor) == N + 1 ||
			throw(
				ArgumentError(
					"ndims(coeff_tensor) must be length(ls)+1; got $(ndims(coeff_tensor)), N=$N",
				),
			)

		return new(ls, Lf, Lseq, atoms, coeff_tensor)
	end
end

function Base.show(io::IO, lc::CoupledBasis)
	print(io, "CoupledBasis(")
	print(io, "ls=$(lc.ls), ")
	print(io, "Lf=$(lc.Lf), ")
	print(io, "Lseq=$(lc.Lseq), ")
	print(io, "atoms=$(lc.atoms), ")
	print(io, "coeff_tensor=$(size(lc.coeff_tensor))")
	print(io, ")")
end

function Base.isless(
	lc1::CoupledBasis,
	lc2::CoupledBasis,
)
	# First compare by number of sites
	length(lc1.ls) != length(lc2.ls) && return length(lc1.ls) < length(lc2.ls)

	# Compare ls
	lc1.ls != lc2.ls && return lc1.ls < lc2.ls

	# Compare Lf
	lc1.Lf != lc2.Lf && return lc1.Lf < lc2.Lf

	# Compare Lseq
	lc1.Lseq != lc2.Lseq && return lc1.Lseq < lc2.Lseq

	# Compare atoms
	lc1.atoms != lc2.atoms && return lc1.atoms < lc2.atoms

	# If all fields are equal, compare coeff_tensor sizes (should be same, but for completeness)
	size(lc1.coeff_tensor) != size(lc2.coeff_tensor) &&
		return size(lc1.coeff_tensor) < size(lc2.coeff_tensor)
	return false
end


"""
	reorder_atoms(cb::CoupledBasis, new_atoms::AbstractVector{<:Integer})

Reorder the atom indices in the coupled angular momentum basis and reorder the corresponding
tensor dimensions accordingly.

This function takes a new set of atom indices `new_atoms` and:
1. Sorts `new_atoms` to maintain the sorted order requirement
2. Permutes the orbital angular momenta `ls` to match the sorted atom order
3. Permutes the first N dimensions of `coeff_tensor` to correspond to the reordered sites,
   while keeping the last dimension (Mf) unchanged

**Arguments:**
- `cb`        : The `CoupledBasis` to modify
- `new_atoms` : New atom indices (length must equal the number of sites N)

**Returns:**
A new `CoupledBasis` with updated atom indices and reordered tensor dimensions.

**Example:**
```julia
# If new_atoms = [5, 1, 2], it will be sorted to [1, 2, 5]
# The permutation p = [2, 3, 1] is applied to reorder ls and coeff_tensor dimensions
```
"""
function reorder_atoms(cb::CoupledBasis, new_atoms::AbstractVector{<:Integer})
	N = length(cb.ls)
	length(new_atoms) == N ||
		throw(ArgumentError("length(new_atoms) must be $N, got $(length(new_atoms))"))

	nd = ndims(cb.coeff_tensor)
	nd == N + 1 ||
		throw(ArgumentError("coeff_tensor must have N+1 dims, got $nd (N=$N)"))

	# Find permutation p that sorts new_atoms
	# Example: new_atoms = [5,1,2] -> p = [2,3,1], new_atoms[p] = [1,2,5]
	p = sortperm(new_atoms)
	atoms_sorted = Int.(new_atoms[p])

	# Permute ls to match the sorted atom order (ls is associated with sites)
	ls_sorted = cb.ls[p]

	# Permute the first N dimensions of coeff_tensor according to p,
	# keeping the last dimension (Mf) unchanged
	dims_perm = vcat(p, nd)  # [p..., N+1]
	coeff_perm = permutedims(cb.coeff_tensor, dims_perm)

	return CoupledBasis(ls_sorted, cb.Lf, cb.Lseq, atoms_sorted, coeff_perm)
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


# Cache for angular momentum coupling results
# Key: (ls::Vector{Int}, normalize::Symbol, isotropy::Bool)
# Value: (bases_by_L, paths_by_L)
const _angular_momentum_cache = Dict{
	Tuple{Vector{Int}, Symbol, Bool},
	Tuple{Dict{Int, Vector{Array{Float64}}}, Dict{Int, Vector{Vector{Int}}}}
}()

"""
	tesseral_coupled_bases_from_tesseral_bases(ls, atoms; normalize=:none, isotropy::Bool=false)

Construct a flat list of `CoupledBasis` objects from the output of
`AngularMomentumCoupling.build_all_real_bases`.

The input tensors are already in real (tesseral) basis for all sites and the final multiplet,
so no additional transformation is needed. Each coupled basis tensor (for a given coupling path
and final `Lf`) is wrapped into a `CoupledBasis` with a single scalar coefficient `1.0` for
each Mf tesseral index.

This function caches the angular momentum coupling results (`bases_by_L` and `paths_by_L`)
based on `ls`, `normalize`, and `isotropy` parameters, so that repeated calls with the same
`ls` values but different `atoms` can reuse the cached results for better performance.
"""
function tesseral_coupled_bases_from_tesseral_bases(
	ls::AbstractVector{<:Integer},
	atoms::AbstractVector{<:Integer};
	normalize::Symbol = :none,
	isotropy::Bool = false,
)
	# Create cache key from ls (as Vector for hashing), normalize, and isotropy
	# Note: ls order matters for angular momentum coupling, so we don't sort it
	ls_vec = collect(Int.(ls))
	cache_key = (ls_vec, normalize, isotropy)
	
	# Check cache or compute
	if haskey(_angular_momentum_cache, cache_key)
		bases_by_L, paths_by_L = _angular_momentum_cache[cache_key]
	else
		bases_by_L, paths_by_L =
			build_all_real_bases(
				collect(Int.(ls));
				normalize = normalize,
				isotropy = isotropy,
			)
		# Cache the results
		_angular_momentum_cache[cache_key] = (bases_by_L, paths_by_L)
	end

	coupled_basis_list = CoupledBasis[]

	for Lf in sort(collect(keys(bases_by_L)))
		tensors::Vector{Array{Float64}} = bases_by_L[Lf]
		Lseqs::Vector{Vector{Int}} = paths_by_L[Lf]

		for (tensor, Lseq) in zip(tensors, Lseqs)
			nd = ndims(tensor)                  # should be N+1
			nd == length(ls) + 1 ||
				throw(ArgumentError("tensor rank must be N+1; got $nd, expected $(length(ls) + 1)"))

			push!(coupled_basis_list, CoupledBasis(ls, Lf, Lseq, atoms, tensor))
		end
	end

	return coupled_basis_list
end


"""
	CoupledBasis_with_coefficient

Container type for coupled angular momentum bases with Mf-dependent coefficients.

- `ls`          : orbital angular momenta for each site (length N)
- `Lf`          : total angular momentum of the final multiplet
- `Lseq`        : intermediate L values along a left-coupling tree (length N-1)
- `atoms`       : atom indices associated with each site (stored as a resizable vector)
- `coeff_tensor`: N-way tensor over m₁,…,m_N (one tensor per fixed final M_f)
- `coefficient` : coefficients for each Mf value (length must match the last dimension of coeff_tensor)
"""
struct CoupledBasis_with_coefficient
	ls::Vector{Int}
	Lf::Int
	Lseq::Vector{Int}
	atoms::Vector{Int}
	coeff_tensor::AbstractArray
	coefficient::Vector{Float64}
	multiplicity::Int

	function CoupledBasis_with_coefficient(
		ls::AbstractVector{<:Integer},
		Lf::Integer,
		Lseq::AbstractVector{<:Integer},
		atoms::AbstractVector{<:Integer},
		coeff_tensor::AbstractArray{<:Number},
		coefficient::AbstractVector{<:Number},
		multiplicity::Int,
	)
		N = length(ls)
		nd = ndims(coeff_tensor)
		nd == N + 1 ||
			throw(
				ArgumentError(
					"ndims(coeff_tensor) must be length(ls)+1; got $nd, expected $(N+1)",
				),
			)

		Mf_size = size(coeff_tensor, nd)
		length(coefficient) == Mf_size ||
			throw(
				ArgumentError(
					"length(coefficient) must match the last dimension of coeff_tensor; got $(length(coefficient)), expected $Mf_size",
				),
			)

		return new(
			collect(Int.(ls)),
			Int(Lf),
			collect(Int.(Lseq)),
			collect(Int.(atoms)),
			coeff_tensor,
			collect(Float64.(coefficient)),
			multiplicity,
		)
	end
end

function CoupledBasis_with_coefficient(
	cb::CoupledBasis,
	coefficient::AbstractVector{<:Number},
	multiplicity::Int,
)
	return CoupledBasis_with_coefficient(
		cb.ls,
		cb.Lf,
		cb.Lseq,
		cb.atoms,
		cb.coeff_tensor,
		coefficient,
		multiplicity,
	)
end

function Base.show(io::IO, cbc::CoupledBasis_with_coefficient)
	print(io, "CoupledBasis_with_coefficient(")
	print(io, "ls=$(cbc.ls), ")
	print(io, "Lf=$(cbc.Lf), ")
	print(io, "Lseq=$(cbc.Lseq), ")
	print(io, "atoms=$(cbc.atoms), ")
	print(io, "coeff_tensor=$(size(cbc.coeff_tensor)), ")
	print(io, "coefficient=$(cbc.coefficient), ")
	print(io, "multiplicity=$(cbc.multiplicity), ")
	print(io, ")")
end

end # module Basis
