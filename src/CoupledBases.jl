module CoupledBases

using ..AngularMomentumCoupling
using LinearAlgebra
using DataStructures

export CoupledBasis,
    CoupledBasis_with_coefficient,
    AngularMomentumCouplingResult,
    reorder_atoms,
    tesseral_coupled_bases_from_tesseral_bases,
    convert_to_coupled_basis,
    enumerate_orbit_clusters

"""
    CoupledBasis{R}

Container type for coupled angular momentum bases for N-body (site) angular-momentum couplings.

The type parameter `R` is the tensor rank `ndims(coeff_tensor) = length(ls) + 1`
(N dimensions for sites + 1 dimension for the final Mf).

- `ls`          : orbital angular momenta for each site (length N = R - 1)
- `Lf`          : total angular momentum of the final multiplet
- `Lseq`        : intermediate L values along a left-coupling tree
                  (length `max(0, N-2)`: empty for N≤2, N-2 entries for N≥3,
                  since `Lf` is stored separately)
- `atoms`       : atom indices associated with each site (length N)
- `coeff_tensor`: rank-`R` tensor of `Float64` over m₁,…,m_N and final M_f
"""
struct CoupledBasis{R}
    ls::Vector{Int}
    Lf::Int
    Lseq::Vector{Int}
    atoms::Vector{Int}
    coeff_tensor::Array{Float64, R}

    """
            CoupledBasis(ls, Lf, Lseq, atoms, coeff_tensor)

    User-facing constructor. The tensor rank `R` is inferred from
    `ndims(coeff_tensor)` and the input tensor is converted to
    `Array{Float64, R}` (a no-op when it is already that type).

- `ls`           : `AbstractVector{<:Integer}` (length N)
- `Lf`           : `Integer`
- `Lseq`         : `AbstractVector{<:Integer}` (length `max(0, N-2)`)
- `atoms`        : `AbstractVector{<:Integer}` (length N), stored internally as `Vector{Int}`
- `coeff_tensor` : `AbstractArray{<:Number}` with `ndims(coeff_tensor) == N+1`
                   (N dimensions for sites, 1 dimension for final Mf)
    """
    function CoupledBasis(
        ls::AbstractVector{<:Integer},
        Lf::Integer,
        Lseq::AbstractVector{<:Integer},
        atoms::AbstractVector{<:Integer},
        coeff_tensor::AbstractArray{<:Number},
    )
        N = length(ls)
        expected_Lseq_len = max(0, N - 2)

        length(Lseq) == expected_Lseq_len || throw(ArgumentError(
            "length(Lseq) must be max(0, N-2) = $expected_Lseq_len; " *
            "got $(length(Lseq)) (N=$N)",
        ))

        length(atoms) == N || throw(ArgumentError(
            "length(atoms) must equal length(ls) = $N; " *
            "got $(length(atoms)) (N=$N)",
        ))

        R = ndims(coeff_tensor)
        R == N + 1 || throw(ArgumentError(
            "ndims(coeff_tensor) must be length(ls)+1 = $(N + 1); " *
            "got $R (N=$N)",
        ))

        tensor_concrete = convert(Array{Float64, R}, coeff_tensor)

        return new{R}(
            collect(Int.(ls)),
            Int(Lf),
            collect(Int.(Lseq)),
            collect(Int.(atoms)),
            tensor_concrete,
        )
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

# Arguments
- `cb`        : The `CoupledBasis` to modify
- `new_atoms` : New atom indices (length must equal the number of sites N)

# Returns
A new `CoupledBasis` with updated atom indices and reordered tensor dimensions.

# Example
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


"""
    AngularMomentumCouplingResult{R}

Stores angular momentum coupling results without atom index information.
The type parameter `R` is the tensor rank `ndims(coeff_tensor)`.

# Fields
- `ls::Vector{Int}`: Angular momenta for each site (e.g., [1, 2])
- `Lseq::Vector{Int}`: Intermediate L values sequence (e.g., [3])
- `Lf::Int`: Final total angular momentum
- `coeff_tensor::Array{Float64, R}`: Coefficient tensor
"""
struct AngularMomentumCouplingResult{R}
    ls::Vector{Int}
    Lseq::Vector{Int}
    Lf::Int
    coeff_tensor::Array{Float64, R}

    function AngularMomentumCouplingResult(
        ls::AbstractVector{<:Integer},
        Lseq::AbstractVector{<:Integer},
        Lf::Integer,
        coeff_tensor::AbstractArray{<:Number},
    )
        R = ndims(coeff_tensor)
        tensor_concrete = convert(Array{Float64, R}, coeff_tensor)
        return new{R}(
            collect(Int.(ls)),
            collect(Int.(Lseq)),
            Int(Lf),
            tensor_concrete,
        )
    end
end

# Memoization of angular-momentum coupling results. Building the coupled
# tensors for a given key is the expensive part of basis construction, and the
# same key recurs across many atom tuples within one build, so a cache shared
# across those calls avoids the rebuild. The value depends only on the key (the
# computation is deterministic). A cache is created per basis construction and
# threaded explicitly through the build, rather than held as module state, so
# that concurrent constructions never share it.
#
# Key:   (ls::Vector{Int}, normalize::Symbol, isotropy::Bool)
# Value: (bases_by_L, paths_by_L)
const CouplingCache = Dict{
    Tuple{Vector{Int}, Symbol, Bool},
    Tuple{Dict{Int, Vector{Array{Float64}}}, Dict{Int, Vector{Vector{Int}}}},
}

"""
    cached_coupling_results(cache, ls, isotropy::Bool)

Look up previously-built angular-momentum coupling results for the per-site
angular momenta `ls` (with the canonical `normalize = :none`) in `cache`,
without triggering a build. Returns the `(bases_by_L, paths_by_L)` tuple if the
combination is present, or `nothing` otherwise, where:
- `bases_by_L`: coupled-tensor arrays keyed by the total angular momentum `Lf`.
- `paths_by_L`: the corresponding coupling-path (`Lseq`) sequences, keyed by `Lf`.

`cache` is the same object passed to
[`tesseral_coupled_bases_from_tesseral_bases`](@ref) during the build, so this
collects coupling results for `ls` combinations already built without depending
on the internal cache representation. Not exported; called across the module
boundary by `SALCBases`.
"""
function cached_coupling_results(
    cache::CouplingCache,
    ls::AbstractVector{<:Integer},
    isotropy::Bool,
)::Union{
    Nothing,
    Tuple{Dict{Int, Vector{Array{Float64}}}, Dict{Int, Vector{Vector{Int}}}},
}
    return get(cache, (collect(Int.(ls)), :none, isotropy), nothing)
end

"""
    tesseral_coupled_bases_from_tesseral_bases(ls, atoms; normalize=:none, isotropy=false, cache=CouplingCache())

Construct a flat list of `CoupledBasis` objects from the output of
`AngularMomentumCoupling.build_all_real_bases`.

The input tensors are already in real (tesseral) basis for all sites and the final multiplet,
so no additional transformation is needed. Each coupled basis tensor (for a given coupling path
and final `Lf`) is wrapped into a `CoupledBasis` with a single scalar coefficient `1.0` for
each Mf tesseral index.

The coupling results (`bases_by_L` and `paths_by_L`) are memoized in `cache`, keyed by
`(ls, normalize, isotropy)`, so repeated calls with the same `ls` but different `atoms` reuse
the build. Pass a shared `cache` across one basis construction to get that reuse; the default
builds a fresh, single-use cache.
"""
function tesseral_coupled_bases_from_tesseral_bases(
    ls::AbstractVector{<:Integer},
    atoms::AbstractVector{<:Integer};
    normalize::Symbol = :none,
    isotropy::Bool = false,
    cache::CouplingCache = CouplingCache(),
)
    # ls order matters for angular momentum coupling, so it is not sorted.
    ls_vec = collect(Int.(ls))
    bases_by_L, paths_by_L = get!(cache, (ls_vec, normalize, isotropy)) do
        build_all_real_bases(ls_vec; normalize = normalize, isotropy = isotropy)
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
    enumerate_orbit_clusters(atoms, map_sym, symnum_translation) -> Vector{Vector{Int}}

Enumerate the distinct symmetry-equivalent clusters generated from a seed
atom-tuple under the supercell's pure translations.

Each entry of the returned vector is an `N`-element atom index list in the
same site order as `atoms` (so it matches the orbital-angular-momentum
order on the coupled-basis tensor). Two translations that map the seed
to the same multiset of atoms are deduplicated by sorting and comparing,
keeping the first occurrence.

This is build-time (not hot-path) work: the result is cached on
`CoupledBasis_with_coefficient.clusters` so the design-matrix kernels
can iterate the orbit directly without re-running the translation /
hash / dedup machinery per matrix element.

# Arguments
- `atoms::AbstractVector{<:Integer}`: seed atom indices (length `N`)
- `map_sym::AbstractMatrix{<:Integer}`: `symmetry.map_sym`
- `symnum_translation::AbstractVector{<:Integer}`: indices of pure-translation operations

# Returns
- `Vector{Vector{Int}}`: distinct translated clusters, each in site order

# Examples
```julia
# Three atoms, three pure translations: identity (col 1) and two cyclic
# shifts (cols 2, 3). `map_sym[a, t]` is the image of atom `a` under
# translation `t`.
map_sym = [1 2 3;
           2 3 1;
           3 1 2]
symnum_translation = [1, 2, 3]
enumerate_orbit_clusters([1, 2], map_sym, symnum_translation)
# returns Vector{Int}[[1, 2], [2, 3], [3, 1]]
# — three distinct pair clusters, each with the same site order as `atoms`.
```
"""
function enumerate_orbit_clusters(
    atoms::AbstractVector{<:Integer},
    map_sym::AbstractMatrix{<:Integer},
    symnum_translation::AbstractVector{<:Integer},
)::Vector{Vector{Int}}
    N = length(atoms)
    seen_sorted = Set{Vector{Int}}()
    clusters = Vector{Vector{Int}}()
    for itrans in symnum_translation
        translated = Vector{Int}(undef, N)
        @inbounds for site_idx = 1:N
            translated[site_idx] = Int(map_sym[atoms[site_idx], itrans])
        end
        sorted_key = sort(translated)
        if sorted_key in seen_sorted
            continue
        end
        push!(seen_sorted, sorted_key)
        push!(clusters, translated)
    end
    return clusters
end


"""
    CoupledBasis_with_coefficient{R, N}

Container type for coupled angular momentum bases with Mf-dependent coefficients.

The type parameter `R` is the tensor rank `ndims(coeff_tensor) = length(ls) + 1`.
The second parameter `N = R - 1` is the rank of `folded_tensor`; it is a
separate parameter because Julia does not allow arithmetic on a TypeVar
inside a struct field declaration.

- `ls`            : orbital angular momenta for each site (length N = R - 1)
- `Lf`            : total angular momentum of the final multiplet
- `Lseq`          : intermediate L values along a left-coupling tree
                    (length `max(0, N-2)`: empty for N≤2, N-2 entries for N≥3,
                    since `Lf` is stored separately)
- `atoms`         : atom indices associated with each site (length N)
- `coeff_tensor`  : rank-`R` tensor of `Float64` over m₁,…,m_N and final M_f
- `coefficient`   : coefficients for each Mf value (length must match the last dimension of coeff_tensor)
- `multiplicity`  : symmetry-orbit multiplicity
- `folded_tensor` : rank-`R-1` tensor obtained by contracting `coeff_tensor`
                    against `coefficient` along the Mf axis,
                    `folded_tensor[m₁,…,m_N] = Σ_Mf coefficient[Mf] · coeff_tensor[m₁,…,m_N, Mf]`,
                    precomputed at construction time. Hot-path consumers
                    (the design-matrix kernels) read this directly to avoid
                    iterating the Mf axis per element.
- `clusters`      : distinct translated cluster atom-tuples in this orbit,
                    each an `N`-element index list in the same site order as
                    `atoms`. Built via [`enumerate_orbit_clusters`](@ref) at
                    SALC construction (or XML load) time and read directly
                    by the design-matrix kernels in place of the previous
                    per-element translation / dedup loop.
"""
struct CoupledBasis_with_coefficient{R, N}
    ls::Vector{Int}
    Lf::Int
    Lseq::Vector{Int}
    atoms::Vector{Int}
    coeff_tensor::Array{Float64, R}
    coefficient::Vector{Float64}
    multiplicity::Int
    folded_tensor::Array{Float64, N}
    clusters::Vector{Vector{Int}}

    function CoupledBasis_with_coefficient(
        ls::AbstractVector{<:Integer},
        Lf::Integer,
        Lseq::AbstractVector{<:Integer},
        atoms::AbstractVector{<:Integer},
        coeff_tensor::AbstractArray{<:Number},
        coefficient::AbstractVector{<:Number},
        multiplicity::Int,
        clusters::AbstractVector{<:AbstractVector{<:Integer}},
    )
        N = length(ls)
        expected_Lseq_len = max(0, N - 2)

        length(Lseq) == expected_Lseq_len || throw(ArgumentError(
            "length(Lseq) must be max(0, N-2) = $expected_Lseq_len; " *
            "got $(length(Lseq)) (N=$N)",
        ))

        length(atoms) == N || throw(ArgumentError(
            "length(atoms) must equal length(ls) = $N; " *
            "got $(length(atoms)) (N=$N)",
        ))

        R = ndims(coeff_tensor)
        R == N + 1 || throw(ArgumentError(
            "ndims(coeff_tensor) must be length(ls)+1 = $(N + 1); " *
            "got $R (N=$N)",
        ))

        Mf_size = size(coeff_tensor, R)
        length(coefficient) == Mf_size || throw(ArgumentError(
            "length(coefficient) must match the last dimension of coeff_tensor = $Mf_size; " *
            "got $(length(coefficient))",
        ))

        for (i, c) in enumerate(clusters)
            length(c) == N || throw(ArgumentError(
                "clusters[$i] must have length N = $N (same as length(ls)); " *
                "got $(length(c))",
            ))
        end

        tensor_concrete = convert(Array{Float64, R}, coeff_tensor)
        coefficient_concrete = collect(Float64.(coefficient))
        folded = _fold_mf(tensor_concrete, coefficient_concrete)
        clusters_concrete = Vector{Int}[collect(Int, c) for c in clusters]

        return new{R, R - 1}(
            collect(Int.(ls)),
            Int(Lf),
            collect(Int.(Lseq)),
            collect(Int.(atoms)),
            tensor_concrete,
            coefficient_concrete,
            multiplicity,
            folded,
            clusters_concrete,
        )
    end
end

# Contract `coeff_tensor` against `coefficient` along the trailing Mf axis,
# producing the (R-1)-rank tensor consumed by the hot-path design-matrix
# kernels. Accumulating in increasing Mf order matches the previous
# in-kernel summation order, so the change in reduction order from
# distributing the Mf weight inwards is the only source of floating-point
# drift (always within a few ULPs in practice).
function _fold_mf(
    coeff_tensor::Array{Float64, R},
    coefficient::Vector{Float64},
) where {R}
    Mf_size = size(coeff_tensor, R)
    out_size = ntuple(i -> size(coeff_tensor, i), Val(R - 1))
    folded = zeros(Float64, out_size...)
    for mf = 1:Mf_size
        c = coefficient[mf]
        slice = selectdim(coeff_tensor, R, mf)
        @inbounds for I in eachindex(folded, slice)
            folded[I] += c * slice[I]
        end
    end
    return folded
end

function CoupledBasis_with_coefficient(
    cb::CoupledBasis,
    coefficient::AbstractVector{<:Number},
    multiplicity::Int,
    clusters::AbstractVector{<:AbstractVector{<:Integer}},
)
    return CoupledBasis_with_coefficient(
        cb.ls,
        cb.Lf,
        cb.Lseq,
        cb.atoms,
        cb.coeff_tensor,
        coefficient,
        multiplicity,
        clusters,
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
    print(io, "folded_tensor=$(size(cbc.folded_tensor)), ")
    print(io, "clusters=$(length(cbc.clusters))")
    print(io, ")")
end

"""
    reorder_atoms(cbc::CoupledBasis_with_coefficient, new_atoms::AbstractVector{<:Integer})

Reorder the atom indices in the coupled angular momentum basis with coefficients and reorder 
the corresponding tensor dimensions accordingly.

This function takes a new set of atom indices `new_atoms` and:
1. Sorts `new_atoms` to maintain the sorted order requirement
2. Permutes the orbital angular momenta `ls` to match the sorted atom order
3. Permutes the first N dimensions of `coeff_tensor` to correspond to the reordered sites,
   while keeping the last dimension (Mf) unchanged
4. Keeps the `coefficient` vector unchanged (since it corresponds to the Mf dimension which is not permuted)

# Arguments
- `cbc`       : The `CoupledBasis_with_coefficient` to modify
- `new_atoms` : New atom indices (length must equal the number of sites N)

# Returns
A new `CoupledBasis_with_coefficient` with updated atom indices and reordered tensor dimensions.
The `coefficient` vector is preserved as-is since it corresponds to the Mf dimension.

# Example
```julia
# If new_atoms = [5, 1, 2], it will be sorted to [1, 2, 5]
# The permutation p = [2, 3, 1] is applied to reorder ls and coeff_tensor dimensions
# The coefficient vector remains unchanged
```
"""
function reorder_atoms(cbc::CoupledBasis_with_coefficient, new_atoms::AbstractVector{<:Integer})
    N = length(cbc.ls)
    length(new_atoms) == N ||
        throw(ArgumentError("length(new_atoms) must be $N, got $(length(new_atoms))"))

    nd = ndims(cbc.coeff_tensor)
    nd == N + 1 ||
        throw(ArgumentError("coeff_tensor must have N+1 dims, got $nd (N=$N)"))

    # Find permutation p that sorts new_atoms
    # Example: new_atoms = [5,1,2] -> p = [2,3,1], new_atoms[p] = [1,2,5]
    p = sortperm(new_atoms)
    atoms_sorted = Int.(new_atoms[p])

    # Permute ls to match the sorted atom order (ls is associated with sites)
    ls_sorted = cbc.ls[p]

    # Permute the first N dimensions of coeff_tensor according to p,
    # keeping the last dimension (Mf) unchanged
    dims_perm = vcat(p, nd)  # [p..., N+1]
    coeff_perm = permutedims(cbc.coeff_tensor, dims_perm)

    # coefficient vector corresponds to the Mf dimension (last dimension),
    # which is not permuted, so we keep it unchanged.
    #
    # Each cluster's atom indices are listed in the same site order as
    # `cbc.atoms`, so permuting `atoms` by `p` requires permuting each
    # cluster by the same `p`.
    clusters_perm = Vector{Int}[c[p] for c in cbc.clusters]
    return CoupledBasis_with_coefficient(
        ls_sorted,
        cbc.Lf,
        cbc.Lseq,
        atoms_sorted,
        coeff_perm,
        cbc.coefficient,  # unchanged - corresponds to Mf dimension
        cbc.multiplicity,
        clusters_perm,
    )
end

"""
    convert_to_coupled_basis(cbc::CoupledBasis_with_coefficient)::CoupledBasis

Convert a `CoupledBasis_with_coefficient` to a `CoupledBasis` by removing the coefficients.
"""
function convert_to_coupled_basis(cbc::CoupledBasis_with_coefficient)::CoupledBasis
    return CoupledBasis(
        cbc.ls,
        cbc.Lf,
        cbc.Lseq,
        cbc.atoms,
        cbc.coeff_tensor,
    )
end

end # module CoupledBases
