"""
Symmetry-Adapted Linear Combination (SALC) basis construction.

Builds, from a crystal structure, its space-group symmetry, and a cluster
expansion, the symmetry-adapted basis used to assemble the spin-cluster
expansion design matrix. The entry point is the `SALCBasis` constructor; the
heavy step is the per-key-group projection-matrix diagonalization.

The order of the key groups in `SALCBasis.salc_list` is load-bearing: it fixes
the design-matrix column order and therefore the physical meaning of the fitted
coefficients, and it must match the order used by the XML I/O. Do not reorder it
without updating those sites together.
"""
module SALCBases

using Base.Threads
using Combinat
using DataStructures
using LinearAlgebra
using OffsetArrays
using Printf

using ..SortedCounters: SortedCounter
using ..InputSpecs: InteractionSpec, SymmetryOptions
using ..Structures
using ..Symmetries
using ..Clusters
using ..RotationMatrix
using ..CoupledBases
using ..ProgressReporting: with_progress, tick!

export SALCBasis

# Default for the `check_irrep_unitary` keyword of
# `_projection_matrix_coupled_basis`. Read once at load time so the value is
# fixed for the process instead of re-read from the environment on every call.
# Set `MAGESTY_CHECK_IRREP_UNITARY=1` in the environment before loading the
# package to enable the per-operation unitarity check by default.
const CHECK_IRREP_UNITARY_DEFAULT = get(ENV, "MAGESTY_CHECK_IRREP_UNITARY", "0") == "1"


"""
    struct SALCBasis

Represents a set of basis functions for atomic interactions in a crystal structure.
This structure is used to store and manage basis functions that are adapted to the symmetry of the crystal.

# Fields
- `coupled_basislist::SortedCounter{CoupledBases.CoupledBasis}`: List of coupled angular momentum basis functions with their multiplicities
- `salc_list::Vector{Vector{CoupledBases.CoupledBasis_with_coefficient}}`: List of symmetry-adapted linear combinations (SALCs), where each element is a vector of coupled basis functions belonging to the same key group

# Constructors
    SALCBasis(structure, symmetry, cluster, body1_lmax, bodyn_lsum, nbody; isotropy=false, verbosity=true)
    SALCBasis(structure, symmetry, cluster, config; verbosity=true)

Constructs a new `SALCBasis` instance for atomic interactions in a crystal structure.

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information for atomic interactions
- `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions for each atomic species
- `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
- `nbody::Integer`: Maximum number of bodies in interactions
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`
- `verbosity::Bool`: Whether to print progress information (default: `true`)

# Returns
- `SALCBasis`: A new basis set instance containing coupled basis functions and symmetry-adapted linear combinations

# Examples
```julia
# Create a basis set using explicit parameters
body1_lmax = [2, 3]  # lmax for each atomic species
bodyn_lsum = OffsetArray([0, 0, 4, 6], 0:3)  # lsum for each body
basis = SALCBasis(structure, symmetry, cluster, body1_lmax, bodyn_lsum, 3)

# Create a basis set using configuration
basis = SALCBasis(structure, symmetry, cluster, config)
```

# Note
The constructor performs the following steps:
1. Constructs and classifies coupled basis functions using orbit information
2. Constructs projection matrices for each symmetry label
3. Generates symmetry-adapted linear combinations (SALCs) of `CoupledBasis_with_coefficient` objects
"""
struct SALCBasis
    coupled_basislist::SortedCounter{CoupledBases.CoupledBasis}
    salc_list::Vector{Vector{CoupledBases.CoupledBasis_with_coefficient}}
    angular_momentum_couplings::Vector{CoupledBases.AngularMomentumCouplingResult}
end

"""
    salc_fingerprint(b::SALCBasis) -> UInt64

Compute a stable structural fingerprint of `b.salc_list`. Two `SALCBasis`
values that share the same key-group ordering and the same per-SALC
structural identifiers produce the same fingerprint, even after a round
trip through `Magesty.save` / `Magesty.load`. Used by the basis-identity
check between an `SCEModel` / `SCEFit` and an `SCEDataset`.

# Recipe

For each `coupled::CoupledBasis_with_coefficient` visited in the order
`(outer key-group index, inner SALC index)`, the tuple

    (coupled.ls, coupled.Lf, coupled.Lseq, coupled.atoms,
     coupled.multiplicity)

is mixed into a running hash seeded with
`hash(:Magesty_SALC_fingerprint_v1)`. The versioned seed lets a future
recipe revision (`v2`, ...) coexist with values produced by the current
implementation.

# Fields included

All structural identifiers are integer-valued and bit-identical across
`save` / `load` (XML parses them via `parse.(Int, ...)`):

- `ls::Vector{Int}` — per-site orbital angular momenta
- `Lf::Int` — final coupled L
- `Lseq::Vector{Int}` — intermediate L sequence
- `atoms::Vector{Int}` — atom indices
- `multiplicity::Int` — SALC multiplicity

# Fields deliberately excluded

The floating-point payload is **not** mixed in:

- `coefficient::Vector{Float64}` — serialized to XML as decimal text;
  reloaded values can differ by a few ULPs.
- `coeff_tensor::Array{Float64,R}` — derived from upstream tensor
  algebra; same round-off concern.

Including either would make the fingerprint flip on every save/load
cycle and defeat the purpose of disk-aware basis matching.

# Returns
- `UInt64` — usable as an O(1) basis-identity check.
"""
function salc_fingerprint(b::SALCBasis)::UInt64
    h = hash(:Magesty_SALC_fingerprint_v1)
    for group in b.salc_list, coupled in group
        h = hash(
            (coupled.ls, coupled.Lf, coupled.Lseq, coupled.atoms,
             coupled.multiplicity),
            h,
        )
    end
    return h
end

"""
    _canonicalize_sign!(v::AbstractVector{<:Real}, tol::Real = 1e-8) -> v

Apply a deterministic sign convention to `v` in place. Scans from index
1 and flips the entire vector iff the first entry with `abs(x) > tol`
is negative. Entries with `abs(x) ≤ tol` — including `+0.0` and `-0.0`
— are skipped, so the choice is robust to round-off noise and to the
IEEE-754 negative-zero / positive-zero distinction.

Pivoting on the first significant entry (rather than the sign of
`sum(v)`) is what makes the convention deterministic: a `sum`-based
rule cannot disambiguate vectors with `sum(v) ≈ 0` — e.g. `(0, a, -a)`
and `(0, -a, a)` — and would keep them as distinct, breaking
cross-platform reproducibility.
"""
function _canonicalize_sign!(
    v::AbstractVector{<:Real},
    tol::Real = 1e-8,
)::AbstractVector{<:Real}
    @inbounds for x in v
        if abs(x) > tol
            x < 0 && (v .= -v)
            return v
        end
    end
    return v
end

"""
    _canonicalize_eigenspace(V::AbstractMatrix{<:Real}; tol::Real = 1e-8) -> Matrix{Float64}

Given `V` whose `d` columns are an orthonormal basis of a subspace
`S ⊂ ℝⁿ`, return a deterministic orthonormal basis `W` of the same
subspace. The construction projects standard axis vectors `eⱼ` for
`j = 1, 2, ..., n` (fixed order) onto `S` via the basis-invariant
projector `P := V Vᵀ`, then keeps each successive projection that is
linearly independent of the previously accepted ones (modified
Gram-Schmidt with axis-order pivoting).

Any orthonormal `V` spanning the same `S` produces the same `W` (up
to round-off below `tol`), so the result is LAPACK-implementation
independent and lifts the gauge ambiguity within degenerate
eigenspaces of the SALC projection matrix.
"""
function _canonicalize_eigenspace(
    V::AbstractMatrix{<:Real};
    tol::Real = 1e-8,
)::Matrix{Float64}
    n, d = size(V)
    P = V * V'                                # projector onto span(V); basis-invariant
    W = Matrix{Float64}(undef, n, d)
    k = 0
    for j = 1:n
        u = P[:, j]
        for i = 1:k
            u .-= dot(view(W, :, i), u) .* view(W, :, i)
        end
        nu = norm(u)
        if nu > tol
            k += 1
            @views W[:, k] .= u ./ nu
        end
        k == d && break
    end
    k == d || error(
        "Canonical eigenspace extraction failed: kept $k of $d basis vectors " *
        "(tol = $tol).",
    )
    return W
end

function _compute_salc_groups(
    coupled_basislist::SortedCounter{CoupledBases.CoupledBasis},
    symmetry::Symmetry,
)::Vector{Vector{CoupledBases.CoupledBasis_with_coefficient}}
    projection_mat = _projection_matrix_coupled_basis(coupled_basislist, symmetry)
    h_projection = Hermitian(projection_mat)
    eigenvals, eigenvecs = eigen!(h_projection)
    eigenvals = real.(round.(eigenvals, digits = 6))
    eigenvecs = round.(eigenvecs .* (abs.(eigenvecs) .≥ 1e-8), digits = 10)
    if !_is_proper_eigenvals(eigenvals)
        # Projection eigenvalues must be exactly 0 or 1 (the projector is
        # idempotent). A drift here means the coupled-basis tensor or the
        # symmetry input is inconsistent — downstream SALCs would be wrong,
        # so stop instead of silently continuing.
        error("SALC projection eigenvalues must be 0 or 1, got: $eigenvals")
    end
    Lf = coupled_basislist[1].Lf
    submatrix_dim = 2 * Lf + 1

    # LAPACK eigensolvers do not promise a platform-stable basis within a
    # degenerate eigenspace, so replace the raw eigenvectors with a
    # canonical, axis-ordered orthonormal basis of the same subspace. This
    # makes SALC coefficients reproducible across BLAS/LAPACK
    # implementations.
    eigval1_indices = findall(x -> isapprox(x, 1.0, atol = 1e-8), eigenvals)
    isempty(eigval1_indices) &&
        return Vector{Vector{CoupledBases.CoupledBasis_with_coefficient}}()
    canonical_basis = _canonicalize_eigenspace(eigenvecs[:, eigval1_indices])

    # Precompute the symmetry-orbit cluster list once per coupled basis.
    # The same `cb.atoms` reappears across every Mf channel of the projection
    # (the outer `col` loop), so caching avoids redoing the translation /
    # dedup work per channel. Multiple `CoupledBasis_with_coefficient`
    # instances built from the same `cb` share the same cluster vector by
    # reference, which is harmless since `clusters` is read-only after
    # construction.
    clusters_per_basis = Vector{Vector{Vector{Int}}}(undef, length(coupled_basislist))
    for (idx_basis, cb) in enumerate(coupled_basislist)
        clusters_per_basis[idx_basis] = CoupledBases.enumerate_orbit_clusters(
            cb.atoms, symmetry.map_sym, symmetry.symnum_translation,
        )
    end

    key_salc_groups = Vector{Vector{CoupledBases.CoupledBasis_with_coefficient}}()
    for col in axes(canonical_basis, 2)
        eigenvec = canonical_basis[:, col]
        eigenvec = round.(eigenvec .* (abs.(eigenvec) .≥ 1e-8), digits = 10)
        # IEEE 754: -0.0 + +0.0 == +0.0, so this normalizes any negative
        # zeros that survived rounding into canonical positive zeros.
        eigenvec .+= 0.0
        eigenvec = eigenvec / norm(eigenvec)
        _canonicalize_sign!(eigenvec)
        salc_group = Vector{CoupledBases.CoupledBasis_with_coefficient}()
        for (idx_basis, cb) in enumerate(coupled_basislist)
            coeff_start = (idx_basis - 1) * submatrix_dim + 1
            coeff_end = idx_basis * submatrix_dim
            coefficient = eigenvec[coeff_start:coeff_end]
            if isapprox(norm(coefficient), 0.0, atol = 1e-10)
                continue
            end
            multiplicity = coupled_basislist.counts[cb]
            cbc = CoupledBases.CoupledBasis_with_coefficient(
                cb, coefficient, multiplicity, clusters_per_basis[idx_basis],
            )
            push!(salc_group, cbc)
        end
        if !isempty(salc_group)
            push!(key_salc_groups, salc_group)
        end
    end
    return key_salc_groups
end

function SALCBasis(
    structure::Structure,
    symmetry::Symmetry,
    cluster::Cluster,
    body1_lmax::Vector{Int},
    bodyn_lsum::OffsetArray{Int, 1},
    nbody::Integer,
    ;
    isotropy::Bool = false,
    verbosity::Bool = true,
)
    # Start timing
    start_time = time_ns()

    if verbosity
        println(
            """

            BASIS SET
            =========
            """,
        )
    end
    if verbosity
        print("Constructing and classifying coupled basis list...")
    end
    # One cache per construction, threaded through the build below and read
    # back when collecting `angular_momentum_couplings`. Keeping it local (not
    # module state) means concurrent constructions never share it.
    coupling_cache = CoupledBases.CouplingCache()
    classified_coupled_basisdict::Dict{Int, SortedCounter{CoupledBases.CoupledBasis}} =
        _construct_and_classify_coupled_basislist(
            structure,
            symmetry,
            cluster,
            body1_lmax,
            bodyn_lsum,
            nbody,
            isotropy = isotropy,
            cache = coupling_cache,
        )
    classified_coupled_basisdict = _filter_basisdict(classified_coupled_basisdict, symmetry)


    if verbosity
        println(" Done.")
    end

    keys_list = sort(collect(keys(classified_coupled_basisdict)))
    num_keys = length(keys_list)
    coupled_basislists::Vector{SortedCounter{CoupledBases.CoupledBasis}} =
        [classified_coupled_basisdict[k] for k in keys_list]
    # Pre-allocate array to store results for each key (preserving order)
    # Each key can have multiple SALC groups (one per eigenvector)
    _SALCGroup = Vector{CoupledBases.CoupledBasis_with_coefficient}
    salc_list_per_key = Vector{Vector{_SALCGroup}}(undef, num_keys)

    if verbosity
        println(@sprintf(
            "Threading: %d key groups across %d thread(s).",
            num_keys, nthreads(),
        ))
    end

    with_progress(num_keys, "Constructing projection matrix"; verbosity = verbosity) do prog
        @threads for idx = 1:num_keys
            salc_list_per_key[idx] = _compute_salc_groups(coupled_basislists[idx], symmetry)
            tick!(prog)
        end
    end

    # Collect all SALC groups in order
    salc_list = Vector{Vector{CoupledBases.CoupledBasis_with_coefficient}}()
    for key_salc_groups in salc_list_per_key
        for salc_group in key_salc_groups
            push!(salc_list, salc_group)
        end
    end
    if verbosity
        elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds
        println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
        println("-------------------------------------------------------------------")
    end




    # Reconstruct basislist from classified dictionary for SALCBasis storage
    result_basislist = SortedCounter{CoupledBases.CoupledBasis}()
    for (key, classified_basislist) in classified_coupled_basisdict
        for cb in classified_basislist
            count = classified_basislist.counts[cb]
            push!(result_basislist, cb, count)
        end
    end

    # Collect angular momentum coupling results from cache
    angular_momentum_couplings = Vector{CoupledBases.AngularMomentumCouplingResult}()
    # Collect unique ls combinations from coupled_basislist (preserve original order)
    ls_combinations_set = Set{Vector{Int}}()
    for cb in result_basislist
        push!(ls_combinations_set, cb.ls)
    end

    # Extract angular momentum coupling results from cache.
    # The cache key encodes the isotropy flag, so look it up with the same
    # isotropy value that drove the basis construction above. When
    # `isotropy=true`, the cached `bases_by_L` only contains `Lf == 0`, so
    # the resulting `angular_momentum_couplings` is naturally restricted to
    # the isotropic (scalar) sector.
    for ls_vec in ls_combinations_set
        results = CoupledBases.cached_coupling_results(coupling_cache, ls_vec, isotropy)
        results === nothing && continue
        bases_by_L, paths_by_L = results
        for Lf in sort(collect(keys(bases_by_L)))
            tensors = bases_by_L[Lf]
            Lseqs = paths_by_L[Lf]
            for (tensor, Lseq) in zip(tensors, Lseqs)
                push!(
                    angular_momentum_couplings,
                    CoupledBases.AngularMomentumCouplingResult(
                        ls_vec,
                        Lseq,
                        Lf,
                        copy(tensor),
                    ),
                )
            end
        end
    end

    return SALCBasis(result_basislist, salc_list, angular_momentum_couplings)
end

function SALCBasis(
    structure::Structure,
    symmetry::Symmetry,
    cluster::Cluster,
    interaction::InteractionSpec,
    options::SymmetryOptions,
    ;
    verbosity::Bool = true,
)
    return SALCBasis(
        structure,
        symmetry,
        cluster,
        interaction.body1_lmax,
        interaction.bodyn_lsum,
        interaction.nbody,
        isotropy = options.isotropy,
        verbosity = verbosity,
    )
end

"""
    _construct_and_classify_coupled_basislist(
        structure::Structure,
        symmetry::Symmetry,
        cluster::Cluster,
        body1_lmax::Vector{Int},
        bodyn_lsum::OffsetArray{Int, 1},
        nbody::Integer;
        isotropy::Bool = false,
    ) -> Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}

Construct coupled basis functions and classify them simultaneously using orbit information.

Basis functions are generated orbit-by-orbit and classified on-the-fly, which is
memory-efficient and keeps the resulting basis organized by symmetry label.

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information containing orbit classification
- `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions
- `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
- `nbody::Integer`: Maximum number of bodies in interactions
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`

# Returns
- `Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}`: Dictionary keyed by classification labels
  (based on `(nbody, Lf, sum(ls), Tuple(sort(ls)...))`), containing classified basis functions
"""
function _construct_and_classify_coupled_basislist(
    structure::Structure,
    symmetry::Symmetry,
    cluster::Cluster,
    body1_lmax::Vector{Int},
    bodyn_lsum::OffsetArray{Int, 1},
    nbody::Integer;
    isotropy::Bool = false,
    cache::CoupledBases.CouplingCache = CoupledBases.CouplingCache(),
)::Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}
    # Result dictionary for classified basis functions
    classified_dict = OrderedDict{Int, SortedCounter{CoupledBases.CoupledBasis}}()
    label_map = Dict{Any, Int}()
    next_label = 0

    irreducible_cluster_dict::Dict{Int, SortedCounter{Vector{Int}}} =
        cluster.irreducible_cluster_dict
    cluster_orbits_dict::Dict{Int, Dict{Int, Vector{Vector{Int}}}} =
        cluster.cluster_orbits_dict

    # Handle 1-body case
    for iat in symmetry.atoms_in_prim
        lmax = body1_lmax[structure.supercell.kd_int_list[iat]]
        for l = 2:lmax # skip l = 1 because it is prohibited by the time-reversal symmetry
            if l % 2 == 1 # skip odd l cases due to the time-reversal symmetry
                continue
            end
            # For 1-body case, create CoupledBasis with single atom
            cb_list = tesseral_coupled_bases_from_tesseral_bases(
                [l],
                [iat];
                isotropy = isotropy,
                cache = cache,
            )
            for cb::CoupledBases.CoupledBasis in cb_list
                # Classify on-the-fly: key is (nbody, Lf, sum(ls), Tuple(sort(ls)...))
                nbody_val = length(cb.ls)
                ls_sorted = Tuple(sort(cb.ls))
                key = (nbody_val, cb.Lf, sum(cb.ls), ls_sorted)

                label = get(label_map, key, 0)
                if label == 0
                    next_label += 1
                    label = next_label
                    label_map[key] = label
                    classified_dict[label] = SortedCounter{CoupledBases.CoupledBasis}()
                end
                push!(classified_dict[label], cb, 1)
            end
        end
    end

    # Process multi-body cases using cluster_orbits_dict for efficient processing
    # Group by orbit first, then classify within each orbit
    for body = 2:nbody
        # Use cluster_orbits_dict to process clusters grouped by symmetry orbits
        if haskey(cluster_orbits_dict, body)
            for (orbit_index, orbit_clusters) in cluster_orbits_dict[body]
                # Process all clusters in this orbit
                # Clusters in the same orbit have the same projection matrix structure
                # Collect all basis functions from this orbit first
                orbit_basis_list = Vector{CoupledBases.CoupledBasis}()
                # Use IdDict instead of Dict since CoupledBasis doesn't implement hash
                orbit_basis_counts = IdDict{CoupledBases.CoupledBasis, Int}()

                for atom_list::Vector{Int} in orbit_clusters
                    # Get multiplicity from irreducible_cluster_dict
                    count = irreducible_cluster_dict[body].counts[atom_list]
                    sorted_atom_list = sort(atom_list)
                    cb_list::Vector{CoupledBases.CoupledBasis} = _listup_coupled_basislist(
                        sorted_atom_list,
                        bodyn_lsum[body];
                        isotropy = isotropy,
                        cache = cache,
                    )

                    # Collect basis functions from this cluster
                    for cb::CoupledBases.CoupledBasis in cb_list
                        # Check for translationally equivalent within orbit
                        found_equivalent = false
                        for existing_cb in orbit_basis_list
                            if _is_translationally_equivalent_coupled_basis(
                                cb,
                                existing_cb,
                                symmetry,
                            )
                                found_equivalent = true
                                # Safe access: get existing count or 0, then add
                                orbit_basis_counts[existing_cb] =
                                    get(orbit_basis_counts, existing_cb, 0) + count
                                break
                            end
                        end
                        if !found_equivalent
                            push!(orbit_basis_list, cb)
                            orbit_basis_counts[cb] = count
                        end
                    end
                end

                # Classify basis functions from this orbit
                # Use orbit index in classification key to group by orbit
                for cb::CoupledBases.CoupledBasis in orbit_basis_list
                    ls_sorted = Tuple(sort(cb.ls))
                    # Include orbit_index in classification key to utilize orbit information
                    key = (body, orbit_index, cb.Lf, sum(cb.ls), ls_sorted)

                    label = get(label_map, key, 0)
                    if label == 0
                        next_label += 1
                        label = next_label
                        label_map[key] = label
                        classified_dict[label] = SortedCounter{CoupledBases.CoupledBasis}()
                    end

                    count = orbit_basis_counts[cb]
                    push!(classified_dict[label], cb, count)
                end
            end
        end
    end

    return classified_dict
end

"""
    _listup_coupled_basislist(atom_list, lsum; isotropy=false) -> Vector{CoupledBases.CoupledBasis}

List up all coupled angular momentum basis functions for a given atom list and maximum angular momentum sum.

# Arguments
- `atom_list::Vector{<:Integer}`: List of atom indices
- `lsum::Integer`: Maximum sum of angular momentum values
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`
- `cache::CoupledBases.CouplingCache`: Memoization cache for coupling results,
  shared across calls within one basis construction; default: a fresh cache

# Returns
- `Vector{CoupledBases.CoupledBasis}`: List of coupled basis functions

# Note
- Skips l values less than the number of atoms or odd l values
- Generates all possible combinations of angular momenta that sum to `lsum` or less
"""
function _listup_coupled_basislist(
    atom_list::Vector{<:Integer},
    lsum::Integer;
    isotropy::Bool = false,
    cache::CoupledBases.CouplingCache = CoupledBases.CouplingCache(),
)::Vector{CoupledBases.CoupledBasis}
    result_basislist = Vector{CoupledBases.CoupledBasis}()
    for l = 2:lsum
        # Skip if l is less than number of atoms (each atom needs at least l=1)
        # or if l is odd (prohibited by time-reversal symmetry)
        if l < length(atom_list) || isodd(l)
            continue
        end
        l_list = Combinat.compositions(l, length(atom_list); min = 1)
        for l_vec::Vector{Int} in l_list
            cb_list = tesseral_coupled_bases_from_tesseral_bases(
                l_vec, atom_list; isotropy = isotropy, cache = cache,
            )
            append!(result_basislist, cb_list)
        end
    end
    return result_basislist
end

"""
Check if the product of two `CoupledBasis` objects is obviously zero.
This is a quick check for the step constructing projection matrix.
"""
function _is_obviously_zero_coupled_basis_product(
    cb1::CoupledBases.CoupledBasis,
    cb2::CoupledBases.CoupledBasis,
)::Bool
    if cb1.Lf != cb2.Lf
        return true
    end
    if cb1.ls != cb2.ls
        return true
    end
    if cb1.atoms != cb2.atoms
        return true
    end
    return false
end

"""
    _tensor_inner_product(tensor1::AbstractArray{T, N}, tensor2::AbstractArray{T, N}) where {T, N} -> Float64

Compute the inner product of two tensors.

Both tensors must have the same element type `T` and the same number of dimensions `N`.

# Arguments
- `tensor1::AbstractArray{T, N}`: First tensor
- `tensor2::AbstractArray{T, N}`: Second tensor (must have same element type and dimensions as tensor1)

# Returns
- `Float64`: Inner product of the two tensors

# Examples
```julia
t1 = [1.0 2.0; 3.0 4.0]
t2 = [5.0 6.0; 7.0 8.0]
_tensor_inner_product(t1, t2)  # 70.0
```
"""
function _tensor_inner_product(
    tensor1::AbstractArray{T, N},
    tensor2::AbstractArray{T, N},
) where {T, N}
    size(tensor1) == size(tensor2) || throw(DimensionMismatch(
        "tensor sizes must match: got $(size(tensor1)) and $(size(tensor2))",
    ))
    # Scalar accumulation (no `@simd`) gives a reduction order that does
    # not depend on the BLAS / LAPACK SIMD lane width, so the SALC
    # projection matrix is bit-identical across x86 and arm64 builds.
    acc = zero(T)
    @inbounds for i in eachindex(tensor1, tensor2)
        acc += conj(tensor1[i]) * tensor2[i]
    end
    return acc
end


"""
    _collect_cluster_atoms(coupled_basislist) -> Set{Vector{Int}}

Collect all unique atom lists from a list of coupled basis functions.

# Arguments
- `coupled_basislist::SortedCounter{CoupledBases.CoupledBasis}`: List of coupled basis functions

# Returns
- `Set{Vector{Int}}`: Set of unique atom lists (as vectors)
"""
function _collect_cluster_atoms(
    coupled_basislist::SortedCounter{CoupledBases.CoupledBasis},
)::Set{Vector{Int}}
    result_set = Set{Vector{Int}}()
    for cb in coupled_basislist
        push!(result_set, cb.atoms)
    end
    return result_set
end

"""
    _find_translation_atoms(atom_list, cluster_atoms, symmetry) -> Vector{Int}

Fold `atom_list` (same site order as a `CoupledBasis` before sorting) onto one of the
clusters in `cluster_atoms` using a supercell→primitive translation recipe built from
`map_s2p`, `symnum_translation`, and `map_sym_inv`.

The recipe yields a group action whose averaged projector is symmetric and idempotent;
naively double-looping over translations can pick inconsistent images at cell boundaries.

# Arguments
- `atom_list::Vector{<:Integer}`: Atom indices after a point-group map (`map_sym`), one per site
- `cluster_atoms::Set{Vector{Int}}`: Sorted atom lists appearing in the coupled-basis orbit
- `symmetry::Symmetry`: Symmetry data (must include `map_s2p`, `atoms_in_prim`)

# Returns
- `Vector{Int}`: Folded atom indices in **site order parallel to** `atom_list` (for `reorder_atoms`)

When several reference-site folds land on the same sorted cluster, the implementation picks one
whose **sort permutation matches** `sortperm(atom_list)`. Otherwise the first successful fold
(e.g. `[17, 1]` → `[1, 17]` in site order) would make `reorder_atoms` apply the identity
permutation to `ls`, incorrectly turning a site-swapping translation into the identity on the
coupled basis.

# Throws
- `ErrorException`: If no folding matches any cluster in `cluster_atoms`
"""
function _find_translation_atoms(
    atom_list::AbstractVector{<:Integer},
    cluster_atoms::Set{Vector{Int}},
    symmetry::Symmetry,
)::Vector{Int}
    # Read-only over atom_list, so skip the upfront defensive copy.
    atoms_in_prim = symmetry.atoms_in_prim
    map_s2p = symmetry.map_s2p
    symnum_translation = symmetry.symnum_translation
    map_sym_inv = symmetry.map_sym_inv

    p_target = sortperm(atom_list)
    fallback_candidate::Union{Nothing, Vector{Int}} = nothing

    # Reference-site recipe: collect all folds that hit `cluster_atoms`.
    for i in eachindex(atom_list)
        ref_atom = atom_list[i]
        ref_atom_in_prim = atoms_in_prim[map_s2p[ref_atom].atom]
        corresponding_translation = symnum_translation[map_s2p[ref_atom].translation]
        atom_translated = Vector{Int}(undef, length(atom_list))
        atom_translated[i] = ref_atom_in_prim
        for (n, new_atom) in enumerate(atom_list)
            if n == i
                continue
            end
            atom_translated[n] = map_sym_inv[new_atom, corresponding_translation]
        end
        if sort(atom_translated) ∈ cluster_atoms
            # Prefer a fold that preserves sort-order permutation.
            if sortperm(atom_translated) == p_target
                return atom_translated
            end
            if fallback_candidate === nothing
                fallback_candidate = atom_translated
            end
        end
    end

    if fallback_candidate === nothing
        error(
            "Failed to find translation atoms (map_s2p fold): atom_list=$atom_list, cluster_atoms=$cluster_atoms",
        )
    end
    return fallback_candidate
end


"""
    _projection_matrix_coupled_basis(coupled_basislist, symmetry) -> Matrix{Float64}

Construct the projection matrix for coupled basis functions under symmetry operations.

The projection matrix is computed by averaging over all symmetry operations (including time-reversal)
and projects onto the subspace of basis functions that are invariant under the symmetry group.

# Arguments
- `coupled_basislist::SortedCounter{CoupledBases.CoupledBasis}`: List of coupled basis functions
- `symmetry::Symmetry`: Symmetry information containing all symmetry operations

# Keywords
- `check_irrep_unitary::Bool`: If `true`, verify that each per-operation
  representation matrix `D(g)` (built before the average) is unitary, as
  required for a unitary irrep. Defaults to `CHECK_IRREP_UNITARY_DEFAULT`,
  which is read once at load time from the environment variable
  `MAGESTY_CHECK_IRREP_UNITARY` (`"1"` to enable, anything else off).

# Returns
- `Matrix{Float64}`: Projection matrix of size (nbasis * submatrix_dim) × (nbasis * submatrix_dim),
  where submatrix_dim = 2*Lf + 1 for the total angular momentum Lf

# Note
- The matrix is averaged over all symmetry operations and time-reversal symmetry
- Each basis function contributes a submatrix of size (2*Lf + 1) × (2*Lf + 1)
"""
function _projection_matrix_coupled_basis(
    coupled_basislist::SortedCounter{CoupledBases.CoupledBasis},
    symmetry::Symmetry,
    ;
    check_irrep_unitary::Bool = CHECK_IRREP_UNITARY_DEFAULT,
)::Matrix{Float64}
    Lf = coupled_basislist[1].Lf
    submatrix_dim = 2 * Lf + 1
    cluster_atoms::Set{Vector{Int}} = _collect_cluster_atoms(coupled_basislist)
    nbasis = length(coupled_basislist)
    full_matrix_dim = nbasis * submatrix_dim

    # Precompute Wigner D matrices Δl(Lf, …) once per symmetry op. base_rot_mat depends
    # only on (Lf, symop), not on time_rev_sym or the basis indices, so it is hoisted
    # out of the per-basis / per-time-reversal inner loops.
    base_rot_mats = Vector{Matrix{Float64}}(undef, length(symmetry.symdata))
    for (n, symop) in enumerate(symmetry.symdata)
        rotmat = symop.is_proper ? symop.rotation_cart : -1 * symop.rotation_cart
        α, β, γ = rotmat2euler(rotmat)
        base_rot_mats[n] = Δl(Lf, α, β, γ)
    end

    projection_mat = zeros(Float64, full_matrix_dim, full_matrix_dim)
    # Single-operation representation matrix D(g) for the current (symop,
    # time_rev_sym). Reused as a scratch buffer to avoid `2*nsym` heap
    # allocations of `full_matrix_dim^2` Float64 matrices; the accumulated
    # projector is the average of these `D(g)`s.
    representation_mat = zeros(Float64, full_matrix_dim, full_matrix_dim)
    # All cbs in this list share the same cluster size, so allocate once and refill.
    N_atoms = length(coupled_basislist[1].atoms)
    atoms_shifted_list = Vector{Int}(undef, N_atoms)
    for (n, symop) in enumerate(symmetry.symdata), time_rev_sym in (false, true)
        fill!(representation_mat, 0.0)
        base_rot_mat = base_rot_mats[n]

        for (i, cb1) in enumerate(coupled_basislist)
            @inbounds for k = 1:N_atoms
                atoms_shifted_list[k] = symmetry.map_sym[cb1.atoms[k], n]
            end
            primitive_atoms = _find_translation_atoms(atoms_shifted_list, cluster_atoms, symmetry)
            reordered_cb = reorder_atoms(cb1, primitive_atoms)
            multiplier = time_rev_sym ? (-1)^sum(reordered_cb.ls) : 1
            col_range = ((i-1)*submatrix_dim+1):(i*submatrix_dim)
            for (j, cb2) in enumerate(coupled_basislist)
                if _is_obviously_zero_coupled_basis_product(reordered_cb, cb2)
                    continue
                end
                phase = _tensor_inner_product(cb2.coeff_tensor, reordered_cb.coeff_tensor) / (2*Lf+1)
                row_range = ((j-1)*submatrix_dim+1):(j*submatrix_dim)
                # Scale the base rotation by the phase and the time-reversal sign
                # (`multiplier == 1` when `!time_rev_sym`) and write in place,
                # avoiding a temporary matrix per (symop, basis pair).
                factor = multiplier * phase
                @views representation_mat[row_range, col_range] .= base_rot_mat .* factor
            end
        end

        if check_irrep_unitary && !_is_unitary(representation_mat, tol = 1e-8)
            # A unitary irrep requires each per-operation D(g) to be unitary.
            # The accumulated projector itself is not unitary (it is
            # idempotent), so the check must run before the average.
            error(
                "Representation matrix D(g) is not unitary " *
                "(symmetry operation index n = $n). D(g)'D(g) = " *
                "$(representation_mat' * representation_mat)",
            )
        end

        projection_mat .+= representation_mat
    end

    # Average over all symmetry operations (2 for time reversal)
    return projection_mat ./ (2 * symmetry.nsym)
end


"""
    _is_translationally_equivalent_coupled_basis(
        cb1::CoupledBases.CoupledBasis,
        cb2::CoupledBases.CoupledBasis,
        symmetry::Symmetry,
    ) -> Bool

Check if two `CoupledBasis` objects are translationally equivalent.

Two `CoupledBasis` objects are translationally equivalent if:
- They are physically equivalent (same `Lf`, `Lseq`, and `(atom, l)` pairs)
- Their atom lists are related by a translation operation in the supercell

This function checks if the atom lists can be mapped to each other via translation operations
defined in `symmetry.symnum_translation`.

# Arguments
- `cb1::CoupledBases.CoupledBasis`: First `CoupledBasis` to compare
- `cb2::CoupledBases.CoupledBasis`: Second `CoupledBasis` to compare
- `symmetry::Symmetry`: Symmetry information containing translation mappings

# Returns
- `Bool`: `true` if the `CoupledBasis` objects are translationally equivalent, `false` otherwise
"""
function _is_translationally_equivalent_coupled_basis(
    cb1::CoupledBases.CoupledBasis,
    cb2::CoupledBases.CoupledBasis,
    symmetry::Symmetry,
)::Bool
    # Different number of sites
    length(cb1.ls) != length(cb2.ls) && return false

    # Different Lf
    cb1.Lf != cb2.Lf && return false

    # Different Lseq
    cb1.Lseq != cb2.Lseq && return false

    # Different ls values (as multisets) means different basis functions
    ls1_sorted = sort(cb1.ls)
    ls2_sorted = sort(cb2.ls)
    if ls1_sorted != ls2_sorted
        return false
    end

    # Check if (atom, l) pairs match as multisets (required for physical equivalence)
    atom_l_pairs1 = sort!(collect(zip(cb1.atoms, cb1.ls)))
    atom_l_pairs2 = sort!(collect(zip(cb2.atoms, cb2.ls)))
    if atom_l_pairs1 != atom_l_pairs2
        return false
    end

    # Check if coeff_tensor has the same size and values
    if size(cb1.coeff_tensor) != size(cb2.coeff_tensor)
        return false
    end
    if !isapprox(cb1.coeff_tensor, cb2.coeff_tensor, atol = 1e-10)
        return false
    end

    # Check if atom lists are translationally equivalent
    atom_list1 = cb1.atoms
    atom_list2 = cb2.atoms

    # Early return if atom lists are the same
    if atom_list1 == atom_list2
        return false
    end

    # Early return if first atom is the same
    # because this function is intended to be used for different first atoms but translationally equivalent clusters
    if atom_list1[1] == atom_list2[1]
        return false
    end

    # Check translation operations. `atom_list2` is fixed across the loop, so
    # sort it once; the shifted image of `atom_list1` is rebuilt in a reused
    # buffer (atom lists are tiny) instead of allocating a fresh vector per
    # itran. Both comparisons are multiset comparisons via sorted vectors.
    sorted_list2 = sort(atom_list2)
    shifted = Vector{Int}(undef, length(atom_list1))
    for itran in symmetry.symnum_translation
        # Method 1: forward translation (map_sym)
        for k in eachindex(atom_list1)
            shifted[k] = symmetry.map_sym[atom_list1[k], itran]
        end
        sort!(shifted) == sorted_list2 && return true

        # Method 2: inverse translation (map_sym_inv)
        for k in eachindex(atom_list1)
            shifted[k] = symmetry.map_sym_inv[atom_list1[k], itran]
        end
        sort!(shifted) == sorted_list2 && return true
    end

    return false
end


"""
    _is_proper_eigenvals(eigenval::AbstractVector; tol = 1e-8)::Bool

Checks whether all eigenvalues in the given vector are approximately 0 or 1.

# Arguments
- `eigenval::AbstractVector`: Vector of eigenvalues to check
- `tol::Real = 1e-8`: Tolerance for floating-point comparisons

# Returns
- `Bool`: `true` if all eigenvalues are approximately 0 or 1, `false` otherwise

# Examples
```julia
# Check eigenvalues from a projection matrix
eigenvals = [0.0, 1.0, 0.0, 1.0]
is_valid = _is_proper_eigenvals(eigenvals)  # true

# With some tolerance
eigenvals = [0.0, 1.0 + 1e-9, 0.0, 1.0]
is_valid = _is_proper_eigenvals(eigenvals)  # true

# Invalid eigenvalues
eigenvals = [0.0, 1.0, 0.5, 1.0]
is_valid = _is_proper_eigenvals(eigenvals)  # false
```
"""
function _is_proper_eigenvals(eigenval::AbstractVector; tol = 1e-8)::Bool
    return all(x -> isapprox(x, 0, atol = tol) || isapprox(x, 1, atol = tol), eigenval)
end





"""
    _is_unitary(mat::AbstractMatrix{<:Real}; tol::Float64 = 1e-10) -> Bool

Check if a matrix is unitary (i.e., UᵀU ≈ I and UUᵀ ≈ I within tolerance).

# Arguments
- `mat::AbstractMatrix{<:Real}`: Matrix to check
- `tol::Float64`: Tolerance for floating-point comparison (default: 1e-10)

# Returns
- `Bool`: `true` if the matrix is unitary within the tolerance, `false` otherwise
"""
function _is_unitary(mat::AbstractMatrix{<:Real}; tol::Float64 = 1e-10)::Bool
    size(mat, 1) == size(mat, 2) || return false
    return isapprox(mat' * mat, I, atol = tol) && isapprox(mat * mat', I, atol = tol)
end


"""
    _filter_basisdict(
        basisdict::Dict{Int, SortedCounter{CoupledBases.CoupledBasis}},
        symmetry::Symmetry,
    ) -> Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}

Filter out basis entries that correspond to clusters where all atoms map to the same
primitive atom, have odd total angular momentum Lf, and have multiplicity greater than 1.

# Description
For clusters consisting of atom sites located at cell boundaries that are connected by
pure translations, the symmetry-adapted linear combinations (SALCs) vanish when the total
angular momentum Lf is odd. This function filters out such basis entries beforehand to
avoid unnecessary computation in the SALC construction process.

A cluster is filtered if all of the following conditions are satisfied:
1. All atoms in the cluster map to the same atom in the primitive cell (i.e., they are
   translationally equivalent)
2. The total angular momentum Lf of the coupled basis is odd
3. The multiplicity (count) of the basis is greater than 1

The third condition ensures that only clusters at cell boundaries are filtered, as
translationally equivalent sites at cell boundaries always have multiplicity > 1 due to
periodic boundary conditions.

# Arguments
- `basisdict::Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}`: Dictionary of
  classified coupled basis sets, where keys are labels and values are lists of basis
  functions belonging to the same symmetry class. The `counts` field of each
  `SortedCounter` stores the multiplicity of each basis.
- `symmetry::Symmetry`: Symmetry information containing the mapping from supercell atoms
  to primitive cell atoms (`symmetry.map_s2p`)

# Returns
- `Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}`: Filtered dictionary with
  renumbered labels. Entries that satisfy all filtering conditions are removed, and
  remaining entries are assigned new consecutive labels starting from 1.

# Notes
- The function checks only the first basis in each `SortedCounter` to determine
  the filtering condition, as all bases in the same list share the same Lf value.
- The multiplicity is accessed via `basislist.counts[first_basis]`, which counts how many
  times the same basis appears due to translational symmetry.
- The input dictionary keys are sorted before processing to ensure deterministic output.
"""
function _filter_basisdict(
    basisdict::Dict{Int, SortedCounter{CoupledBases.CoupledBasis}},
    symmetry::Symmetry,
)::Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}
    result_basisdict = Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}()
    new_label = 1
    for label in sort!(collect(keys(basisdict)))
        basislist::SortedCounter{CoupledBases.CoupledBasis} = basisdict[label]
        first_basis = basislist[begin]
        atom_list = first_basis.atoms
        atom_list_in_prim = [symmetry.map_s2p[atom].atom for atom in atom_list]
        unique_atom_list_in_prim = unique(atom_list_in_prim)
        if length(unique_atom_list_in_prim) == 1 && isodd(first_basis.Lf) && basislist.counts[first_basis] > 1
            continue
        end

        # basislist is treated as immutable in this flow; avoid costly deep copy.
        result_basisdict[new_label] = basislist
        new_label += 1
    end
    return result_basisdict
end


end # module SALCBases
