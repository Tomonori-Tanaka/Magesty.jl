using Magesty.InputSpecs: SystemSpec, InteractionSpec, SymmetryOptions
using OffsetArrays
using Test

# Reaffirmation tests for the interaction-cluster cutoff criterion.
#
# The cluster definition implemented in `src/Clusters.jl` is uniform across
# every body order n >= 2: an n-body cluster is accepted iff every one of
# its C(n, 2) pairwise distances is within the body- and pair-specific
# cutoff `bodyn_cutoff[body, kd_i, kd_j]`. These tests pin that rule by
# constructing a fixture whose pairwise distances are known and asserting
# which clusters appear / disappear as we vary the cutoff matrix.
#
# Same rule, same test pattern, at n = 2, n = 3, n = 4. Do not read these
# tests as if n = 3 is the canonical case — the all-pairs criterion is not
# specific to any particular body order.

@testset "Cluster cutoff criterion (all-pairs rule, uniform for n >= 2)" begin
    # Fixture: 4 atoms of 4 distinct species along the x axis. Distinct
    # species disable symmetry reduction so atoms_in_prim = [1, 2, 3, 4]
    # and cluster_dict[n][prim] enumerates every accepted cluster starting
    # from each primitive atom.
    #
    # Lattice 20 x 20 x 20 (Angstrom). Fractional positions place atoms
    # at cartesian x = 1, 3, 5, 7 (others at y = z = 10). In-cell pairwise
    # cartesian distances:
    #
    #   pair  distance
    #   1-2     2
    #   1-3     4
    #   1-4     6
    #   2-3     2
    #   2-4     4
    #   3-4     2
    #
    # Periodic images live in cells displaced by +/- 20 A; the minimum
    # image is always the in-cell one for every pair, so the periodic-
    # image disambiguation branch of `is_within_cutoff` is a no-op here.
    kd_name = ["A", "B", "C", "D"]
    nkd = length(kd_name)
    kd_int_list = [1, 2, 3, 4]
    lattice_vectors = [20.0 0.0 0.0; 0.0 20.0 0.0; 0.0 0.0 20.0]
    x_fractional = hcat(
        [0.05, 0.5, 0.5],
        [0.15, 0.5, 0.5],
        [0.25, 0.5, 0.5],
        [0.35, 0.5, 0.5],
    )
    system = SystemSpec(
        name = "cutoff-criterion-fixture",
        num_atoms = 4,
        kd_name = kd_name,
        kd_int_list = kd_int_list,
        lattice_vectors = lattice_vectors,
        x_fractional = x_fractional,
    )
    structure = Magesty.Structure(system; verbosity = false)
    options = SymmetryOptions()
    symmetry = Magesty.Symmetry(structure, options; verbosity = false)

    # All four atoms should be primitive (no symmetry relates distinct species).
    @test sort(symmetry.atoms_in_prim) == [1, 2, 3, 4]

    # Helpers ----------------------------------------------------------
    function build_interaction(nbody::Int, cutoff_matrix_per_body::Dict{Int, Matrix{Float64}})
        # Body-1 lmax is irrelevant for cluster generation; supply zeros.
        body1_lmax = zeros(Int, nkd)
        if nbody == 1
            bodyn_lsum = OffsetArray(Int[], 2:1)
            bodyn_cutoff = OffsetArray(zeros(Float64, 0, nkd, nkd), 2:1, 1:nkd, 1:nkd)
        else
            bodyn_lsum = OffsetArray(fill(0, nbody - 1), 2:nbody)
            bodyn_cutoff = OffsetArray(zeros(Float64, nbody - 1, nkd, nkd),
                                       2:nbody, 1:nkd, 1:nkd)
            for n in 2:nbody
                M = cutoff_matrix_per_body[n]
                for i in 1:nkd, j in 1:nkd
                    bodyn_cutoff[n, i, j] = M[i, j]
                end
            end
        end
        return InteractionSpec(
            nbody = nbody,
            body1_lmax = body1_lmax,
            bodyn_lsum = bodyn_lsum,
            bodyn_cutoff = bodyn_cutoff,
            kd_name = kd_name,
        )
    end

    # cluster_dict[body][prim] :: OrderedDict{Vector{Int}, Int}. Keys are
    # atom-index lists of length `body`, starting with `prim` followed by
    # the partners in increasing order. `keys_for` returns those keys as a
    # sorted Vector of Vector{Int} for stable comparison.
    function keys_for(cluster, body, prim)
        haskey(cluster.cluster_dict, body) || return Vector{Vector{Int}}()
        haskey(cluster.cluster_dict[body], prim) || return Vector{Vector{Int}}()
        return sort(collect(keys(cluster.cluster_dict[body][prim])))
    end

    # Body-2: rule says one pair distance, must be within cutoff. ------
    @testset "body 2: pair distance within cutoff" begin
        # Cutoff matrix for body 2: uniform 3.0. Pairs at distance 2 pass;
        # pairs at distance 4 or 6 fail.
        rc2_tight = fill(3.0, nkd, nkd)
        inter = build_interaction(2, Dict(2 => rc2_tight))
        cluster = Magesty.Cluster(structure, symmetry, inter; verbosity = false)
        @test keys_for(cluster, 2, 1) == [[1, 2]]
        @test keys_for(cluster, 2, 2) == [[2, 1], [2, 3]]
        @test keys_for(cluster, 2, 3) == [[3, 2], [3, 4]]
        @test keys_for(cluster, 2, 4) == [[4, 3]]

        # Wider cutoff = 5.0: pairs at distance <= 5 pass (excludes only 1-4 at dist 6).
        rc2_wide = fill(5.0, nkd, nkd)
        inter = build_interaction(2, Dict(2 => rc2_wide))
        cluster = Magesty.Cluster(structure, symmetry, inter; verbosity = false)
        @test keys_for(cluster, 2, 1) == [[1, 2], [1, 3]]
        @test keys_for(cluster, 2, 4) == [[4, 2], [4, 3]]

        # Per-pair cutoff: only the (A, C) pair (atoms 1-3, distance 4) has a tight cutoff.
        rc2_mixed = fill(7.0, nkd, nkd)
        rc2_mixed[1, 3] = 3.0
        rc2_mixed[3, 1] = 3.0
        inter = build_interaction(2, Dict(2 => rc2_mixed))
        cluster = Magesty.Cluster(structure, symmetry, inter; verbosity = false)
        # Atom 1's pairs at 7.0: with atom 2 (d=2) yes, atom 3 (d=4) — but pair AC cutoff
        # is 3.0 < 4, so atom 3 rejected; atom 4 (d=6) yes (AD cutoff 7).
        @test keys_for(cluster, 2, 1) == [[1, 2], [1, 4]]
    end

    # Body-3: rule says all three pairs must be within cutoff. ---------
    @testset "body 3: every pair within cutoff" begin
        # Wide cutoff = 7.0 across the board: every triplet accepted.
        rc3_wide = fill(7.0, nkd, nkd)
        inter = build_interaction(3, Dict(2 => rc3_wide, 3 => rc3_wide))
        cluster = Magesty.Cluster(structure, symmetry, inter; verbosity = false)
        # From atom 1, all C(3,2) = 3 triplets accepted.
        @test keys_for(cluster, 3, 1) == [[1, 2, 3], [1, 2, 4], [1, 3, 4]]

        # Same wide cutoff for body 3, but tighten the (B, D) = (2, 4) pair
        # to 3.0. Then the triplet (1, 2, 4) loses pair (2, 4) = 4 > 3 and is
        # rejected, even though the prim-pair distances 1-2 and 1-4 are both
        # within the wide cutoff. The all-pairs rule is what excludes it.
        rc3_break_bd = copy(rc3_wide)
        rc3_break_bd[2, 4] = 3.0
        rc3_break_bd[4, 2] = 3.0
        inter = build_interaction(3, Dict(2 => rc3_wide, 3 => rc3_break_bd))
        cluster = Magesty.Cluster(structure, symmetry, inter; verbosity = false)
        # (1,2,4) gone; (1,2,3) and (1,3,4) remain.
        @test keys_for(cluster, 3, 1) == [[1, 2, 3], [1, 3, 4]]
    end

    # Body-4: same rule again, this time at n = 4. ---------------------
    @testset "body 4: every pair within cutoff" begin
        # Wide cutoff = 7.0: all C(4,2) = 6 pairs in the quartet (1,2,3,4) are
        # within cutoff, so the quartet is accepted (one cluster from atom 1).
        rc_wide = fill(7.0, nkd, nkd)
        inter = build_interaction(4,
            Dict(2 => rc_wide, 3 => rc_wide, 4 => rc_wide))
        cluster = Magesty.Cluster(structure, symmetry, inter; verbosity = false)
        @test keys_for(cluster, 4, 1) == [[1, 2, 3, 4]]

        # Tighten a single non-primitive pair: (C, D) = (3, 4), distance 2,
        # cutoff lowered to 1.5. Pair (3, 4) now fails — even though every
        # pair involving the prim atom 1 passes the wide cutoff. The quartet
        # is rejected.
        rc_break_cd = copy(rc_wide)
        rc_break_cd[3, 4] = 1.5
        rc_break_cd[4, 3] = 1.5
        inter = build_interaction(4,
            Dict(2 => rc_wide, 3 => rc_wide, 4 => rc_break_cd))
        cluster = Magesty.Cluster(structure, symmetry, inter; verbosity = false)
        @test keys_for(cluster, 4, 1) == Vector{Vector{Int}}()
    end
end
