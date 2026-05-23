#!/usr/bin/env julia

using Profile
using TOML
using Magesty

"""
Per-stage benchmark and line-profile for `Cluster` construction.

`Cluster` construction (`src/Clusters.jl`) is a fixed four-stage pipeline
(in execution order):

    set_mindist_pairs -> generate_clusters -> irreducible_clusters -> cluster_orbits

This script times each stage individually with `@timed` and then runs a
single `Profile`-sampled `Cluster(...)` call for a flat/tree breakdown, so
the stage that dominates three-body cluster generation can be identified
with numbers. The printed stage table preserves the execution order.

`@elapsed` / `@timed` is used instead of `BenchmarkTools.@benchmark`: the
heavy three-body case takes thousands of seconds and cannot be sampled
repeatedly. One warmup call pays the JIT cost; one measured pass per stage
follows.

`--input` may be passed more than once; each fixture is benchmarked in
turn so several systems can be compared in one invocation.

Usage:
  julia bench/benchmark_cluster.jl
  julia bench/benchmark_cluster.jl --input test/integration/fege_2x2x2/input.toml
  julia bench/benchmark_cluster.jl --input a/input.toml --input b/input.toml
  julia bench/benchmark_cluster.jl --profile-delay 1.0
  julia bench/benchmark_cluster.jl --no-profile
"""
function parse_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :inputs => String[],
        :profile_delay => 0.001,
        :profile => true,
        :verbosity => false,
    )

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--input"
            push!(cfg[:inputs], args[i + 1]); i += 2
        elseif a == "--profile-delay"
            cfg[:profile_delay] = parse(Float64, args[i + 1]); i += 2
        elseif a == "--no-profile"
            cfg[:profile] = false; i += 1
        elseif a == "--verbose"
            cfg[:verbosity] = true; i += 1
        else
            error("Unknown argument: $a")
        end
    end
    # Default: all three-body fixtures shipped under bench/fixtures/,
    # in increasing cost order.
    if isempty(cfg[:inputs])
        cfg[:inputs] = [
            joinpath(@__DIR__, "fixtures", "fege_2x2x2_3body_light",
                     "input.toml"),
            joinpath(@__DIR__, "fixtures", "fege_2x2x2_3body_fefe_open",
                     "input.toml"),
            joinpath(@__DIR__, "fixtures", "fege_2x2x2_3body_all_open",
                     "input.toml"),
        ]
    end
    return cfg
end

# Counts emitted alongside the stage timings, so a slow stage can be read
# against the size of the data it processes.
function count_clusters(cluster_dict)
    counts = Dict{Int, Int}()
    for (body, per_prim) in cluster_dict
        counts[body] = sum(length(d) for d in values(per_prim); init = 0)
    end
    return counts
end

count_irreducible(irr) = Dict(body => length(sc) for (body, sc) in irr)

count_orbits(orb) = Dict(body => length(orbit_map) for (body, orbit_map) in orb)

function load_context(cfg::Dict{Symbol, Any}, input_path::AbstractString)
    input_path = abspath(input_path)
    input = TOML.parsefile(input_path)
    workdir = dirname(input_path)

    system_spec, interaction, options = Magesty.parse_toml_inputs(input)
    structure, symmetry = cd(workdir) do
        structure = Magesty.Structures.Structure(
            system_spec; verbosity = cfg[:verbosity],
        )
        symmetry = Magesty.Symmetries.Symmetry(
            structure, options; verbosity = cfg[:verbosity],
        )
        return structure, symmetry
    end

    return (
        input_path = input_path,
        structure = structure,
        symmetry = symmetry,
        interaction = interaction,
    )
end

function run_bench(cfg::Dict{Symbol, Any}, input_path::AbstractString)
    ctx = load_context(cfg, input_path)
    Clusters = Magesty.Clusters

    structure = ctx.structure
    symmetry = ctx.symmetry
    interaction = ctx.interaction
    cutoff_radii = interaction.bodyn_cutoff
    nbody = interaction.nbody

    println("Input TOML: ", ctx.input_path)
    println("nbody: ", nbody, ", num_atoms: ", structure.supercell.num_atoms)
    println("cutoff_radii (body axis ", axes(cutoff_radii, 1), "):")
    for body in axes(cutoff_radii, 1)
        println("  body ", body, ": ", collect(cutoff_radii[body, :, :]))
    end

    # Warmup: one full construction to pay the JIT cost before measuring.
    # With --verbose this also prints the constructor's own `Time Elapsed`,
    # which can be cross-checked against the per-stage total below.
    println("\nWarming up (one full Cluster construction)...")
    Clusters.Cluster(structure, symmetry, interaction;
                     verbosity = cfg[:verbosity])

    # --- Per-stage wall-time measurement -------------------------------
    # The stage functions are called in the same order as the `Cluster`
    # constructor, reproducing its execution path.
    println("\n=== Per-stage timing ===")

    # `min_distance_pairs` is computed once and threaded into both
    # `generate_clusters` and the `Cluster` constructor; the standalone
    # measurement below reproduces that single call.
    mdp = @timed Clusters.set_mindist_pairs(
        structure.supercell.num_atoms,
        structure.x_image_cart,
        structure.exist_image;
        tol = symmetry.tol,
    )
    min_distance_pairs = mdp.value

    gen = @timed Clusters.generate_clusters(
        structure, symmetry, cutoff_radii, nbody, min_distance_pairs,
    )
    cluster_dict = gen.value

    irr = @timed Clusters.irreducible_clusters(cluster_dict, symmetry)
    irreducible_cluster_dict = irr.value

    orb = @timed Clusters.cluster_orbits(irreducible_cluster_dict, symmetry)
    cluster_orbits_dict = orb.value

    stages = [
        ("set_mindist_pairs", mdp.time, mdp.bytes),
        ("generate_clusters", gen.time, gen.bytes),
        ("irreducible_clusters", irr.time, irr.bytes),
        ("cluster_orbits", orb.time, orb.bytes),
    ]
    total_time = sum(s[2] for s in stages)

    println(rpad("stage", 22), rpad("time [s]", 16), rpad("alloc [MiB]", 14),
            "share")
    for (name, t, bytes) in stages
        share = total_time > 0 ? 100 * t / total_time : 0.0
        println(
            rpad(name, 22),
            rpad(string(round(t; digits = 4)), 16),
            rpad(string(round(bytes / 2^20; digits = 1)), 14),
            string(round(share; digits = 1), "%"),
        )
    end
    println(rpad("TOTAL", 22), round(total_time; digits = 4), " s")

    # --- Output-size metrics -------------------------------------------
    println("\n=== Cluster counts ===")
    cdc = count_clusters(cluster_dict)
    irc = count_irreducible(irreducible_cluster_dict)
    orc = count_orbits(cluster_orbits_dict)
    for body in sort(collect(keys(cdc)))
        println(
            "  body ", body,
            ": raw=", cdc[body],
            ", irreducible=", get(irc, body, 0),
            ", orbits=", get(orc, body, 0),
        )
    end

    # `set_mindist_pairs` is computed once and reused by both
    # `generate_clusters` and the `Cluster` constructor; the timing above
    # reflects that single shared call.
    dominant = stages[argmax([s[2] for s in stages])]
    println("\n=== Conclusion ===")
    println(
        "Dominant stage: ", dominant[1],
        " (", round(dominant[2]; digits = 4), " s, ",
        round(100 * dominant[2] / total_time; digits = 1), "% of total)",
    )

    # --- Profile a single full construction ----------------------------
    if cfg[:profile]
        println("\n=== Profile (full Cluster construction) ===")
        Profile.clear()
        Profile.init(; delay = cfg[:profile_delay])
        full = @timed Profile.@profile Clusters.Cluster(
            structure, symmetry, interaction; verbosity = false,
        )
        # Cross-check: a compiled full construction should be close to the
        # per-stage total (the stage functions reproduce the constructor
        # path). Profiling adds minor sampling overhead.
        println(
            "full construction (compiled): ", round(full.time; digits = 4),
            " s  vs per-stage TOTAL ", round(total_time; digits = 4), " s",
        )
        println("\nprofile (flat, mincount=10)")
        Profile.print(; format = :flat, sortedby = :count, mincount = 10)
        println("\nprofile (tree, mincount=10)")
        Profile.print(; format = :tree, mincount = 10)
    end

    return nothing
end

function main()
    cfg = parse_args(ARGS)
    for (idx, input_path) in enumerate(cfg[:inputs])
        println("\n", "#"^70)
        println("# Fixture ", idx, " / ", length(cfg[:inputs]), ": ",
                input_path)
        println("#"^70)
        run_bench(cfg, input_path)
    end
end

main()
