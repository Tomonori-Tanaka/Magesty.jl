#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using Profile
using TOML
using Magesty

"""
Micro-benchmark and line-profile for SALCBases.jl hotspots.

Usage:
  julia test/benchmark_basisset_hotspots.jl
  julia test/benchmark_basisset_hotspots.jl --input test/examples/fege_2x2x2/input.toml
  julia test/benchmark_basisset_hotspots.jl --samples 10 --evals 1 --profile-iters-basisset 2
"""
function parse_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :input => joinpath(@__DIR__, "examples", "fege_2x2x2", "input.toml"),
        :samples => 10,
        :evals => 1,
        :profile_iters_listup => 100,
        :profile_iters_projection => 5,
        :profile_iters_basisset => 2,
        :verbosity => false,
    )

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--input"
            cfg[:input] = args[i + 1]; i += 2
        elseif a == "--samples"
            cfg[:samples] = parse(Int, args[i + 1]); i += 2
        elseif a == "--evals"
            cfg[:evals] = parse(Int, args[i + 1]); i += 2
        elseif a == "--profile-iters"
            v = parse(Int, args[i + 1])
            cfg[:profile_iters_listup] = v
            cfg[:profile_iters_projection] = v
            cfg[:profile_iters_basisset] = v
            i += 2
        elseif a == "--profile-iters-listup"
            cfg[:profile_iters_listup] = parse(Int, args[i + 1]); i += 2
        elseif a == "--profile-iters-projection"
            cfg[:profile_iters_projection] = parse(Int, args[i + 1]); i += 2
        elseif a == "--profile-iters-basisset"
            cfg[:profile_iters_basisset] = parse(Int, args[i + 1]); i += 2
        elseif a == "--verbose"
            cfg[:verbosity] = true; i += 1
        else
            error("Unknown argument: $a")
        end
    end
    return cfg
end

function pick_atom_list_and_lsum(cluster, interaction)
    # Prefer a 2-body orbit example for listup benchmark.
    if haskey(cluster.cluster_orbits_dict, 2)
        orbit_map = cluster.cluster_orbits_dict[2]
        first_orbit_clusters = first(values(orbit_map))
        atom_list = sort(first(first_orbit_clusters))
        return atom_list, interaction.bodyn_lsum[2]
    end
    # Fallback to a tiny synthetic 2-site list if no 2-body data exists.
    return [1, 2], max(2, interaction.bodyn_lsum[min(2, interaction.nbody)])
end

# Simplified classifier for `CoupledBasis` lists used in this benchmark.
# Groups basis functions by `(nbody, Lf, sum(ls), Tuple(sort(ls)))`, ignoring
# spatial symmetry — this gives a deterministic, fixture-independent grouping.
function _classify_coupled_basislist(coupled_basislist)
    CoupledBases = Magesty.Basis
    SortedCounter = Magesty.SortedCounters.SortedCounter

    if isempty(coupled_basislist)
        return Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}()
    end

    label_map = Dict{Any, Int}()
    label_list = Vector{Int}(undef, length(coupled_basislist))
    next_label = 0
    for (idx, cb) in enumerate(coupled_basislist)
        key = (length(cb.atoms), cb.Lf, sum(cb.ls), Tuple(sort(cb.ls)))
        label = get(label_map, key, 0)
        if label == 0
            next_label += 1
            label = next_label
            label_map[key] = label
        end
        label_list[idx] = label
    end

    dict = Dict{Int, SortedCounter{CoupledBases.CoupledBasis}}()
    for label in 1:next_label
        dict[label] = SortedCounter{CoupledBases.CoupledBasis}()
    end
    for (cb, label) in zip(coupled_basislist, label_list)
        count_val =
            coupled_basislist isa SortedCounter ? get(coupled_basislist.counts, cb, 1) : 1
        push!(dict[label], cb, count_val)
    end
    return dict
end

function load_context(cfg::Dict{Symbol, Any})
    input_path = abspath(cfg[:input])
    input = TOML.parsefile(input_path)
    workdir = dirname(input_path)

    system_spec, interaction, options = Magesty.parse_toml_inputs(input)
    structure, symmetry, cluster = cd(workdir) do
        Magesty._build_structure_skeleton(system_spec, options, interaction;
                                          verbosity = cfg[:verbosity])
    end
    salcbasis = cd(workdir) do
        Magesty.SALCBases.SALCBasis(structure, symmetry, cluster, interaction, options;
                                    verbosity = cfg[:verbosity])
    end

    classified = _classify_coupled_basislist(salcbasis.coupled_basislist)
    keys_sorted = sort(collect(keys(classified)))
    key = keys_sorted[1]
    representative_group = classified[key]
    atom_list, lsum = pick_atom_list_and_lsum(cluster, interaction)

    return (
        input_path = input_path,
        structure = structure,
        symmetry = symmetry,
        cluster = cluster,
        interaction = interaction,
        options = options,
        atom_list = atom_list,
        lsum = lsum,
        representative_group = representative_group,
        isotropy = options.isotropy,
    )
end

function print_profile_for(f::Function, label::AbstractString)
    Profile.clear()
    @profile f()
    println("\n[$label] profile (flat, mincount=10)")
    Profile.print(; format = :flat, sortedby = :count, mincount = 10)
    println("\n[$label] profile (tree, mincount=10)")
    Profile.print(; format = :tree, mincount = 10)
end

function run_bench(cfg::Dict{Symbol, Any})
    ctx = load_context(cfg)

    println("Input TOML: ", ctx.input_path)
    println("listup atom_list: ", ctx.atom_list, ", lsum: ", ctx.lsum)
    println("representative coupled_basis count: ", length(ctx.representative_group))
    println("benchmark samples: ", cfg[:samples], ", evals: ", cfg[:evals])
    println(
        "profile iters: listup=",
        cfg[:profile_iters_listup],
        ", projection=",
        cfg[:profile_iters_projection],
        ", basisset=",
        cfg[:profile_iters_basisset],
    )

    println("\nWarming up target functions...")
    Magesty.SALCBases.listup_coupled_basislist(
        ctx.atom_list,
        ctx.lsum;
        isotropy = ctx.isotropy,
    )
    Magesty.SALCBases.projection_matrix_coupled_basis(ctx.representative_group, ctx.symmetry)
    Magesty.SALCBasis(
        ctx.structure,
        ctx.symmetry,
        ctx.cluster,
        ctx.interaction, ctx.options;
        verbosity = false,
    )

    println("\n=== BenchmarkTools results ===")
    nsamples = cfg[:samples]
    nevals = cfg[:evals]

    bench_listup = @benchmark Magesty.SALCBases.listup_coupled_basislist(
        $(ctx.atom_list),
        $(ctx.lsum);
        isotropy = $(ctx.isotropy),
    ) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_listup)
    println()

    bench_projection = @benchmark Magesty.SALCBases.projection_matrix_coupled_basis(
        $(ctx.representative_group),
        $(ctx.symmetry),
    ) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_projection)
    println()

    bench_basisset = @benchmark Magesty.SALCBasis(
        $(ctx.structure),
        $(ctx.symmetry),
        $(ctx.cluster),
        $(ctx.interaction), $(ctx.options);
        verbosity = false,
    ) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_basisset)
    println()

    print_profile_for(() -> begin
        for _ in 1:cfg[:profile_iters_listup]
            Magesty.SALCBases.listup_coupled_basislist(
                ctx.atom_list,
                ctx.lsum;
                isotropy = ctx.isotropy,
            )
        end
    end, "listup_coupled_basislist")

    print_profile_for(() -> begin
        for _ in 1:cfg[:profile_iters_projection]
            Magesty.SALCBases.projection_matrix_coupled_basis(
                ctx.representative_group,
                ctx.symmetry,
            )
        end
    end, "projection_matrix_coupled_basis")

    print_profile_for(() -> begin
        for _ in 1:cfg[:profile_iters_basisset]
            Magesty.SALCBasis(
                ctx.structure,
                ctx.symmetry,
                ctx.cluster,
                ctx.interaction, ctx.options;
                verbosity = false,
            )
        end
    end, "SALCBasis constructor")

    return nothing
end

function main()
    cfg = parse_args(ARGS)
    run_bench(cfg)
end

main()
