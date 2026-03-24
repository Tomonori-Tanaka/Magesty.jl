#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using Profile
using TOML
using Magesty

"""
Micro-benchmark and line-profile for Optimize.jl hotspots.

Usage:
  julia test/benchmark_optimize_hotspots.jl
  julia test/benchmark_optimize_hotspots.jl --input test/examples/fege_2x2x2/input.toml
  julia test/benchmark_optimize_hotspots.jl --samples 20 --evals 1 --profile-iters 40
"""
function parse_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :input => joinpath(@__DIR__, "examples", "fege_2x2x2", "input.toml"),
        :samples => 15,
        :evals => 1,
        :profile_iters => 30,
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
            cfg[:profile_iters] = parse(Int, args[i + 1]); i += 2
        elseif a == "--verbose"
            cfg[:verbosity] = true; i += 1
        else
            error("Unknown argument: $a")
        end
    end
    return cfg
end

function load_context(cfg::Dict{Symbol, Any})
    input_path = abspath(cfg[:input])
    input = TOML.parsefile(input_path)
    workdir = dirname(input_path)

    system = cd(workdir) do
        Magesty.System(input; verbosity = cfg[:verbosity])
    end
    spincluster = cd(workdir) do
        Magesty.SpinCluster(system, input; verbosity = cfg[:verbosity])
    end

    embset_relpath = input["regression"]["datafile"]
    spinconfig_list = cd(workdir) do
        Magesty.read_embset(embset_relpath)
    end

    salc_list = spincluster.basisset.salc_list
    symmetry = spincluster.symmetry
    num_atoms = spincluster.structure.supercell.num_atoms

    sc = spinconfig_list[1]
    key_group = salc_list[1]
    cbc = key_group[1]
    iatom = cbc.atoms[1]

    return (
        input_path = input_path,
        spincluster = spincluster,
        spinconfig_list = spinconfig_list,
        salc_list = salc_list,
        symmetry = symmetry,
        num_atoms = num_atoms,
        sc = sc,
        key_group = key_group,
        cbc = cbc,
        iatom = iatom,
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
    println("spinconfigs: ", length(ctx.spinconfig_list))
    println("salc key groups: ", length(ctx.salc_list))
    println("num atoms: ", ctx.num_atoms)
    println("benchmark samples: ", cfg[:samples], ", evals: ", cfg[:evals])
    println("profile iters: ", cfg[:profile_iters])

    println("\nWarming up target functions...")
    Magesty.Optimize.calc_∇ₑu(
        ctx.cbc,
        ctx.iatom,
        ctx.sc.spin_directions,
        ctx.symmetry,
        ctx.key_group,
    )
    Magesty.Optimize.build_design_matrix_torque(
        ctx.salc_list,
        [ctx.sc],
        ctx.num_atoms,
        ctx.symmetry,
    )

    println("\n=== BenchmarkTools results ===")
    nsamples = cfg[:samples]
    nevals = cfg[:evals]

    bench_grad = @benchmark Magesty.Optimize.calc_∇ₑu(
        $(ctx.cbc),
        $(ctx.iatom),
        $(ctx.sc.spin_directions),
        $(ctx.symmetry),
        $(ctx.key_group),
    ) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_grad)
    println()

    bench_torque = @benchmark Magesty.Optimize.build_design_matrix_torque(
        $(ctx.salc_list),
        [$(ctx.sc)],
        $(ctx.num_atoms),
        $(ctx.symmetry),
    ) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_torque)
    println()

    iters = cfg[:profile_iters]

    print_profile_for(() -> begin
        for _ in 1:iters
            Magesty.Optimize.calc_∇ₑu(
                ctx.cbc,
                ctx.iatom,
                ctx.sc.spin_directions,
                ctx.symmetry,
                ctx.key_group,
            )
        end
    end, "calc_∇ₑu")

    print_profile_for(() -> begin
        for _ in 1:iters
            Magesty.Optimize.build_design_matrix_torque(
                ctx.salc_list,
                [ctx.sc],
                ctx.num_atoms,
                ctx.symmetry,
            )
        end
    end, "build_design_matrix_torque(1 spinconfig)")

    return nothing
end

function main()
    cfg = parse_args(ARGS)
    run_bench(cfg)
end

main()
