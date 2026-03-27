#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using TOML
using Magesty

function parse_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :input => joinpath(@__DIR__, "examples", "fege_2x2x2", "input.toml"),
        :samples => 10,
        :evals => 1,
        :verbosity => false,
        :bench_fit => false,
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
        elseif a == "--verbose"
            cfg[:verbosity] = true; i += 1
        elseif a == "--with-fit"
            cfg[:bench_fit] = true; i += 1
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
    spinconfig_list = spincluster.optimize.spinconfig_list

    salc_list = spincluster.basisset.salc_list
    symmetry = spincluster.symmetry
    num_atoms = spincluster.structure.supercell.num_atoms
    sc = spinconfig_list[1]
    key_group = salc_list[1]
    cbc = key_group[1]
    iatom = cbc.atoms[1]

    return (
        input_path = input_path,
        system = system,
        spincluster = spincluster,
        spinconfig_list = spinconfig_list,
        salc_list = salc_list,
        symmetry = symmetry,
        num_atoms = num_atoms,
        sc = sc,
        cbc = cbc,
        iatom = iatom,
    )
end

function run_bench(cfg::Dict{Symbol, Any})
    ctx = load_context(cfg)
    regression = TOML.parsefile(ctx.input_path)["regression"]
    estimator = Magesty.ElasticNet(alpha = regression["alpha"], lambda = regression["lambda"])
    weight = regression["weight"]

    println("Input TOML: ", ctx.input_path)
    println("spinconfigs: ", length(ctx.spinconfig_list))
    println("salc key groups: ", length(ctx.salc_list))
    println("num atoms: ", ctx.num_atoms)
    println("benchmark samples: ", cfg[:samples], ", evals: ", cfg[:evals])

    nsamples = cfg[:samples]
    nevals = cfg[:evals]

    println("\n=== calc_∇ₑu ===")
    bench_grad = @benchmark Magesty.Optimize.calc_∇ₑu(
        $(ctx.cbc),
        $(ctx.iatom),
        $(ctx.sc.spin_directions),
        $(ctx.symmetry),
    ) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_grad)
    println()

    println("\n=== build_design_matrix_torque (1 spinconfig) ===")
    bench_torque = @benchmark Magesty.Optimize.build_design_matrix_torque(
        $(ctx.salc_list),
        [$(ctx.sc)],
        $(ctx.num_atoms),
        $(ctx.symmetry),
    ) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_torque)
    println()

    if cfg[:bench_fit]
        println("\n=== fit_sce_model (full dataset) ===")
        bench_fit = @benchmark Magesty.fit_sce_model(
            $(ctx.system),
            $(ctx.spinconfig_list),
            $(estimator),
            $(weight),
            verbosity = false,
        ) samples=nsamples evals=nevals
        show(stdout, MIME"text/plain"(), bench_fit)
        println()
    end
end

function main()
    cfg = parse_args(ARGS)
    run_bench(cfg)
end

main()
