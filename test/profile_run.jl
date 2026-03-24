#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Profile
using TOML
using Logging
using Magesty

"""
Run a lightweight profile for Magesty model construction.

Usage examples:
  julia test/profile_run.jl
  julia test/profile_run.jl --input test/examples/fege_2x2x2/input.toml --build-iters 10
"""
function parse_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :input => joinpath(@__DIR__, "examples", "fege_2x2x2", "input.toml"),
        :build_iters => 5,
        :verbosity => false,
    )

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--input"
            cfg[:input] = args[i + 1]; i += 2
        elseif a == "--build-iters"
            cfg[:build_iters] = parse(Int, args[i + 1]); i += 2
        elseif a == "--verbose"
            cfg[:verbosity] = true; i += 1
        else
            error("Unknown argument: $a")
        end
    end
    return cfg
end

function run_profile(cfg::Dict{Symbol, Any})
    input_path = abspath(cfg[:input])
    @info "Loading input TOML" input_path
    input = TOML.parsefile(input_path)
    workdir = dirname(input_path)

    println("Input TOML: ", input_path)
    println("build iters: ", cfg[:build_iters])
    build_log_step = max(1, cfg[:build_iters] ÷ 10)

    Profile.clear()
    @info "Build-only profile started" total = cfg[:build_iters] step = build_log_step
    @profile begin
        for i in 1:cfg[:build_iters]
            system = cd(workdir) do
                Magesty.System(input; verbosity = cfg[:verbosity])
            end
            cd(workdir) do
                Magesty.SpinCluster(system, input; verbosity = cfg[:verbosity])
            end
            if i == 1 || i % build_log_step == 0 || i == cfg[:build_iters]
                @info "Build profile progress" done = i total = cfg[:build_iters] percent = round(100 * i / cfg[:build_iters]; digits = 1)
            end
        end
    end
    @info "Build-only profile finished"

    println("\nTop profile results (flat):")
    Profile.print(; format = :flat, sortedby = :count, mincount = 20)
    println("\nTop profile results (tree):")
    Profile.print(; format = :tree, mincount = 20)
    return nothing
end

function main()
    cfg = parse_args(ARGS)
    run_profile(cfg)
end

main()
