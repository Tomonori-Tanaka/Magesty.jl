#!/usr/bin/env julia
# Benchmark thread scaling of build_design_matrix_torque.
# Run with: julia -t N test/benchmark_threads.jl  (N = 1, 2, 4)

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using TOML
using Magesty
using Printf
using Statistics

const INPUT_PATH = joinpath(@__DIR__, "examples", "fege_2x2x2", "input.toml")
const NREPS = 1  # timing repetitions (after 1 warmup)

function median_elapsed(f, n)
    f()  # warmup
    times = [(@elapsed f()) for _ in 1:n]
    return median(times)
end

function fmt(t_s)
    t_ms = t_s * 1e3
    t_ms >= 1000 ? @sprintf("%.3f s ", t_ms / 1000) : @sprintf("%.3f ms", t_ms)
end

function main()
    nthreads = Threads.nthreads()
    println("nthreads = $nthreads")
    println("Loading context...")
    flush(stdout)

    input    = TOML.parsefile(INPUT_PATH)
    workdir  = dirname(INPUT_PATH)
    system   = cd(workdir) do; Magesty.System(input; verbosity=false); end
    spincluster = cd(workdir) do; Magesty.SpinCluster(system, input; verbosity=false); end

    spinconfig_list = spincluster.optimize.spinconfig_list
    salc_list       = spincluster.basisset.salc_list
    symmetry        = spincluster.symmetry
    num_atoms       = spincluster.structure.supercell.num_atoms

    println("  spinconfigs : $(length(spinconfig_list))")
    println("  salc groups : $(length(salc_list))")
    println("  num atoms   : $num_atoms")
    println("Benchmarking ($NREPS reps + 1 warmup)...")
    flush(stdout)

    # --- build_design_matrix_energy (sequential, no @threads) ---
    t_energy = median_elapsed(NREPS) do
        Magesty.Optimize.build_design_matrix_energy(salc_list, spinconfig_list, symmetry)
    end

    # --- build_design_matrix_torque (@threads over spinconfigs) ---
    t_torque = median_elapsed(NREPS) do
        Magesty.Optimize.build_design_matrix_torque(salc_list, spinconfig_list, num_atoms, symmetry)
    end

    println()
    println("=" ^ 55)
    @printf("  %-38s  %s\n", "build_design_matrix_energy (serial)", fmt(t_energy))
    @printf("  %-38s  %s\n", "build_design_matrix_torque (@threads)", fmt(t_torque))
    println("=" ^ 55)
    flush(stdout)
end

main()
