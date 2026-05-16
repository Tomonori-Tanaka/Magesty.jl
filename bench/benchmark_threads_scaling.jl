#!/usr/bin/env julia
# Single-thread-count datapoint for thread scaling of
# `build_design_matrix_energy` and `build_design_matrix_torque`.
#
# Usage:
#   julia --project=bench -t N bench/benchmark_threads_scaling.jl
#
# ENV overrides:
#   BENCH_INPUT     — input.toml path (default: fept_tetragonal_2x2x2 fixture)
#   BENCH_SAMPLES   — BenchmarkTools samples (default: 5)
#   BENCH_SECONDS   — per-trial time budget (default: 10.0)
#
# Output (one block per kernel, parseable by run_threads_scaling.sh):
#   THREADS=4
#   build_design_matrix_energy: median=12.300 ms minimum=11.800 ms allocs=9876
#   build_design_matrix_torque: median=42.100 ms minimum=40.500 ms allocs=112345

using BenchmarkTools
using Magesty
using Magesty.Fitting: build_design_matrix_energy, build_design_matrix_torque
using Printf
using TOML

const DEFAULT_INPUT = joinpath(
    @__DIR__, "..", "test", "examples", "fept_tetragonal_2x2x2", "input.toml",
)
const INPUT_TOML = get(ENV, "BENCH_INPUT", DEFAULT_INPUT)
const SAMPLES = parse(Int, get(ENV, "BENCH_SAMPLES", "5"))
const SECONDS = parse(Float64, get(ENV, "BENCH_SECONDS", "10.0"))

function load_context()
    input = TOML.parsefile(INPUT_TOML)
    basis = SCEBasis(input; verbosity = false)
    embset = joinpath(dirname(INPUT_TOML), "EMBSET.dat")
    spinconfigs = Magesty.SpinConfigs.read_embset(embset)
    num_atoms = basis.structure.supercell.num_atoms
    return basis, spinconfigs, num_atoms
end

function format_trial(label::AbstractString, t::BenchmarkTools.Trial)
    @printf(
        "%s: median=%.3f ms minimum=%.3f ms allocs=%d\n",
        label,
        median(t).time / 1e6,
        minimum(t).time / 1e6,
        median(t).allocs,
    )
end

function main()
    basis, spinconfigs, num_atoms = load_context()
    salc_list = basis.salcbasis.salc_list
    symmetry = basis.symmetry

    # Warm up so JIT/precompile cost is not folded into the timing.
    build_design_matrix_energy(salc_list, spinconfigs, symmetry)
    build_design_matrix_torque(salc_list, spinconfigs, num_atoms, symmetry)

    println("THREADS=", Threads.nthreads())
    println("INPUT=", INPUT_TOML)
    println("n_spinconfigs=", length(spinconfigs),
            "  n_salcs=", length(salc_list),
            "  n_atoms=", num_atoms)

    t_e = @benchmark build_design_matrix_energy($salc_list, $spinconfigs, $symmetry) samples =
        SAMPLES seconds = SECONDS evals = 1
    format_trial("build_design_matrix_energy", t_e)

    t_t = @benchmark build_design_matrix_torque(
        $salc_list, $spinconfigs, $num_atoms, $symmetry,
    ) samples = SAMPLES seconds = SECONDS evals = 1
    format_trial("build_design_matrix_torque", t_t)
end

main()
