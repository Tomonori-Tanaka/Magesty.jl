#!/usr/bin/env julia
# Fitting design-matrix benchmark.
# Times `build_design_matrix_energy` and `build_design_matrix_torque`
# on the fept_tetragonal_2x2x2 example.
#
# Usage:
#   julia --project=bench bench/bench_b1_design_matrix.jl

using Magesty
using Magesty.Fitting: build_design_matrix_energy, build_design_matrix_torque
using Magesty.SpinConfigs: read_embset
using Statistics
using Printf

const EXAMPLE_DIR = joinpath(@__DIR__, "..", "test", "integration", "fept_tetragonal_2x2x2")
const SYSTEM_XML  = joinpath(EXAMPLE_DIR, "system.xml")
const EMBSET_PATH = joinpath(EXAMPLE_DIR, "EMBSET")
const NTRIALS     = 5

function bench_one(label::String, f::Function)
    f()  # warm-up
    times_s    = Float64[]
    allocs_b   = Float64[]
    alloc_cnt  = Float64[]
    for _ in 1:NTRIALS
        stats = @timed f()
        push!(times_s, stats.time)
        push!(allocs_b, stats.bytes)
        push!(alloc_cnt, stats.gctime)  # we ignore gctime; reuse field
    end
    # @timed gives time, bytes, gctime; allocation count not directly given.
    # Use @allocations for count.
    counts = Float64[]
    for _ in 1:NTRIALS
        push!(counts, Float64(@allocations f()))
    end
    @printf("%-32s  min=%.4fs  med=%.4fs  bytes(med)=%.2e  allocs(med)=%d\n",
            label,
            minimum(times_s),
            median(times_s),
            median(allocs_b),
            Int(median(counts)))
end

function main()
    println("Loading SCEBasis from ", SYSTEM_XML)
    basis = Magesty.load(SCEBasis, SYSTEM_XML)
    println("Loading EMBSET from ", EMBSET_PATH)
    spinconfigs = read_embset(EMBSET_PATH)
    println("num_spinconfigs = ", length(spinconfigs))
    println("num_salcs       = ", length(basis.salcbasis.salc_list))
    println("num_atoms       = ", basis.structure.supercell.num_atoms)
    println()

    println("=== Fitting design-matrix benchmark (5 trials each) ===")
    bench_one("build_design_matrix_energy", () ->
        build_design_matrix_energy(
            basis.salcbasis.salc_list,
            spinconfigs,
            basis.symmetry,
        ))
    bench_one("build_design_matrix_torque", () ->
        build_design_matrix_torque(
            basis.salcbasis.salc_list,
            spinconfigs,
            basis.structure.supercell.num_atoms,
            basis.symmetry,
        ))
end

main()
