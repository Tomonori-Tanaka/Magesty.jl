#!/usr/bin/env julia
# Bottleneck profiling for build_design_matrix_energy / build_design_matrix_torque
# on the 3-body FeGe 2x2x2 fixtures.
#
# Usage:
#   julia --project=bench bench/bench_b1_3body_fege.jl
#   julia --project=bench bench/bench_b1_3body_fege.jl --fixture fege_open
#
# Flags:
#   --fixture <name>   Choose a fixture (default: light). One of: fege_open, light.
#   --ntrials <n>      Number of timed trials per matrix (default: 3).
#   --no-profile       Skip the Profile dump after timing.

using Magesty
using Magesty.Fitting: build_design_matrix_energy, build_design_matrix_torque
using Magesty.SpinConfigs: read_embset
using Profile
using Statistics
using Printf

const REPO_ROOT   = joinpath(@__DIR__, "..")
const EMBSET_PATH = joinpath(REPO_ROOT, "test", "integration", "fege_2x2x2", "EMBSET")

const FIXTURES = Dict(
    "light"     => joinpath(@__DIR__, "fixtures", "fege_2x2x2_3body_light",    "input.toml"),
    "fege_open" => joinpath(@__DIR__, "fixtures", "fege_2x2x2_3body_fefe_open", "input.toml"),
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
function parse_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :fixture     => "light",
        :ntrials     => 3,
        :profile     => true,
        :profile_delay => 0.001,
        :mincount    => 8,
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--fixture"
            cfg[:fixture] = args[i + 1]; i += 2
        elseif a == "--no-profile"
            cfg[:profile] = false; i += 1
        elseif a == "--ntrials"
            cfg[:ntrials] = parse(Int, args[i + 1]); i += 2
        else
            error("Unknown argument: $a")
        end
    end
    return cfg
end

# ---------------------------------------------------------------------------
# Simple timing helper (no BenchmarkTools needed)
# ---------------------------------------------------------------------------
function bench_one(label::String, f::Function, ntrials::Int)
    f()  # warm-up
    times_s  = Float64[]
    allocs_b = Float64[]
    counts   = Int[]
    for _ in 1:ntrials
        stats = @timed f()
        push!(times_s,  stats.time)
        push!(allocs_b, stats.bytes)
        push!(counts,   @allocations f())
    end
    @printf("  %-38s  min=%7.3fs  med=%7.3fs  allocs(med)=%d  bytes(med)=%.2e\n",
            label,
            minimum(times_s),
            median(times_s),
            Int(median(counts)),
            median(allocs_b))
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    cfg  = parse_args(ARGS)
    fixture_name = cfg[:fixture]
    fixture_path = get(FIXTURES, fixture_name, nothing)
    fixture_path === nothing && error("Unknown fixture: $fixture_name. Choose from: $(sort(collect(keys(FIXTURES))))")

    println("="^70)
    println("Fixture : ", fixture_name, "  (", fixture_path, ")")
    println("EMBSET  : ", EMBSET_PATH)
    println("Threads : ", Threads.nthreads())
    println("="^70)

    # ----- build SCEBasis (includes Cluster + SALC construction) -----------
    println("\nBuilding SCEBasis (Cluster + SALC construction)...")
    t_basis = @elapsed basis = SCEBasis(fixture_path; verbosity = false)
    @printf("  SCEBasis build: %.3f s\n", t_basis)

    println("\nLoading EMBSET...")
    spinconfigs = read_embset(EMBSET_PATH)

    num_spinconfigs = length(spinconfigs)
    num_salcs       = length(basis.salcbasis.salc_list)
    num_atoms       = basis.structure.supercell.num_atoms
    salc_list       = basis.salcbasis.salc_list
    symmetry        = basis.symmetry

    # Report SALC body distribution
    nbody_counts = Dict{Int, Int}()
    for kg in salc_list
        n = length(kg[1].atoms)
        nbody_counts[n] = get(nbody_counts, n, 0) + 1
    end

    println()
    println("num_spinconfigs = ", num_spinconfigs)
    println("num_salcs       = ", num_salcs)
    println("  by nbody      = ", sort(collect(nbody_counts)))
    println("num_atoms       = ", num_atoms)
    println()

    # ----- Timing -----------------------------------------------------------
    ntrials = cfg[:ntrials]
    println("=== Timing ($(ntrials) trials each) ===")
    bench_one("build_design_matrix_energy", () ->
        build_design_matrix_energy(salc_list, spinconfigs, symmetry),
        ntrials)
    bench_one("build_design_matrix_torque", () ->
        build_design_matrix_torque(salc_list, spinconfigs, num_atoms, symmetry),
        ntrials)

    # ----- Profile ----------------------------------------------------------
    if cfg[:profile]
        println()
        println("=== Profile: build_design_matrix_energy ===")
        Profile.clear()
        Profile.init(; delay = cfg[:profile_delay])
        build_design_matrix_energy(salc_list, spinconfigs, symmetry)  # ensure compiled
        Profile.@profile build_design_matrix_energy(salc_list, spinconfigs, symmetry)
        println("\n--- flat, sortedby=:count, mincount=$(cfg[:mincount]) ---")
        Profile.print(; format = :flat, sortedby = :count, mincount = cfg[:mincount])

        println()
        println("=== Profile: build_design_matrix_torque ===")
        Profile.clear()
        Profile.init(; delay = cfg[:profile_delay])
        build_design_matrix_torque(salc_list, spinconfigs, num_atoms, symmetry)
        Profile.@profile build_design_matrix_torque(salc_list, spinconfigs, num_atoms, symmetry)
        println("\n--- flat, sortedby=:count, mincount=$(cfg[:mincount]) ---")
        Profile.print(; format = :flat, sortedby = :count, mincount = cfg[:mincount])
    end
end

main()
