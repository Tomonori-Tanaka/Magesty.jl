using Test
using Serialization

# Run the worker once per thread count and compare the serialized results.
# The worker (`_thread_safety_worker.jl`) is invoked through `julia -t N`
# subprocesses so that `@threads` loops actually execute with `N` threads.
# This catches data races as well as any accidental dependence on the
# scheduling order across the `@threads` partitions used by
# `build_design_matrix_energy`, `build_design_matrix_torque`,
# `Symmetries`, and `SALCBases.SALCBasis`.
@testset "thread safety" begin
    worker = joinpath(@__DIR__, "_thread_safety_worker.jl")
    project = dirname(dirname(@__DIR__))  # repo root (contains Project.toml)

    candidate_thread_counts = [1, 2, 4]
    max_cpu = max(1, Sys.CPU_THREADS)
    thread_counts = unique(min.(candidate_thread_counts, max_cpu))

    payloads = Dict{Int, Any}()
    mktempdir() do tmp
        for n in thread_counts
            out = joinpath(tmp, "result_t$(n).jls")
            julia_cmd = `$(Base.julia_cmd()) --project=$(project) -t $(n) $(worker) $(out)`
            run(pipeline(julia_cmd; stdout = devnull, stderr = stderr))
            payloads[n] = open(deserialize, out)
        end
    end

    @testset "worker actually used $(n) thread(s)" for n in thread_counts
        @test payloads[n].nthreads == n
    end

    ref_n = first(thread_counts)
    ref = payloads[ref_n]
    others = filter(!=(ref_n), thread_counts)

    @testset "design matrices invariant ($(ref_n) vs $(n) threads)" for n in others
        p = payloads[n]
        @test size(p.X_E) == size(ref.X_E)
        @test size(p.X_T) == size(ref.X_T)
        @test isapprox(p.X_E, ref.X_E; rtol = 1e-12)
        @test isapprox(p.X_T, ref.X_T; rtol = 1e-12)
    end

    @testset "fit coefficients invariant ($(ref_n) vs $(n) threads)" for n in others
        p = payloads[n]
        @test isapprox(p.j0, ref.j0; rtol = 1e-12)
        @test length(p.jphi) == length(ref.jphi)
        @test isapprox(p.jphi, ref.jphi; rtol = 1e-12)
    end

    @testset "predictions invariant ($(ref_n) vs $(n) threads)" for n in others
        p = payloads[n]
        @test isapprox(p.e_pred, ref.e_pred; rtol = 1e-12)
        @test size(p.t_pred) == size(ref.t_pred)
        @test isapprox(p.t_pred, ref.t_pred; rtol = 1e-12)
    end
end
