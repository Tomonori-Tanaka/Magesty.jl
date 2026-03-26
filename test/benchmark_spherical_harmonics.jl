#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using Profile
using LinearAlgebra
using Magesty
using Magesty.MySphericalHarmonics: P̄ₗₘ, dP̄ₗₘ, dP̄ₗₘ_unsafe, Yₗₘ, Yₗₘ_unsafe, Zₗₘ, Zₗₘ_unsafe, ∂ᵢZlm, ∂ᵢZlm_unsafe

"""
Micro-benchmark and line-profile for MySphericalHarmonics.jl.
Compares safe APIs (with `validate_lm` / `validate_uvec`) vs `*_unsafe` where applicable.
Note: `P̄ₗₘ` has no validation in the library, so only `dP̄ₗₘ` is shown for Legendre-derivative safe/unsafe.

Usage:
  julia test/benchmark_spherical_harmonics.jl
  julia test/benchmark_spherical_harmonics.jl --lmax 4
  julia test/benchmark_spherical_harmonics.jl --samples 20 --evals 3 --profile-iters 500
"""
function parse_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :lmax => 4,
        :samples => 15,
        :evals => 3,
        :profile_iters_plm => 2000,
        :profile_iters_zlm => 1000,
        :profile_iters_grad => 500,
    )

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--lmax"
            cfg[:lmax] = parse(Int, args[i + 1]); i += 2
        elseif a == "--samples"
            cfg[:samples] = parse(Int, args[i + 1]); i += 2
        elseif a == "--evals"
            cfg[:evals] = parse(Int, args[i + 1]); i += 2
        elseif a == "--profile-iters"
            v = parse(Int, args[i + 1])
            cfg[:profile_iters_plm] = v
            cfg[:profile_iters_zlm] = v
            cfg[:profile_iters_grad] = v
            i += 2
        elseif a == "--profile-iters-plm"
            cfg[:profile_iters_plm] = parse(Int, args[i + 1]); i += 2
        elseif a == "--profile-iters-zlm"
            cfg[:profile_iters_zlm] = parse(Int, args[i + 1]); i += 2
        elseif a == "--profile-iters-grad"
            cfg[:profile_iters_grad] = parse(Int, args[i + 1]); i += 2
        else
            error("Unknown argument: $a")
        end
    end
    return cfg
end

# All (l, m) pairs up to lmax
function lm_pairs(lmax::Int)
    return [(l, m) for l in 0:lmax for m in -l:l]
end

# Sweep over all (l, m) pairs for a given function and uvec
function sweep_zlm(lmax::Int, uvec::Vector{Float64})
    s = 0.0
    for l in 0:lmax, m in -l:l
        s += Zₗₘ(l, m, uvec)
    end
    return s
end

function sweep_zlm_unsafe(lmax::Int, uvec::Vector{Float64})
    s = 0.0
    for l in 0:lmax, m in -l:l
        s += Zₗₘ_unsafe(l, m, uvec)
    end
    return s
end

function sweep_plm(lmax::Int, z::Float64)
    s = 0.0
    for l in 0:lmax, m in 0:l
        s += P̄ₗₘ(l, m, z)
    end
    return s
end

function sweep_grad_zlm(lmax::Int, uvec::Vector{Float64})
    s = zeros(3)
    for l in 0:lmax, m in -l:l
        s .+= ∂ᵢZlm(l, m, uvec)
    end
    return s
end

function sweep_grad_zlm_unsafe(lmax::Int, uvec::Vector{Float64})
    s = zeros(3)
    for l in 0:lmax, m in -l:l
        s .+= ∂ᵢZlm_unsafe(l, m, uvec)
    end
    return s
end

"""Print median times and speedup factor (how many × faster the second benchmark is)."""
function print_speedup(label::AbstractString, bench_safe, bench_unsafe)
    t_s = median(bench_safe).time
    t_u = median(bench_unsafe).time
    ratio = t_s / t_u
    println("\n[$label] safe vs unsafe (median time per sample)")
    println("  safe:   $(t_s) ns")
    println("  unsafe: $(t_u) ns")
    println("  unsafe is ~$(round(ratio, digits=2))× faster than safe")
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
    lmax = cfg[:lmax]
    nsamples = cfg[:samples]
    nevals = cfg[:evals]

    # Representative inputs
    uvec = normalize([1.0, 2.0, 3.0])
    z = uvec[3]
    pairs = lm_pairs(lmax)

    println("=== MySphericalHarmonics benchmark ===")
    println("lmax: $lmax  ($(length(pairs)) (l,m) pairs)")
    println("uvec: $uvec")
    println("samples: $nsamples,  evals: $nevals")
    println()

    # --- Warm-up ---
    println("Warming up...")
    sweep_plm(lmax, z)
    sweep_zlm(lmax, uvec)
    sweep_zlm_unsafe(lmax, uvec)
    sweep_grad_zlm(lmax, uvec)
    sweep_grad_zlm_unsafe(lmax, uvec)
    println()

    # --- BenchmarkTools ---
    println("=== BenchmarkTools results ===")

    println("\n--- P̄ₗₘ sweep (l=0..$(lmax), m=0..l) ---")
    bench_plm = @benchmark sweep_plm($lmax, $z) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_plm)
    println()

    println("\n--- Zₗₘ sweep (l=0..$(lmax), all m) ---")
    bench_zlm = @benchmark sweep_zlm($lmax, $uvec) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_zlm)
    println()
    bench_zlm_u = @benchmark sweep_zlm_unsafe($lmax, $uvec) samples=nsamples evals=nevals
    println("\n--- Zₗₘ_unsafe sweep (same work) ---")
    show(stdout, MIME"text/plain"(), bench_zlm_u)
    println()
    print_speedup("Zₗₘ full sweep", bench_zlm, bench_zlm_u)

    println("\n--- ∂ᵢZlm sweep (l=0..$(lmax), all m) ---")
    bench_grad = @benchmark sweep_grad_zlm($lmax, $uvec) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_grad)
    println()
    bench_grad_u = @benchmark sweep_grad_zlm_unsafe($lmax, $uvec) samples=nsamples evals=nevals
    println("\n--- ∂ᵢZlm_unsafe sweep (same work) ---")
    show(stdout, MIME"text/plain"(), bench_grad_u)
    println()
    print_speedup("∂ᵢZlm full sweep", bench_grad, bench_grad_u)

    # Single-call benchmarks for a representative (l, m) = (lmax, lmax-1)
    l_rep = lmax
    m_rep = lmax - 1
    println("\n--- Single call: P̄ₗₘ(l=$l_rep, m=$(abs(m_rep)), z) ---")
    bench_plm_single = @benchmark P̄ₗₘ($l_rep, $(abs(m_rep)), $z) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_plm_single)
    println()

    println("\n--- Single call: Yₗₘ(l=$l_rep, m=$m_rep, uvec) ---")
    bench_ylm_single = @benchmark Yₗₘ($l_rep, $m_rep, $uvec) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_ylm_single)
    println()
    bench_ylm_single_u = @benchmark Yₗₘ_unsafe($l_rep, $m_rep, $uvec) samples=nsamples evals=nevals
    println("\n--- Single call: Yₗₘ_unsafe ---")
    show(stdout, MIME"text/plain"(), bench_ylm_single_u)
    println()
    print_speedup("Yₗₘ single call", bench_ylm_single, bench_ylm_single_u)

    println("\n--- Single call: Zₗₘ(l=$l_rep, m=$m_rep, uvec) ---")
    bench_zlm_single = @benchmark Zₗₘ($l_rep, $m_rep, $uvec) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_zlm_single)
    println()
    bench_zlm_single_u = @benchmark Zₗₘ_unsafe($l_rep, $m_rep, $uvec) samples=nsamples evals=nevals
    println("\n--- Single call: Zₗₘ_unsafe (same l, m, uvec) ---")
    show(stdout, MIME"text/plain"(), bench_zlm_single_u)
    println()
    print_speedup("Zₗₘ single call", bench_zlm_single, bench_zlm_single_u)

    println("\n--- Single call: ∂ᵢZlm(l=$l_rep, m=$m_rep, uvec) ---")
    bench_grad_single = @benchmark ∂ᵢZlm($l_rep, $m_rep, $uvec) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_grad_single)
    println()
    bench_grad_single_u = @benchmark ∂ᵢZlm_unsafe($l_rep, $m_rep, $uvec) samples=nsamples evals=nevals
    println("\n--- Single call: ∂ᵢZlm_unsafe ---")
    show(stdout, MIME"text/plain"(), bench_grad_single_u)
    println()
    print_speedup("∂ᵢZlm single call", bench_grad_single, bench_grad_single_u)

    println("\n--- Single call: dP̄ₗₘ (still validates l,m,r̂z) vs dP̄ₗₘ_unsafe ---")
    bench_dp = @benchmark dP̄ₗₘ($l_rep, $(abs(m_rep)), $z) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_dp)
    println()
    bench_dp_u = @benchmark dP̄ₗₘ_unsafe($l_rep, $(abs(m_rep)), $z) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_dp_u)
    println()
    print_speedup("dP̄ₗₘ single call", bench_dp, bench_dp_u)

    # --- Profile ---
    print_profile_for(() -> begin
        for _ in 1:cfg[:profile_iters_plm]
            sweep_plm(lmax, z)
        end
    end, "P̄ₗₘ sweep")

    print_profile_for(() -> begin
        for _ in 1:cfg[:profile_iters_zlm]
            sweep_zlm(lmax, uvec)
        end
    end, "Zₗₘ sweep")

    print_profile_for(() -> begin
        for _ in 1:cfg[:profile_iters_grad]
            sweep_grad_zlm(lmax, uvec)
        end
    end, "∂ᵢZlm sweep")

    return nothing
end

function main()
    cfg = parse_args(ARGS)
    run_bench(cfg)
end

main()
