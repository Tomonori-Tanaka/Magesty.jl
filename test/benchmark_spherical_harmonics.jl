#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using Profile
using LinearAlgebra
using Magesty
using Magesty: MySphericalHarmonics
using Magesty.MySphericalHarmonics: P̄ₗₘ, dP̄ₗₘ, Yₗₘ, Zₗₘ,
    ∂Zₗₘ_∂r̂x, ∂Zₗₘ_∂r̂y, ∂Zₗₘ_∂r̂z, zzₗₘ,
    ∂Zₗₘ_∂x, ∂Zₗₘ_∂y, ∂Zₗₘ_∂z, ∂ᵢZlm

"""
Micro-benchmark and line-profile for MySphericalHarmonics.jl.

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
    sweep_grad_zlm(lmax, uvec)
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

    println("\n--- ∂ᵢZlm sweep (l=0..$(lmax), all m) ---")
    bench_grad = @benchmark sweep_grad_zlm($lmax, $uvec) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_grad)
    println()

    # Single-call benchmarks for a representative (l, m) = (lmax, lmax-1)
    l_rep = lmax
    m_rep = lmax - 1
    println("\n--- Single call: P̄ₗₘ(l=$l_rep, m=$(abs(m_rep)), z) ---")
    bench_plm_single = @benchmark P̄ₗₘ($l_rep, $(abs(m_rep)), $z) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_plm_single)
    println()

    println("\n--- Single call: Zₗₘ(l=$l_rep, m=$m_rep, uvec) ---")
    bench_zlm_single = @benchmark Zₗₘ($l_rep, $m_rep, $uvec) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_zlm_single)
    println()

    println("\n--- Single call: ∂ᵢZlm(l=$l_rep, m=$m_rep, uvec) ---")
    bench_grad_single = @benchmark ∂ᵢZlm($l_rep, $m_rep, $uvec) samples=nsamples evals=nevals
    show(stdout, MIME"text/plain"(), bench_grad_single)
    println()

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
