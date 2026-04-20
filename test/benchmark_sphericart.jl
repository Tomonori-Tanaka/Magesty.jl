#!/usr/bin/env julia
# Benchmark: MySphericalHarmonics (Zₗₘ / ∂ᵢZlm) vs SpheriCart
#
# MySphericalHarmonics computes one (l,m) at a time.
# SpheriCart computes all harmonics up to lmax in a single batched call.
#
# Run via Makefile (uses the test environment where SpheriCart is available):
#   make bench-sphericart
#
# ENV overrides:
#   BENCH_LMAX=6 BENCH_NPOINTS=500 make bench-sphericart

using BenchmarkTools
using LinearAlgebra
using Printf
using StaticArrays
using Magesty
using Magesty.MySphericalHarmonics: Zₗₘ_unsafe, ∂ᵢZlm_unsafe
using SpheriCart

# ── config from ENV ───────────────────────────────────────────────────────────

const LMAX    = parse(Int, get(ENV, "BENCH_LMAX",    "4"))
const NPOINTS = parse(Int, get(ENV, "BENCH_NPOINTS", "100"))
const SAMPLES = parse(Int, get(ENV, "BENCH_SAMPLES", "10"))
const EVALS   = parse(Int, get(ENV, "BENCH_EVALS",   "3"))

# ── sweep functions (MySphericalHarmonics) ────────────────────────────────────

function my_sweep_zlm(lmax::Int, uvec::Vector{Float64})
    s = 0.0
    for l in 0:lmax, m in -l:l
        s += Zₗₘ_unsafe(l, m, uvec)
    end
    return s
end

function my_sweep_zlm_npoints(lmax::Int, pts::Vector{Vector{Float64}})
    s = 0.0
    for uvec in pts
        for l in 0:lmax, m in -l:l
            s += Zₗₘ_unsafe(l, m, uvec)
        end
    end
    return s
end

function my_sweep_grad(lmax::Int, uvec::Vector{Float64})
    s = SVector(0.0, 0.0, 0.0)
    for l in 0:lmax, m in -l:l
        s = s + ∂ᵢZlm_unsafe(l, m, uvec)
    end
    return s
end

function my_sweep_grad_npoints(lmax::Int, pts::Vector{Vector{Float64}})
    s = SVector(0.0, 0.0, 0.0)
    for uvec in pts
        for l in 0:lmax, m in -l:l
            s = s + ∂ᵢZlm_unsafe(l, m, uvec)
        end
    end
    return s
end

# ── display ───────────────────────────────────────────────────────────────────

function print_comparison(label::String, bench_mine, bench_sph)
    t_mine = median(bench_mine).time
    t_sph  = median(bench_sph).time
    ratio  = t_mine / t_sph
    faster = ratio >= 1 ? "SpheriCart          " : "MySphericalHarmonics"
    println("  [$label]")
    @printf "    MySphericalHarmonics : %8.2f μs\n" t_mine / 1e3
    @printf "    SpheriCart           : %8.2f μs\n" t_sph  / 1e3
    @printf "    → %s is ~%.2f× faster\n\n" faster max(ratio, 1 / ratio)
end

# ── main ──────────────────────────────────────────────────────────────────────

function run_bench()
    lmax  = LMAX
    npts  = NPOINTS
    ns    = SAMPLES
    ne    = EVALS
    nharm = (lmax + 1)^2

    println("=" ^ 64)
    println("  MySphericalHarmonics vs SpheriCart")
    println("  lmax = $lmax  →  $nharm harmonics per point")
    println("  N = $npts points,  samples = $ns,  evals = $ne")
    println("=" ^ 64)

    # Inputs: single point
    uvec_vec = normalize([1.0, 2.0, 3.0])
    uvec_sv  = SVector{3,Float64}(uvec_vec...)

    # Inputs: N points
    rng_pts  = [normalize(randn(3)) for _ in 1:npts]
    pts_vec  = rng_pts                              # Vector{Vector{Float64}}
    pts_sv   = SVector{3,Float64}.(rng_pts)         # Vector{SVector{3,Float64}}

    sph = SphericalHarmonics(lmax)

    # warm-up
    print("Warming up ... ")
    my_sweep_zlm(lmax, uvec_vec);           my_sweep_grad(lmax, uvec_vec)
    my_sweep_zlm_npoints(lmax, pts_vec);    my_sweep_grad_npoints(lmax, pts_vec)
    SpheriCart.compute(sph, uvec_sv);       SpheriCart.compute(sph, pts_sv)
    SpheriCart.compute_with_gradients(sph, uvec_sv)
    SpheriCart.compute_with_gradients(sph, pts_sv)
    println("done.\n")

    # ── 1 point: Zₗₘ values ──────────────────────────────────────────────────
    println("─── Single point: Zₗₘ values ───")

    b_my_1   = @benchmark my_sweep_zlm($lmax, $uvec_vec) samples=ns evals=ne
    b_sph_1  = @benchmark SpheriCart.compute($sph, $uvec_sv) samples=ns evals=ne

    println("MySphericalHarmonics:")
    show(stdout, MIME"text/plain"(), b_my_1);  println()
    println("SpheriCart:")
    show(stdout, MIME"text/plain"(), b_sph_1); println()
    print_comparison("1 point, values", b_my_1, b_sph_1)

    # ── N points: Zₗₘ values ─────────────────────────────────────────────────
    println("─── $npts points: Zₗₘ values ───")

    b_my_n   = @benchmark my_sweep_zlm_npoints($lmax, $pts_vec) samples=ns evals=ne
    b_sph_n  = @benchmark SpheriCart.compute($sph, $pts_sv)     samples=ns evals=ne

    println("MySphericalHarmonics:")
    show(stdout, MIME"text/plain"(), b_my_n);  println()
    println("SpheriCart:")
    show(stdout, MIME"text/plain"(), b_sph_n); println()
    print_comparison("$npts points, values", b_my_n, b_sph_n)

    # ── 1 point: gradients ────────────────────────────────────────────────────
    println("─── Single point: ∂ᵢZlm gradients ───")

    b_my_g1  = @benchmark my_sweep_grad($lmax, $uvec_vec) samples=ns evals=ne
    b_sph_g1 = @benchmark SpheriCart.compute_with_gradients($sph, $uvec_sv) samples=ns evals=ne

    println("MySphericalHarmonics:")
    show(stdout, MIME"text/plain"(), b_my_g1);  println()
    println("SpheriCart:")
    show(stdout, MIME"text/plain"(), b_sph_g1); println()
    print_comparison("1 point, gradients", b_my_g1, b_sph_g1)

    # ── N points: gradients ───────────────────────────────────────────────────
    println("─── $npts points: ∂ᵢZlm gradients ───")

    b_my_gn  = @benchmark my_sweep_grad_npoints($lmax, $pts_vec) samples=ns evals=ne
    b_sph_gn = @benchmark SpheriCart.compute_with_gradients($sph, $pts_sv) samples=ns evals=ne

    println("MySphericalHarmonics:")
    show(stdout, MIME"text/plain"(), b_my_gn);  println()
    println("SpheriCart:")
    show(stdout, MIME"text/plain"(), b_sph_gn); println()
    print_comparison("$npts points, gradients", b_my_gn, b_sph_gn)

    println("=" ^ 64)
end

run_bench()
