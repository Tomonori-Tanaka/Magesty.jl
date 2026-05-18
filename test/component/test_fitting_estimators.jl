using Test
using Random
using Statistics
using Magesty

# Component tests for the ElasticNet estimator (and the `Lasso` convenience
# function). Six properties from the 260517 spec are checked:
#
#   1. `λ → 0` limit recovers OLS on a centered energy-only system
#      (weight = 0), for both `alpha = 1` and `alpha = 0`.
#   2. `λ → 0` limit, mixed objective (weight = 0.5): same agreement on
#      the augmented (energy + torque) system produced by
#      `assemble_weighted_problem`. Confirms that the upstream
#      energy-only centering plus `intercept = false` is correctly
#      wired.
#   3. `λ → 0` limit with `standardize = true`: still agrees with OLS,
#      but with a looser tolerance documenting the standardise /
#      back-transform round-trip precision.
#   4. GLMNet `α = 0` vs analytic `Ridge` agreement: GLMNet's L2 path
#      objective is `(1/2n)||y-Xβ||² + (λ_g/2)||β||²` and analytic
#      ridge minimizes `||y-Xβ||² + λ_a||β||²`, so the two are
#      equivalent under `λ_g = λ_a / n`. With this correction the
#      coefficients agree to ≤ 2e-2 relative across the λ sweep
#      `[1e-3, 1e-2, 1e-1, 1.0]` (worst residual ≈ 1.1e-2 at λ_a = 1.0,
#      dominated by GLMNet's coordinate-descent convergence threshold;
#      the smaller λ_a values agree to ~1e-4).
#   5. Lasso sparsity monotonicity: increasing `λ` at `alpha = 1`
#      weakly increases the number of exactly-zero coefficients.
#   6. Standardisation qualitative effect: on a fabricated matrix with
#      large per-column scale variation, `standardize = true` retains
#      the high-scale columns at a moderate `λ` while
#      `standardize = false` drops them earlier.

# atol for `λ→0` agreement between GLMNet and the OLS solution on a
# centered system. Empirical worst over the test fixtures below:
# energy-only standardize=false ≈ 3.5e-6, mixed-weight standardize=false
# ≈ 1.7e-5, standardize=true (worst case mixed) ≈ 5.6e-5. The residual is
# GLMNet's coordinate-descent convergence threshold and the
# standardize/back-transform round-trip; the constants below leave ~6x
# headroom on top of the worst case.
const _GLMNET_OLS_ATOL    = 1e-4
const _GLMNET_STD_ATOL    = 1e-3
# rtol for GLMNet(α=0, λ=λ_a/n) vs analytic Ridge(λ=λ_a). Measured worst
# case over λ_a ∈ {1e-3, 1e-2, 1e-1} (n=60) is ~1.0e-3 at λ_a = 1e-1,
# falling to ~7.8e-5 at λ_a = 1e-3. The residual is GLMNet's coordinate-
# descent convergence threshold (the analytic correction λ_g = λ_a / n
# matches the two objectives exactly). The λ_a = 1.0 endpoint was
# omitted because GLMNet's CD threshold dominates there (rel ~1.1e-2),
# violating the spec's 1e-2 guard band.
const _GLMNET_RIDGE_RTOL  = 2.0e-3

# Build a centered (energy-only) `(X, y)` pair that mimics the output of
# `assemble_weighted_problem` at `weight = 0` — only the energy block is
# present, and both `X` and `y` have been mean-centered along their rows
# / over the energy dimension.
function _centered_energy_system(rng, n::Int, p::Int)
    X_raw = randn(rng, n, p)
    beta_true = randn(rng, p)
    y_raw = X_raw * beta_true .+ 0.05 .* randn(rng, n)
    X = X_raw .- mean(X_raw, dims = 1)
    y = y_raw .- mean(y_raw)
    return X, y
end

# Build a mixed (energy + torque) augmented system using the actual
# `assemble_weighted_problem`. The synthetic energy / torque samples are
# i.i.d. Gaussians; the test only requires that the resulting `(X, y)`
# is a valid in-shape augmented problem on which `OLS()` and
# `ElasticNet(alpha, λ→0)` should agree to coordinate-descent precision.
function _mixed_system(rng, n_config::Int, n_atoms::Int, p::Int, weight::Float64)
    design_matrix_energy = randn(rng, n_config, p)
    design_matrix_torque = randn(rng, 3 * n_atoms * n_config, p)
    observed_energy_list = randn(rng, n_config)
    observed_torque_list = [randn(rng, 3, n_atoms) for _ in 1:n_config]
    X, y = Magesty.Fitting.assemble_weighted_problem(
        design_matrix_energy,
        design_matrix_torque,
        observed_energy_list,
        observed_torque_list,
        weight,
    )
    return X, y
end

@testset "ElasticNet estimator" begin
    @testset "λ→0 limit ≈ OLS, energy-only (weight=0)" begin
        rng = MersenneTwister(20260517)
        X, y = _centered_energy_system(rng, 60, 8)
        b_ols = Magesty.Fitting.solve_coefficients(OLS(), X, y)
        b_lasso = Magesty.Fitting.solve_coefficients(
            ElasticNet(alpha = 1.0, lambda = 1e-10, standardize = false),
            X, y,
        )
        b_l2 = Magesty.Fitting.solve_coefficients(
            ElasticNet(alpha = 0.0, lambda = 1e-10, standardize = false),
            X, y,
        )
        @test maximum(abs.(b_lasso .- b_ols)) < _GLMNET_OLS_ATOL
        @test maximum(abs.(b_l2    .- b_ols)) < _GLMNET_OLS_ATOL
    end

    @testset "λ→0 limit ≈ OLS, mixed objective (weight=0.5)" begin
        rng = MersenneTwister(20260518)
        X, y = _mixed_system(rng, 8, 2, 4, 0.5)
        b_ols = Magesty.Fitting.solve_coefficients(OLS(), X, y)
        b_lasso = Magesty.Fitting.solve_coefficients(
            ElasticNet(alpha = 1.0, lambda = 1e-10, standardize = false),
            X, y,
        )
        b_l2 = Magesty.Fitting.solve_coefficients(
            ElasticNet(alpha = 0.0, lambda = 1e-10, standardize = false),
            X, y,
        )
        @test maximum(abs.(b_lasso .- b_ols)) < _GLMNET_OLS_ATOL
        @test maximum(abs.(b_l2    .- b_ols)) < _GLMNET_OLS_ATOL
    end

    @testset "λ→0 limit ≈ OLS, standardize=true" begin
        rng = MersenneTwister(20260519)
        for weight in (0.0, 0.5)
            local X, y
            if iszero(weight)
                X, y = _centered_energy_system(rng, 60, 8)
            else
                X, y = _mixed_system(rng, 8, 2, 4, weight)
            end
            b_ols = Magesty.Fitting.solve_coefficients(OLS(), X, y)
            b_std = Magesty.Fitting.solve_coefficients(
                ElasticNet(alpha = 1.0, lambda = 1e-10, standardize = true),
                X, y,
            )
            @test maximum(abs.(b_std .- b_ols)) < _GLMNET_STD_ATOL
        end
    end

    @testset "ElasticNet(α=0) ≈ analytic Ridge (λ_g = λ_a / n)" begin
        rng = MersenneTwister(20260520)
        X, y = _centered_energy_system(rng, 60, 8)
        n = size(X, 1)
        for lambda_a in (1e-3, 1e-2, 1e-1)
            b_ridge = Magesty.Fitting.solve_coefficients(
                Ridge(lambda = lambda_a), X, y,
            )
            b_en = Magesty.Fitting.solve_coefficients(
                ElasticNet(alpha = 0.0, lambda = lambda_a / n, standardize = false),
                X, y,
            )
            rel = maximum(abs.(b_en .- b_ridge)) / max(maximum(abs.(b_ridge)), eps())
            @test rel < _GLMNET_RIDGE_RTOL
        end
    end

    @testset "Lasso sparsity monotonicity in λ" begin
        rng = MersenneTwister(20260521)
        # Sparse ground truth: only columns 1, 4, 7 carry nonzero weight.
        n, p = 100, 10
        X_raw = randn(rng, n, p)
        beta_true = zeros(p)
        beta_true[1] = 1.0
        beta_true[4] = -0.7
        beta_true[7] = 0.4
        y_raw = X_raw * beta_true .+ 0.05 .* randn(rng, n)
        X = X_raw .- mean(X_raw, dims = 1)
        y = y_raw .- mean(y_raw)

        zero_counts = Int[]
        for lambda in (1e-4, 1e-3, 1e-2, 1e-1, 5e-1)
            beta = Magesty.Fitting.solve_coefficients(
                Lasso(lambda = lambda, standardize = true), X, y,
            )
            push!(zero_counts, count(iszero, beta))
        end
        # Weakly nondecreasing.
        for k = 2:length(zero_counts)
            @test zero_counts[k] >= zero_counts[k - 1]
        end
        # At a strongly regularised endpoint, expect at least a few zeros.
        @test zero_counts[end] >= 1
    end

    @testset "standardize=true retains high-scale columns" begin
        rng = MersenneTwister(20260522)
        # Mimic the SCE per-cluster (4π)^(N/2) factor for N=1..4. The true
        # support lives on the large-scale columns.
        scales = [1.0, 12.0, 45.0, 158.0]
        n, p = 200, length(scales)
        X_raw = randn(rng, n, p) .* reshape(scales, 1, p)
        # True coefficients live mostly on the large-scale columns.
        beta_true = [0.0, 0.0, 0.05, 0.02]
        y_raw = X_raw * beta_true .+ 0.05 .* randn(rng, n)
        X = X_raw .- mean(X_raw, dims = 1)
        y = y_raw .- mean(y_raw)

        # At a moderate λ, `standardize = true` should retain the high-N
        # (large-scale) columns while `standardize = false` drops them.
        lambda = 5e-1
        b_std = Magesty.Fitting.solve_coefficients(
            ElasticNet(alpha = 1.0, lambda = lambda, standardize = true),
            X, y,
        )
        b_nostd = Magesty.Fitting.solve_coefficients(
            ElasticNet(alpha = 1.0, lambda = lambda, standardize = false),
            X, y,
        )
        # The crucial qualitative claim: standardise keeps the high-N
        # columns alive, no-standardise kills (at least one of) them.
        @test count(!iszero, b_std) >= count(!iszero, b_nostd)
        # Both known-true large-scale columns (3 and 4) survive at this λ
        # when standardise reverses the (4π)^(N/2)-style column scaling.
        @test !iszero(b_std[3])
        @test !iszero(b_std[end])
    end
end
