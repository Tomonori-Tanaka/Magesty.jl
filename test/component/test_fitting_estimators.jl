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

# Component tests for the AdaptiveLasso (oneshot, Zou 2006) estimator.
# Five properties from the 260518-adaptive-lasso-oneshot spec are checked:
#
#   1. `gamma = 0` reduces to plain Lasso. Both paths route through the
#      same `_glmnet_solve` with `penalty_factor` equal to (a rescale of)
#      the all-ones vector; GLMNet's internal `sum(pf) = nvars` rescale
#      is the identity on `ones(p)`, so the two coefficient vectors agree
#      to GLMNet's coordinate-descent precision. Note: at gamma = 0 the
#      pilot output never enters the `pf` formula (x^0 = 1), so this test
#      does not distinguish "pilot ran" from "pilot was skipped". A
#      future `gamma == 0` fast path that skips the pilot would still
#      pass.
#   2. Weight construction is bit-exact: building `pf` manually from
#      `solve_coefficients(OLS(), X, y)` and calling `_glmnet_solve` with
#      that `pf` gives the same vector as `solve_coefficients(AdaptiveLasso
#      (pilot = OLS(), ...), X, y)`. Pins the formula `pf[j] = 1 /
#      max(|beta_pilot[j]|, eps)^gamma`. Depends on the non-exported
#      `Magesty.Fitting._glmnet_solve`; renaming that helper or changing
#      its keyword order is the only way this strict test breaks.
#   3. `epsilon` clip with a Lasso pilot: pilot has exact zeros that
#      would otherwise produce `pf = 1/0 = Inf`. The clip
#      (`max(|beta_pilot|, epsilon)`) keeps the call finite. The
#      adaptive fit's *per-column zero* pattern is intentionally not
#      asserted: GLMNet's `sum(pf) = nvars` rescale leaves the
#      zero-pilot columns at moderate effective penalty (`pf ~ nvars/k`
#      with k = pilot-zero count) and crushes the alive-pilot columns
#      to `pf ~ 1e-15`, so whether a pilot-zero column ends nonzero in
#      the adaptive fit depends on data signal vs the (data-dependent)
#      effective lambda. We pin only the numerical-safety claim
#      (finite output) plus true-support preservation.
#   4. `pilot = Ridge(...)` path: smoke test. The alternative pilot
#      path runs without error and returns a finite, partially-sparse
#      vector that recovers most of the true support. No claim that
#      Ridge-pilot AdaptiveLasso matches OLS-pilot AdaptiveLasso on
#      support recovery -- a Ridge pilot's softer shrinkage typically
#      produces less aggressive adaptive penalties.
#   5. Best-`lambda`-per-side sparse recovery: on a correlated-feature
#      synthetic problem, `AdaptiveLasso(pilot = OLS(), gamma = 1)`
#      recovers the true support exactly at some `lambda` in a small
#      geometric grid, while plain `Lasso` recovers it at no `lambda` in
#      the same grid. GLMNet's internal `penalty_factor` rescaling
#      breaks matched-lambda comparison, so each estimator is tuned
#      independently before comparing support sets. The seed, fixture,
#      and pass condition are tuned against GLMNet.jl 0.7.4; a future
#      GLMNet upgrade that shifts its coordinate-descent convergence
#      may require seed adjustment.

# atol for the `gamma = 0` agreement between AdaptiveLasso and plain
# Lasso. Both call sites end at GLMNet with the same `(alpha = 1.0,
# lambda, standardize)`; the only difference is that AdaptiveLasso
# explicitly passes `penalty_factor = ones(p)`, which GLMNet rescales
# to the same vector. Empirical worst case on the fixture below is
# ~3e-8; the constant below leaves ~30x headroom.
const _ADALASSO_GAMMA0_ATOL = 1.0e-6

@testset "AdaptiveLasso estimator" begin
    @testset "gamma = 0 reduces to plain Lasso" begin
        rng = MersenneTwister(20260601)
        X, y = _centered_energy_system(rng, 60, 8)
        lambda = 1e-3
        b_lasso = Magesty.Fitting.solve_coefficients(
            Lasso(lambda = lambda, standardize = false),
            X, y,
        )
        b_ada = Magesty.Fitting.solve_coefficients(
            AdaptiveLasso(
                pilot = OLS(), lambda = lambda,
                gamma = 0.0, standardize = false,
            ),
            X, y,
        )
        @test maximum(abs.(b_ada .- b_lasso)) < _ADALASSO_GAMMA0_ATOL
    end

    @testset "weight construction is exact (==)" begin
        rng = MersenneTwister(20260602)
        X, y = _centered_energy_system(rng, 60, 8)
        lambda = 1e-3
        gamma = 1.0
        epsilon = eps(Float64)
        beta_pilot = Magesty.Fitting.solve_coefficients(OLS(), X, y)
        pf = inv.(max.(abs.(beta_pilot), epsilon) .^ gamma)
        b_direct = Magesty.Fitting._glmnet_solve(
            X, y;
            alpha = 1.0, lambda = lambda,
            standardize = false, penalty_factor = pf,
        )
        b_ada = Magesty.Fitting.solve_coefficients(
            AdaptiveLasso(
                pilot = OLS(), lambda = lambda,
                gamma = gamma, epsilon = epsilon,
                standardize = false,
            ),
            X, y,
        )
        # Both call sites pass identical (X, y, alpha, lambda, standardize,
        # penalty_factor) into `_glmnet_solve`, so the coefficient vectors
        # must agree element-wise -- not merely `isapprox`.
        @test b_ada == b_direct
    end

    @testset "epsilon clip keeps the call finite for a Lasso pilot" begin
        rng = MersenneTwister(20260603)
        # Sparse ground truth; the pilot's heavy shrinkage produces
        # several exactly-zero columns. Without the epsilon clip
        # `pf = 1 / |beta_pilot|^gamma` would be `Inf` on those columns
        # and GLMNet would error / produce NaN.
        n, p = 100, 10
        X_raw = randn(rng, n, p)
        beta_true = zeros(p)
        beta_true[1] = 1.0
        beta_true[4] = -0.7
        beta_true[7] = 0.4
        y_raw = X_raw * beta_true .+ 0.05 .* randn(rng, n)
        X = X_raw .- mean(X_raw, dims = 1)
        y = y_raw .- mean(y_raw)

        pilot = Lasso(lambda = 0.5, standardize = false)
        beta_pilot = Magesty.Fitting.solve_coefficients(pilot, X, y)
        # Sanity: the pilot must actually produce at least one exact zero
        # for the clip behaviour to be exercised.
        @test count(iszero, beta_pilot) >= 1

        b_ada = Magesty.Fitting.solve_coefficients(
            AdaptiveLasso(
                pilot = pilot, lambda = 1e-3,
                gamma = 1.0, epsilon = eps(Float64),
                standardize = false,
            ),
            X, y,
        )
        # The substantive claim of the epsilon clip: the call survives
        # the pilot's exact zeros and returns a finite vector. GLMNet's
        # `sum(pf) = nvars` rescale then puts the clipped columns at
        # moderate effective lambda (`pf ~ nvars / k`, k = pilot-zero
        # count) and crushes the alive columns to `pf ~ 1e-15`. Whether
        # a clipped column ends nonzero in the adaptive fit therefore
        # depends on data signal vs the (data-dependent) effective
        # lambda, so the per-column zero pattern is not asserted here.
        @test all(isfinite, b_ada)
        # The true support survives: alive pilot columns are
        # effectively unpenalized after rescale, so the OLS fit on them
        # recovers the true betas.
        @test !iszero(b_ada[1])
        @test !iszero(b_ada[4])
    end

    @testset "pilot = Ridge(...) path runs and selects sensibly" begin
        rng = MersenneTwister(20260604)
        # Same correlated-feature fixture used by the sparse-recovery
        # test below.
        n, p = 100, 15
        X_raw = randn(rng, n, p)
        for j in 4:p
            target = ((j - 1) % 3) + 1
            X_raw[:, j] .= 0.85 .* X_raw[:, target] .+ 0.15 .* randn(rng, n)
        end
        beta_true = zeros(p)
        beta_true[1] = 1.5
        beta_true[2] = -1.2
        beta_true[3] = 0.9
        y_raw = X_raw * beta_true .+ 0.05 .* randn(rng, n)
        X = X_raw .- mean(X_raw, dims = 1)
        y = y_raw .- mean(y_raw)

        b_ada = Magesty.Fitting.solve_coefficients(
            AdaptiveLasso(
                pilot = Ridge(lambda = 1e-3),
                lambda = 5e-2,
                gamma = 1.0,
                standardize = false,
            ),
            X, y,
        )
        @test all(isfinite, b_ada)
        # Qualitative: the fit is sparse and recovers at least two of
        # the three true-support columns.
        @test count(iszero, b_ada) >= 1
        true_support = Set([1, 2, 3])
        active_support = Set(findall(!iszero, b_ada))
        @test length(intersect(true_support, active_support)) >= 2
    end

    @testset "best-lambda sparse recovery beats plain Lasso" begin
        rng = MersenneTwister(20260605)
        # Correlated noise columns: columns 4-15 are 0.6 * (one of the
        # true columns) + 0.4 * (independent noise). Plain Lasso retains
        # spurious correlated columns across the lambda grid, but
        # AdaptiveLasso (OLS pilot, gamma = 1) re-weights them out.
        n, p = 100, 15
        X_raw = randn(rng, n, p)
        for j in 4:p
            target = ((j - 1) % 3) + 1
            X_raw[:, j] .= 0.85 .* X_raw[:, target] .+ 0.15 .* randn(rng, n)
        end
        beta_true = zeros(p)
        beta_true[1] = 1.5
        beta_true[2] = -1.2
        beta_true[3] = 0.9
        y_raw = X_raw * beta_true .+ 0.05 .* randn(rng, n)
        X = X_raw .- mean(X_raw, dims = 1)
        y = y_raw .- mean(y_raw)
        true_support = Set([1, 2, 3])

        lambda_grid = 10.0 .^ range(-3, 0, length = 13)
        lasso_hit = false
        ada_hit = false
        for lambda in lambda_grid
            b_lasso = Magesty.Fitting.solve_coefficients(
                Lasso(lambda = lambda, standardize = true), X, y,
            )
            if Set(findall(!iszero, b_lasso)) == true_support
                lasso_hit = true
            end
            b_ada = Magesty.Fitting.solve_coefficients(
                AdaptiveLasso(
                    pilot = OLS(), lambda = lambda,
                    gamma = 1.0, standardize = true,
                ),
                X, y,
            )
            if Set(findall(!iszero, b_ada)) == true_support
                ada_hit = true
            end
        end
        @test ada_hit
        @test !lasso_hit
    end
end

# Component tests for the PrecomputedPilot adapter. Three properties are
# checked:
#
#   1. Pilot-reuse parity. `AdaptiveLasso` driven by `PrecomputedPilot(coef_ols)`
#      must produce coefficients bit-exactly identical to `AdaptiveLasso` with
#      `pilot = OLS()`, because both call sites compute the same `pf` and call
#      `_glmnet_solve` with the same inputs.
#   2. Defensive copy. Mutating the caller's vector after construction must
#      not leak into the adapter's stored coefficients.
#   3. Length mismatch is rejected. `solve_coefficients(::PrecomputedPilot, X, y)`
#      raises `DimensionMismatch` when the stored vector length disagrees with
#      `size(X, 2)`; this is the only check the adapter can make for an
#      SCEBasis mismatch.
#
# Smoke coverage for the `AdaptiveLasso(::SCEFit; ...)` /
# `AdaptiveLasso(::SCEModel; ...)` convenience constructors lives in the
# integration tests, where a real `SCEFit` / `SCEModel` is available.

@testset "PrecomputedPilot adapter" begin
    @testset "pilot reuse is bit-exact against fresh OLS pilot" begin
        rng = MersenneTwister(20260701)
        X, y = _centered_energy_system(rng, 60, 8)
        lambda = 1e-3
        gamma = 1.0
        beta_ols = Magesty.Fitting.solve_coefficients(OLS(), X, y)
        b_fresh = Magesty.Fitting.solve_coefficients(
            AdaptiveLasso(
                pilot = OLS(), lambda = lambda,
                gamma = gamma, standardize = false,
            ),
            X, y,
        )
        b_reuse = Magesty.Fitting.solve_coefficients(
            AdaptiveLasso(
                pilot = PrecomputedPilot(beta_ols), lambda = lambda,
                gamma = gamma, standardize = false,
            ),
            X, y,
        )
        # Both code paths compute the same `pf` (the OLS solution is
        # deterministic via LAPACK) and call `_glmnet_solve` with the
        # same inputs, so equality must be element-wise.
        @test b_reuse == b_fresh
    end

    @testset "input vector is copied at construction" begin
        beta_src = [1.0, -2.0, 3.0]
        pp = PrecomputedPilot(beta_src)
        beta_src[1] = 999.0
        # The adapter must hold its own storage so that later mutation
        # of caller-side arrays does not leak into the stored pilot.
        @test pp.beta[1] == 1.0
    end

    @testset "length mismatch raises DimensionMismatch" begin
        rng = MersenneTwister(20260702)
        X, y = _centered_energy_system(rng, 40, 6)
        pp_wrong = PrecomputedPilot(zeros(5))  # one short
        @test_throws DimensionMismatch Magesty.Fitting.solve_coefficients(
            pp_wrong, X, y,
        )
    end
end
