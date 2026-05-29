using Test
using Statistics
using LinearAlgebra: norm
using Magesty
import TOML

# Component tests for `refit` -- post-selection refit on the basis
# support of an `SCEFit`.
#
# Properties checked:
#
#   1. OLS idempotence. `refit(fit_ols, OLS(); threshold = 0)` reproduces
#      `coef(fit_ols)` to round-off. With OLS there are no exact-zero
#      coefficients on the test fixture, so the support is the full set
#      and the refit objective is identical to the original.
#
#   2. Scaled-magnitude support selection on a manually zeroed `SCEFit`.
#      Constructing an `SCEFit` whose `jphi` has been masked to a
#      caller-chosen index set, refitting with OLS reproduces the
#      independent OLS solve restricted to those columns of the
#      assembled `(X, y)`. This is the load-bearing correctness check;
#      no Lasso fixture is required.
#
#   3. L1 fit at the default `threshold = 0.0`. A `Lasso` fit's
#      exact-zero coefficients fall out at the default threshold; refit
#      support equals `findall(!iszero, coef(lasso_fit))`, with dropped
#      basis coefficients exactly `0.0`.
#
#   4. Empty support. A `threshold` above every scaled magnitude warns,
#      returns an all-zero `jphi`, and yields `j0 == mean(y_E)`.
#
#   5. `torque_weight` reuse. `refit` of a fit made with
#      `torque_weight = w` carries `w` through to the new fit.
#
#   6. `PrecomputedPilot`-backed estimator rejection (both the
#      `PrecomputedPilot` itself and an `AdaptiveLasso` whose pilot is
#      one, including the `AdaptiveLasso(::SCEFit)` form).
#
#   7. Negative `threshold` rejection.

const DIMER_TOML_REFIT = joinpath(@__DIR__, "..", "integration", "dimer", "input.toml")

# Build a non-isotropic dimer basis with multiple SCE coefficients so the
# refit tests exercise partial-support behavior. The shipped dimer TOML
# has `isotropy = true` plus a single `(l1, l2, Lf)` triple, which yields
# `num_salcs = 1` -- too small for any non-trivial support test.
function _refit_basis()
    input = TOML.parsefile(DIMER_TOML_REFIT)
    input["symmetry"]["isotropy"] = false
    input["interaction"]["body1"]["lmax"]["X"] = 1
    input["interaction"]["body2"]["lsum"] = 4
    return SCEBasis(input; verbosity = false)
end

# Spin configurations on the same 2-atom dimer. The set is chosen with
# enough orientational variety that the augmented design matrix is full
# column rank (rank == num_salcs); the four high-symmetry configs above
# alone leave rank deficient by ~3 against the (lmax=1, lsum=4,
# isotropy=false) basis, so this set adds five generic off-axis
# orientations to break the residual symmetry-induced dependencies.
function _refit_configs()
    SC = Magesty.SpinConfigs.SpinConfig
    # Helper: column-normalize a 3xN spin-direction matrix to unit length.
    _col_normalize(v) = v ./ sqrt.(sum(abs2, v; dims = 1))
    cfgs = [
        SC(-1.0, [1.0, 1.0],
           [0.0 0.0; 0.0 0.0; 1.0 1.0],
           [0.1 0.0; 0.0 -0.1; 0.0 0.0]),
        SC( 1.0, [1.0, 1.0],
           [0.0 0.0; 0.0 0.0; -1.0 1.0],
           [-0.1 0.0; 0.0 0.1; 0.0 0.0]),
        SC( 0.3, [1.0, 1.0],
           [1.0 0.0; 0.0 1.0; 0.0 0.0],
           [0.0 0.2; -0.2 0.0; 0.0 0.0]),
        SC(-0.2, [1.0, 1.0],
           [1.0/sqrt(3) 1.0/sqrt(3); 1.0/sqrt(3) -1.0/sqrt(3); 1.0/sqrt(3) 1.0/sqrt(3)],
           [0.05 -0.05; -0.05 0.05; 0.0 0.0]),
        SC( 0.7, [1.0, 1.0],
           _col_normalize([0.2 0.6; 0.7 -0.3; 0.5 0.4]),
           [0.03 -0.01; -0.02 0.04; 0.01 -0.03]),
        SC(-0.4, [1.0, 1.0],
           _col_normalize([0.5 -0.4; -0.2 0.8; 0.6 -0.3]),
           [-0.02 0.03; 0.04 -0.01; -0.03 0.02]),
        SC( 0.1, [1.0, 1.0],
           _col_normalize([0.3 0.5; -0.4 0.2; 0.7 -0.6]),
           [0.04 -0.02; -0.01 0.03; 0.02 -0.04]),
        SC(-0.6, [1.0, 1.0],
           _col_normalize([-0.5 0.7; 0.6 -0.2; -0.3 0.4]),
           [0.01 -0.04; 0.03 0.02; -0.02 0.01]),
        SC( 0.5, [1.0, 1.0],
           _col_normalize([0.4 -0.3; -0.5 0.6; 0.2 -0.7]),
           [-0.03 0.01; 0.02 -0.04; 0.04 -0.02]),
        SC(-0.1, [1.0, 1.0],
           _col_normalize([0.6 -0.5; 0.3 0.4; -0.5 0.7]),
           [0.02 -0.03; -0.04 0.01; 0.03 -0.02]),
        SC( 0.4, [1.0, 1.0],
           _col_normalize([-0.2 0.4; -0.7 -0.5; 0.6 0.3]),
           [-0.01 0.04; 0.03 -0.02; -0.02 0.01]),
    ]
    return cfgs
end

@testset "refit" begin
    basis = _refit_basis()
    configs = _refit_configs()
    dataset = SCEDataset(basis, configs)
    num_salcs = length(basis.salcbasis.salc_list)
    @assert num_salcs >= 2 "refit tests require >= 2 SCE basis functions"

    # Reference weighted problem -- the same call `fit` (and therefore
    # `refit`) makes when `torque_weight = 0.3`.
    torque_weight = 0.3
    X, y = Magesty.Fitting.assemble_weighted_problem(
        dataset.X_E,
        dataset.X_T,
        dataset.y_E,
        dataset.y_T,
        torque_weight,
    )

    @testset "OLS idempotence at threshold = 0" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = torque_weight, verbosity = false)
        # OLS produces no exact-zero coefficients here, so the support
        # is the full set and `refit` should reproduce the analytic OLS
        # solve. Compare against the independent `X \ y` (not `coef(f)`)
        # so a latent bug in either `fit` or `refit` is caught.
        @assert all(!iszero, coef(f)) "OLS fit produced a zero coefficient; pick a different fixture"
        jphi_expected = X \ y
        j0_expected = mean(dataset.y_E .- dataset.X_E * jphi_expected)
        rf = refit(f, OLS(); verbosity = false)
        @test rf isa SCEFit
        @test rf.torque_weight == f.torque_weight
        @test rf.estimator isa OLS
        @test length(coef(rf)) == length(coef(f))
        @test coef(rf) ≈ jphi_expected atol = 1e-10
        @test intercept(rf) ≈ j0_expected atol = 1e-10
        # The two paths also agree with each other.
        @test coef(rf) ≈ coef(f) atol = 1e-10
        @test intercept(rf) ≈ intercept(f) atol = 1e-10
        @test rf.dataset === f.dataset
        @test length(rf.residuals) == size(dataset.X_E, 1) + size(dataset.X_T, 1)
    end

    @testset "support equals an independent OLS solve on the selected columns" begin
        f_full = fit(SCEFit, dataset, OLS(); torque_weight = torque_weight, verbosity = false)
        # Mask `jphi` to a caller-chosen support; at `threshold = 0`,
        # `refit` should pick exactly that set.
        keep = [1, num_salcs]
        @assert length(unique(keep)) == length(keep) "keep indices must be distinct"
        jphi_masked = zeros(Float64, num_salcs)
        jphi_masked[keep] .= coef(f_full)[keep]
        f_masked = Magesty.SCEFit(
            dataset,
            f_full.j0,
            jphi_masked,
            OLS(),
            f_full.torque_weight,
            f_full.residuals,
        )

        rf = refit(f_masked, OLS(); verbosity = false)

        # Analytic expected jphi: OLS on `X[:, keep]` (same weighted
        # problem the refit assembles), then `extract_j0_jphi` recovers
        # `j0 = mean(y_E - X_E * jphi)`.
        jphi_expected = zeros(Float64, num_salcs)
        jphi_expected[keep] .= X[:, keep] \ y
        j0_expected = mean(dataset.y_E .- dataset.X_E * jphi_expected)

        @test coef(rf) ≈ jphi_expected atol = 1e-10
        @test intercept(rf) ≈ j0_expected atol = 1e-10
        for j = 1:num_salcs
            j in keep && continue
            @test coef(rf)[j] == 0.0
        end
    end

    @testset "L1 exact-zero support at threshold = 0.0" begin
        # Sweep lambda until the Lasso drops at least one coefficient
        # but keeps at least one, so this exercises the "exact-zero is
        # dropped" path rather than degrading to OLS idempotence.
        f_lasso = nothing
        for lambda in (1e-3, 1e-2, 1e-1, 0.5, 1.0, 5.0)
            cand = fit(
                SCEFit,
                dataset,
                Lasso(lambda = lambda, standardize = false);
                torque_weight = torque_weight,
                verbosity = false,
            )
            if any(iszero, coef(cand)) && any(!iszero, coef(cand))
                f_lasso = cand
                break
            end
        end
        if f_lasso === nothing
            @info "skipping L1 exact-zero subtest: Lasso did not produce a partial support on this fixture"
        else
            expected_support = findall(!iszero, coef(f_lasso))
            rf = refit(f_lasso, OLS(); verbosity = false)
            @test findall(!iszero, coef(rf)) == expected_support
            jphi_expected = zeros(Float64, length(coef(rf)))
            jphi_expected[expected_support] .= X[:, expected_support] \ y
            @test coef(rf) ≈ jphi_expected atol = 1e-9
        end
    end

    @testset "AdaptiveRidge thresholded support" begin
        # `AdaptiveRidge` produces continuously small but nonzero
        # coefficients (no exact zeros), so the default threshold = 0
        # is a no-op and the caller must pass a positive threshold to
        # drop near-negligible bases.
        f_ar = fit(
            SCEFit, dataset, AdaptiveRidge(lambda = 1.0);
            torque_weight = torque_weight, verbosity = false,
        )
        # threshold = 0 should keep everything.
        rf_noop = refit(f_ar, OLS(); threshold = 0.0, verbosity = false)
        @test count(!iszero, coef(rf_noop)) == num_salcs

        # Pick a threshold partway through the scaled-magnitude
        # distribution so the support strictly shrinks but is non-empty.
        scaled = [abs(coef(f_ar)[j]) * norm(@view X[:, j]) for j in eachindex(coef(f_ar))]
        sorted = sort(scaled)
        cutoff = (sorted[end - 1] + sorted[end]) / 2
        rf = refit(f_ar, OLS(); threshold = cutoff, verbosity = false)

        expected_support = findall(>(cutoff), scaled)
        @test 0 < length(expected_support) < num_salcs
        @test findall(!iszero, coef(rf)) == expected_support
        # Dropped bases stay at exactly 0.0.
        for j = 1:num_salcs
            j in expected_support && continue
            @test coef(rf)[j] == 0.0
        end
        # Refit support equals an independent OLS solve on the
        # selected columns.
        jphi_expected = zeros(Float64, num_salcs)
        jphi_expected[expected_support] .= X[:, expected_support] \ y
        @test coef(rf) ≈ jphi_expected atol = 1e-9
    end

    @testset "empty support warns and returns all-zero jphi" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = torque_weight, verbosity = false)
        max_scaled = maximum(
            j -> abs(coef(f)[j]) * norm(@view X[:, j]),
            eachindex(coef(f)),
        )
        rf = @test_logs (:warn, r"empty support") refit(
            f, OLS();
            threshold = max_scaled + 1.0,
            verbosity = false,
        )
        @test all(iszero, coef(rf))
        @test intercept(rf) ≈ mean(dataset.y_E)
    end

    @testset "torque_weight is reused" begin
        for w in (0.0, 0.25, 0.5, 0.75, 1.0)
            f = fit(SCEFit, dataset, OLS(); torque_weight = w, verbosity = false)
            rf = refit(f, OLS(); verbosity = false)
            @test rf.torque_weight == w
        end
    end

    @testset "rejects PrecomputedPilot-backed estimators" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = torque_weight, verbosity = false)
        # Bare `PrecomputedPilot`: fixed vector has the original column
        # count, not the support length.
        @test_throws ArgumentError refit(f, PrecomputedPilot(coef(f)); verbosity = false)
        # `AdaptiveLasso(::SCEFit; ...)` builds an AdaptiveLasso whose
        # pilot is a `PrecomputedPilot`; reject that too.
        est_from_fit = AdaptiveLasso(f; lambda = 1e-3)
        @test_throws ArgumentError refit(f, est_from_fit; verbosity = false)
    end

    @testset "rejects negative threshold" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = torque_weight, verbosity = false)
        @test_throws ArgumentError refit(f, OLS(); threshold = -1e-12, verbosity = false)
    end
end
