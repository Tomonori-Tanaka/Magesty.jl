using Test
using Magesty
using LinearAlgebra
import Random

# GCV diagnostics. The core-numerics tests are analytic cross-checks: the
# SVD-based effective dof and lambda-path residuals are validated against an
# independent dense recomputation of the same quantity (matrix inverse / an
# explicit-intercept hat matrix), never against the implementation's own output.

const F = Magesty.Fitting
const GCV_FIXTURE = joinpath(@__DIR__, "..", "integration", "fept_tetragonal_2x2x2")

# Non-comment, non-empty lines of a written table.
_gcv_data_lines(text::AbstractString) =
    [ln for ln in split(text, '\n')
     if !isempty(strip(ln)) && !startswith(strip(ln), "#")]

@testset "GCV diagnostics" begin

    @testset "core numerics (synthetic, analytic cross-checks)" begin
        Random.seed!(20260610)
        n_samples, p = 60, 6
        design_raw = randn(n_samples, p)
        # Nonzero intercept so the +1 dof term is genuinely exercised.
        observed = design_raw * randn(p) .+ 0.5 .+ 0.1 .* randn(n_samples)
        # Mirror the fit pipeline: energy-block centering eliminates j0.
        design_centered = design_raw .- sum(design_raw; dims = 1) ./ n_samples
        observed_centered = observed .- sum(observed) / n_samples

        @testset "OLS GCV equals explicit-intercept hat matrix" begin
            # Independent reference: OLS with an explicit constant column.
            design_aug = hcat(ones(n_samples), design_raw)
            beta_aug = design_aug \ observed
            resid_aug = observed - design_aug * beta_aug
            hat_aug = design_aug * (Symmetric(design_aug' * design_aug) \ Matrix(design_aug'))
            dof_ref = tr(hat_aug)                       # = p + 1 at full rank
            gcv_ref = (sum(abs2, resid_aug) / n_samples) / (1 - dof_ref / n_samples)^2

            jphi = design_centered \ observed_centered
            resid = observed_centered - design_centered * jphi
            # All rows live (no zero-weight block), intercept present.
            score, dof = F._gcv_single(
                design_centered, observed_centered, resid, F.OLS(), jphi, n_samples, 1)
            @test dof ≈ dof_ref atol = 1e-9
            @test dof ≈ p + 1 atol = 1e-9
            @test score ≈ gcv_ref rtol = 1e-9
        end

        @testset "Ridge dof and GCV match dense recomputation" begin
            for lam in (1e-3, 0.1, 1.0, 5.0)
                jphi = (design_centered' * design_centered + lam * I) \
                       (design_centered' * observed_centered)
                resid = observed_centered - design_centered * jphi
                dof_ref = 1 + tr(design_centered *
                    ((design_centered' * design_centered + lam * I) \ design_centered'))
                gcv_ref = (sum(abs2, resid) / n_samples) / (1 - dof_ref / n_samples)^2
                score, dof = F._gcv_single(
                    design_centered, observed_centered, resid,
                    F.Ridge(lambda = lam), jphi, n_samples, 1)
                @test dof ≈ dof_ref rtol = 1e-10
                @test score ≈ gcv_ref rtol = 1e-10
            end
        end

        @testset "lambda path matches per-lambda re-solve" begin
            lambdas = [1e-4, 1e-2, 1.0, 10.0]
            gcvs, dofs = F._gcv_lambda_path(
                design_centered, observed_centered, lambdas, n_samples, 1)
            @test length(gcvs) == length(lambdas)
            @test length(dofs) == length(lambdas)
            for (k, lam) in enumerate(lambdas)
                jphi = (design_centered' * design_centered + lam * I) \
                       (design_centered' * observed_centered)
                resid = observed_centered - design_centered * jphi
                dof_ref = 1 + tr(design_centered *
                    ((design_centered' * design_centered + lam * I) \ design_centered'))
                gcv_ref = (sum(abs2, resid) / n_samples) / (1 - dof_ref / n_samples)^2
                @test dofs[k] ≈ dof_ref rtol = 1e-9
                @test gcvs[k] ≈ gcv_ref rtol = 1e-8
            end
        end

        @testset "AdaptiveRidge conditional dof matches dense" begin
            est = F.AdaptiveRidge(lambda = 0.3, epsilon = 1e-8)
            beta = (design_centered' * design_centered + est.lambda * I) \
                   (design_centered' * observed_centered)
            w = inv.(abs2.(beta) .+ est.epsilon)
            dof_ref = 1 + tr(
                (design_centered' * design_centered + est.lambda * Diagonal(w)) \
                (design_centered' * design_centered))
            @test F._effective_dof(design_centered, est, beta, 1) ≈ dof_ref rtol = 1e-9
        end

        @testset "_gcv_value saturation -> NaN" begin
            @test isnan(F._gcv_value(1.0, 10, 10))      # dof == N
            @test isnan(F._gcv_value(1.0, 11, 10))      # dof > N
            @test F._gcv_value(2.0, 5, 10) ≈ (2.0 / 10) / (1 - 5 / 10)^2
        end

        @testset "_gcv_sample_count excludes the zero-weight block" begin
            n_E, n_T = 30, 720
            # 0 < w < 1: both blocks live, intercept present.
            @test F._gcv_sample_count(n_E, n_T, 0.5) == (n_E + n_T, 1)
            # w == 1: energy block zeroed -> only torque rows, no intercept.
            @test F._gcv_sample_count(n_E, n_T, 1.0) == (n_T, 0)
            # w == 0: torque block zeroed -> only energy rows, intercept present.
            @test F._gcv_sample_count(n_E, n_T, 0.0) == (n_E, 1)
        end

        @testset "estimator gating" begin
            @test F._is_linear_estimator(OLS())
            @test F._is_linear_estimator(Ridge(lambda = 1.0))
            @test F._is_linear_estimator(AdaptiveRidge(lambda = 1.0))
            @test !F._is_linear_estimator(Lasso(lambda = 1.0))
            @test !F._is_linear_estimator(ElasticNet(alpha = 0.5, lambda = 1.0))
            @test !F._is_linear_estimator(AdaptiveLasso(lambda = 1.0))
            @test F._require_linear_estimator(OLS()) === nothing
            @test_throws ArgumentError F._require_linear_estimator(Lasso(lambda = 1.0))
        end
    end

    @testset "public API (FePt fixture)" begin
        basis = SCEBasis(joinpath(GCV_FIXTURE, "input.toml"); verbosity = false)
        dataset = SCEDataset(basis, joinpath(GCV_FIXTURE, "EMBSET"); verbosity = false)
        n_total = length(dataset)

        @testset "gcv(f) for linear estimators" begin
            for est in (OLS(), Ridge(lambda = 1e-3), AdaptiveRidge(lambda = 1e-3))
                f = fit(SCEFit, dataset, est; verbosity = false)
                g = gcv(f)
                @test isfinite(g)
                @test g > 0
            end
        end

        @testset "gcv rejects a non-linear fit" begin
            fl = fit(SCEFit, dataset, Lasso(lambda = 1e-3); verbosity = false)
            @test_throws ArgumentError gcv(fl)
        end

        @testset "gcv_lambda consistency with gcv(fit)" begin
            lam = 1e-2
            path = gcv_lambda(dataset, [lam])
            fr = fit(SCEFit, dataset, Ridge(lambda = lam); verbosity = false)
            @test path.gcv_scores[1] ≈ gcv(fr) rtol = 1e-8
            @test path.torque_weight == 1.0

            lambdas = 10.0 .^ (-6:1.0:1)
            sweep = gcv_lambda(dataset, lambdas)
            @test length(sweep.gcv_scores) == length(lambdas)
            @test length(sweep.dof) == length(lambdas)
            @test sweep.lambda_best == lambdas[argmin(sweep.gcv_scores)]
            # dof decreases monotonically with penalty.
            @test issorted(sweep.dof; rev = true)
        end

        @testset "gcv_lambda input validation" begin
            @test_throws ArgumentError gcv_lambda(dataset, Float64[])
            @test_throws ArgumentError gcv_lambda(dataset, [-1.0, 1.0])
            @test_throws ArgumentError gcv_lambda(dataset, [1.0]; torque_weight = 1.5)
            @test_throws ArgumentError gcv_lambda(dataset, [1.0]; torque_weight = -0.1)
        end

        @testset "two-block weighting: dead rows excluded from N and dof" begin
            # A block zeroed by the weighting (energy at w=1, torque at w=0)
            # must not inflate N, and the j0 intercept dof counts only when the
            # energy block is live. Cross-check gcv(f) against a dense reference
            # computed over the live block alone.
            n_E = length(dataset.y_E)
            n_T = length(dataset.y_T)
            lam = 1e-3
            for w in (1.0, 0.0)
                f = fit(SCEFit, dataset, Ridge(lambda = lam);
                    torque_weight = w, verbosity = false)
                X, y = Magesty.Fitting.assemble_weighted_problem(
                    dataset.X_E, dataset.X_T, dataset.y_E, dataset.y_T, w)
                # Live rows: the block whose whitening scale is nonzero.
                rows = w == 1.0 ? ((n_E + 1):(n_E + n_T)) : (1:n_E)
                intercept = w < 1.0 ? 1 : 0
                Xl = X[rows, :]
                yl = y[rows]
                resid = yl - Xl * coef(f)
                n_live = length(yl)
                dof_ref = intercept + tr(Xl * ((Xl' * Xl + lam * I) \ Xl'))
                gcv_ref = (sum(abs2, resid) / n_live) / (1 - dof_ref / n_live)^2
                @test gcv(f) ≈ gcv_ref rtol = 1e-7
            end
        end

        @testset "gcv_learning_curve reproducibility and shape" begin
            sizes = [10, 20, 30]
            c1 = gcv_learning_curve(dataset, Ridge(lambda = 1e-3); sizes = sizes, repeats = 4, seed = 7)
            c2 = gcv_learning_curve(dataset, Ridge(lambda = 1e-3); sizes = sizes, repeats = 4, seed = 7)
            @test c1.sizes == sizes
            @test c1.gcv_mean == c2.gcv_mean          # exact reproducibility under seed
            @test c1.gcv_std == c2.gcv_std
            @test all(isfinite, c1.gcv_mean)
            @test c1.repeats == 4
            @test c1.seed == 7
            @test c1.estimator == Ridge(lambda = 1e-3)

            c3 = gcv_learning_curve(dataset, Ridge(lambda = 1e-3); sizes = sizes, repeats = 4, seed = 99)
            @test c1.gcv_mean != c3.gcv_mean          # a different seed draws differently

            # A full-size draw repeats the identical subset, so its spread is ~0.
            cfull = gcv_learning_curve(dataset, Ridge(lambda = 1e-3);
                sizes = [n_total], repeats = 3, seed = 1)
            @test cfull.gcv_std[1] < 1e-12
        end

        @testset "gcv_learning_curve statistics match a manual recomputation" begin
            sizes = [12, 24]
            repeats = 5
            seed = 3
            curve = gcv_learning_curve(dataset, Ridge(lambda = 1e-3);
                sizes = sizes, repeats = repeats, seed = seed)
            # Reproduce the seeded draws and score them through the public gcv(f).
            rng = Random.MersenneTwister(seed)
            for (i, n) in enumerate(sizes)
                scores = Float64[]
                for _ in 1:repeats
                    idx = Random.randperm(rng, n_total)[1:n]
                    f = fit(SCEFit, dataset[idx], Ridge(lambda = 1e-3); verbosity = false)
                    push!(scores, gcv(f))
                end
                @test curve.gcv_mean[i] ≈ sum(scores) / repeats rtol = 1e-9
            end
        end

        @testset "gcv_learning_curve validation" begin
            @test_throws ArgumentError gcv_learning_curve(dataset, Lasso(lambda = 1e-3); sizes = [10, 20])
            @test_throws ArgumentError gcv_learning_curve(dataset, OLS(); sizes = [10], repeats = 0)
            @test_throws ArgumentError gcv_learning_curve(dataset, OLS(); sizes = [0])
            @test_throws ArgumentError gcv_learning_curve(dataset, OLS(); sizes = [n_total + 1])
            @test_throws ArgumentError gcv_learning_curve(dataset, OLS();
                sizes = [10], torque_weight = 1.5)
        end

        @testset "writers" begin
            path = gcv_lambda(dataset, 10.0 .^ (-4:1.0:0))
            curve = gcv_learning_curve(dataset, Ridge(lambda = 1e-3); sizes = [15, 30], repeats = 3, seed = 0)
            dir = mktempdir()

            @testset "write_gcv_lambda format" begin
                fp = joinpath(dir, "gcv_lambda.txt")
                @test write_gcv_lambda(path, fp) === nothing
                text = read(fp, String)
                @test occursin("gcv_lambda", text)
                @test occursin("lambda_best", text)
                rows = _gcv_data_lines(text)
                @test length(rows) == length(path.lambdas)
                for (i, row) in enumerate(rows)
                    cols = split(row)
                    @test length(cols) == 3
                    @test parse(Float64, cols[1]) ≈ path.lambdas[i] rtol = 1e-9
                    @test parse(Float64, cols[2]) ≈ path.gcv_scores[i] rtol = 1e-7
                    @test parse(Float64, cols[3]) ≈ path.dof[i] rtol = 1e-7
                end
            end

            @testset "write_gcv_learning_curve format" begin
                fp = joinpath(dir, "gcv_learning_curve.txt")
                @test write_gcv_learning_curve(curve, fp) === nothing
                text = read(fp, String)
                @test occursin("gcv_learning_curve", text)
                rows = _gcv_data_lines(text)
                @test length(rows) == length(curve.sizes)
                for (i, row) in enumerate(rows)
                    cols = split(row)
                    @test length(cols) == 3
                    @test parse(Int, cols[1]) == curve.sizes[i]
                    @test parse(Float64, cols[2]) ≈ curve.gcv_mean[i] rtol = 1e-7
                    @test isapprox(parse(Float64, cols[3]), curve.gcv_std[i];
                        atol = 1e-12, rtol = 1e-7)
                end
            end

            @testset "default filenames" begin
                cd(mktempdir()) do
                    write_gcv_lambda(path)
                    write_gcv_learning_curve(curve)
                    @test isfile("gcv_lambda.txt")
                    @test isfile("gcv_learning_curve.txt")
                end
            end
        end
    end
end
