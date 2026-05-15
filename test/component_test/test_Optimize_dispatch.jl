using Test
using Random
using Magesty

# Regression test for the 3-layer estimator dispatch in src/Optimize.jl.
#
# Golden (j0, jphi) values were recaptured on 2026-05-14 after
# `assemble_weighted_problem` switched to per-sample MSE normalization
# (energy rows scaled by √((1-w)/n_E), torque rows by √(w/n_T)). The
# fixture below uses `Random.seed!(20260513)`, WEIGHT = 0.3. Any change
# here that alters numerical output trips this test.
#
# See: docs/specs/260513-estimator-dispatch/

@testset "Optimize estimator dispatch" begin
    Random.seed!(20260513)

    n_config = 4
    n_atoms  = 2
    n_basis  = 2

    design_matrix_energy = hcat(ones(n_config), randn(n_config, n_basis))
    design_matrix_torque = randn(3 * n_atoms * n_config, n_basis)
    observed_energy_list = randn(n_config)
    observed_torque_list = [randn(3, n_atoms) for _ in 1:n_config]
    weight = 0.3

    @testset "assemble_weighted_problem" begin
        X, y, bias_col = Magesty.Optimize.assemble_weighted_problem(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            weight,
        )
        @test bias_col == 1
        # Energy rows come first; bias column is reset to 1.0 after scaling.
        @test all(X[1:n_config, 1] .== 1.0)
        # Torque block has a leading zero bias column.
        @test all(X[(n_config + 1):end, 1] .== 0.0)
        @test size(X, 1) == n_config + 3 * n_atoms * n_config
        @test size(X, 2) == n_basis + 1
        @test length(y) == size(X, 1)

        # Per-sample MSE normalization: energy block scaled by √((1-w)/n_E),
        # torque block by √(w/n_T) with n_T = 3 * n_atoms * n_config. This
        # makes `weight` a convex combination of the two MSEs independent
        # of each block's row count.
        n_T = 3 * n_atoms * n_config
        scale_e = sqrt((1 - weight) / n_config)
        scale_m = sqrt(weight / n_T)
        @test X[1:n_config, 2:end] ≈ design_matrix_energy[:, 2:end] .* scale_e
        @test y[1:n_config] ≈ observed_energy_list .* scale_e
        @test X[(n_config + 1):end, 2:end] ≈ design_matrix_torque .* scale_m
    end

    @testset "solve_coefficients(::OLS) matches X \\ y" begin
        X, y, bias_col = Magesty.Optimize.assemble_weighted_problem(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            weight,
        )
        j_ols = Magesty.Optimize.solve_coefficients(
            Magesty.OLS(), X, y; bias_col = bias_col,
        )
        @test j_ols ≈ X \ y
    end

    # Helper: run the three-stage estimator dispatch (assemble -> solve ->
    # extract) end-to-end and return (j0, jphi). Mirrors the kernel that
    # `fit(SCEFit, dataset, est; torque_weight)` runs under the hood.
    function _fit_kernel(estimator)
        X, y, bias_col = Magesty.Optimize.assemble_weighted_problem(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            weight,
        )
        j_values = Magesty.Optimize.solve_coefficients(
            estimator, X, y; bias_col = bias_col,
        )
        return Magesty.Optimize.extract_j0_jphi(
            j_values, design_matrix_energy, observed_energy_list,
        )
    end

    @testset "golden (j0, jphi) — OLS" begin
        j0, jphi = _fit_kernel(Magesty.OLS())
        # Golden values recaptured 2026-05-14 after per-sample MSE
        # normalization (was: 0.77052.., 0.15051.., -0.31073..).
        @test isapprox(j0, 0.73174756944306574; atol = 1e-12, rtol = 1e-10)
        @test isapprox(jphi[1],  0.15794098335672488; atol = 1e-12, rtol = 1e-10)
        @test isapprox(jphi[2], -0.24272817525866192; atol = 1e-12, rtol = 1e-10)
    end

    @testset "golden (j0, jphi) — Ridge (lambda=0.1)" begin
        j0, jphi = _fit_kernel(Magesty.Ridge(lambda = 0.1))
        # Golden values recaptured 2026-05-14 after per-sample MSE
        # normalization (was: 0.76894.., 0.14995.., -0.30766..).
        @test isapprox(j0, 0.71912905195775911; atol = 1e-12, rtol = 1e-10)
        @test isapprox(jphi[1],  0.15469233524793524; atol = 1e-12, rtol = 1e-10)
        @test isapprox(jphi[2], -0.2185533681028447; atol = 1e-12, rtol = 1e-10)
    end

    @testset "extract_j0_jphi splits augmented coefficients" begin
        X, y, bias_col = Magesty.Optimize.assemble_weighted_problem(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            weight,
        )
        j_values = Magesty.Optimize.solve_coefficients(
            Magesty.OLS(), X, y; bias_col = bias_col,
        )
        j0, jphi = Magesty.Optimize.extract_j0_jphi(
            j_values, design_matrix_energy, observed_energy_list,
        )
        @test length(jphi) == n_basis
        # j0 is recomputed from the unscaled energy residual, not lifted from j_values[1].
        @test j0 ≈ sum(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi) / n_config
    end
end
