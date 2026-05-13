using Test
using Random
using Magesty

# Regression test for the 3-layer estimator dispatch in src/Optimize.jl.
#
# Golden (j0, jphi) values were captured from the pre-refactor
# `elastic_net_regression` implementation on the fixture below
# (`Random.seed!(20260513)`, WEIGHT = 0.3). Any change here that alters
# numerical output trips this test.
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
        # Energy rows come first; bias column is reset to 1.0 after √weight scaling.
        @test all(X[1:n_config, 1] .== 1.0)
        # Torque block has a leading zero bias column.
        @test all(X[(n_config + 1):end, 1] .== 0.0)
        @test size(X, 1) == n_config + 3 * n_atoms * n_config
        @test size(X, 2) == n_basis + 1
        @test length(y) == size(X, 1)
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

    @testset "golden (j0, jphi) — OLS-equivalent (alpha=0, lambda=0)" begin
        j0, jphi = Magesty.Optimize._fit_sce_model_internal(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            Magesty.OLS(),
            weight,
        )
        # Golden values captured 2026-05-13 from pre-refactor
        # `elastic_net_regression(..., 0.0, 0.0, 0.3)`.
        @test isapprox(j0, 0.77052178555705997; atol = 1e-12, rtol = 1e-10)
        @test isapprox(jphi[1],  0.15051107443703846; atol = 1e-12, rtol = 1e-10)
        @test isapprox(jphi[2], -0.3107340994547646;  atol = 1e-12, rtol = 1e-10)
    end

    @testset "golden (j0, jphi) — Ridge (lambda=0.1)" begin
        j0, jphi = Magesty.Optimize._fit_sce_model_internal(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            Magesty.Ridge(lambda = 0.1),
            weight,
        )
        # Golden values captured 2026-05-13 from the pre-refactor
        # `elastic_net_regression(..., 0.0, 0.1, 0.3)` call. The current
        # Ridge(lambda=0.1) path is numerically identical because the
        # historical `alpha` field was always ignored.
        @test isapprox(j0, 0.76894461111080004; atol = 1e-12, rtol = 1e-10)
        @test isapprox(jphi[1],  0.14994907999948034; atol = 1e-12, rtol = 1e-10)
        @test isapprox(jphi[2], -0.30765628547929808; atol = 1e-12, rtol = 1e-10)
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
