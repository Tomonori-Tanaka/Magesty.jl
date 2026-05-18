using Test
using Random
using Statistics
using Magesty

# Regression test for the 3-layer estimator dispatch in src/Fitting.jl.
#
# Golden (j0, jphi) values are mathematically equivalent to the
# pre-refactor bias-in-X formulation up to floating-point rounding
# noise: `assemble_weighted_problem` now mean-centers the energy block
# and drops the bias column, and `extract_j0_jphi` recovers j0 via
# `mean(y_E - X_E * jphi)` (block elimination of j0 from the joint
# energy/torque objective). The fixture below uses
# `Random.seed!(20260513)`, WEIGHT = 0.3. Any change here that alters
# numerical output trips this test.
#
# See: docs/specs/260518-energy-centered-design-matrix/

@testset "Fitting estimator dispatch" begin
    Random.seed!(20260513)

    n_config = 4
    n_atoms  = 2
    n_basis  = 2

    design_matrix_energy = randn(n_config, n_basis)
    design_matrix_torque = randn(3 * n_atoms * n_config, n_basis)
    observed_energy_list = randn(n_config)
    observed_torque_list = [randn(3, n_atoms) for _ in 1:n_config]
    weight = 0.3

    @testset "assemble_weighted_problem" begin
        X, y = Magesty.Fitting.assemble_weighted_problem(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            weight,
        )
        @test size(X, 1) == n_config + 3 * n_atoms * n_config
        @test size(X, 2) == n_basis             # no bias column
        @test length(y) == size(X, 1)

        # Per-sample MSE normalization: energy block scaled by √((1-w)/n_E)
        # *after* mean-centering, torque block by √(w/n_T) with
        # n_T = 3 * n_atoms * n_config. This makes `weight` a convex
        # combination of the two MSEs independent of each block's row
        # count and removes j0 from the system via block elimination
        # (mean(y_E - X_E * jphi) = j0_star).
        n_T = 3 * n_atoms * n_config
        scale_e = sqrt((1 - weight) / n_config)
        scale_m = sqrt(weight / n_T)
        centered_X_E = design_matrix_energy .- mean(design_matrix_energy, dims = 1)
        centered_y_E = observed_energy_list .- mean(observed_energy_list)
        @test X[1:n_config, :]        ≈ centered_X_E .* scale_e rtol = 1e-12
        @test y[1:n_config]           ≈ centered_y_E .* scale_e rtol = 1e-12
        @test X[(n_config + 1):end, :] ≈ design_matrix_torque .* scale_m rtol = 1e-12
        @test y[(n_config + 1):end]    ≈ vec(reduce(hcat, observed_torque_list)) .* scale_m rtol = 1e-12
    end

    @testset "solve_coefficients(::OLS) matches X \\ y" begin
        X, y = Magesty.Fitting.assemble_weighted_problem(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            weight,
        )
        j_ols = Magesty.Fitting.solve_coefficients(Magesty.OLS(), X, y)
        @test j_ols ≈ X \ y
    end

    # Helper: run the three-stage estimator dispatch (assemble -> solve ->
    # extract) end-to-end and return (j0, jphi). Mirrors the kernel that
    # `fit(SCEFit, dataset, est; torque_weight)` runs under the hood.
    function _fit_kernel(estimator)
        X, y = Magesty.Fitting.assemble_weighted_problem(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            weight,
        )
        j_values = Magesty.Fitting.solve_coefficients(estimator, X, y)
        return Magesty.Fitting.extract_j0_jphi(
            j_values, design_matrix_energy, observed_energy_list,
        )
    end

    @testset "golden (j0, jphi) — OLS" begin
        j0, jphi = _fit_kernel(Magesty.OLS())
        # Pre-refactor (bias-in-X) values from 2026-05-14:
        #   j0 = 0.73174756944306574
        #   jphi = [0.15794098335672488, -0.24272817525866192]
        # The centered/eliminated formulation must agree within
        # floating-point rounding noise.
        @test isapprox(j0,        0.73174756944306574; rtol = 1e-9)
        @test isapprox(jphi[1],   0.15794098335672488; rtol = 1e-9)
        @test isapprox(jphi[2],  -0.24272817525866192; rtol = 1e-9)
    end

    @testset "golden (j0, jphi) — Ridge (lambda=0.1)" begin
        j0, jphi = _fit_kernel(Magesty.Ridge(lambda = 0.1))
        # Pre-refactor (bias-in-X) values from 2026-05-14:
        #   j0 = 0.71912905195775911
        #   jphi = [0.15469233524793524, -0.2185533681028447]
        @test isapprox(j0,        0.71912905195775911; rtol = 1e-9)
        @test isapprox(jphi[1],   0.15469233524793524; rtol = 1e-9)
        @test isapprox(jphi[2],  -0.2185533681028447;  rtol = 1e-9)
    end

    @testset "extract_j0_jphi recovers j0 from unscaled residual" begin
        X, y = Magesty.Fitting.assemble_weighted_problem(
            design_matrix_energy,
            design_matrix_torque,
            observed_energy_list,
            observed_torque_list,
            weight,
        )
        j_values = Magesty.Fitting.solve_coefficients(Magesty.OLS(), X, y)
        j0, jphi = Magesty.Fitting.extract_j0_jphi(
            j_values, design_matrix_energy, observed_energy_list,
        )
        @test length(jphi) == n_basis
        @test jphi == j_values                 # no bias slot to discard
        # j0 is the closed-form solution `mean(y_E - X_E * jphi)`.
        @test j0 ≈ sum(observed_energy_list .- design_matrix_energy * jphi) / n_config
    end
end
