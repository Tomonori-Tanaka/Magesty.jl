using Test
using Magesty
import TOML

# Fixtures -------------------------------------------------------------

const DIMER_TOML_FIT = joinpath(@__DIR__, "..", "integration", "dimer", "input.toml")

# FM / AFM spin configurations for the 2-atom dimer.
function _dimer_configs_fit()
    SC = Magesty.SpinConfigs.SpinConfig
    fm = SC(-1.0, [1.0, 1.0],
            [0.0 0.0; 0.0 0.0; 1.0 1.0],
            [0.1 0.0; 0.0 -0.1; 0.0 0.0])
    afm = SC(1.0, [1.0, 1.0],
             [0.0 0.0; 0.0 0.0; -1.0 1.0],
             [-0.1 0.0; 0.0 0.1; 0.0 0.0])
    return [fm, afm]
end

# Tests ----------------------------------------------------------------

@testset "SCEFit" begin
    input = TOML.parsefile(DIMER_TOML_FIT)
    basis = SCEBasis(input; verbosity = false)
    configs = _dimer_configs_fit()
    dataset = SCEDataset(basis, configs)

    @testset "fit returns a populated SCEFit" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)
        @test f isa SCEFit
        @test f isa Magesty.StatsAPI.RegressionModel
        @test f.dataset === dataset
        @test f.estimator isa OLS
        @test f.torque_weight == 0.3
        # residuals span the augmented system: n_E + n_T rows
        @test length(f.residuals) == size(dataset.X_E, 1) + size(dataset.X_T, 1)
    end

    @testset "StatsAPI verbs" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)
        @test coef(f) === f.jphi
        @test intercept(f) === f.j0
        @test nobs(f) == 2                       # n_configs
        @test dof(f) == length(coef(f)) + 1

        # coef / intercept also work on a bare SCEModel
        m = SCEModel(basis, intercept(f), coef(f))
        @test coef(m) == f.jphi
        @test intercept(m) == f.j0
    end

    @testset "default torque_weight" begin
        f = fit(SCEFit, dataset, OLS())
        @test f.torque_weight == 0.5
    end

    @testset "Ridge estimator" begin
        f = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.3)
        @test f isa SCEFit
        @test f.estimator isa Ridge
    end

    @testset "SCEModel conversion" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)
        m = SCEModel(f)
        @test m isa SCEModel
        @test m.basis === basis
        @test m.j0 == f.j0
        @test m.jphi === f.jphi
    end

    @testset "predict_energy / predict_torque" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)
        m = SCEModel(f)
        n_atoms = basis.structure.supercell.num_atoms

        # single-config forms: SCEModel and SCEFit agree
        E1 = predict_energy(m, configs[1])
        @test E1 isa Float64
        @test predict_energy(f, configs[1]) == E1
        @test predict_energy(f, configs[1].spin_directions) ≈ E1

        T1 = predict_torque(m, configs[1])
        @test size(T1) == (3, n_atoms)
        @test predict_torque(f, configs[1]) == T1
        @test predict_torque(f, configs[1].spin_directions) ≈ T1

        # dataset-batch forms match the per-config forms and the design matrices
        E_batch = predict_energy(m, dataset)
        @test E_batch isa Vector{Float64}
        @test length(E_batch) == 2
        @test E_batch ≈ [predict_energy(m, sc) for sc in configs]
        @test E_batch ≈ dataset.X_E[:, 2:end] * f.jphi .+ f.j0

        T_batch = predict_torque(m, dataset)
        @test T_batch isa Vector{Matrix{Float64}}
        @test length(T_batch) == 2
        @test all(size(t) == (3, n_atoms) for t in T_batch)
        @test reduce(vcat, vec.(T_batch)) ≈ dataset.X_T * f.jphi
        for (i, sc) in enumerate(configs)
            @test T_batch[i] ≈ predict_torque(m, sc)
        end
    end

    @testset "evaluation verbs are self-consistent" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)

        # residuals / rss are self-consistent
        re = residuals_energy(f)
        @test length(re) == 2
        @test re == dataset.y_E .- predict_energy(f, dataset)
        @test rss_energy(f) ≈ sum(abs2, re)

        rt = residuals_torque(f)
        @test length(rt) == 3 * basis.structure.supercell.num_atoms * 2
        @test rss_torque(f) ≈ sum(abs2, rt)
        @test rmse_energy(f) ≈ sqrt(rss_energy(f) / length(re))
        @test rmse_torque(f) ≈ sqrt(rss_torque(f) / length(rt))
    end

    @testset "evaluation verbs accept configs and datasets" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)
        m = SCEModel(f)

        # (model, dataset) is the base method; (fit, dataset) delegates
        @test r2_energy(m, dataset) == r2_energy(f, dataset)
        @test r2_energy(f, dataset) == r2_energy(f)

        # raw configs path agrees with the dataset path
        @test r2_energy(f, configs) ≈ r2_energy(f, dataset)
        @test rmse_torque(m, configs) ≈ rmse_torque(m, dataset)
        @test residuals_energy(m, configs) ≈ residuals_energy(m, dataset)

        # out-of-sample evaluation on a slice runs and is finite
        test_slice = dataset[2:2]
        @test isfinite(rmse_energy(f, test_slice))
        @test length(residuals_torque(f, test_slice)) ==
              3 * basis.structure.supercell.num_atoms
    end

    @testset "basis-identity check" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)
        m = SCEModel(f)
        # a dataset built from a different SCEBasis is rejected
        other_basis = SCEBasis(input; verbosity = false)
        other_dataset = SCEDataset(other_basis, configs)
        @test_throws ArgumentError r2_energy(f, other_dataset)
        @test_throws ArgumentError r2_energy(m, other_dataset)
        @test_throws ArgumentError predict_energy(m, other_dataset)
        @test_throws ArgumentError predict_torque(f, other_dataset)
    end
end
