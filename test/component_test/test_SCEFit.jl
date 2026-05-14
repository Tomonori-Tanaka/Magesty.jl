using Test
using Magesty
import TOML

# Fixtures -------------------------------------------------------------

const DIMER_TOML_FIT = joinpath(@__DIR__, "..", "examples", "dimer", "input.toml")

# FM / AFM spin configurations for the 2-atom dimer.
function _dimer_configs_fit()
    SC = Magesty.SpinConfigs.SpinConfig
    fm = SC(-1.0, [1.0, 1.0],
            [0.0 0.0; 0.0 0.0; 1.0 1.0],
            [0.0 0.0; 0.0 0.0; 0.0 0.0])
    afm = SC(1.0, [1.0, 1.0],
             [0.0 0.0; 0.0 0.0; -1.0 1.0],
             [0.0 0.0; 0.0 0.0; 0.0 0.0])
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
        @test Set(keys(f.metrics)) ==
              Set([:rmse_energy, :rmse_torque, :r2score_energy, :r2score_torque])
    end

    @testset "StatsAPI verbs" begin
        f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)
        @test coef(f) === f.jphi
        @test intercept(f) === f.j0
        @test nobs(f) == 2                       # n_configs
        @test dof(f) == length(coef(f)) + 1

        # coef / intercept also work on a bare SCEModel
        m = SCEModel(intercept(f), coef(f), basis.salcbasis, basis.symmetry,
                     basis.structure.supercell.num_atoms)
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

    @testset "golden (j0, jphi) match the legacy fit_sce_model path" begin
        system = System(input; verbosity = false)
        for (est, w) in ((OLS(), 0.3), (Ridge(lambda = 1e-4), 0.3),
                         (OLS(), 0.0), (OLS(), 1.0))
            f = fit(SCEFit, dataset, est; torque_weight = w)
            opt = fit_sce_model(system, configs, est, w; verbosity = false)
            @test intercept(f) ≈ opt.reference_energy
            @test coef(f) ≈ opt.SCE
        end
    end
end
