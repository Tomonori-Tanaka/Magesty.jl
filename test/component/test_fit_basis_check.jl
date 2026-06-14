using Test
using Magesty
import TOML

const DIMER_TOML_BC = joinpath(@__DIR__, "..", "integration", "dimer", "input.toml")

function _dimer_configs_bc()
    SC = Magesty.SpinConfigs.SpinConfig
    fm = SC(-1.0, [1.0, 1.0],
            [0.0 0.0; 0.0 0.0; 1.0 1.0],
            [0.1 0.0; 0.0 -0.1; 0.0 0.0])
    afm = SC(1.0, [1.0, 1.0],
             [0.0 0.0; 0.0 0.0; -1.0 1.0],
             [-0.1 0.0; 0.0 0.1; 0.0 0.0])
    return [fm, afm]
end

@testset "basis fingerprint mismatch is rejected" begin
    input = TOML.parsefile(DIMER_TOML_BC)
    basis = SCEBasis(input; verbosity = false)
    configs = _dimer_configs_bc()
    dataset = SCEDataset(basis, configs)
    f = fit(SCEFit, dataset, OLS(); torque_weight = 0.3, verbosity = false)
    m = SCEModel(f)

    # The 5-arg default constructor lets a test inject a deliberately
    # corrupted fingerprint while keeping every other structural field
    # identical. This simulates a SALCBasis whose key-group ordering
    # diverged from the predictor's basis even though the column count
    # happens to match.
    bad_fp = basis.salc_fingerprint + UInt64(1)
    bad_basis = SCEBasis(basis.structure, basis.symmetry, basis.salcbasis,
                         basis.isotropy, bad_fp)
    bad_dataset = SCEDataset(bad_basis, configs)

    @test bad_dataset.basis.salc_fingerprint != m.basis.salc_fingerprint

    # Every entry point that runs `_check_basis` rejects the mismatch.
    @test_throws ArgumentError predict_energy(m, bad_dataset)
    @test_throws ArgumentError predict_energy(f, bad_dataset)
    @test_throws ArgumentError predict_torque(m, bad_dataset)
    @test_throws ArgumentError predict_torque(f, bad_dataset)
    @test_throws ArgumentError r2_energy(m, bad_dataset)
    @test_throws ArgumentError r2_energy(f, bad_dataset)
    @test_throws ArgumentError residuals_energy(f, bad_dataset)
    @test_throws ArgumentError rmse_torque(m, bad_dataset)
end
