# Cache the SALC basis to XML, then reuse it later for fitting.
#
# Run from the repository root:
#     julia --project=. examples/03_save_load.jl
#
# Uses the FePt L1_0 fixture from the test tree. The SALC construction is
# the expensive step; caching the basis avoids recomputing it across
# sessions or estimator sweeps.

using Magesty

const FIXTURE = joinpath(@__DIR__, "..", "test", "integration", "fept_tetragonal_2x2x2")

function main()
    mktempdir() do dir
        basis_path = joinpath(dir, "basis.xml")
        model_path = joinpath(dir, "model.xml")

        # --- Session 1: build and cache the basis, fit and save a model.
        basis = SCEBasis(joinpath(FIXTURE, "input.toml"); verbosity = false)
        Magesty.save(basis, basis_path)

        dataset = SCEDataset(basis, joinpath(FIXTURE, "EMBSET.dat"))
        f = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
        Magesty.save(SCEModel(f), model_path)

        println("Saved basis to ", basis_path)
        println("Saved model to ", model_path)
        println("In-sample RMSE energy: ", rmse_energy(f))

        # --- Session 2: reload basis and model without rebuilding SALCs.
        reloaded_basis = Magesty.load(SCEBasis, basis_path)
        reloaded_model = Magesty.load(SCEModel, model_path)
        println("Reloaded basis: ",
                length(reloaded_basis.salcbasis.salc_list), " SALCs")
        println("Reloaded model j0: ", reloaded_model.j0)
        println("Reloaded model jphi[1:3]: ",
                reloaded_model.jphi[1:min(3, length(reloaded_model.jphi))])

        # The reloaded predictor produces identical predictions.
        sc1 = dataset.spinconfigs[1]
        println("Predicted energy from reloaded model: ",
                predict_energy(reloaded_model, sc1))
    end
    return nothing
end

main()
