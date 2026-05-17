# Basic flow: TOML -> SCEBasis -> SCEDataset -> fit -> predict.
#
# Run from the repository root:
#     julia --project=. examples/01_basic_flow.jl
#
# Uses the FePt L1_0 fixture from the test tree.

using Magesty
using LinearAlgebra: norm

const FIXTURE = joinpath(@__DIR__, "..", "test", "integration", "fept_tetragonal_2x2x2")

function main()
    # 1. Build the SCE basis from a TOML template.
    basis = SCEBasis(joinpath(FIXTURE, "input.toml"); verbosity = false)
    println("Number of SALCs: ", length(basis.salcbasis.salc_list))

    # 2. Build a dataset from the basis and an EMBSET file.
    dataset = SCEDataset(basis, joinpath(FIXTURE, "EMBSET.dat"))
    println("Number of training configurations: ", length(dataset))

    # 3. Fit with Ridge regression and inspect in-sample metrics.
    #    `torque_weight` defaults to 1.0 (torque-only fit); pass
    #    `torque_weight = 0.5` for a balanced joint fit.
    f = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
    println("Reference energy j0:  ", intercept(f))
    println("RMSE energy (in-sample): ", rmse_energy(f))
    println("R^2  energy (in-sample): ", r2_energy(f))
    println("RMSE torque (in-sample): ", rmse_torque(f))

    # 4. Predict on the first training configuration.
    sc1 = dataset.spinconfigs[1]
    println("Predicted energy for config 1: ", predict_energy(f, sc1))
    println("Observed  energy for config 1: ", sc1.energy)

    # 5. Predict per-atom torques and report the largest in magnitude.
    #    T is a 3 x num_atoms matrix; columns are per-atom torque vectors.
    T = predict_torque(f, sc1)
    max_atom = argmax(norm.(eachcol(T)))
    println("Largest predicted |torque| at atom $(max_atom): ", norm(T[:, max_atom]))

    return nothing
end

main()
