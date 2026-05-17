#!/usr/bin/env julia
# Worker for test_thread_safety.jl.
#
# Builds the SCE design matrices, fits coefficients, and evaluates
# predictions on the fept_tetragonal_2x2x2 fixture, then serializes the
# results to ARGS[1] so the parent process can compare across thread
# counts.

using Magesty
using Serialization
using TOML

const EXAMPLE_DIR = joinpath(@__DIR__, "..", "integration", "fept_tetragonal_2x2x2")

function compute_payload()
    input = TOML.parsefile(joinpath(EXAMPLE_DIR, "input.toml"))
    basis = SCEBasis(input; verbosity = false)
    dataset = SCEDataset(basis, joinpath(EXAMPLE_DIR, "EMBSET.dat"))
    fitted = fit(SCEFit, dataset, Ridge(lambda = 0.0); torque_weight = 0.5, verbosity = false)
    model = SCEModel(fitted)

    spinconfigs = dataset.spinconfigs
    e_pred = predict_energy(model, spinconfigs[1])
    t_pred = predict_torque(model, spinconfigs[1])

    return (
        nthreads = Threads.nthreads(),
        X_E = dataset.X_E,
        X_T = dataset.X_T,
        j0 = intercept(fitted),
        jphi = coef(fitted),
        e_pred = e_pred,
        t_pred = t_pred,
    )
end

length(ARGS) == 1 || error("usage: _thread_safety_worker.jl <output-path>")
serialize(ARGS[1], compute_payload())
