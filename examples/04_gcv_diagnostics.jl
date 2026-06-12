# GCV diagnostics: choose a ridge penalty and check data sufficiency.
#
# Run from the repository root:
#     julia --project=. examples/04_gcv_diagnostics.jl
#
# Generalized cross-validation (GCV) estimates out-of-sample prediction error
# from a single fit, on the same combined energy+torque weighted objective that
# `fit` minimizes. Two uses are shown:
#   (a) gcv_lambda  -- sweep the ridge penalty and pick the GCV minimizer.
#   (b) gcv_learning_curve -- sweep the training-set size and watch for a plateau.
#
# Uses the FePt L1_0 fixture from the test tree.

using Magesty

const FIXTURE = joinpath(@__DIR__, "..", "test", "integration", "fept_tetragonal_2x2x2")

function main()
    basis = SCEBasis(joinpath(FIXTURE, "input.toml"); verbosity = false)
    dataset = SCEDataset(basis, joinpath(FIXTURE, "EMBSET"); verbosity = false)
    println("Number of training configurations: ", length(dataset))

    # (a) Ridge penalty selection. A single SVD serves the whole path, so a
    #     fine lambda grid is cheap.
    lambdas = 10.0 .^ (-6:0.5:0)
    path = gcv_lambda(dataset, lambdas)
    # The raw GCV score is in the weighted-objective unit and is only meaningful
    # in relative comparison; the companion predictive R2 (gcv_r2 = 1 - GCV/MSY)
    # reads on a fixed scale (1 perfect, 0 matches the null model).
    println("\nGCV penalty sweep:")
    for (l, g, r, d) in zip(path.lambdas, path.gcv_scores, path.gcv_r2, path.dof)
        println("  lambda = $(l)\tGCV = $(g)\tR2 = $(r)\tdof = $(d)")
    end
    println("Selected lambda_best = ", path.lambda_best)
    write_gcv_lambda(path, "gcv_lambda.txt")
    # Plot:  python tools/FitCheck_gcv_lambda.py gcv_lambda.txt --r2

    # A single GCV score (and its predictive R2) for the fit at the selected penalty.
    f = fit(SCEFit, dataset, Ridge(lambda = path.lambda_best); verbosity = false)
    println("\nGCV of the fit at lambda_best: ", gcv(f), "  (R2 = ", gcv_r2(f), ")")

    # (b) Data-sufficiency sweep. At each size, several random subsets are fit
    #     and their GCV scores averaged; a flattening curve means enough data.
    curve = gcv_learning_curve(dataset, Ridge(lambda = path.lambda_best); repeats = 8, seed = 0)
    println("\nData-sufficiency sweep:")
    for (s, m, sd, r) in zip(curve.sizes, curve.gcv_mean, curve.gcv_std, curve.gcv_r2_mean)
        println("  size = $(s)\tGCV = $(m) +/- $(sd)\tR2 = $(r)")
    end
    write_gcv_learning_curve(curve, "gcv_learning_curve.txt")
    # Plot:  python tools/FitCheck_gcv_learning_curve.py gcv_learning_curve.txt --r2

    return nothing
end

main()
