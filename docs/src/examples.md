# Examples

This page demonstrates the SCE public API:
`SCEBasis` → `SCEDataset` → `fit(SCEFit, ...)` → `SCEModel` →
`Magesty.save` / `Magesty.load`. Runnable scripts for several of these
examples live in the `examples/` directory of the repository.

`save` and `load` are intentionally **not exported** to avoid clashing
with the generic names from JLD2, FileIO, CSV.jl, etc.; call them as
`Magesty.save(...)` / `Magesty.load(...)`.

## Example 1: Basic flow

Load a TOML template, build the dataset, fit, and inspect in-sample
metrics.

```julia
using Magesty

basis   = SCEBasis("input.toml")
dataset = SCEDataset(basis, "EMBSET")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4))

println("RMSE energy: ", rmse_energy(f) * 1000, " meV")
println("R^2  energy: ", r2_energy(f))
println("RMSE torque: ", rmse_torque(f) * 1000, " meV")
println("j0: ", intercept(f), "  number of SCE coefficients: ", length(coef(f)))

Magesty.save(SCEModel(f), "model.xml")
```

## Example 2: Building from a CIF file (AtomsIO)

`AtomsIO` is not a Magesty dependency — install it into the active
environment when you want to load structure files. The `interaction`
keyword mirrors the TOML `[interaction]` table.

```julia
using Magesty
using AtomsIO   # ] add AtomsIO

system  = load_system("FePt.cif")
basis   = SCEBasis(
    system;
    interaction = (
        body1 = (lmax = Dict(:Fe => 2, :Pt => 2),),
        body2 = (lsum = 2,
                 cutoff = Dict((:Fe, :Fe) => -1.0,
                               (:Fe, :Pt) => -1.0,
                               (:Pt, :Pt) => -1.0)),
    ),
    tolerance_sym = 1e-5,
    isotropy = false,
)
dataset = SCEDataset(basis, "EMBSET")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
```

## Example 3: Estimator comparison (reuse one dataset)

`SCEDataset` stores the unweighted energy and torque design matrices, so
you can sweep estimators and torque weights without rebuilding it.

```julia
using Magesty

dataset = SCEDataset(SCEBasis("input.toml"), "EMBSET")
for est in (OLS(), Ridge(lambda = 1e-4), Ridge(lambda = 1e-2))
    f = fit(SCEFit, dataset, est)
    println(est, ": RMSE energy = ", rmse_energy(f))
end
```

## Example 4: Train / test split

`SCEDataset` supports vector indexing; the slice is a fresh dataset that
shares the same `SCEBasis`.

```julia
using Magesty

dataset = SCEDataset(SCEBasis("input.toml"), "EMBSET")
n_train = round(Int, 0.8 * length(dataset))
train   = dataset[1:n_train]
test    = dataset[(n_train + 1):end]

f = fit(SCEFit, train, Ridge(lambda = 1e-4))
println("in-sample  RMSE energy: ", rmse_energy(f))
println("out-sample RMSE energy: ", rmse_energy(f, test))
```

## Example 5: Cache the SALC basis

Building an `SCEBasis` runs the heavy SALC construction. Persist it to
XML and reload it later to skip that step.

```julia
using Magesty

# Session 1: cache the basis.
Magesty.save(SCEBasis("input.toml"), "basis.xml")

# Session 2: reload, fit, save the model.
basis   = Magesty.load(SCEBasis, "basis.xml")
dataset = SCEDataset(basis, "EMBSET")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
Magesty.save(SCEModel(f), "model.xml")
```

## Example 6: Predict and evaluate on custom configurations

```julia
using Magesty
using LinearAlgebra

basis   = SCEBasis("input.toml")
dataset = SCEDataset(basis, "EMBSET")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
model   = SCEModel(f)

num_atoms = basis.structure.supercell.num_atoms

# Ferromagnetic spin configuration along +z.
spins_fm = zeros(3, num_atoms)
spins_fm[3, :] .= 1.0
println("FM energy: ", predict_energy(model, spins_fm), " eV")

# Spiral configuration in the xy-plane.
spins_spiral = zeros(3, num_atoms)
for i in 1:num_atoms
    angle = 2π * (i - 1) / num_atoms
    spins_spiral[1, i] = cos(angle)
    spins_spiral[2, i] = sin(angle)
end
T = predict_torque(model, spins_spiral)
println("Max torque magnitude: ", maximum(norm.(eachcol(T))))
```

## Example 7: Evaluate a saved model on new data

`SCEDataset` accepts a fitted `SCEModel` or `SCEFit` in place of an
`SCEBasis`: it reuses the basis embedded in it (`model.basis` /
`fit.dataset.basis`) without rebuilding the SALCs. This is convenient
when reloading a saved model to score it on a fresh EMBSET file.

```julia
using Magesty

# Reload a model saved earlier and score it on new reference data.
model    = Magesty.load(SCEModel, "model.xml")
new_data = SCEDataset(model, "EMBSET_new")   # basis taken from the model

ŷ_E = predict_energy(model, new_data)
ŷ_T = predict_torque(model, new_data)
```

## Example 8: Inspect symmetry, structure, and SALCs

```julia
using Magesty

basis = SCEBasis("input.toml")

sym = basis.symmetry
println("Space group: ", sym.international_symbol,
        " (#", sym.spacegroup_number, ")")
println("Number of symmetry operations: ", sym.nsym)

cell = basis.structure.supercell
println("Number of atoms: ", cell.num_atoms)
println("Lattice (columns):")
display(Matrix(cell.lattice_vectors))

println("Number of SALCs: ", length(basis.salcbasis.salc_list))
```

## Example 9: GCV diagnostics (penalty and data sufficiency)

Generalized cross-validation scores the combined energy+torque objective
from a single fit. Use `gcv_lambda` to pick a ridge penalty and
`gcv_learning_curve` to check whether there is enough training data (a flattening
curve). Both work only for the linear estimators (`OLS`, `Ridge`,
`AdaptiveRidge`). The raw GCV score is in the weighted-objective unit and only
meaningful in relative comparison; `gcv_r2` (and the `gcv_r2` columns on the
sweep results) gives the companion predictive R² on a fixed scale (`1` perfect,
`0` matches the null model). See [Cross-validation diagnostics](@ref) for the
formula.

```julia
using Magesty

dataset = SCEDataset(SCEBasis("input.toml"), "EMBSET")

# Penalty selection: one SVD serves the whole lambda path.
path = gcv_lambda(dataset, 10.0 .^ (-6:0.5:0))
println("lambda_best = ", path.lambda_best)
f = fit(SCEFit, dataset, Ridge(lambda = path.lambda_best))
println("GCV at lambda_best: ", gcv(f), "  (R2 = ", gcv_r2(f), ")")
write_gcv_lambda(path, "gcv_lambda.txt")
# python tools/FitCheck_gcv_lambda.py gcv_lambda.txt --r2

# Data sufficiency: average GCV over random subsets at each size.
curve = gcv_learning_curve(dataset, Ridge(lambda = path.lambda_best); repeats = 8)
write_gcv_learning_curve(curve, "gcv_learning_curve.txt")
# python tools/FitCheck_gcv_learning_curve.py gcv_learning_curve.txt --r2
```
