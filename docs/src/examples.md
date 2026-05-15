# Examples

This page demonstrates the new SCE public API:
`SCEBasis` → `SCEDataset` → `fit(SCEFit, ...)` → `SCEModel` →
`save` / `load`. Runnable versions of these examples live in the
`examples/` directory of the repository.

## Example 1: Basic flow

Load a TOML template, build the dataset, fit, and inspect in-sample
metrics.

```julia
using Magesty

basis   = SCEBasis("input.toml")
dataset = SCEDataset(basis, "EMBSET.dat")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)

println("RMSE energy: ", rmse_energy(f) * 1000, " meV")
println("R^2  energy: ", r2_energy(f))
println("RMSE torque: ", rmse_torque(f) * 1000, " meV")
println("j0: ", intercept(f), "  number of SCE coefficients: ", length(coef(f)))

save(SCEModel(f), "model.xml")
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
dataset = SCEDataset(basis, "EMBSET.dat")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)
```

## Example 3: Estimator comparison (reuse one dataset)

`SCEDataset` stores the unweighted energy and torque design matrices, so
you can sweep estimators and torque weights without rebuilding it.

```julia
using Magesty

dataset = SCEDataset(SCEBasis("input.toml"), "EMBSET.dat")
for est in (OLS(), Ridge(lambda = 1e-4), Ridge(lambda = 1e-2))
    f = fit(SCEFit, dataset, est; torque_weight = 0.5)
    println(est, ": RMSE energy = ", rmse_energy(f))
end
```

## Example 4: Train / test split

`SCEDataset` supports vector indexing; the slice is a fresh dataset that
shares the same `SCEBasis`.

```julia
using Magesty

dataset = SCEDataset(SCEBasis("input.toml"), "EMBSET.dat")
n_train = round(Int, 0.8 * length(dataset))
train   = dataset[1:n_train]
test    = dataset[(n_train + 1):end]

f = fit(SCEFit, train, Ridge(lambda = 1e-4); torque_weight = 0.5)
println("in-sample  RMSE energy: ", rmse_energy(f))
println("out-sample RMSE energy: ", rmse_energy(f, test))
```

## Example 5: Cache the SALC basis

Building an `SCEBasis` runs the heavy SALC construction. Persist it to
XML and reload it later to skip that step.

```julia
using Magesty

# Session 1: cache the basis.
save(SCEBasis("input.toml"), "basis.xml")

# Session 2: reload, fit, save the model.
basis   = load(SCEBasis, "basis.xml")
dataset = SCEDataset(basis, "EMBSET.dat")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)
save(SCEModel(f), "model.xml")
```

## Example 6: Predict and evaluate on custom configurations

```julia
using Magesty
using LinearAlgebra

basis   = SCEBasis("input.toml")
dataset = SCEDataset(basis, "EMBSET.dat")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)
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

## Example 7: Inspect symmetry, structure, and SALCs

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
