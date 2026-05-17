# Tutorial

This tutorial walks through the basic usage of Magesty.jl for spin-cluster
expansion (SCE) analysis.

## Configuration file format

Magesty.jl reads structure and interaction settings from a TOML file.
Fit parameters (estimator, regularization, torque weight) are passed in
Julia at `fit` time and **do not** live in the TOML.

```toml:input.toml
[general]
name = "bccfe"
kd   = ["Fe"]           # list of element names
nat  = 16               # total number of atoms in the supercell
periodicity = [true, true, true]  # apply periodic boundary to all (x, y, z) directions

[symmetry]
tolerance = 1e-5        # symmetry detection tolerance
isotropy = true         # restrict to Lf = 0 (isotropic exchange) terms

[interaction]
nbody = 2               # maximum interaction body order
[interaction.body1]
lmax.Fe = 0             # on-site (1-body) maximum angular momentum per element
[interaction.body2]
lsum = 2                # cutoff sum of l values for 2-body basis functions
cutoff."Fe-Fe" = -1     # pairwise cutoff radius in Å (-1 = include all pairs)

[structure]
kd_list  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # element index per atom
lattice  = [
  [5.66, 0.0, 0.0],
  [0.0, 5.66, 0.0],
  [0.0, 0.0, 5.66],
]
position = [
  [0.00, 0.00, 0.00],
  [0.25, 0.25, 0.25],
  # ... (16 atoms total)
]
```

For a full key reference see [Input Keys](input_keys.md).

## Basic workflow

```julia
using Magesty

# 1. Build the SCE basis (heavy: runs SALC construction).
basis = SCEBasis("input.toml")

# 2. Pair the basis with training data to form a dataset.
dataset = SCEDataset(basis, "EMBSET.dat")

# 3. Fit. `torque_weight` is the convex combination of per-sample
#    energy MSE and torque MSE that the augmented least-squares problem
#    minimizes — `0.5` weights both equally.
f = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)

println("Reference energy j0: ", intercept(f), " eV")
println("Number of SCE coefficients: ", length(coef(f)))
println("RMSE energy: ", rmse_energy(f) * 1000, " meV")
println("R²    energy: ", r2_energy(f))
```

## Exporting and reloading results

`Magesty.save` / `Magesty.load` round-trip a basis or model via XML. They
are deliberately **not exported** so that the generic `save` / `load`
exported by JLD2, FileIO, CSV.jl, etc. can coexist in the same script;
call them through the module:

```julia
# Persist a fitted model.
Magesty.save(SCEModel(f), "model.xml")

# Persist just the basis (skip the expensive SALC build next time).
Magesty.save(basis, "basis.xml")

# Later session: reload and refit without recomputing SALCs.
basis   = Magesty.load(SCEBasis, "basis.xml")
dataset = SCEDataset(basis, "EMBSET.dat")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)
```

The XML schema is shared with the Monte Carlo package `SpinClusterMC.jl`.

## Predicting on new configurations

```julia
model = SCEModel(f)

# Single configuration: a 3 × num_atoms matrix of unit spin directions
# (column = atom; rows are x, y, z).
num_atoms = basis.structure.supercell.num_atoms
spin_directions = zeros(3, num_atoms)
spin_directions[3, :] .= 1.0   # all spins along +z (ferromagnetic)

E = predict_energy(model, spin_directions)
T = predict_torque(model, spin_directions)   # 3 × num_atoms matrix
```

## Evaluation on held-out data

`SCEDataset` supports indexing, so a train / test split is one line:

```julia
n_train = round(Int, 0.8 * length(dataset))
train, test = dataset[1:n_train], dataset[(n_train + 1):end]

f = fit(SCEFit, train, Ridge(lambda = 1e-4); torque_weight = 0.5)
println("in-sample  RMSE energy: ", rmse_energy(f))
println("out-sample RMSE energy: ", rmse_energy(f, test))
```

## Inspecting the basis

```julia
sym = basis.symmetry
println("Space group: ", sym.international_symbol,
        " (#", sym.spacegroup_number, ")")
println("Number of symmetry operations: ", sym.nsym)

cell = basis.structure.supercell
println("Number of atoms: ", cell.num_atoms)
println("Number of SALCs: ", length(basis.salcbasis.salc_list))
```

For detailed function documentation see the [API Reference](@ref). More
complex use cases live in [Examples](examples.md).
