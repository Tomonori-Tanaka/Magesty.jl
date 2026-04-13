# Tutorial

This tutorial will guide you through the basic usage of Magesty.jl for magnetic structure analysis and spin cluster expansion calculations.

## Configuration File Format

Magesty.jl uses TOML configuration files. Below is an annotated example for a BCC Fe supercell:

```toml:input.toml
[general]
name = "bccfe"
kd   = ["Fe"]           # list of element names
nat  = 16               # total number of atoms in the supercell
periodicity = [true, true, true]  # apply periodic boundary to all (x, y, z) directions

[symmetry]
tolerance = 1e-5        # symmetry detection tolerance (optional, default 1e-3)
isotropy = true

[interaction]
nbody = 2               # maximum interaction body
[interaction.body1]
lmax.Fe = 0             # 1-body maximum angular momentum per element, i.e. this represents on-site anisotropy
[interaction.body2]
lsum = 2                # cutoff summation of l values for basis functions
cutoff."Fe-Fe" = -1     # pairwise cutoff radius in Å (-1 uses all possible pairs)

[regression]
datafile = "EMBSET" # path to training data
weight   = 0.5          # 0 = torque only, 1 = energy only, 0.5 = balanced
alpha    = 0.0          # elastic-net mixing (0 = ridge)
lambda   = 0.0          # regularization strength (0 = no regularization)

[structure]
kd_list  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # element index per atom
lattice  = [            # lattice parameters
  [5.66, 0.0, 0.0],
  [0.0, 5.66, 0.0],
  [0.0, 0.0, 5.66],
]
position = [            # fractional coordinates
  [0.00, 0.00, 0.00],
  [0.25, 0.25, 0.25],
  [0.00, 0.00, 0.50],
  [0.50, 0.00, 0.00],
  [0.00, 0.50, 0.00],
  [0.25, 0.25, 0.75],
  [0.75, 0.25, 0.25],
  [0.25, 0.75, 0.25],
  [0.00, 0.50, 0.50],
  [0.50, 0.00, 0.50],
  [0.50, 0.50, 0.00],
  [0.25, 0.75, 0.75],
  [0.75, 0.25, 0.75],
  [0.75, 0.75, 0.25],
  [0.50, 0.50, 0.50],
  [0.75, 0.75, 0.75],
]
```

For a full key reference see [Input Keys](input_keys.md).

## Basic Workflow

### 1. Creating a SpinCluster (All-in-One)

The simplest approach reads the TOML file, builds the basis, loads the training set, and fits the SCE model in one call:

```julia
using Magesty

sc = SpinCluster("config.toml")
```

### 2. Accessing Fitted Results

```julia
# Reference energy (eV)
j0 = Magesty.get_j0(sc)

# SCE coefficients
jphi = Magesty.get_jphi(sc)

# Both at once
j0, jphi = Magesty.get_j0_jphi(sc)
```

### 3. Exporting Results

```julia
# Write SCE coefficients to XML
# This XML file can be used for the Monte Carlo package, `SpinClusterMC.jl`.
write_xml(sc, "results.xml")

# Write energy and torque comparison files
Magesty.write_energies(sc, "energy_list.txt")
Magesty.write_torques(sc, "torque_list.txt")
```

## Programmatic Fitting with `fit_sce_model`

For more control—e.g., cross-validation or custom training sets—use the
lower-level `fit_sce_model` API:

```julia
using Magesty, TOML

input  = TOML.parsefile("config.toml")
system = build_sce_basis(input)

# Load training data
spinconfigs = read_embset("EMBSET.dat")

# OLS fit
optimizer = fit_sce_model(system, spinconfigs)

# Elastic-Net fit
estimator = ElasticNet(lambda = 1e-4)
optimizer = fit_sce_model(system, spinconfigs, estimator, 0.5)

j0   = optimizer.reference_energy
jphi = optimizer.SCE
println("RMSE energy: ", optimizer.metrics[:rmse_energy] * 1000, " meV")
println("RMSE torque: ", optimizer.metrics[:rmse_torque] * 1000, " meV")
```

### Creating a SpinCluster from an Existing System

```julia
# Reuse a System built above
sc = SpinCluster(system, input)
```

### Creating a SpinCluster with a Pre-loaded Training Set

```julia
spinconfigs = read_embset("EMBSET.dat")
sc = SpinCluster(system, input, spinconfigs)
```

## Building Basis from XML

If the basis set has already been saved to XML, load it directly to skip the
expensive SALC computation:

```julia
using Magesty, TOML

input  = TOML.parsefile("config.toml")
system = build_sce_basis_from_xml(input, "scecoeffs.xml")
```

## Advanced Usage

### Symmetry Information

```julia
sym = sc.symmetry
println("Space group: ", sym.international_symbol, " (#", sym.spacegroup_number, ")")
println("Number of symmetry operations: ", sym.nsym)

for (i, symop) in enumerate(sym.symdata)
    println("Operation $i:")
    println("  Rotation (fractional): ", symop.rotation_frac)
    println("  Translation (fractional): ", symop.translation_frac)
    println("  Is proper: ", symop.is_proper)
end
```

### Structure Information

```julia
cell = sc.structure.supercell
println("Number of atoms: ", cell.num_atoms)
println("Lattice vectors:")
for i in 1:3
    println("  a$i: ", cell.lattice_vectors[:, i])
end
println("Atomic positions (fractional):")
for i in 1:cell.num_atoms
    elem = sc.structure.kd_name[cell.kd_int_list[i]]
    println("  $elem: ", cell.x_frac[:, i])
end
```

### Basis Set Information

```julia
basisset = sc.basisset
println("Number of SALCs: ", length(basisset.salc_list))
```

For detailed function documentation see the [API Reference](@ref).
For more complex use cases see [Examples](examples.md).
