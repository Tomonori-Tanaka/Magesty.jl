# Magesty.jl

```@meta
CurrentModule = Magesty
```

Magesty.jl is a Julia package for magnetic structure analysis and spin cluster expansion calculations. It provides comprehensive tools for:

- Magnetic structure setup and management
- Symmetry analysis and operations
- Cluster expansion calculations
- Basis function generation
- Spin configuration optimization

## Quick Start

### Installation

```julia
using Pkg
Pkg.add("Magesty")
```

### Basic Usage

```julia
using Magesty

# Fit SCE model from a TOML configuration file
sc = SpinCluster("config.toml")

# Write results to XML
write_xml(sc, "results.xml")

# Write energy / torque comparison files
Magesty.write_energies(sc, "energy_list.txt")
Magesty.write_torques(sc, "torque_list.txt")

# Get fitted coefficients
j0   = Magesty.get_j0(sc)    # reference energy (eV)
jphi = Magesty.get_jphi(sc)  # SCE coefficients
```

### Programmatic Fitting with `fit_sce_model`

```julia
using Magesty, TOML

input  = TOML.parsefile("config.toml")
system = build_sce_basis(input)

spinconfigs = read_embset("EMBSET.dat")
optimizer   = fit_sce_model(system, spinconfigs)   # OLS by default

j0   = optimizer.reference_energy
jphi = optimizer.SCE

# Export
write_xml(system.structure, system.symmetry, system.basisset, optimizer, "results.xml")
```

### Elastic-Net Regularization

```julia
estimator = ElasticNet(lambda = 1e-4)
optimizer = fit_sce_model(system, spinconfigs, estimator)
```

### Load Basis Set from XML

If the basis set has already been computed and saved, skip the expensive SALC
construction:

```julia
using Magesty, TOML

input  = TOML.parsefile("config.toml")
system = build_sce_basis_from_xml(input, "scecoeffs.xml")
```

### Configuration File Format

Magesty.jl uses TOML files with four required sections:

```toml
[general]
name = "bccfe"
kd   = ["Fe"]          # element names
nat  = 16              # total number of atoms in supercell
periodicity = [true, true, true]

[symmetry]
tolerance = 1e-5

[interaction]
nbody = 2
[interaction.body1]
lmax.Fe = 0            # on-site angular momentum per element
[interaction.body2]
lsum = 2               # max L for two-body basis
cutoff."Fe-Fe" = 5.66  # cutoff radius in bohr (-1 = all pairs)

[regression]
datafile = "EMBSET.dat"
weight   = 1.0         # 0 = torque only, 1 = energy only
alpha    = 0.0
lambda   = 0.0

[structure]
kd_list  = [1, 1, 1, 1, ...]   # element index per atom (1-based)
lattice  = [[5.66, 0.0, 0.0], [0.0, 5.66, 0.0], [0.0, 0.0, 5.66]]
position = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], ...]  # fractional
```

## Main Components

### System
A collection of structure, symmetry, cluster, and basis set information.

### SpinCluster
An extension of `System` with fitted SCE coefficients for energy/torque evaluation.

### Submodules
- `Structures`: Crystal structure processing
- `Symmetries`: Symmetry operations processing
- `Clusters`: Cluster expansion processing
- `BasisSets`: Basis function generation
- `Optimize`: SCE coefficient fitting

## Features

- **Magnetic Structure Analysis**: Process and analyze magnetic crystal structures
- **Symmetry Operations**: Apply and analyze symmetry operations on magnetic structures
- **Cluster Expansion**: Perform spin cluster expansion calculations
- **Basis Functions**: Generate symmetry-adapted linear combinations (SALCs)
- **Optimization**: Fit SCE coefficients using OLS or Elastic-Net regression
- **XML Output**: Export results in SCE (Spin Cluster Expansion) format

## Examples

See the [API Reference](@ref) for detailed documentation of all functions and types.

## Tools

Magesty.jl includes a comprehensive set of utility tools in the `tools/` directory for:

- **Data Processing**: Convert between different file formats (VASP, TOML, XML)
- **Analysis**: Compare energies, perform cross-validation, analyze magnetic moments
- **Visualization**: Create scatter plots for energy/torque fit quality
- **Sampling**: Generate spin configurations using Mean-Field Approximation
- **Advanced Analysis**: Calculate micromagnetics parameters

See the [Tools](tools.md) page for detailed documentation of all available tools.
