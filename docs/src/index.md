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

# Load configuration from a TOML file
sc = SpinCluster("config.toml")

# Calculate energy for a spin configuration
spin_config = rand(3, sc.structure.supercell.num_atoms)
energy = calc_energy(sc, spin_config)

# Get optimization results
j0 = get_j0(sc)  # Reference energy
jphi = get_jphi(sc)  # Spin-cluster coefficients

# Write results to files
write_xml(sc, "results.xml")
write_energies(sc, "energy_list.txt")
write_torques(sc, "torque_list.txt")
```

### Creating a System

```julia
# From a TOML configuration file
system = System("input.toml")

# From a dictionary
input_dict = Dict(
    "structure" => Dict(
        "lattice" => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "atoms" => [["Fe", [0.0, 0.0, 0.0]]]
    )
)
system = System(input_dict)
```

## Main Components

### System
A collection of structure, symmetry, cluster, and basis set information.

### SpinCluster
An extension of System with optimization capabilities for spin configurations.

### Submodules
- `Structures`: Crystal structure processing
- `Symmetries`: Symmetry operations processing  
- `Clusters`: Cluster expansion processing
- `BasisSets`: Basis function generation
- `Optimize`: Spin configuration optimization

## Features

- **Magnetic Structure Analysis**: Process and analyze magnetic crystal structures
- **Symmetry Operations**: Apply and analyze symmetry operations on magnetic structures
- **Cluster Expansion**: Perform spin cluster expansion calculations
- **Basis Functions**: Generate symmetry-adapted basis functions
- **Optimization**: Optimize spin configurations using various algorithms
- **XML Output**: Export results in SCE (Spin Cluster Expansion) format

## Examples

See the [API Reference](@ref) for detailed documentation of all functions and types.

## Tools

Magesty.jl includes a comprehensive set of utility tools in the `tools/` directory for:

- **Data Processing**: Convert between different file formats (VASP, TOML, XML)
- **Analysis**: Compare energies, perform cross-validation, analyze magnetic moments
- **Visualization**: Create plots for Bader analysis, charge density, and magnetic structures
- **Sampling**: Generate spin configurations using various statistical methods
- **Advanced Analysis**: Calculate micromagnetics parameters and spin wave dispersion

See the [Tools](tools.md) page for detailed documentation of all available tools.