# Magesty.jl

Magesty is a Julia package for constructing effective spin models in magnetic materials within the Spin-Cluster Expansion (SCE) formalism [1].

## Features

- Generation of SCE models considering crystal symmetry
- Optimization of SCE coefficients based on energy and local magnetic field
- High-performance parallel computation
- Rich analysis tools and visualization capabilities
- Automatic generation of symmetry-adapted linear combinations (SALCs)
- Energy calculation for spin configurations
- XML format output for results

## Installation

You can install the package using Julia's package manager. If the package is not yet registered, add it by URL:

```julia
using Pkg
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl")
```

## Usage

```julia
using Magesty

# Quick start from a TOML configuration
sc = SpinCluster("input.toml")

# Optionally compute energy/torque for a given spin configuration
# spins must be a 3 x N matrix (N: number of atoms in supercell)
# energy = calc_energy(sc, spins)
# torque  = calc_torque(sc, spins)

# Write results
write_xml(sc, "jphi.xml")
write_energies(sc, "energy_list.txt")
write_torques(sc, "torque_list.txt")
```

For detailed examples, please refer to the `test/examples` directory. Additional runnable inputs are also under `test/develop_tmp`.

## Main Types and Functions

### Main Types
- `System`: Collection of structure, symmetry, cluster, and basis set
- `SpinCluster`: Extension of System with optimization capabilities
- `Structure`: Crystal structure information
- `Symmetry`: Symmetry operations
- `Cluster`: Cluster information
- `BasisSet`: Basis set information
- `Optimizer`: Optimizer instance

### Main Functions
- `System(input_dict | toml_file)`: Create a system
- `SpinCluster(input_dict | toml_file | system[, input_dict])`: Create a spin cluster
- `calc_energy(spin_cluster, spin_config)`: Calculate energy for spin configuration
- `calc_torque(spin_cluster, spin_config)`: Calculate torques for spin configuration
- `write_xml(spin_cluster, filename="jphi.xml"; write_jphi=true)`: Output structure, basis, optimizer in XML
- `write_energies(spin_cluster, filename="energy_list.txt")`: Output energy lists
- `write_torques(spin_cluster, filename="torque_list.txt")`: Output torque lists
- `get_j0(spin_cluster)`: Get reference energy
- `get_jphi(spin_cluster)`: Get SCE coefficients

## Core Modules

- `Structures`: Crystal structure definition and processing
- `Symmetries`: Space group analysis and symmetry operation generation
- `Clusters`: Cluster generation and analysis
- `BasisSets`: Generation of symmetry-adapted SCE basis functions
- `Optimize`: Estimation of SCE coefficients
- `SpinConfigs`: Spin configuration management
- `SALCs`: Symmetry-adapted linear combinations

## Dependencies

- Julia 1.11 or higher
- TOML.jl: Configuration file parsing
- Spglib.jl: Crystal symmetry analysis
- StaticArrays.jl: Static arrays
- LinearAlgebra.jl, SparseArrays.jl, Statistics.jl, Random.jl, Dates.jl, Printf.jl
- Combinatorics.jl, DataStructures.jl
- EzXML.jl: XML processing
- LegendrePolynomials.jl, WignerD.jl, Rotations.jl: Angular functions and rotations
- MultivariateStats.jl
- ArgParse.jl, JLD2.jl (tools and I/O utilities)

For exact versions, see `Project.toml`.

## Input File Format

The package uses TOML format for input files. A typical input file includes:

### General Settings
```toml
[general]
name = "system_name"
kd = ["Fe", "Pt"]  # Element names
nat = 16           # Number of atoms
periodicity = [true, true, true]  # Periodicity
```

### Symmetry Settings
```toml
[symmetry]
tolerance = 1e-5   # Tolerance for symmetry determination
```

### Interaction Settings
```toml
[interaction]
nbody = 2          # Maximum order of many-body interactions
[interaction.lmax]
Fe = [0, 1]        # Maximum angular momentum for each element
[interaction.cutoff]
Fe-Fe = [-1, -1]   # Cutoff distances (in Bohr units)
```

### Regression Settings
```toml
[regression]
datafile = "EMBSET.dat"  # Data file
weight = 1.0             # Weight
alpha = 0                # Regularization parameter
lambda = 0.0             # Regularization parameter
```

### Structure Settings
```toml
[structure]
kd_list = [1, 1, 1, 1, ...]  # List of atomic species
lattice = [                   # Lattice vectors
  [5.66, 0.0, 0.0],
  [0.0, 5.66, 0.0],
  [0.0, 0.0, 5.66],
]
position = [                  # Atomic positions (fractional coordinates)
  [0.00, 0.00, 0.00],
  [0.25, 0.25, 0.25],
  # ...
]
```

Example input files can be found in the `test/examples` directory.

## Analysis Tools

The `tools` directory contains the following analysis utilities (non-exhaustive):

- `CrossValidation.jl`: Cross-validation
- `FitCheck_energy.jl`, `FitCheck_torque.jl`: Fit quality checks
- `micromagnetics.jl`: Micromagnetic analysis
- `sampling_mfa.jl`: Sampling using mean field approximation
- `pos2toml.jl`: Generate TOML configuration files from position files
- `convert2tensor.jl`: Tensor conversion
- `extract.jl`: Data extraction helpers
- `histogram_magmom.jl`: Magnetic moment histograms
- `vasptools.jl`: VASP-related tools

Personal scripts are under `tools/personal/`.

## Performance

The package is optimized for:
- Parallel computation using multiple threads
- Optimized cluster generation
- Fast computation using static arrays
- Efficient basis function generation utilizing symmetry

## License

This project is licensed under the MIT License.

## Author

- Tomonori Tanaka ([@TomonoriTanaka](https://github.com/Tomonori-Tanaka)) ([@ORCID](https://orcid.org/0000-0001-7306-6770))

## Citation

If you use this software in your research, please cite:

[Citation information to be added]

## Contributing

We welcome contributions! Please feel free to:
- Report bugs and request features through [Issues](https://github.com/Tomonori-Tanaka/Magesty.jl/issues)
- Submit pull requests
- Improve documentation
- Share usage examples

## Documentation

Detailed documentation is available at `docs/build/index.html`.

## Support

For questions and support:
1. Check the documentation
2. Search existing issues
3. Open a new issue if needed

## References

1. R. Drautz and M. Fähnle, "Spin-cluster expansion: Parametrization of the general adiabatic magnetic energy surface with ab initio accuracy", Phys. Rev. B 69, 104404 (2004). DOI: [10.1103/PhysRevB.69.104404](https://doi.org/10.1103/PhysRevB.69.104404)



