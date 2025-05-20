# Magesty.jl

Magesty is a Julia package for constructing effective spin models in magnetic materials within the Spin-Cluster Expansion (SCE) formalism [1].

## Features

- Generation of SCE model considering crystal symmetry
- Optimization of SCE coefficients based on energy and local magnetic field 
- High-performance parallel computation
- Rich analysis tools and visualization capabilities

## Installation

You can install the package using Julia's package manager:

```julia
using Pkg
Pkg.add("Magesty")
```

## Usage

Here's a basic example:

```julia
using Magesty
using TOML

# Load system data
input = TOML.parse("input.tom")
system = System(input)

# estimate SCE coefficients
sce = SpinCluster(system, input)
```

For more detailed examples, please refer to the `test/examples` directory.

## Core Modules

- `Structure`: Definition of crystal structures
- `Symmetries`: Space group analysis and symmetry operation generation
- `Cluster`: Cluster generation and analysis
- `BasisSet`: Generation of symmetry-datapted SCE basis functions
- `Optimize`: Estimation of SCE coefficients

## Dependencies

- Julia 1.11 or higher
- Spglib.jl: Crystal symmetry analysis
- Optim.jl: Optimization algorithms
- For other dependencies, see `Project.toml`

## Input File Format

The package uses TOML format for input files. A typical input file includes:
- Crystal structure parameters
- Interaction parameters
- Symmetry tolerance settings
- Optimization parameters

Example input files can be found in the `test/examples` directory.

## Performance

The package is optimized for:
- Parallel computation using multiple threads
- Optimized cluster generation

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

Detailed documentation is available at [Documentation Link].

## Support

For questions and support:
1. Check the documentation
2. Search existing issues
3. Open a new issue if needed

## References

1. R. Drautz and M. Fähnle, "Spin-cluster expansion: Parametrization of the general adiabatic magnetic energy surface with ab initio accuracy", Phys. Rev. B 69, 104404 (2004). DOI: [10.1103/PhysRevB.69.104404](https://doi.org/10.1103/PhysRevB.69.104404)



