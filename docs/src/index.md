# Magesty.jl

```@meta
CurrentModule = Magesty
```

Magesty.jl (MAGnetic model ESTimator) is a Julia package for construction of Spin-Cluster Expansion (SCE) model. It provides comprehensive tools for:

- Construction of symmetry-adapted SCE basis set
- Derivation of SCE coefficients
- Tools for spin configuration sampling

## Installation

```julia
using Pkg
Pkg.add("https://github.com/Tomonori-Tanaka/Magesty.jl")
```

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
