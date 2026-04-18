# Magesty.jl

```@meta
CurrentModule = Magesty
```

Magesty.jl (MAGnetic model ESTimator) is a Julia package for construction of Spin-Cluster Expansion (SCE) [1] model. It provides comprehensive tools for:

- Construction of symmetry-adapted SCE basis set
- Derivation of SCE coefficients
- Tools for spin configuration sampling

## Installation

```julia
using Pkg
Pkg.add("https://github.com/Tomonori-Tanaka/Magesty.jl")
```

## Documentation

| Page | Description |
|------|-------------|
| [Tutorial](tutorial.md) | Step-by-step guide for first-time users |
| [Examples](examples.md) | Code examples for common tasks |
| [Input Keys](input_keys.md) | Full reference for TOML configuration keys |
| [API Reference](api.md) | Detailed documentation of all exported functions and types |
| [Tools](tools.md) | Utility scripts in the `tools/` directory |
| [Technical Notes](technical_notes.md) | Theory behind the SCE formalism |
| [Tips](tips/index.md) | Practical tips (RWIGS dependence, penalty parameter tuning) |

## Tools

Magesty.jl includes a comprehensive set of utility tools in the `tools/` directory for:

- **Data Processing**: Convert between different file formats (VASP, TOML, XML)
- **Analysis**: Compare energies, perform cross-validation, analyze magnetic moments
- **Visualization**: Create scatter plots for energy/torque fit quality
- **Sampling**: Generate spin configurations using Mean-Field Approximation
- **Advanced Analysis**: Calculate micromagnetics parameters

See the [Tools](tools.md) page for detailed documentation of all available tools.

## Citation

If you use Magesty.jl in your research, please cite:

> Tomonori Tanaka and Yoshihiro Gohda, "General spin models from noncollinear spin density functional theory and spin-cluster expansion", [arXiv:2512.04458](https://arxiv.org/abs/2512.04458)

## References

1. R. Drautz and M. Fähnle, "Spin-cluster expansion: Parametrization of the general adiabatic magnetic energy surface with ab initio accuracy", *Phys. Rev. B* **69**, 104404 (2004). DOI: [10.1103/PhysRevB.69.104404](https://doi.org/10.1103/PhysRevB.69.104404)
