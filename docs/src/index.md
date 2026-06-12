# Magesty.jl

```@meta
CurrentModule = Magesty
```

Magesty.jl (MAGnetic model ESTimator) is a Julia package for construction of Spin-Cluster Expansion (SCE) [1] model. It provides comprehensive tools for:

- Construction of symmetry-adapted SCE basis set
- Derivation of SCE coefficients
- Tools for spin configuration sampling

## Documentation

| Page | Description |
|------|-------------|
| [Installation](installation.md) | Package installation and how to run the CLI tools |
| [Tutorial](tutorial.md) | Step-by-step guide for first-time users |
| [Examples](examples/index.md) | Worked examples applying the SCE workflow to specific systems |
| [API Examples](api_examples.md) | Short, feature-by-feature snippets of the public API |
| [Input Keys](input_keys.md) | Full reference for TOML configuration keys |
| [API Reference](api.md) | Detailed documentation of all exported functions and types |
| [Internal API](api_internal.md) | Lower-level building blocks, not covered by the stability guarantee |
| [Tools](tools.md) | The `magesty` command-line interface |
| [Theoretical Background](theory/index.md) | The SCE theory, end to end, linked to the implementation |
| [Tips](tips/index.md) | Practical tips|

## Tools

Magesty.jl provides a `magesty` command-line interface for converting VASP
output into the formats used by the spin-cluster expansion workflow:

- `magesty vasp extxyz` — convert a VASP run to extended XYZ
- `magesty vasp toml` — convert a VASP POSCAR to a Magesty input TOML configuration
- `magesty vasp embset` — convert VASP OSZICAR files to the EMBSET training-data format

See the [Tools](tools.md) page for details.

## Citation

If you use Magesty.jl in your research, please cite:

> Tomonori Tanaka and Yoshihiro Gohda, "General spin models from noncollinear spin density functional theory and spin-cluster expansion", [arXiv:2512.04458](https://arxiv.org/abs/2512.04458)

## References

1. R. Drautz and M. Fähnle, "Spin-cluster expansion: Parametrization of the general adiabatic magnetic energy surface with ab initio accuracy", *Phys. Rev. B* **69**, 104404 (2004). DOI: [10.1103/PhysRevB.69.104404](https://doi.org/10.1103/PhysRevB.69.104404)
