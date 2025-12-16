# Magesty.jl

Julia package for constructing effective spin models in magnetic materials using Spin-Cluster Expansion (SCE) formalism.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl")
```

## Quick Start

```julia
using Magesty
using TOML

# Build SCE basis set
input = TOML.parsefile("input.toml")
system = build_sce_basis(input, verbosity = false)

# Read spin configurations and fit coefficients
spinconfig_list = read_embset("EMBSET.dat")
optimizer = fit_sce_model(system, spinconfig_list)

# Calculate energy and torque
using Magesty.EnergyTorque
num_atoms = system.structure.supercell.num_atoms
spin_config = rand(3, num_atoms)
for i in 1:num_atoms
    spin_config[:, i] ./= norm(spin_config[:, i])
end

energy = EnergyTorque.calc_energy(
    system.basisset.salc_list,
    spin_config,
    system.symmetry,
    optimizer
)

# Access coefficients
j0 = optimizer.reference_energy
jphi = optimizer.SCE

# Write results
write_xml(system.structure, system.symmetry, system.basisset, optimizer, "jphi.xml")
```

## Main Functions

- `build_sce_basis(input_dict; verbosity=false)`: Build SCE basis set (returns `System`)
- `fit_sce_model(system, spinconfig_list, estimator=OLS(), weight=0.5)`: Fit SCE coefficients (returns `Optimizer`)
- `read_embset(filename)`: Read spin configuration data from EMBSET format file
- `write_xml(structure, symmetry, basis_set, optimizer, filename)`: Write results to XML file

## Examples

See `test/examples` directory for complete examples.

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



