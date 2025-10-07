# Tutorial

This tutorial will guide you through the basic usage of Magesty.jl for magnetic structure analysis and spin cluster expansion calculations.

## Installation

First, install Magesty.jl:

```julia
using Pkg
Pkg.add("Magesty")
```

## Basic Workflow

### 1. Setting up a Configuration

Magesty.jl uses TOML configuration files to specify input parameters. Here's a basic example:

```toml
[structure]
lattice = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
atoms = [["Fe", [0.0, 0.0, 0.0]], ["Fe", [0.5, 0.5, 0.5]]]

[optimization]
method = "LBFGS"
max_iterations = 1000
tolerance = 1e-6
```

### 2. Creating a System

```julia
using Magesty

# Create a system from a TOML file
system = System("config.toml")

# Or create from a dictionary
input_dict = Dict(
    "structure" => Dict(
        "lattice" => [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
        "atoms" => [["Fe", [0.0, 0.0, 0.0]], ["Fe", [0.5, 0.5, 0.5]]]
    )
)
system = System(input_dict)
```

### 3. Creating a SpinCluster for Optimization

```julia
# Create a SpinCluster for optimization
sc = SpinCluster("config.toml")

# Or create from an existing System
sc = SpinCluster(system, input_dict)
```

### 4. Calculating Energies

```julia
# Generate a random spin configuration
num_atoms = sc.structure.supercell.num_atoms
spin_config = rand(3, num_atoms)

# Normalize spins to unit length
for i in 1:num_atoms
    spin_config[:, i] ./= norm(spin_config[:, i])
end

# Calculate energy
energy = calc_energy(sc, spin_config)
println("Energy: ", energy)
```

### 5. Accessing Results

```julia
# Get reference energy
j0 = get_j0(sc)

# Get spin-cluster coefficients
jphi = get_jphi(sc)

# Get both
j0, jphi = get_j0_jphi(sc)
```

### 6. Exporting Results

```julia
# Write results to XML file
write_xml(sc, "results.xml")

# Write without J_ij parameters
write_xml(sc, "structure_only.xml", write_jphi=false)
```

## Advanced Usage

### Working with Spin Configurations

```julia
# Create specific spin configurations
spin_config = zeros(3, num_atoms)
spin_config[1, :] .= 1.0  # All spins along x-axis

# Calculate energy and torque
energy = calc_energy(sc, spin_config)
torque = calc_torque(sc, spin_config)
```

### Symmetry Analysis

```julia
# Access symmetry information
symmetry = sc.symmetry
println("Number of symmetry operations: ", length(symmetry.symmetry_operations))

# Access structure information
structure = sc.structure
println("Number of atoms: ", structure.supercell.num_atoms)
println("Lattice vectors: ", structure.supercell.lattice_vectors)
```

### Cluster Information

```julia
# Access cluster information
cluster = sc.cluster
println("Number of clusters: ", length(cluster.interaction_clusters))

# Access basis set information
basisset = sc.basisset
println("Number of basis functions: ", length(basisset.salc_list))
```

## Configuration Options

### Structure Parameters

- `lattice`: 3x3 matrix of lattice vectors
- `atoms`: List of [element, position] pairs
- `supercell`: Supercell dimensions (optional)

### Optimization Parameters

- `method`: Optimization method ("LBFGS", "CG", etc.)
- `max_iterations`: Maximum number of iterations
- `tolerance`: Convergence tolerance
- `alpha`: Regularization parameter
- `lambda`: Weight parameter

### Symmetry Parameters

- `tolerance`: Symmetry detection tolerance
- `use_spglib`: Whether to use Spglib for symmetry detection

## Troubleshooting

### Common Issues

1. **File not found**: Make sure the TOML file path is correct
2. **Invalid configuration**: Check that all required parameters are present
3. **Memory issues**: For large systems, consider reducing the cluster cutoff radius
4. **Convergence problems**: Try adjusting the tolerance or maximum iterations

### Getting Help

- Check the [API Reference](@ref) for detailed function documentation
- Look at the [Examples](@ref) for more complex use cases
- Report issues on the GitHub repository
