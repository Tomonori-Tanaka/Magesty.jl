# Examples

This page provides comprehensive examples of using Magesty.jl for various magnetic structure analysis tasks.

## Example 1: Simple BCC Iron

This example demonstrates the basic workflow for a simple BCC iron structure.

```julia
using Magesty

# Define input parameters
input_dict = Dict(
    "structure" => Dict(
        "lattice" => [[2.87, 0.0, 0.0], [0.0, 2.87, 0.0], [0.0, 0.0, 2.87]],
        "atoms" => [["Fe", [0.0, 0.0, 0.0]]]
    ),
    "optimization" => Dict(
        "method" => "LBFGS",
        "max_iterations" => 1000,
        "tolerance" => 1e-6
    )
)

# Create system and spin cluster
system = System(input_dict)
sc = SpinCluster(system, input_dict)

# Generate random spin configuration
num_atoms = sc.structure.supercell.num_atoms
spin_config = rand(3, num_atoms)
for i in 1:num_atoms
    spin_config[:, i] ./= norm(spin_config[:, i])
end

# Calculate energy
energy = calc_energy(sc, spin_config)
println("Energy: ", energy)

# Export results
write_xml(sc, "bcc_fe_results.xml")
```

## Example 2: Ferromagnetic vs Antiferromagnetic

This example compares ferromagnetic and antiferromagnetic configurations.

```julia
using Magesty

# Load configuration
sc = SpinCluster("config.toml")
num_atoms = sc.structure.supercell.num_atoms

# Ferromagnetic configuration (all spins aligned)
fm_config = ones(3, num_atoms)
for i in 1:num_atoms
    fm_config[:, i] ./= norm(fm_config[:, i])
end

# Antiferromagnetic configuration (alternating spins)
afm_config = zeros(3, num_atoms)
for i in 1:num_atoms
    if i % 2 == 1
        afm_config[1, i] = 1.0  # +x direction
    else
        afm_config[1, i] = -1.0  # -x direction
    end
end

# Calculate energies
fm_energy = calc_energy(sc, fm_config)
afm_energy = calc_energy(sc, afm_config)

println("Ferromagnetic energy: ", fm_energy)
println("Antiferromagnetic energy: ", afm_energy)
println("Energy difference: ", afm_energy - fm_energy)
```

## Example 3: Spin Wave Analysis

This example demonstrates how to analyze spin wave excitations.

```julia
using Magesty
using LinearAlgebra

# Create spin cluster
sc = SpinCluster("config.toml")
num_atoms = sc.structure.supercell.num_atoms

# Create ground state configuration
ground_state = zeros(3, num_atoms)
ground_state[1, :] .= 1.0  # All spins along x-axis

# Calculate ground state energy
ground_energy = calc_energy(sc, ground_state)

# Create small perturbation
perturbation = 0.1 * randn(3, num_atoms)
excited_state = ground_state + perturbation

# Normalize spins
for i in 1:num_atoms
    excited_state[:, i] ./= norm(excited_state[:, i])
end

# Calculate excited state energy
excited_energy = calc_energy(sc, excited_state)

# Calculate excitation energy
excitation_energy = excited_energy - ground_energy
println("Excitation energy: ", excitation_energy)
```

## Example 4: Symmetry Analysis

This example demonstrates how to analyze the symmetry properties of a magnetic structure.

```julia
using Magesty

# Create system
system = System("config.toml")

# Access symmetry information
symmetry = system.symmetry
println("Number of symmetry operations: ", length(symmetry.symmetry_operations))

# Print symmetry operations
for (i, symop) in enumerate(symmetry.symmetry_operations)
    println("Symmetry operation $i:")
    println("  Rotation (fractional): ", symop.rotation_frac)
    println("  Translation (fractional): ", symop.translation_frac)
    println("  Is proper: ", symop.is_proper)
    println()
end

# Access structure information
structure = system.structure
println("Lattice vectors:")
for i in 1:3
    println("  a$i: ", structure.supercell.lattice_vectors[i, :])
end

println("Atomic positions (fractional):")
for i in 1:structure.supercell.num_atoms
    element = structure.kd_name[structure.supercell.kd_int_list[i]]
    pos = structure.supercell.x_frac[:, i]
    println("  $element: ", pos)
end
```

## Example 5: Cluster Analysis

This example demonstrates how to analyze the cluster expansion.

```julia
using Magesty

# Create spin cluster
sc = SpinCluster("config.toml")

# Access cluster information
cluster = sc.cluster
println("Number of interaction clusters: ", length(cluster.interaction_clusters))

# Print cluster information
for (i, interaction_cluster) in enumerate(cluster.interaction_clusters)
    println("Cluster $i:")
    println("  Number of atoms: ", length(interaction_cluster.atom_indices))
    println("  Atom indices: ", interaction_cluster.atom_indices)
    println("  Maximum distance: ", interaction_cluster.distmax)
    println()
end

# Access basis set information
basisset = sc.basisset
println("Number of basis functions: ", length(basisset.salc_list))

# Print basis function information
for (i, salc) in enumerate(basisset.salc_list)
    println("Basis function $i:")
    println("  Cluster index: ", salc.cluster_index)
    println("  Angular momentum: ", salc.l)
    println("  Magnetic quantum number: ", salc.m)
    println()
end
```

## Example 6: Optimization Results

This example demonstrates how to access and analyze optimization results.

```julia
using Magesty

# Create spin cluster with optimization
sc = SpinCluster("config.toml")

# Get optimization results
j0 = get_j0(sc)
jphi = get_jphi(sc)

println("Reference energy (J0): ", j0)
println("Number of spin-cluster coefficients: ", length(jphi))
println("First 10 coefficients: ", jphi[1:min(10, length(jphi))])

# Analyze coefficient magnitudes
max_coeff = maximum(abs.(jphi))
min_coeff = minimum(abs.(jphi))
println("Maximum coefficient magnitude: ", max_coeff)
println("Minimum coefficient magnitude: ", min_coeff)

# Find most important coefficients
important_indices = findall(x -> abs(x) > 0.1 * max_coeff, jphi)
println("Important coefficient indices: ", important_indices)
```

## Example 7: Custom Spin Configurations

This example demonstrates how to work with custom spin configurations.

```julia
using Magesty
using LinearAlgebra

# Create spin cluster
sc = SpinCluster("config.toml")
num_atoms = sc.structure.supercell.num_atoms

# Create custom spin configuration
spin_config = zeros(3, num_atoms)

# Set spins in a specific pattern
for i in 1:num_atoms
    # Create a spiral pattern
    angle = 2π * i / num_atoms
    spin_config[1, i] = cos(angle)
    spin_config[2, i] = sin(angle)
    spin_config[3, i] = 0.0
end

# Calculate energy and torque
energy = calc_energy(sc, spin_config)
torque = calc_torque(sc, spin_config)

println("Spiral configuration energy: ", energy)
println("Maximum torque magnitude: ", maximum(norm.(eachcol(torque))))

# Create another configuration (random)
random_config = randn(3, num_atoms)
for i in 1:num_atoms
    random_config[:, i] ./= norm(random_config[:, i])
end

random_energy = calc_energy(sc, random_config)
println("Random configuration energy: ", random_energy)
```

## Example 8: Batch Processing

This example demonstrates how to process multiple configurations efficiently.

```julia
using Magesty

# Create spin cluster
sc = SpinCluster("config.toml")
num_atoms = sc.structure.supercell.num_atoms

# Generate multiple random configurations
n_configs = 100
energies = Float64[]
configurations = []

for i in 1:n_configs
    # Generate random configuration
    config = randn(3, num_atoms)
    for j in 1:num_atoms
        config[:, j] ./= norm(config[:, j])
    end
    
    # Calculate energy
    energy = calc_energy(sc, config)
    
    push!(energies, energy)
    push!(configurations, config)
end

# Analyze results
println("Number of configurations: ", length(energies))
println("Minimum energy: ", minimum(energies))
println("Maximum energy: ", maximum(energies))
println("Average energy: ", mean(energies))
println("Standard deviation: ", std(energies))

# Find lowest energy configuration
min_index = argmin(energies)
println("Lowest energy configuration index: ", min_index)
println("Lowest energy: ", energies[min_index])
```

These examples demonstrate the versatility and power of Magesty.jl for magnetic structure analysis and spin cluster expansion calculations. Each example can be adapted and extended for specific research needs.
