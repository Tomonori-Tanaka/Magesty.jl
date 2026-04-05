# Examples

This page provides examples of using Magesty.jl for various magnetic structure analysis tasks.

## Example 1: BCC Iron — Full Workflow

Load a TOML configuration, fit the SCE model, and export results.

```julia
using Magesty

sc = SpinCluster("input.toml")

write_xml(sc, "bcc_fe_results.xml")
Magesty.write_energies(sc, "energy_list.txt")
Magesty.write_torques(sc, "torque_list.txt")

j0   = Magesty.get_j0(sc)
jphi = Magesty.get_jphi(sc)
println("Reference energy: ", j0, " eV")
println("Number of SCE coefficients: ", length(jphi))
```

## Example 2: Programmatic Fitting with `fit_sce_model`

Use the lower-level API for more control over the fitting process.

```julia
using Magesty, TOML

input  = TOML.parsefile("input.toml")
system = build_sce_basis(input)

spinconfigs = read_embset("EMBSET.dat")

# OLS fit
optimizer = fit_sce_model(system, spinconfigs)
println("RMSE energy: ", optimizer.metrics[:rmse_energy] * 1000, " meV")
println("RMSE torque: ", optimizer.metrics[:rmse_torque] * 1000, " meV")

# Elastic-Net fit
estimator = ElasticNet(lambda = 1e-4)
optimizer_reg = fit_sce_model(system, spinconfigs, estimator)
```

## Example 3: Ferromagnetic vs Antiferromagnetic Energies

```julia
using Magesty, LinearAlgebra

sc        = SpinCluster("input.toml")
num_atoms = sc.structure.supercell.num_atoms

# Ferromagnetic: all spins along +z
fm_config       = zeros(3, num_atoms)
fm_config[3, :] .= 1.0

# Antiferromagnetic: alternating +z/-z
afm_config = zeros(3, num_atoms)
for i in 1:num_atoms
    afm_config[3, i] = iseven(i) ? -1.0 : 1.0
end

fm_energy  = Magesty.calc_energy(sc, fm_config)
afm_energy = Magesty.calc_energy(sc, afm_config)

println("FM  energy: ", fm_energy,  " eV")
println("AFM energy: ", afm_energy, " eV")
println("ΔE (AFM-FM): ", (afm_energy - fm_energy) * 1000, " meV")
```

## Example 4: Symmetry Analysis

```julia
using Magesty

system = build_sce_basis(TOML.parsefile("input.toml"))
sym    = system.symmetry

println("Space group: ", sym.international_symbol, " (#", sym.spacegroup_number, ")")
println("Number of symmetry operations: ", sym.nsym)
println("Number of pure translations: ", sym.ntran)
println("Atoms in primitive cell: ", sym.atoms_in_prim)

for (i, symop) in enumerate(sym.symdata)
    println("Operation $i: proper=", symop.is_proper,
            "  translation=", symop.translation_frac)
end
```

## Example 5: Structure Information

```julia
using Magesty

system = build_sce_basis(TOML.parsefile("input.toml"))
cell   = system.structure.supercell

println("Number of atoms: ", cell.num_atoms)
println("Lattice vectors (columns):")
display(Matrix(cell.lattice_vectors))

println("Atomic positions (fractional):")
for i in 1:cell.num_atoms
    elem = system.structure.kd_name[cell.kd_int_list[i]]
    println("  $elem: ", cell.x_frac[:, i])
end
```

## Example 6: Basis Set Information

```julia
using Magesty

system = build_sce_basis(TOML.parsefile("input.toml"))
bs     = system.basisset

println("Number of SALCs: ", length(bs.salc_list))
println("Number of coupled basis functions: ", length(bs.coupled_basislist))
```

## Example 7: Custom Spin Configurations

Evaluate energy and torque for arbitrary spin configurations.

```julia
using Magesty, LinearAlgebra

sc        = SpinCluster("input.toml")
num_atoms = sc.structure.supercell.num_atoms

# Spiral configuration in the xy-plane
spiral = zeros(3, num_atoms)
for i in 1:num_atoms
    angle = 2π * (i - 1) / num_atoms
    spiral[1, i] = cos(angle)
    spiral[2, i] = sin(angle)
end

energy = Magesty.calc_energy(sc, spiral)
torque = Magesty.calc_torque(sc, spiral)
println("Spiral energy: ", energy, " eV")
println("Max torque magnitude: ", maximum(norm.(eachcol(torque))), " eV")
```

## Example 8: Batch Energy Evaluation

```julia
using Magesty, LinearAlgebra, Statistics

sc        = SpinCluster("input.toml")
num_atoms = sc.structure.supercell.num_atoms

n_configs = 100
energies  = Vector{Float64}(undef, n_configs)

for i in 1:n_configs
    config = randn(3, num_atoms)
    for j in 1:num_atoms
        config[:, j] ./= norm(config[:, j])
    end
    energies[i] = Magesty.calc_energy(sc, config)
end

println("Min energy: ",  minimum(energies) * 1000, " meV")
println("Max energy: ",  maximum(energies) * 1000, " meV")
println("Mean energy: ", mean(energies)    * 1000, " meV")
println("Std dev: ",     std(energies)     * 1000, " meV")
```

## Example 9: Load Basis from XML and Refit

Reuse a previously computed basis set to refit with different training data or
regularization without recomputing SALCs:

```julia
using Magesty, TOML

input  = TOML.parsefile("input.toml")
system = build_sce_basis_from_xml(input, "scecoeffs.xml")

spinconfigs = read_embset("new_EMBSET.dat")
estimator   = ElasticNet(lambda = 1e-3)
optimizer   = fit_sce_model(system, spinconfigs, estimator, 0.5)

write_xml(system.structure, system.symmetry, system.basisset, optimizer, "new_results.xml")
```
