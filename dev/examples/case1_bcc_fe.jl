using Magesty

basis = SCEBasis("input.toml")

sym = basis.symmetry
println("Space group: ", sym.international_symbol, " (#", sym.spacegroup_number, ")")
println("Symmetry operations: ", sym.nsym)
println("Atoms in the supercell: ", basis.structure.supercell.num_atoms)
println("Number of SALCs: ", length(basis.salcbasis.salc_list))

dataset = SCEDataset(basis, "EMBSET"; verbosity = false)
f = fit(SCEFit, dataset, OLS(); verbosity = false)
println("Fitted ", length(coef(f)), " coefficients on ", length(dataset),
        " configurations.")

using Printf

@printf("RMSE energy : %.3f meV\n", rmse_energy(f) * 1000)
@printf("R^2  energy : %.4f\n", r2_energy(f))
@printf("RMSE torque : %.3f meV\n", rmse_torque(f) * 1000)
@printf("j0 (reference energy) : %.4f eV\n", intercept(f))

write_energies(f, "energy_list.txt")
write_torques(f, "torque_list.txt");

Magesty.save(SCEModel(f), "model.xml");

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
