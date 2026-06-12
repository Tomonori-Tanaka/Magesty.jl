# ## Building the basis
#
# Build the symmetry-adapted SCE basis from the input TOML. This runs the SALC
# construction, the heavy step of the workflow.

using Magesty

basis = SCEBasis("input.toml")

# Inspect the detected symmetry and the size of the basis:

sym = basis.symmetry
println("Space group: ", sym.international_symbol, " (#", sym.spacegroup_number, ")")
println("Symmetry operations: ", sym.nsym)
println("Atoms in the supercell: ", basis.structure.supercell.num_atoms)
println("Number of SALCs: ", length(basis.salcbasis.salc_list))

# ## Fitting
#
# Pair the basis with the `EMBSET` reference data to form a dataset, then fit the
# SCE coefficients with ordinary least squares.

dataset = SCEDataset(basis, "EMBSET"; verbosity = false)
f = fit(SCEFit, dataset, OLS(); verbosity = false)
println("Fitted ", length(coef(f)), " coefficients on ", length(dataset),
        " configurations.")

# ## Validation
#
# In-sample fit quality: the root-mean-square errors and the energy coefficient
# of determination.

using Printf

@printf("RMSE energy : %.3f meV\n", rmse_energy(f) * 1000)
@printf("R^2  energy : %.4f\n", r2_energy(f))
@printf("RMSE torque : %.3f meV\n", rmse_torque(f) * 1000)
@printf("j0 (reference energy) : %.4f eV\n", intercept(f))

# ### Parity plots
#
# For a visual check, write the observed-vs-predicted energies and torques, then
# render parity plots with the `FitCheck` helper scripts under `tools/`.

write_energies(f, "energy_list.txt")
write_torques(f, "torque_list.txt");

# Render the parity plots with the `FitCheck` scripts —
# [`FitCheck_energy.py`](https://github.com/Tomonori-Tanaka/Magesty.jl/blob/main/tools/FitCheck_energy.py)
# and
# [`FitCheck_torque.py`](https://github.com/Tomonori-Tanaka/Magesty.jl/blob/main/tools/FitCheck_torque.py)
# from the repository's `tools/` directory. (When Magesty is installed with
# `] add Magesty`, this directory lives under a version-specific path in the
# package depot, so the links above are the easiest way to reach the scripts.)
#
# ```bash
# python FitCheck_energy.py energy_list.txt -o energy_parity.png
# python FitCheck_torque.py torque_list.txt -o torque_parity.png
# ```
#
# ![Energy parity plot for bcc Fe](case1_inputs/energy_parity.png)
#
# ![Torque parity plot for bcc Fe](case1_inputs/torque_parity.png)

# ## Physical interpretation
#
# (to be written)

# ## Spin-wave dispersion with Sunny.jl
#
# The lowest-order fitted terms map onto a conventional spin Hamiltonian, which
# lets us compute a linear spin-wave-theory (LSWT) magnon dispersion with
# [Sunny.jl](https://github.com/SunnySuite/Sunny.jl). Save the fitted model,
# then generate a runnable Sunny script from it with the `magesty sunny script`
# command.

Magesty.save(SCEModel(f), "model.xml");

# ```bash
# magesty sunny script model.xml -o sunny.jl
# ```
#
# The generated script builds the spin system from the fitted exchange constants
# and evaluates the dispersion along the standard bcc high-symmetry path
# Γ–H–N–Γ–P–H. Magesty itself takes on no Sunny dependency — run the script in an
# environment that has Sunny installed. The script generated for this example is
# [`sunny.jl`](case1_inputs/sunny.jl).
#
# Running it produces the magnon dispersion below:
#
# ![Spin-wave dispersion of bcc Fe along the Γ–H–N–Γ–P–H path](case1_inputs/dispersion.png)
#
# As a check against experiment, the SCE dispersion along the Γ–N direction
# closely tracks inelastic neutron-scattering measurements of bcc Fe[^expt]:
#
# ![SCE (Sunny) magnon dispersion along Γ–N compared with experiment](case1_inputs/dispersion_exp.png)
#
# [^expt]: Experimental spin-wave energies from C.-K. Loong, J. M. Carpenter,
#     J. W. Lynn, R. A. Robinson, and H. A. Mook, "Neutron scattering study of the
#     magnetic excitations in ferromagnetic iron at high energy transfers",
#     J. Appl. Phys. **55**, 1895–1897 (1984).
#     DOI: [10.1063/1.333511](https://doi.org/10.1063/1.333511).
