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
# python FitCheck_energy.py energy_list.txt --output energy_parity.png
# python FitCheck_torque.py torque_list.txt --output torque_parity.png
# ```
#
# ![Energy parity plot for bcc Fe](case1_inputs/energy_parity.png)
#
# ![Torque parity plot for bcc Fe](case1_inputs/torque_parity.png)

# ## Exchange interactions
#
# The fitted couplings can be read back out as a conventional Heisenberg-style
# exchange. Save the model, then plot the isotropic ``J_{ij}`` against pair
# distance with the
# [`plot_jij.jl`](https://github.com/Tomonori-Tanaka/Magesty.jl/blob/main/tools/plot_jij.jl)
# script from the repository's `tools/` directory:

Magesty.save(SCEModel(f), "model.xml");

# ```bash
# julia plot_jij.jl model.xml -i -H
# ```
#
# The two flags adapt Magesty's conventions to the form of the Heisenberg model
# most common in the literature,
#
# ```math
# \mathcal{H} = -\sum_{i \neq j} J_{ij}\,\hat{\boldsymbol{e}}_i \cdot \hat{\boldsymbol{e}}_j,
# ```
#
# where ``\hat{\boldsymbol{e}}_i`` is the unit vector along the spin on site
# ``i`` and a positive ``J_{ij}`` favors ferromagnetic alignment:
#
# - `-i` inverts the sign. Magesty writes the pair energy as
#   ``+J_{ij}\,\hat{\boldsymbol{e}}_i \cdot \hat{\boldsymbol{e}}_j``, whereas the
#   Hamiltonian above carries an overall minus sign.
# - `-H` halves the values. By default — and in the SCE coefficients
#   themselves — every interaction term is counted **once** (each pair appears
#   once, and likewise a three-body or higher cluster contributes a single term;
#   no double counting). The sum over ``i \neq j`` above instead visits every
#   pair twice, so the count-once ``J_{ij}`` must be divided by two to match it.
#
# See the [technical notes](../technical_notes.md) for the full mapping from SCE
# coefficients to conventional spin-model parameters.
#
# ![Isotropic Jij versus pair distance for bcc Fe](case1_inputs/jij.png)
#
# The nearest-neighbor coupling dominates and the interaction decays, changing
# sign, at larger distances — consistent with the ferromagnetic ground state of
# bcc Fe.

# ## Spin-wave dispersion with Sunny.jl
#
# The lowest-order fitted terms map onto a conventional spin Hamiltonian, which
# lets us compute a linear spin-wave-theory (LSWT) magnon dispersion with
# [Sunny.jl](https://github.com/SunnySuite/Sunny.jl). Reusing the model XML saved
# above, generate a runnable Sunny script from it with the `magesty sunny script`
# command:
#
# ```bash
# magesty sunny script model.xml --output sunny.jl
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
