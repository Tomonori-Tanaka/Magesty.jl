# Tools

Magesty.jl provides a `magesty` command-line interface for converting VASP
output into the formats used by the spin-cluster expansion workflow, and for
exporting a fitted model to other spin-dynamics packages. Install the `magesty`
command as described in [Installation](installation.md). Each operation is also
available programmatically as an exported function (see the
[API Reference](@ref)).

## `magesty vasp extxyz`

Convert a single VASP calculation directory (`vasprun.xml` + optional
`OSZICAR`) to an [extended XYZ (extxyz)](https://github.com/libAtoms/extxyz)
file — the standard training-data format for machine-learning interatomic
potentials.

**Usage:**
```bash
# Structure, forces, stress, and energies only
magesty vasp extxyz vasprun.xml

# Add per-atom magnetic moments and constraint field
magesty vasp extxyz vasprun.xml --oszicar OSZICAR --output out.xyz
```

**Arguments:**
- `vasprun` (positional): Path to `vasprun.xml` (required)
- `--oszicar`: Path to `OSZICAR` — enables per-atom magnetic moment (`magmom_smoothed`, `magmom_raw`) and constraint field (`constr_field`) columns
- `--output`, `-o`: Output filename (appends `.extxyz` automatically if omitted; default: stdout)

The same conversion is available programmatically as the exported `vasp_to_extxyz` function.

**Output columns (extxyz Properties field):**

| Column | Type | Description |
|--------|------|-------------|
| `species` | S:1 | Element symbol |
| `pos` | R:3 | Cartesian positions (Å) |
| `forces` | R:3 | Forces (eV/Å) |
| `rwigs` | R:1 | Per-atom RWIGS radius (Å); present when `RWIGS` is set in INCAR |
| `magmom_smoothed` | R:3 | Smooth-window integrated magnetic moments `MW_int` (μB; integrated with a weight that decays to zero at the Wigner-Seitz sphere boundary — the form used by VASP `I_CONSTRAINED_M`); present with `--oszicar` |
| `magmom_raw` | R:3 | Directly integrated magnetic moments `M_int` (μB); present with `--oszicar` |
| `constr_field` | R:3 | Penalty constraint field λ·MW_perp (eV/μB); present with `--oszicar` |

**Global header keys:** `energy_free` (eV), `energy_zero` (eV), `stress` (eV/Å³, Voigt notation), `comment` (VASP version, ENCUT, KPOINTS, ICONST, LAMBDA).

## `magesty vasp toml`

Convert a VASP POSCAR structure file to a Magesty input TOML configuration.

**Usage:**
```bash
# Print the TOML configuration to stdout
magesty vasp toml POSCAR

# Write it to a file
magesty vasp toml POSCAR --output input.toml
```

**Arguments:**
- `poscar` (positional): Path to a POSCAR structure file (required)
- `--output`, `-o`: Output filename (appends `.toml` automatically if omitted; default: stdout)

The same conversion is available programmatically as the exported `poscar_to_toml` function.

The generated configuration fills `[general]`, `[symmetry]`, `[interaction]`, and `[structure]` from the POSCAR. POSCAR direct and Cartesian coordinates are both accepted (Cartesian is converted to direct). The interaction settings are placeholders (`lmax = 0`, `cutoff = -1`) and should be edited before use.

## `magesty vasp embset`

Extract the energy, magnetic moments, and local magnetic field from one or more VASP `OUTCAR` files and write them in the EMBSET training-data format.

**Usage:**
```bash
# Print the EMBSET text to stdout
magesty vasp embset OUTCAR1 OUTCAR2

# Write it to a file, with a non-default quantization axis
magesty vasp embset OUTCAR1 OUTCAR2 --saxis "1 0 0" --output EMBSET
```

**Arguments:**
- `outcars` (positional): Paths to one or more `OUTCAR` files; each becomes one configuration block, numbered in the given order

**Options:**
- `--saxis`: Quantization axis as three numbers in one quoted argument (default: `"0.0 0.0 1.0"`)
- `--energy-kind`: `f` (free energy) or `e0` (energy σ→0) (default: `f`)
- `--output`, `-o`: Output filename (default: stdout)

**Flags:**
- `--mint`: Extract the magnetic moment from the `M_int` columns instead of `MW_int`

The same conversion is available programmatically as the exported `outcar_to_embset` function.

## `magesty sunny script`

Export a fitted SCE model (a saved `SCEModel` XML file) to a runnable
[Sunny.jl](https://github.com/SunnySuite/Sunny.jl) script that computes a linear
spin-wave-theory (LSWT) magnon dispersion. Magesty itself takes on no Sunny
dependency — the command only emits text; you run the resulting script in an
environment that has Sunny installed.

**Usage:**
```bash
# Print the Sunny.jl script to stdout
magesty sunny script model.xml

# Write it to a file
magesty sunny script model.xml --output lswt.jl
```

**Arguments:**
- `model` (positional): Path to a saved `SCEModel` XML file (from `Magesty.save`)

**Options:**
- `--placement`: cell route, `auto` (default), `primitive`, or `explicit` (see below)
- `--output`, `-o`: Output filename (appends `.jl` automatically if omitted; default: stdout)

The same export is available programmatically as the exported `sce_to_sunny` function.

### What is converted

The lowest-order symmetry-adapted basis functions map onto a conventional spin
Hamiltonian (see the
[technical notes](technical_notes.md) for the exact coefficient correspondence):

| SCE term | Spin-model interaction | Sunny call |
|----------|------------------------|------------|
| two-site, `l₁ = l₂ = 1` | bilinear exchange (Heisenberg + Dzyaloshinskii–Moriya + anisotropic symmetric Γ), as a 3×3 matrix | `set_exchange!` |
| single-site, `l = 2` | single-ion anisotropy | `set_onsite_coupling!` |

Higher-order terms (higher-`l` pairs, three-body and beyond) cannot be
represented in Sunny's spin Hamiltonian and are skipped, with a warning listing
what was dropped. Spins use the reduced convention `s = 1`, `g = 2`, so the
dispersion is in the energy unit of the fit (typically eV). The reference energy
`j0` and any spin-independent terms are dropped (Sunny carries no constant energy
term); the dispersion is unaffected.

### Cell route and the magnetic structure

The interaction parameters are independent of the magnetic order, so the
generated script builds the crystal and interactions and then leaves a clearly
marked block for you to set up the magnetic unit cell / propagation vector before
minimizing and computing the dispersion.

- `--placement=primitive` maps the interactions onto the chemical primitive cell
  for an **unfolded** dispersion. This is exact only when the fitted model
  respects the assumption that the interaction range is below half the supercell.
- `--placement=explicit` keeps the training supercell; it is exact for any model
  but the dispersion is **folded** into the supercell Brillouin zone.
- `--placement=auto` (default) chooses the primitive route when the model is
  cleanly unfoldable and otherwise falls back to the explicit route with a
  warning.

The high-symmetry path (`qs`) in the generated script is a placeholder; edit it
for your crystal.
