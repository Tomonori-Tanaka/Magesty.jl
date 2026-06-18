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

**Global header keys:** `energy_free` (eV), `energy_zero` (eV), `stress` (eV/Å³, Voigt notation), `soc` (`T`/`F`, spin-orbit coupling flag), `converged` (`T`/`F`, electronic SCF convergence; omitted when the status is indeterminate, e.g. no `NELM`, single-shot `NELM = 1`, or `EDIFF <= 0`), `comment` (VASP version, ENCUT, KPOINTS, ICONST, LAMBDA). A non-converged run is still written, with `converged=F` and a warning.

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

Extract the energy, magnetic moments, and local magnetic field from one or more VASP `OSZICAR` files and write them in the EMBSET training-data format.

**Usage:**
```bash
# Print the EMBSET text to stdout
magesty vasp embset OSZICAR1 OSZICAR2

# Write it to a file, with a non-default quantization axis
magesty vasp embset OSZICAR1 OSZICAR2 --saxis "1 0 0" --output EMBSET
```

**Arguments:**
- `oszicars` (positional): Paths to one or more `OSZICAR` files; each becomes one configuration block, numbered in the given order

**Options:**
- `--saxis`: Quantization axis as three numbers in one quoted argument (default: `"0.0 0.0 1.0"`)
- `--energy-kind`: `f` (free energy) or `e0` (energy σ→0) (default: `f`)
- `--output`, `-o`: Output filename (default: stdout)

**Flags:**
- `--mint`: Extract the magnetic moment from the `M_int` columns instead of `MW_int`

The same conversion is available programmatically as the exported `oszicar_to_embset` function.

## `magesty vasp mfa`

Sample thermally conditioned spin configurations from a VASP INCAR using the
Mean-Field Approximation (MFA). The initial spins are read from `MAGMOM` (or
`M_CONSTR` if `MAGMOM` is absent), and for each value in an evenly spaced sweep
of the control variable, per-atom directions are drawn from a von Mises-Fisher
distribution whose concentration follows from the MFA self-consistency
equation. Each drawn configuration is written to its own INCAR file (with both
`MAGMOM` and `M_CONSTR` set, all other input keys preserved).

**Usage:**
```bash
# 9 temperatures from 0.1 to 0.9 (T/Tc), one configuration each
magesty vasp mfa INCAR tau --start 0.1 --stop 0.9 --num-points 9

# Sweep magnetization instead, 100 samples per value, into ./samples
magesty vasp mfa INCAR m --start 0.2 --stop 0.8 --num-points 7 \
    --num-samples 100 --outdir samples
```

**Arguments:**
- `incar` (positional): Path to the input INCAR (required)
- `variable` (positional): Control variable — `tau` (reduced temperature `T/Tc`) or `m` (magnetization) (required)

**Options:**
- `--start`: First sweep value (required)
- `--stop`: Last sweep value (required, `≥ start`)
- `--num-points`: Number of evenly spaced sweep values (required, `≥ 1`); the values are `range(start, stop; length = num_points)`
- `--num-samples`: Configurations drawn per sweep value (default: `1`)
- `--fix`: 1-based atom indices kept at their input directions (rotated by the same global rotation when `--randomize`), e.g. `"1-10,12,20-22"`
- `--uniform-atoms`: 1-based atom indices whose direction is redrawn uniformly on the sphere instead of from the vMF distribution (same index syntax). These atoms carry no mean-field alignment: their direction is fully isotropic for every sweep value, independent of `tau`/`m` — in contrast to a default atom (partially aligned by `κ = 3m/τ`) and a `--fix` atom (frozen at its input). Use it for sites outside the ordered network, e.g. a weakly coupled impurity or a sublattice held paramagnetic as a control. Magnitudes are preserved and zero-moment sites are left untouched; if an index is also in `--fix`, `--fix` wins. See [Mean-Field Sampling](tips/mfa_sampling.md) for the underlying limit
- `--outdir`: Output directory (default: `.`, created if needed)

**Flags:**
- `--randomize`: Apply a random global rotation (quantization-axis randomization) to each drawn configuration

Output files are `<outdir>/sample-NN.INCAR`, numbered `(point-1)·num_samples + sample` and zero-padded to the width of `num_points · num_samples`.

The run conditions (input file, control variable, sweep range and point count, samples per point, randomization, fixed/uniform atom indices, output directory, and the number of files written) are echoed to stdout, so a generated sample set is reproducible from its log alone. Redirect stdout to a file alongside the samples to keep that record.

The same sampling is available programmatically as the exported `sample_mfa_incar` function. See [Mean-Field Sampling](tips/mfa_sampling.md) for the underlying theory (the von Mises-Fisher distribution and the MFA self-consistency equation).

## `magesty sunny script`

Export a fitted SCE model (a saved `SCEModel` XML file) to a runnable
[Sunny.jl](https://github.com/SunnySuite/Sunny.jl) script that computes a linear
spin-wave-theory (LSWT) magnon dispersion. Magesty itself takes on no Sunny
dependency — the command only emits text; you run the resulting script in an
environment that has Sunny installed.

**Usage:**
```bash
# MnTe: Mn²⁺ has S = 5/2. --spin is required.
magesty sunny script model.xml --spin 5//2

# Write it to a file
magesty sunny script model.xml --spin 5//2 --output lswt.jl

# Per-species spin (e.g. a ferrimagnet)
magesty sunny script model.xml --spin Mn=5//2,Fe=3//2
```

**Arguments:**
- `model` (positional): Path to a saved `SCEModel` XML file (from `Magesty.save`)

**Options:**
- `--spin` (**required**): physical effective spin length `S_eff = m/(g μ_B)`. A bare
  number applies to all magnetic species (`5//2`, `2`); a comma list
  `Mn=5//2,Fe=3//2` sets it per species. See "Physical spin" below.
- `--g`: `g`-factor, scalar or per-species (default `2`); does not affect the bare dispersion
- `--mode`: Sunny system mode, `auto` (default), `dipole`, or `dipole_uncorrected`
- `--placement`: cell route, `auto` (default), `primitive`, or `explicit` (see below)
- `--output`, `-o`: Output filename (appends `.jl` automatically if omitted; default: stdout)

The same export is available programmatically as the exported `sce_to_sunny` function.

### Physical spin

The SCE couplings are fit with unit spin directions, so they absorb the spin
magnitude (`J_SCE = J_phys·S²`). The classical energy is independent of the spin
length, but the **magnon dispersion scales as `ħω ∝ 1/s`** for a fixed energy
landscape. You must therefore supply the physical effective spin
`S_eff = m/(g μ_B)` (the local-moment magnitude); the placeholder `s = 1` would
inflate the dispersion by a factor `~S` (for MnTe, `S = 5/2` ⇒ ~2.5× too high). At
emission each bilinear bond is rescaled by `1/(s_i s_j)` and each single-ion term
by a mode-dependent factor, so `energy(sys)` still reproduces
`predict_energy(model, …) - j0` while the dispersion is physical.

Sunny requires `s` to be an exact multiple of `1/2`, so `--spin` must be a
half-integer; a non-half-integer (itinerant) `S_eff` is rejected. `--mode=auto`
uses `:dipole` (the quantum single-ion renormalization `s(2s-1)/2`);
`:dipole_uncorrected` uses the classical `s²` large-`s` limit.

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
what was dropped. The dispersion is in the energy unit of the fit (typically eV).
The reference energy `j0` and any spin-independent terms are dropped (Sunny
carries no constant energy term); the dispersion is unaffected.

### Cell route and the magnetic structure

The interaction parameters are independent of the magnetic order, so the
generated script builds the crystal and interactions and then leaves a clearly
marked block for you to set up the magnetic unit cell / propagation vector before
minimizing and computing the dispersion.

- `--placement=primitive` maps the interactions onto the chemical primitive cell
  for an **unfolded** dispersion. Magesty keeps only minimum-distance images when
  fitting (interactions beyond the representable range are zero), so every fitted
  model unfolds exactly — pairs connected by several equal-distance lattice
  vectors (multiplicity > 1) are placed as separate primitive bonds.
- `--placement=explicit` keeps the training supercell; the dispersion is
  **folded** into the supercell Brillouin zone. Equivalent physics, just a
  larger magnetic cell.
- `--placement=auto` (default) uses the primitive route (a rare unresolvable
  geometry falls back to explicit with a warning).

The high-symmetry path (`qs`) in the generated script is a placeholder; edit it
for your crystal.
