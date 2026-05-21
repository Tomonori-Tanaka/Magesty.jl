# Tools

Magesty.jl provides a `magesty` command-line interface for converting VASP
output into the formats used by the spin-cluster expansion workflow. Install
the `magesty` command as described in [Installation](installation.md). Each
conversion is also available programmatically as an exported function (see
the [API Reference](@ref)).

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

The generated configuration fills `[general]`, `[symmetry]`, `[interaction]`, `[regression]`, and `[structure]` from the POSCAR. POSCAR direct and Cartesian coordinates are both accepted (Cartesian is converted to direct). The interaction settings are placeholders (`lmax = 0`, `cutoff = -1`) and should be edited before use.

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
