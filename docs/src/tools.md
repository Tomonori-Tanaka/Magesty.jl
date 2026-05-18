# Tools

This page describes the utility tools available in the `tools/` directory of Magesty.jl.

## ExtXYZ Tools

These tools convert VASP output to the [extended XYZ (extxyz)](https://github.com/libAtoms/extxyz) format, which is the standard training-data format for machine-learning interatomic potentials.

### tools/vasp/vasp2extxyz.jl

Convert a single VASP calculation directory (vasprun.xml + optional OSZICAR) to an extxyz file.

**Usage:**
```bash
# Structure, forces, stress, and energies only
julia tools/vasp/vasp2extxyz.jl --vasprun vasprun.xml

# Add per-atom magnetic moments and constraint field
julia tools/vasp/vasp2extxyz.jl --vasprun vasprun.xml --oszicar OSZICAR --output out.xyz
```

If `~/.julia/bin` is in your PATH (see [Installation](installation.md)), you can also call it as a CLI command:
```bash
vasp2extxyz --vasprun vasprun.xml
vasp2extxyz --vasprun vasprun.xml --oszicar OSZICAR --output out.xyz
```

**Arguments:**
- `--vasprun`, `-v`: Path to `vasprun.xml` (required)
- `--oszicar`, `-s`: Path to `OSZICAR` — enables per-atom magnetic moment (`magmom_smoothed`, `magmom_raw`) and constraint field (`constr_field`) columns
- `--output`, `-o`: Output filename (appends `.extxyz` automatically if omitted; default: stdout)

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

---

### tools/vasp/vasp2extxyz_recursive.jl

Recursively walk a directory tree and convert every subdirectory that contains `vasprun.xml` to extxyz. Useful for batch-processing large training datasets.

**Usage:**
```bash
# Write <dirname>.extxyz in each subdirectory
julia tools/vasp/vasp2extxyz_recursive.jl --root ./calculations

# Collect all frames into a single file
julia tools/vasp/vasp2extxyz_recursive.jl --root ./calculations --mode combined --output dataset.extxyz

# Both individual files and a combined file
julia tools/vasp/vasp2extxyz_recursive.jl --root ./calculations --mode both
```

If `~/.julia/bin` is in your PATH (see [Installation](installation.md)), you can also call it as a CLI command:
```bash
vasp2extxyz_recursive --root ./calculations
vasp2extxyz_recursive --root ./calculations --mode combined --output dataset.extxyz
```

**Arguments:**
- `--root`, `-r`: Root directory to search (default: `.`)
- `--vasprun`, `-v`: `vasprun.xml` filename to look for (default: `vasprun.xml`)
- `--oszicar`, `-s`: `OSZICAR` filename to look for (default: `OSZICAR`)
- `--mode`, `-m`: Output mode — `each`, `combined`, or `both` (default: `each`)
- `--output`, `-o`: Output filename for `combined` mode (default: `combined.extxyz`)

**Output modes:**

| Mode | Behavior |
|------|----------|
| `each` | Write `<dirname>.extxyz` inside each successfully processed directory |
| `combined` | Append all frames to a single file (set with `--output`) |
| `both` | Do both of the above |

**Skip / warning rules:**
- Neither `vasprun.xml` nor `OSZICAR` found → silently skip
- `vasprun.xml` missing but `OSZICAR` present → warn and skip
- `OSZICAR` missing but `vasprun.xml` present → note and proceed (magnetic data omitted)
- Parse error → warn and skip that directory

## Data Processing Tools

### tools/vasp/pos2toml.jl

Convert a VASP structure file (POSCAR) to a Magesty.jl TOML configuration file.

**Usage:**
```bash
julia tools/vasp/pos2toml.jl POSCAR
julia tools/vasp/pos2toml.jl POSCAR --output structure.toml
```

**Arguments:**
- `input` (positional): Input VASP structure file (e.g., POSCAR) (required)
- `--output`, `-o`: Output TOML file (default: `"input.toml"`)

**Features:**
- Parses VASP POSCAR format (direct and Cartesian coordinates)
- Converts to Magesty.jl TOML format
- Handles scaling factors and element symbols

---

### tools/extract.jl
Extract energy, magnetic moments, and local magnetic field from VASP OUTCAR files and write in EMBSET format.

**Usage:**
```bash
julia tools/extract.jl -f OUTCAR1 OUTCAR2 -e f -s 0 0 1
```

**Arguments:**
- `--target_files`, `-f`: Target OUTCAR files (required, multiple files allowed)
- `--energy_kind`, `-e`: Kind of energy: `f` (free energy) or `e0` (energy σ→0) (default: `"f"`)
- `--mint`: Extract magnetic moment from `"M_int"` instead of `"MW_int"` (flag)
- `--saxis`, `-s`: Quantization axis as three real numbers (default: `[0.0, 0.0, 1.0]`)
- `--randomize`: Randomize the order of spin configurations (flag)

**Example:**
```bash
julia tools/extract.jl -f OUTCAR1 OUTCAR2 -e f -s 1 0 0 --randomize
```

---

### tools/vasp/oszicar2magmom.jl

Extract magnetic moments from the last SCF step of a VASP `OSZICAR` file and print them as a single-line INCAR `MAGMOM = ...` entry. Useful for seeding a follow-up VASP run from the magnetic state of a previous calculation.

**Usage:**
```bash
julia tools/vasp/oszicar2magmom.jl OSZICAR
julia tools/vasp/oszicar2magmom.jl OSZICAR --type M_int
julia tools/vasp/oszicar2magmom.jl OSZICAR --output magmom.txt
```

**Arguments:**
- `oszicar_file` (positional): Path to `OSZICAR` (required)
- `--type`, `-t`: Magnetic moment type, `MW_int` (smooth-window integrated moment used by VASP `I_CONSTRAINED_M`; integrated with a weight that vanishes at the Wigner-Seitz sphere boundary) or `M_int` (direct integration over the PAW sphere) (default: `MW_int`)
- `--output`, `-o`: Output filename. If omitted, the `MAGMOM = ...` line is written to stdout (informational messages go to stderr)

## Visualization Tools

### tools/FitCheck_energy.py

Scatter plot comparing observed (DFT) and predicted (SCE) energies. Requires Python with `matplotlib` and `numpy`.

**Usage:**
```bash
python3 tools/FitCheck_energy.py energy_list.txt -o energy.png
python3 tools/FitCheck_energy.py energy_list.txt -o energy.png -l 0.2
```

**Arguments:**
- `files` (positional): Input files; 2 columns (observed predicted) or 3+ columns (index observed predicted)
- `--output`, `-o`: Output filename (format inferred from extension: png, svg, pdf)
- `--lim`, `-l`: Fix X/Y axes to `[-lim, lim]` in **eV**
- `--lim-min`, `--lim-max`: Asymmetric axis limits in eV
- `--zero-min`, `-z`: Shift so the global minimum observed energy is 0
- `--zero-at`: Shift so a specified observed energy (eV) is 0
- `--colored`, `-c`: Highlight specific data indices (e.g. `'5,7-10'`)
- `--marker-size`, `-m`: Marker size (default: 5)
- `--marker-alpha`, `-A`: Marker transparency 0–1 (default: 0.8)
- `--no-legend`, `-L`: Disable legend
- `--tick-interval`, `-T`: Major tick interval in meV (same for X and Y; default: auto)
- `--per`: Divide all energies by this integer (e.g. atoms per formula unit) and report per-unit energies

**Features:**
- eV → meV conversion for display
- Per-file center shift on observed values
- y=x reference line with per-series RMSE in legend

---

### tools/FitCheck_torque.py
Scatter plot comparing observed (DFT) and predicted (SCE) torques. Supports component, magnitude, and direction comparisons, with filtering by atom indices or elements. Requires Python with `matplotlib` and `numpy`.

**Usage:**
```bash
python3 tools/FitCheck_torque.py torque_list.txt -t all -o torque.png
python3 tools/FitCheck_torque.py torque_list.txt -t norm -e Fe
python3 tools/FitCheck_torque.py torque_list.txt -t dir -a 1,2,3 -l 50
```

**Arguments:**
- `files` (positional): Input torque files (format: `index element DFTx DFTy DFTz SCEx SCEy SCEz`)
- `--type`, `-t`: Plot type: `all` (x/y/z components), `norm` (magnitude), `dir` (direction unit vector) (default: `all`)
- `--atom-indices`, `-a`: Atom indices to include (e.g. `'1,3,5-10'`; mutually exclusive with `--elements`)
- `--elements`, `-e`: Element symbols to include (e.g. `'Fe,Co'`; mutually exclusive with `--atom-indices`)
- `--output`, `-o`: Output filename (format inferred from extension)
- `--lim`, `-l`: Fix axes to `[-lim, lim]` in meV (`all`/`norm`) or unitless (`dir`)
- `--lim-min`, `--lim-max`: Asymmetric axis limits in meV
- `--marker-size`, `-m`: Marker size (default: 5)
- `--marker-alpha`, `-A`: Marker transparency 0–1 (default: 0.8)
- `--no-legend`, `-L`: Disable legend
- `--tick-interval`, `-T`: Major tick interval in meV

**Features:**
- eV → meV conversion for `all`/`norm`; unitless for `dir`
- Per-file center shift on observed values
- y=x reference line with per-series RMSE in legend

---

### tools/histogram_magmom.jl
Generate histograms of magnetic moment magnitudes from an EMBSET file.

**Usage:**
```bash
julia tools/histogram_magmom.jl EMBSET
julia tools/histogram_magmom.jl EMBSET --atoms 1,2,3-5 --bin_width 0.1
```

**Arguments:**
- `files` (positional): Input EMBSET files (multiple files allowed)
- `--atoms`, `-a`: Atom indices to analyze; accepts integers, ranges (`1-5`), and comma-separated lists (default: all atoms)
- `--data`, `-d`: SpinConfig indices to analyze (default: all)
- `--min_bound`, `-l`: Lower bound of histogram
- `--max_bound`, `-u`: Upper bound of histogram
- `--bin_width`, `-w`: Bin width
- `--vertical_lines`, `-v`: Values at which to draw vertical reference lines

**Features:**
- Histograms of magnetic moment magnitude distributions from EMBSET data
- Statistical summary (mean, variance, std)

---

### tools/plot_jij_atom.jl

Plot isotropic `Jij` vs interatomic distance for all pairs involving a specified reference atom. Requires Plots.jl. Displays the figure interactively (no `-o` save option).

**Usage:**
```bash
julia tools/plot_jij_atom.jl jphi.xml 1
julia tools/plot_jij_atom.jl jphi.xml 1 --invert-sign --ymin -10 --ymax 10
```

**Arguments:**
- `input` (positional): Input XML file (e.g. `jphi.xml`) (required)
- `reference_atom` (positional): 1-based atom index used as the reference for `Jij` pairs (required)
- `--invert-sign`, `-i`: Invert the sign of `Jij` values (flag)
- `--half-jij`, `-H`: Plot `Jij/2` instead of `Jij` (flag)
- `--ymin`, `--ymax`: Y-axis bounds (default: auto)
- `--markersize`, `-m`: Marker size (default: 5)
- `--no-legend`: Hide the legend (flag)

---

### tools/plot_jij.jl

Plot isotropic `Jij` vs distance for **all pairs** in the system across one or more XML files, with optional element-pair filtering. Multiple input files are drawn as distinct series. Requires Plots.jl.

**Usage:**
```bash
julia tools/plot_jij.jl jphi.xml
julia tools/plot_jij.jl runA.xml runB.xml --label 'A,B'
julia tools/plot_jij.jl jphi.xml --element1 Fe --element2 Co
```

**Arguments:**
- `input` (positional, one or more): Input XML files (required)
- `--label`, `-l`: Comma-separated labels for each input file (default: file basename)
- `--element1`, `-e` / `--element2`, `-E`: Filter to pairs with this element on each side (both required together)
- `--invert-sign`, `-i`: Invert the sign of `Jij` values (flag)
- `--half-jij`, `-H`: Plot `Jij/2` instead of `Jij` (flag)
- `--ymin`, `--ymax`: Y-axis bounds (default: auto)
- `--markersize`, `-m`: Marker size (default: 5)
- `--no-legend`: Hide the legend (flag)

---

### tools/plot_jphi_cluster_distance.jl

Plot each SALC coefficient `jφ` (from `save(model, ...)`) against the maximum over cluster atom pairs of the minimum-image (MIC) distance. Different N-body terms are drawn as separate series. Requires Plots.jl.

**Usage:**
```bash
julia tools/plot_jphi_cluster_distance.jl jphi.xml
julia tools/plot_jphi_cluster_distance.jl jphi.xml --bodies 2,3 -o jphi.png
julia tools/plot_jphi_cluster_distance.jl jphi.xml --per-cluster
```

**Arguments:**
- `input` (positional): Path to a Magesty XML file (e.g. `jphi.xml`) (required)
- `--bodies`, `-b`: Comma-separated N-body values to include (e.g. `2,3`); omit or `all` for all clusters
- `--output`, `-o`: Output path (PNG, etc.). If set, save and exit without interactive display. Sets `GKSwstype=100` for GR when unset (headless-friendly).
- `--title`: Plot title (default: auto)
- `--per-cluster`: Divide each SALC coefficient by `√num_basis` (exact for scalar `Lf=0` SALCs; RMS approximation otherwise) (flag)

## Sampling Tools

### tools/sampling_mfa.jl
Sample spin configurations using the Mean-Field Approximation (MFA) with a von Mises-Fisher distribution.

**Usage:**
```bash
julia tools/sampling_mfa.jl config.toml tau --start 0.1 --end 1.0 --step 0.1 --num_samples 50
julia tools/sampling_mfa.jl config.toml m   --start 0.5 --end 0.9 --step 0.1 --randomize
```

**Arguments:**
- `input` (positional): Path to `input.toml` (required)
- `variable` (positional): Sampling variable: `tau` ($T/T_\mathrm{c}^\mathrm{MFA}$) or `m` (magnetization) (required)
- `--start`, `-s`: Starting value of the sweep (required)
- `--end`, `-e`: Ending value of the sweep (required)
- `--step`, `-w`: Step size (required)
- `--num_samples`, `-n`: Number of samples per step (default: 1)
- `--randomize`, `-r`: Randomize the quantization axis (flag)
- `--fix`: Fix magnetic moments for specified 1-based atom indices (e.g. `"1-10,12"`). Without `--randomize`: keep original values; with `--randomize`: apply the same rotation as other atoms.
- `--uniform-atoms`: Sample specified 1-based atom indices from a uniform distribution on the sphere instead of MFA (e.g. `"1-10,12"`).

**Example:**
```bash
julia tools/sampling_mfa.jl config.toml tau --start 0.1 --end 1.0 --step 0.1 --num_samples 50 --randomize
```

This generates 50 samples at each step for $\tau$ = 0.1, 0.2, …, 1.0.

---

### tools/sampling_mfa_lebedev.jl

Sample spin configurations using the Mean-Field Approximation (MFA) with quantization axes placed at the points of a Lebedev quadrature grid on the sphere. Reuses the MFA core from `sampling_mfa.jl`.

**Usage:**
```bash
julia tools/sampling_mfa_lebedev.jl config.toml tau 0.5 --order 9 --num_samples 5
julia tools/sampling_mfa_lebedev.jl config.toml m 0.8 --lebedev grid.txt
julia tools/sampling_mfa_lebedev.jl config.toml tau 0.5 --order 13 --theta-min 30 --theta-max 90
```

**Arguments:**
- `input` (positional): Path to `input.toml` (required)
- `variable` (positional): Sampling variable: `tau` ($T/T_\mathrm{c}^\mathrm{MFA}$) or `m` (magnetization) (required)
- `value` (positional): Single value of the sampling variable (required, `Float64`)
- `--order`: Lebedev quadrature order (generate grid in-script; e.g. `9` → 38 points). If set (`> 0`), `--lebedev` is ignored. (default: `0`)
- `--num_samples`, `-n`: Number of samples per Lebedev direction (default: `1`)
- `--lebedev`, `-l`: Path to a Lebedev grid file (used only when `--order` is `0`) (default: `grid.txt`)
- `--theta-min` / `--theta-max`: Polar angle bounds in degrees (default: `0`, `180`)
- `--phi-min` / `--phi-max`: Azimuthal angle bounds in degrees (default: `0`, `360`)
- `--output-grid`, `-o`: Write the (filtered) grid points to a file as `x y z [weight]`; no output if omitted
- `--fix`: Fix magnetic moments for 1-based atom indices (e.g. `"1-10,12"`); the same uniform rotation as other atoms is applied

## Advanced Analysis Tools

### tools/micromagnetics.jl

Calculate micromagnetics model parameters (stiffness and spiralization matrices) from SCE coefficients.

**Usage:**
```bash
julia tools/micromagnetics.jl -x jphi.xml
julia tools/micromagnetics.jl -x jphi.xml --cutoff 5.0
```

**Arguments:**
- `--input_xml`, `-x`: Input XML file containing the SCE basis and coefficients (e.g. the output of `save(model, "jphi.xml")`) (required)
- `--cutoff`, `-c`: Cutoff distance for atom pairs in Å (optional)

**Note:** Only the lowest-order two-body interaction terms are included.

---

### tools/convert2tensor.jl
Convert SCE coefficients to a pairwise interaction tensor (isotropic `Jij`, symmetric anisotropy, DM vector) for a single atom pair. The tensor is printed to stdout in human-readable form. Also exports the `ExchangeTensor` module for use from other scripts (`tools/micromagnetics.jl` reuses it).

**Usage:**
```bash
julia tools/convert2tensor.jl jphi.xml --atoms 1 3
```

**Arguments:**
- `input` (positional): Input XML file with SCE basis and coefficients (e.g. `jphi.xml`) (required)
- `--atoms`, `-a`: Two atom indices in the supercell (required, exactly 2)

---

### tools/mfa_analysis.jl

Mean-field analysis of an SCE model. Fourier-transforms the exchange interactions to obtain $J(q)$ on a uniform reciprocal grid, finds the wave vector that minimizes the lowest eigenvalue of $J(q)$ (i.e. the MFA ordering vector), and optionally refines that minimum with gradient descent.

**Usage:**
```bash
julia tools/mfa_analysis.jl -x jphi.xml
julia tools/mfa_analysis.jl -x jphi.xml --nk 40 --spin 2.5
julia tools/mfa_analysis.jl -x jphi.xml --no-refine --eigvec 1,2,3
```

**Arguments:**
- `--xml`, `-x`: Path to `jphi.xml` or `scecoeffs.xml` (required)
- `--nk`, `-n`: Number of q-points per reciprocal direction (default: `20`)
- `--spin`, `-s`: Spin magnitude `S`. Omit for classical spins (`Tc ∝ S^2/3`, `S=1`); provide a value for quantum spins (`Tc ∝ S(S+1)/3`)
- `--no-refine`: Skip gradient-descent refinement of the grid-search minimum (flag)
- `--eigvec`, `-e`: Comma-separated indices of eigenvectors to print at the minimum (e.g. `1,2,3`; `1` is the smallest eigenvalue)
