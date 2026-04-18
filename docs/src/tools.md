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
- `--oszicar`, `-s`: Path to `OSZICAR` — enables per-atom magnetic moment (`MAGMOM_smoothed`, `magmom_raw`) and constraint field (`constr_field`) columns
- `--output`, `-o`: Output filename (appends `.extxyz` automatically if omitted; default: stdout)

**Output columns (extxyz Properties field):**

| Column | Type | Description |
|--------|------|-------------|
| `species` | S:1 | Element symbol |
| `pos` | R:3 | Cartesian positions (Å) |
| `forces` | R:3 | Forces (eV/Å) |
| `rwigs` | R:1 | Per-atom RWIGS radius (Å); present when `RWIGS` is set in INCAR |
| `MAGMOM_smoothed` | R:3 | Wannier-interpolated magnetic moments `MW_int` (μB); present with `--oszicar` |
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

### tools/check_convergence_embset.jl
Scan training-set size `n` along a Julia-style range, compare fitted `jphi` to the all-data reference, and compute energy RMSE over all EMBSET configurations. Writes a summary TSV, `{prefix}_jphi_wide.tsv` and `{prefix}_jphi_long.tsv`, and optionally a two-panel PNG (requires Plots.jl).

**Usage:**
```bash
julia --project=. tools/check_convergence_embset.jl -i input.toml -e EMBSET.dat --n-range 10:10:200 -o convergence_out
julia --project=. tools/check_convergence_embset.jl -i input.toml -s system.jld2 --n-range 5:5:100
```

**Arguments:**
- `--input`, `-i`: Path to `input.toml` (required)
- `--system`, `-s`: Path to `system.jld2` saved with `@save "file.jld2" system` (optional; if omitted, `System` is built from the TOML)
- `--embset`, `-e`: EMBSET path (optional; defaults to `regression.datafile` resolved relative to the TOML directory)
- `--n-range`, `-n`: Training counts to evaluate, e.g. `20:200` or `10:10:200` (required)
- `--out`, `-o`: Output prefix for `.tsv` and `.png` files (default: `convergence_embset`)
- `--shuffle`: Shuffle EMBSET order before taking the first `n` structures for training (flag)
- `--seed`: RNG seed for `--shuffle` (default: `42`)
- `--no-plot`: Skip PNG output (TSV only) (flag)

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
- `--type`, `-t`: Plot type: `all` (x/y/z components), `norm` (magnitude), `dir` (direction unit vector) (required)
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
julia tools/histogram_magmom.jl EMBSET.txt
julia tools/histogram_magmom.jl EMBSET.txt --atoms 1,2,3-5 --bin_width 0.1
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

## Advanced Analysis Tools

### tools/micromagnetics.jl

Calculate micromagnetics model parameters (stiffness and spiralization matrices) from SCE coefficients.

**Usage:**
```bash
julia tools/micromagnetics.jl -x jphi.xml -j system.jld2
julia tools/micromagnetics.jl -x jphi.xml -j system.jld2 --cutoff 5.0
```

**Arguments:**
- `--input_xml`, `-x`: Input XML file containing SCE coefficients (required)
- `--input_jld2`, `-j`: Input JLD2 file containing a `System` struct (required)
- `--cutoff`, `-c`: Cutoff distance for atom pairs in Å (optional)

**Note:** Only the lowest-order two-body interaction terms are included.

---

### tools/convert2tensor.jl
Convert SCE coefficients to pairwise interaction tensor representation.

**Usage:**
```bash
julia tools/convert2tensor.jl jphi.xml --atoms 1 3 --output tensor.dat
```

**Arguments:**
- `input` (positional): Input XML file (required)
- `--atoms`, `-a`: Two atom indices (required)
- `--output`, `-o`: Output file name (default: `"tensor.dat"`)
- `--format`: Output format: `matrix` or `vector` (default: `"matrix"`)
