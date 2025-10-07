# Tools

This page describes the various utility tools available in the `tools/` directory of Magesty.jl. These tools provide additional functionality for data processing, analysis, and visualization.

## Data Processing Tools

### pos2toml.jl
Convert VASP structure files (POSCAR) to TOML configuration files for Magesty.jl.

**Usage:**
```bash
julia pos2toml.jl --input POSCAR --output structure.toml
```

**Arguments:**
- `--input` (positional): Input VASP structure file (e.g., POSCAR) (required)
- `--output`, `-o`: Output TOML file (default: "input.toml")

**Features:**
- Parses VASP POSCAR format
- Converts to Magesty.jl TOML format
- Handles both direct and Cartesian coordinates
- Supports scaling factors and element symbols

**Example:**
```bash
julia pos2toml.jl --input POSCAR --output my_structure.toml
```

### extract.jl
Extract energy, magnetic moments, and local magnetic field from OUTCAR files and print the result in EMBSET format.

**Usage:**
```bash
julia extract.jl --target_files OUTCAR1 OUTCAR2 --energy_kind f --saxis 0 0 1
```

**Arguments:**
- `--target_files`, `-f`: Target OUTCAR files (required, multiple files allowed)
- `--energy_kind`, `-e`: Kind of energy (f: free energy, e0: energy(sigma->0)) (default: "f")
- `--mint`: Flag to extract magnetic moment from "M_int" instead of "MW_int" (usually not needed)
- `--saxis`, `-s`: Quantization axis (three real numbers) (default: [0.0, 0.0, 1.0])
- `--randomize`: Randomize the order of the spin configurations (flag)

**Features:**
- Extracts energy, magnetic moments, and magnetic field from VASP OUTCAR files
- Outputs data in EMBSET format
- Supports multiple input files
- Rotates magnetic moments and fields to specified quantization axis
- Can randomize the order of configurations

**Example:**
```bash
julia extract.jl --target_files OUTCAR1 OUTCAR2 --energy_kind f --saxis 1 0 0 --randomize
```

## Visualization Tools

### histogram_magmom.jl
Generate histograms of magnetic moments.

**Usage:**
```bash
julia histogram_magmom.jl --input EMBSET.txt --n_atoms 8 --atoms 1 2 3 4
```

**Arguments:**
- `--input`, `-i`: Input file (i.e. EMBSET.txt) (required)
- `--n_atoms`, `-n`: Total number of atoms (required)
- `--atoms`, `-a`: Atom indices to analyze (space-separated, default: all atoms)
- `--min_bound`, `-l`: Lower bound of the histogram
- `--max_bound`, `-u`: Upper bound of the histogram
- `--bin_width`, `-w`: Bin width size of the histogram
- `--vertical_lines`, `-v`: Values of vertical lines to add to plot

**Features:**
- Creates histograms of magnetic moment distributions
- Supports EMBSET.txt format
- Customizable bin sizes and ranges
- Can add vertical reference lines
- Statistical analysis (mean, variance, std)

**Example:**
```bash
julia histogram_magmom.jl --input EMBSET.txt --n_atoms 8 --atoms 1 2 3 4 --bin_width 0.1
```

### yyplot.jl
Create YY plots for magnetic structure analysis.

**Usage:**
```bash
julia yyplot.jl file1.dat file2.dat --type e --output yy_plot.png
```

**Arguments:**
- `files` (positional): Input data files (multiple files allowed)
- `--type`, `-t`: Type of data to plot (e: energy, m: magnetic field) (required)
- `--output`, `-o`: Output file name (default: "yy_plot.png")
- `--format`: Output format (png, pdf, svg) (default: "png")
- `--xlabel`: X-axis label (default: "X")
- `--ylabel`: Y-axis label (default: "Y")
- `--title`: Plot title (default: "YY Plot")

**Features:**
- Generates YY plots for magnetic structures
- Supports multiple data formats
- Customizable plot parameters

**Example:**
```bash
julia yyplot.jl data1.dat data2.dat data3.dat --type e --output energy_yy.png
```

## Sampling and Statistical Tools

### sampling_mfa.jl
Sample spin configurations using Mean-Field Approximation (MFA).

**Usage:**
```bash
julia sampling_mfa.jl input.toml tau --start 0.1 --end 2.0 --step 0.1 --num_samples 100
```

**Arguments:**
- `input` (positional): Input file path (required)
- `variable` (positional): Variable name to be sampled (tau or m) (required)
 - `--start`, `-s`: Starting temperature ($T/T_{\rm c}^{\rm MFA}$) (required)
 - `--end`, `-e`: Ending temperature ($T/T_{\rm c}^{\rm MFA}$) (required)
- `--step`, `-w`: Temperature step size (required)
- `--num_samples`, `-n`: Number of samples in each step (default: 1)
- `--randomize`, `-r`: Randomize the quantization axis (flag)

**Features:**
- Samples spin configurations at specified temperatures ($\tau$) or magnetization (m)
- Uses von Mises-Fisher distribution

**Example:**
```bash
julia sampling_mfa.jl config.toml tau --start 0.1 --end 1.0 --step 0.1 --num_samples 50 --randomize
```

This command generates 50 samples at each step for $\tau$ (or $m$) values 0.1, 0.2, ..., 1.0 (i.e., 10 steps × 50 samples).

## Advanced Analysis Tools

### micromagnetics.jl
Calculate micromagnetics model parameters.

**Usage:**
```bash
julia micromagnetics.jl --input_xml jphi.xml --input_jld2 system.jld2
```

**Arguments:**
- `--input_xml`, `-x`: Input XML file (required)
- `--input_jld2`, `-j`: Input JLD2 file for System struct (required)
- `--cutoff`, `-c`: Cutoff distance for atom pairs in Angstrom (optional)

**Features:**
- Derives micromagnetics model parameters
- Calculates stiffness and spiralization matrices
- Supports multi-threading

**Note:**
- Only considers the lowest-order two-body interaction terms in the SCE model

**Example:**
```bash
julia micromagnetics.jl --input_xml jphi.xml --input_jld2 system.jld2 --cutoff 5.0
```

### convert2tensor.jl
Convert between different tensor representations.

**Usage:**
```bash
julia convert2tensor.jl input.xml --atoms 1 2 --output tensor.dat
```

**Arguments:**
- `input` (positional): Input XML file (required)
- `--atoms`, `-a`: Two atom indices (required)
- `--output`, `-o`: Output file name (default: "tensor.dat")
- `--format`: Output format (matrix, vector) (default: "matrix")

**Features:**
- Converts between tensor formats
- Supports various input/output formats
- Handles symmetry operations

**Example:**
```bash
julia convert2tensor.jl jphi.xml --atoms 1 3 --output interaction_tensor.dat
```
