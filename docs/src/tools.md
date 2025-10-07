# Tools

This page describes the various utility tools available in the `tools/` directory of Magesty.jl. These tools provide additional functionality for data processing, analysis, and visualization.

## Data Processing Tools

### pos2toml.jl
Convert VASP structure files (POSCAR) to TOML configuration files for Magesty.jl.

**Usage:**
```bash
julia pos2toml.jl --input POSCAR --output structure.toml
```

**Features:**
- Parses VASP POSCAR format
- Converts to Magesty.jl TOML format
- Handles both direct and Cartesian coordinates
- Supports scaling factors and element symbols

### extract.jl
Extract data from various file formats used in magnetic calculations.

**Usage:**
```bash
julia extract.jl --input data_file --format vasp
```

**Supported formats:**
- VASP output files
- EMBSET files
- XML configuration files

## Visualization Tools

### histogram_magmom.jl
Generate histograms of magnetic moments.

**Usage:**
```bash
julia histogram_magmom.jl --input magmom_data.dat --output histogram.png
```

**Features:**
- Creates histograms of magnetic moment distributions
- Supports multiple data formats
- Customizable bin sizes and ranges

### yyplot.jl
Create YY plots for magnetic structure analysis.

**Usage:**
```bash
julia yyplot.jl --input data.toml --output yy_plot.png
```

**Features:**
- Generates YY plots for magnetic structures
- Supports multiple data formats
- Customizable plot parameters

## Sampling and Statistical Tools

### sampling_mfa.jl
Sample spin configurations using Mean-Field Approximation (MFA).

**Usage:**
```bash
julia sampling_mfa.jl input.toml tau --start 0.1 --end 2.0 --step 0.1 --num_samples 100
```

**Features:**
- Samples spin configurations at specified temperatures
- Uses von Mises-Fisher distribution
- Supports temperature range sampling
- Configurable number of samples

## Advanced Analysis Tools

### micromagnetics.jl
Calculate micromagnetics model parameters.

**Usage:**
```bash
julia micromagnetics.jl --input jphi.xml --output micromagnetics.dat
```

**Features:**
- Derives micromagnetics model parameters
- Calculates stiffness and spiralization matrices
- Supports multi-threading

### convert2tensor.jl
Convert between different tensor representations.

**Usage:**
```julia
using convert2tensor
tensor = convert_to_tensor(data, format="matrix")
```

**Features:**
- Converts between tensor formats
- Supports various input/output formats
- Handles symmetry operations
