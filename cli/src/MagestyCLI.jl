# Command-line interface for Magesty, built with Comonicon.
#
# The `@cast` / `@main` macros define `MagestyCLI.command_main` (the
# `magesty` entry point) and `MagestyCLI.comonicon_install` (writes the
# launcher into ~/.julia/bin); the latter is invoked from `deps/build.jl`.
#
# `@cast` / `@main` are qualified with `Comonicon.` because `@main` would
# otherwise clash with `Base.@main`.

module MagestyCLI

import Comonicon
import Magesty

"""
Print the installed Magesty.jl version.
"""
Comonicon.@cast function version()
	println(pkgversion(Magesty))
	return nothing
end

"""
VASP input/output conversion commands.
"""
Comonicon.@cast module Vasp

import Comonicon
import Magesty
import Magesty: vasp_to_extxyz, poscar_to_toml, oszicar_to_embset

"""
Convert a VASP run to extended XYZ (extxyz) format.

# Introduction

When only `vasprun` is given, the extxyz carries structure, forces, stress,
and energies. Passing `--oszicar` additionally writes per-atom magnetic
moments and the constraint field.

# Args

- `vasprun`: path to `vasprun.xml`.

# Options

- `--oszicar=<path>`: path to `OSZICAR`; adds the magnetic-moment and
  constraint-field columns.
- `-o, --output=<path>`: output file (`.extxyz` is appended if missing).
  Without it, the extxyz text is printed to stdout.
"""
Comonicon.@cast function extxyz(vasprun::String; oszicar::String = "", output::String = "")
	osz = isempty(oszicar) ? nothing : oszicar
	out = isempty(output) ? nothing : output
	text = vasp_to_extxyz(vasprun; oszicar = osz, output = out)
	out === nothing && print(text)
	return nothing
end

"""
Convert a VASP POSCAR structure file to a Magesty input TOML configuration.

# Introduction

The generated configuration is a starting point for an SCE input file;
the placeholder interaction settings (`lmax = 0`, `cutoff = -1`) are meant
to be edited before use.

# Args

- `poscar`: path to a POSCAR structure file.

# Options

- `-o, --output=<path>`: output file (`.toml` is appended if missing).
  Without it, the TOML text is printed to stdout.
"""
Comonicon.@cast function toml(poscar::String; output::String = "")
	out = isempty(output) ? nothing : output
	text = poscar_to_toml(poscar; output = out)
	out === nothing && print(text)
	return nothing
end

"""
Convert one or more VASP OSZICAR files to the EMBSET training-data format.

# Introduction

Each OSZICAR becomes one configuration block (energy, per-atom magnetic
moments, per-atom constraining field), numbered in the given order.

# Args

- `oszicars`: paths to one or more OSZICAR files.

# Options

- `--saxis=<"x y z">`: quantization axis, three numbers in one quoted
  argument (default: `"0.0 0.0 1.0"`).
- `--energy-kind=<f|e0>`: `f` for the free energy, `e0` for
  energy(sigma->0) (default: `f`).
- `-o, --output=<path>`: output file. Without it, the EMBSET text is
  printed to stdout.

# Flags

- `--mint`: read the magnetic moment from the `M_int` columns instead of
  `MW_int`.
"""
Comonicon.@cast function embset(
	oszicars::String...;
	saxis::String = "0.0 0.0 1.0",
	energy_kind::String = "f",
	mint::Bool = false,
	output::String = "",
)
	saxis_vec = parse.(Float64, split(strip(saxis)))
	length(saxis_vec) == 3 ||
		error("--saxis must be three numbers, e.g. \"0 0 1\"")
	isempty(oszicars) && error("at least one OSZICAR file is required")
	out = isempty(output) ? nothing : output
	text = oszicar_to_embset(
		collect(oszicars);
		saxis = saxis_vec,
		energy_kind = energy_kind,
		mint = mint,
		output = out,
	)
	out === nothing && print(text)
	return nothing
end

"""
Sample thermally conditioned spin configurations from a VASP INCAR (MFA).

# Introduction

For each value in an evenly spaced sweep of the control variable, draws
`--num-samples` spin configurations with the Mean-Field Approximation sampler
(per-atom directions from a von Mises-Fisher distribution, magnitudes
preserved) and writes each to its own `<outdir>/sample-NN.INCAR` file. The
initial spins are read from `MAGMOM`, or `M_CONSTR` if `MAGMOM` is absent; each
output sets both `MAGMOM` and `M_CONSTR` and copies all other input keys. The
run conditions (input file, sweep, sample count, options) are echoed to stdout
so a sample set is reproducible from its log alone.

# Args

- `incar`: path to the input INCAR.
- `variable`: control variable, given as the second positional argument
  (immediately after `incar`) — `tau` (scaled temperature T/Tc, in `(0, 1]`)
  or `m` (magnetization, in `[0, 1)`).

# Options

- `--start=<value>`: first sweep value (required).
- `--stop=<value>`: last sweep value (required, `>= start`).
- `--num-points=<n>`: number of evenly spaced sweep values (required, `>= 1`);
  the values are `range(start, stop; length = num_points)`.
- `--num-samples=<n>`: configurations drawn per sweep value (default: `1`).
- `--fix=<spec>`: 1-based atom indices kept at their input directions
  (rotated by the same global rotation when `--randomize`), e.g.
  `"1-10,12,20-22"`.
- `--uniform-atoms=<spec>`: 1-based atom indices whose direction is redrawn
  uniformly on the sphere instead of from the vMF distribution (same index
  syntax). These carry no mean-field alignment — their direction is fully
  isotropic for every sweep value, independent of `tau`/`m` (the disordered
  `κ → 0` limit), unlike a default atom (partially aligned) or a `--fix` atom
  (frozen). Magnitudes are preserved; if an index is also in `--fix`, `--fix`
  wins.
- `--outdir=<path>`: output directory (default: `.`, created if needed).

# Flags

- `--randomize`: apply a random global rotation (quantization-axis
  randomization) to each drawn configuration.
"""
Comonicon.@cast function mfa(
	incar::String,
	variable::String;
	start::Float64 = NaN,
	stop::Float64 = NaN,
	num_points::Int = 0,
	num_samples::Int = 1,
	randomize::Bool = false,
	fix::String = "",
	uniform_atoms::String = "",
	outdir::String = ".",
)
	# `start`/`stop` default to NaN and `num_points` to 0 as "unset" sentinels:
	# Comonicon options always carry a default, so required values are detected
	# here and reported with their CLI flag names.
	variable in ("tau", "m") || error("variable must be \"tau\" or \"m\"; got \"$variable\"")
	isnan(start) && error("--start is required")
	isnan(stop) && error("--stop is required")
	num_points >= 1 || error("--num-points is required and must be >= 1")
	num_samples >= 1 || error("--num-samples must be >= 1; got $num_samples")
	start <= stop || error("--stop ($stop) must be >= --start ($start)")
	paths = Magesty.sample_mfa_incar(
		incar;
		variable = variable,
		start = start,
		stop = stop,
		num_points = num_points,
		num_samples = num_samples,
		randomize = randomize,
		fix = fix,
		uniform_atoms = uniform_atoms,
		outdir = outdir,
	)
	# Echo the run conditions so a sample set is self-documenting and
	# reproducible from its log alone.
	println("MFA spin sampling")
	println("  Input INCAR:       ", incar)
	println("  Variable:          ", variable)
	println("  Sweep:             ", start, " to ", stop, ", ", num_points,
		" point(s) [range(start, stop; length = num_points)]")
	println("  Samples per point: ", num_samples)
	println("  Randomize axis:    ", randomize ? "yes" : "no")
	println("  Fixed atoms:       ", isempty(strip(fix)) ? "(none)" : fix)
	println("  Uniform atoms:     ", isempty(strip(uniform_atoms)) ? "(none)" : uniform_atoms)
	println("  Output directory:  ", outdir)
	println("  Files written:     ", length(paths))
	return nothing
end

end # module Vasp

"""
Sunny.jl interface commands.
"""
Comonicon.@cast module Sunny

import Comonicon
import Magesty
import Magesty: sce_to_sunny

"""
Export a fitted SCE model to a runnable Sunny.jl spin-wave script.

# Introduction

The lowest-order SALCs become a conventional spin Hamiltonian (bilinear
exchange and single-ion anisotropy); higher-order terms are skipped with a
warning. With `--placement=auto` (default) the script uses the chemical
primitive cell when the model is cleanly unfoldable (unfolded dispersion),
otherwise the training supercell (folded dispersion).

# Args

- `model`: path to a saved `SCEModel` XML file.

# Options

- `--placement=<auto|primitive|explicit>`: cell route (default: `auto`).
- `-o, --output=<path>`: output file (`.jl` is appended if missing).
  Without it, the script is printed to stdout.
"""
Comonicon.@cast function script(model::String; placement::String = "auto", output::String = "")
	placement in ("auto", "primitive", "explicit") ||
		error("--placement must be auto, primitive, or explicit; got \"$placement\"")
	out = isempty(output) ? nothing : output
	text = sce_to_sunny(
		Magesty.load(Magesty.SCEModel, model);
		placement = Symbol(placement),
		output = out,
	)
	out === nothing && print(text)
	return nothing
end

end # module Sunny

"""
Magesty command-line interface.

Build general effective spin models from noncollinear spin DFT
calculations using the spin-cluster expansion.
"""
Comonicon.@main

end # module MagestyCLI
