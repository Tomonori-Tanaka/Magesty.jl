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
import Magesty: vasp_to_extxyz, poscar_to_toml, outcar_to_embset

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
Convert one or more VASP OUTCAR files to the EMBSET training-data format.

# Introduction

Each OUTCAR becomes one configuration block (energy, per-atom magnetic
moments, per-atom constraining field), numbered in the given order.

# Args

- `outcars`: paths to one or more OUTCAR files.

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
	outcars::String...;
	saxis::String = "0.0 0.0 1.0",
	energy_kind::String = "f",
	mint::Bool = false,
	output::String = "",
)
	saxis_vec = parse.(Float64, split(strip(saxis)))
	length(saxis_vec) == 3 ||
		error("--saxis must be three numbers, e.g. \"0 0 1\"")
	isempty(outcars) && error("at least one OUTCAR file is required")
	out = isempty(output) ? nothing : output
	text = outcar_to_embset(
		collect(outcars);
		saxis = saxis_vec,
		energy_kind = energy_kind,
		mint = mint,
		output = out,
	)
	out === nothing && print(text)
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
