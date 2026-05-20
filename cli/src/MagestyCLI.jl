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
import Magesty: vasp_to_extxyz

"""
Convert a VASP run to extended XYZ (extxyz) format.

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

end # module Vasp

"""
Magesty command-line interface.

Build general effective spin models from noncollinear spin DFT
calculations using the spin-cluster expansion.
"""
Comonicon.@main

end # module MagestyCLI
