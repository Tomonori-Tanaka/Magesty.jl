# Command-line interface for Magesty, built with Comonicon.
#
# This file is included last in `module Magesty` so the `@cast` / `@main`
# macros see the full package API. They define `Magesty.command_main`
# (the `magesty` entry point) and `Magesty.comonicon_install` (writes the
# launcher into ~/.julia/bin); the latter is invoked from `deps/build.jl`.
#
# `@cast` / `@main` are qualified with `Comonicon.` because `@main` would
# otherwise clash with `Base.@main`.

import Comonicon

# ── VASP → extxyz conversion ────────────────────────────────────────────────

# Expand per-species RWIGS to a 1×N per-atom matrix using atomtype indices.
# atomtype_per_atom[i] is the 1-based species-type index for atom i, matching
# the order of RWIGS values in the INCAR. This correctly handles cases where
# the same element appears as multiple species with distinct RWIGS values.
function _rwigs_per_atom(rwigs::Vector{Float64}, atomtype_per_atom::Vector{Int})::Matrix{Float64}
	out = Matrix{Float64}(undef, 1, length(atomtype_per_atom))
	for (i, t) in enumerate(atomtype_per_atom)
		out[1, i] = rwigs[t]
	end
	return out
end

function _build_vasp_comment(d::VaspIO.VaspRunData)::String
	parts = String[]
	d.encut        !== nothing && push!(parts, "ENCUT=$(Int(round(d.encut)))")
	d.kpoints_mesh !== nothing && push!(parts, "KPOINTS=$(d.kpoints_mesh)")
	d.iconst       !== nothing && push!(parts, "ICONST=$(d.iconst)")
	d.lambda       !== nothing && push!(parts, @sprintf("LAMBDA=%.4g", d.lambda))
	return join(parts, ", ")
end

"""
	vasp_to_extxyz(vasprun; oszicar=nothing, output=nothing) -> String

Convert a VASP run to extended XYZ (extxyz) format and return the extxyz
text.

When only `vasprun` is given, the extxyz carries structure, forces, stress,
and energies. Passing `oszicar` additionally writes per-atom magnetic
moments (`magmom_smoothed`, `magmom_raw`) and the constraint field
(`constr_field`).

# Arguments
- `vasprun::AbstractString`: path to `vasprun.xml`.

# Keyword arguments
- `oszicar::Union{AbstractString, Nothing} = nothing`: path to `OSZICAR`;
  when given, the magnetic-moment and constraint-field columns are added.
- `output::Union{AbstractString, Nothing} = nothing`: when given, the
  extxyz text is also written to this file (`.extxyz` is appended if the
  name does not already end with it).

# Returns
- `String`: the full extxyz text.

# Examples
```julia
text = vasp_to_extxyz("vasprun.xml")
vasp_to_extxyz("vasprun.xml"; oszicar = "OSZICAR", output = "frame.extxyz")
```
"""
function vasp_to_extxyz(
	vasprun::AbstractString;
	oszicar::Union{AbstractString, Nothing} = nothing,
	output::Union{AbstractString, Nothing} = nothing,
)::String
	vd = VaspIO.parse_vasprun(vasprun)

	# Fractional → Cartesian positions (lattice columns = lattice vectors)
	positions_cart = vd.lattice * vd.positions_frac  # 3×N

	# Additional per-atom properties written as extra columns in extxyz.
	# Each entry is "name" => ncols×N matrix (ncols=1 for scalars, 3 for vectors).
	extra_per_atom = Pair{String, Matrix{Float64}}[]

	if vd.rwigs !== nothing
		push!(extra_per_atom, "rwigs" => _rwigs_per_atom(vd.rwigs, vd.atomtype_per_atom))
	end

	if oszicar !== nothing
		md = VaspIO.parse_oszicar_magdata(oszicar, vd.m_constr, vd.num_atoms)
		push!(extra_per_atom, "magmom_smoothed" => md.magmom_smoothed)
		push!(extra_per_atom, "magmom_raw"      => md.magmom_raw)
		push!(extra_per_atom, "constr_field"    => md.constr_field)
	end

	frame = ExtXYZ.AtomFrame(
		num_atoms      = vd.num_atoms,
		lattice        = vd.lattice,
		pbc            = vd.pbc,
		species        = vd.species,
		positions      = positions_cart,
		forces         = vd.forces,
		energy_free    = vd.energy_free,
		energy_zero    = vd.energy_zero,
		stress         = vd.stress,
		soc            = vd.lsorbit,
		code           = "VASP",
		version        = isempty(vd.version) ? nothing : vd.version,
		extra_per_atom = extra_per_atom,
		comment        = _build_vasp_comment(vd),
	)

	buf = IOBuffer()
	ExtXYZ.write_extxyz(buf, frame)
	text = String(take!(buf))

	if output !== nothing
		outfile = endswith(output, ".extxyz") ? output : output * ".extxyz"
		write(outfile, text)
	end
	return text
end

# ── command-line interface ──────────────────────────────────────────────────

"""
Print the installed Magesty.jl version.
"""
Comonicon.@cast function version()
	println(pkgversion(@__MODULE__))
	return nothing
end

"""
VASP input/output conversion commands.
"""
Comonicon.@cast module Vasp

import Comonicon
import ..vasp_to_extxyz

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
