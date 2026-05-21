# VASP-to-extxyz conversion. Orchestrates the `VaspIO` parser and the
# `ExtXYZ` writer into the `vasp_to_extxyz` public API.

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

# ── POSCAR → Magesty input TOML ─────────────────────────────────────────────

# Write a Magesty input TOML configuration built from parsed POSCAR data.
# Ordered dictionaries keep section keys in a deterministic write order.
function _write_input_toml(io::IO, data::VaspIO.PoscarData)
	# [general]
	general = DataStructures.OrderedDict{String, Any}()
	general["name"]        = data.comment
	general["nat"]         = sum(data.numbers)
	general["kd"]          = data.element_symbols
	general["periodicity"] = [true, true, true]

	# [symmetry]
	symmetry = DataStructures.OrderedDict{String, Any}()
	symmetry["tolerance"] = 1e-5
	symmetry["isotropy"]  = true

	# [interaction.body1] lmax: one entry per element symbol
	lmax = DataStructures.OrderedDict{String, Any}()
	for element in data.element_symbols
		lmax[element] = 0
	end

	# [interaction.body2] cutoff: one entry per element pair
	cutoff = DataStructures.OrderedDict{String, Any}()
	for i = 1:length(data.element_symbols)
		for j = i:length(data.element_symbols)
			pair = join(sort([data.element_symbols[i], data.element_symbols[j]]), "-")
			cutoff[pair] = -1
		end
	end

	# [regression]
	regression = DataStructures.OrderedDict{String, Any}()
	regression["datafile"] = "EMBSET"
	regression["weight"]   = 1.0
	regression["alpha"]    = 0
	regression["lambda"]   = 0.0

	# [structure] kd_list / position: expand species counts to per-atom rows
	kd_list = Int[]
	positions = Vector{Vector{Float64}}()
	atom_index = 1
	for (i, count) in enumerate(data.numbers)
		for _ = 1:count
			push!(kd_list, i)
			push!(positions, data.positions[atom_index])
			atom_index += 1
		end
	end

	# [general]
	write(io, "[general]\n")
	for (key, value) in general
		if key == "kd"
			quoted = join(["\"$(v)\"" for v in value], ", ")
			write(io, "kd = [$quoted]\n")
		elseif key == "periodicity"
			joined = join(string.(value), ", ")
			write(io, "periodicity = [$joined]\n")
		elseif key == "name"
			write(io, "name = \"$value\"\n")
		else
			write(io, "$key = $value\n")
		end
	end
	write(io, "\n")

	# [symmetry]
	write(io, "[symmetry]\n")
	for (key, value) in symmetry
		write(io, "$key = $value\n")
	end
	write(io, "\n")

	# [interaction]
	write(io, "[interaction]\n")
	write(io, "nbody = 2\n\n")

	# [interaction.body1]
	write(io, "[interaction.body1]\n")
	for (element, value) in lmax
		write(io, "lmax.$element = $value\n")
	end
	write(io, "\n")

	# [interaction.body2]
	write(io, "[interaction.body2]\n")
	for (pair, value) in cutoff
		write(io, "cutoff.\"$pair\" = $value\n")
	end
	write(io, "lsum = 2\n\n")

	# [regression]
	write(io, "[regression]\n")
	for (key, value) in regression
		if key == "datafile"
			write(io, "$key = \"$value\"\n")
		else
			write(io, "$key = $value\n")
		end
	end
	write(io, "\n")

	# [structure]
	write(io, "[structure]\n")
	write(io, "lattice = [\n")
	for vec in data.lattice_vectors
		write(io, "  [$(join(vec, ", "))],\n")
	end
	write(io, "]\n")

	# kd_list: group consecutive equal values, at most 20 per line
	write(io, "kd_list = [\n")
	current_line = Int[]
	current_value = kd_list[1]
	for value in kd_list
		if value != current_value || length(current_line) >= 20
			if !isempty(current_line)
				write(io, "  $(join(current_line, ", ")),\n")
			end
			current_line = Int[]
			current_value = value
		end
		push!(current_line, value)
	end
	if !isempty(current_line)
		write(io, "  $(join(current_line, ", ")),\n")
	end
	write(io, "]\n")

	write(io, "position = [\n")
	for pos in positions
		formatted_pos = map(x -> @sprintf("%.15f", x), pos)
		write(io, "  [$(join(formatted_pos, ", "))],\n")
	end
	write(io, "]\n")
	return nothing
end

"""
	poscar_to_toml(poscar; output=nothing) -> String

Convert a VASP POSCAR structure file to a Magesty input TOML configuration
and return the TOML text.

The generated configuration is a starting point for an SCE input file: it
fills `[general]`, `[symmetry]`, `[interaction]`, `[regression]`, and
`[structure]` from the POSCAR, with placeholder interaction settings
(`lmax = 0`, `cutoff = -1`) meant to be edited before use.

# Arguments
- `poscar::AbstractString`: path to a POSCAR structure file.

# Keyword arguments
- `output::Union{AbstractString, Nothing} = nothing`: when given, the TOML
  text is also written to this file (`.toml` is appended if the name does
  not already end with it).

# Returns
- `String`: the full TOML text.

# Examples
```julia
text = poscar_to_toml("POSCAR")
poscar_to_toml("POSCAR"; output = "input.toml")
```
"""
function poscar_to_toml(
	poscar::AbstractString;
	output::Union{AbstractString, Nothing} = nothing,
)::String
	data = VaspIO.parse_poscar(poscar)

	buf = IOBuffer()
	_write_input_toml(buf, data)
	text = String(take!(buf))

	if output !== nothing
		outfile = endswith(output, ".toml") ? output : output * ".toml"
		write(outfile, text)
	end
	return text
end
