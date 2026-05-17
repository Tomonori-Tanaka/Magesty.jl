"""
	SpinConfigs

A module for handling spin configurations in magnetic systems.

# Types
- `SpinConfig`: A structure representing a spin configuration with associated properties

# Functions
- `read_embset`: Read spin configurations from an EMBSET file
"""
module SpinConfigs

using Printf
using LinearAlgebra
import Base: show

export SpinConfig, read_embset

"""
	calc_local_magfield_vertical(spin_directions, local_magfield) -> Matrix{Float64}

Calculate the component of the local magnetic field that is perpendicular to the magnetic moments.

# Arguments
- `spin_directions::Matrix{Float64}`: Spin direction vectors [3 × num_atoms]
- `local_magfield::Matrix{Float64}`: Local magnetic field vectors [3 × num_atoms]

# Returns
- `Matrix{Float64}`: Component of the local magnetic field perpendicular to the magnetic moments [3 × num_atoms]

# Throws
- `ArgumentError` if dimensions of input matrices do not match
"""
function calc_local_magfield_vertical(
	spin_directions::Matrix{Float64},
	local_magfield::Matrix{Float64},
)::Matrix{Float64}
	# Validate input dimensions
	if size(spin_directions) != size(local_magfield)
		throw(ArgumentError("Dimensions of spin_directions and local_magfield must match"))
	end

	local_magfield_vertical = zeros(3, size(local_magfield, 2))
	for i = 1:size(local_magfield, 2)
		proj_B = dot(local_magfield[:, i], spin_directions[:, i]) * spin_directions[:, i]
		local_magfield_vertical[:, i] = local_magfield[:, i] - proj_B
	end
	return local_magfield_vertical
end

"""
	SpinConfig

A single noncollinear spin configuration plus the DFT observables that
go with it. Typically produced by `read_embset`, then handed to
`SCEDataset`, `predict_energy`, `predict_torque`, or the evaluation
verbs (`r2_energy`, `rmse_torque`, …).

# Fields
- `energy::Float64`: DFT total energy of the configuration [eV].
- `magmom_size::Vector{Float64}`: Magnitude of the magnetic moment on
  each atom [μ_B]. Length is `num_atoms`.
- `spin_directions::Matrix{Float64}`: Unit-vector spin directions, laid
  out as `3 × num_atoms` (rows = x, y, z; each column has unit norm).
- `local_magfield::Matrix{Float64}`: Local constraining magnetic field
  at each atom [eV / μ_B], same `3 × num_atoms` layout as
  `spin_directions`. The values come from VASP's
  `lambda*MW_perp` block, which has dimensions of energy per magnetic
  moment so that `E = −m · B` is in eV — not the Tesla a conventional
  magnetic field would carry.
- `local_magfield_vertical::Matrix{Float64}`: Component of
  `local_magfield` perpendicular to the spin direction at each atom
  [eV / μ_B]. Computed at construction; the parallel component is
  dropped because classical-spin torque depends only on the
  perpendicular field.
- `torques::Matrix{Float64}`: Per-atom torque vectors
  `τᵢ = −mᵢ (eᵢ × Bᵢ)` (where `mᵢ` is the moment magnitude, `eᵢ` the
  unit-vector spin direction, `Bᵢ` the local field), laid out as
  `3 × num_atoms`. These are the observables compared against
  `predict_torque` during fitting and evaluation.

# Constructors
- `SpinConfig(energy, magmom_size, spin_directions, local_magfield)` —
  `local_magfield_vertical` and `torques` are computed from the other
  four fields.

# Throws
- `ArgumentError` if `spin_directions` and `local_magfield` do not have
  matching `3 × num_atoms` shapes.
- `ArgumentError` if any entry of `magmom_size` is negative.
"""
struct SpinConfig
	energy::Float64
	magmom_size::Vector{Float64}
	spin_directions::Matrix{Float64}
	local_magfield::Matrix{Float64}
	local_magfield_vertical::Matrix{Float64}
	torques::Matrix{Float64}

	function SpinConfig(
		energy::Real,
		magmom_size::AbstractVector{<:Real},
		spin_directions::AbstractMatrix{<:Real},
		local_magfield::AbstractMatrix{<:Real},
	)
		num_atoms = length(magmom_size)

		# Validate dimensions
		if size(spin_directions, 1) != 3
			throw(
				ArgumentError(
					"spin_directions must have 3 rows (x, y, z); got $(size(spin_directions, 1))",
				),
			)
		end
		if size(local_magfield, 1) != 3
			throw(
				ArgumentError(
					"local_magfield must have 3 rows (x, y, z); got $(size(local_magfield, 1))",
				),
			)
		end
		if num_atoms != size(spin_directions, 2)
			throw(
				ArgumentError(
					"Dimension mismatch: spin_directions has $(size(spin_directions, 2)) columns, " *
					"but magmom_size has $num_atoms elements",
				),
			)
		end
		if num_atoms != size(local_magfield, 2)
			throw(
				ArgumentError(
					"Dimension mismatch: local_magfield has $(size(local_magfield, 2)) columns, " *
					"but magmom_size has $num_atoms elements",
				),
			)
		end

		# Validate magnetic moment sizes
		if any(x -> x < 0, magmom_size)
			throw(ArgumentError("Magnetic moment sizes must be non-negative"))
		end

		# Calculate vertical component of local magnetic field
		local_magfield_vertical = calc_local_magfield_vertical(spin_directions, local_magfield)

		new(
			Float64(energy),
			Float64.(magmom_size),
			Float64.(spin_directions),
			Float64.(local_magfield),
			local_magfield_vertical,
			calc_torques(magmom_size, spin_directions, local_magfield),
		)
	end
end

"""
	calc_torques(magmom_size::AbstractVector{<:Real}, spin_directions::AbstractMatrix{<:Real}, local_magfield::AbstractMatrix{<:Real})::Matrix{Float64}

Calculate the torques for each atom in the spin configuration.

# Arguments
- `spin_directions::AbstractMatrix{<:Real}`: Spin direction vectors [3 × num_atoms]
- `local_magfield::AbstractMatrix{<:Real}`: Local magnetic field vectors [3 × num_atoms]

# Returns
- `Matrix{Float64}`: Torques for each atom [3 × num_atoms]
"""
function calc_torques(
	magmom_size::AbstractVector{<:Real},
	spin_directions::AbstractMatrix{<:Real},
	local_magfield::AbstractMatrix{<:Real},
)::Matrix{Float64}
	num_atoms = length(magmom_size)
	torques = zeros(3, num_atoms)
	for i = 1:num_atoms
		torques[:, i] = -1.0 * (magmom_size[i] * spin_directions[:, i]) × local_magfield[:, i]
	end
	return torques
end

"""
	show(io::IO, ::MIME"text/plain", config::SpinConfig)

Display a spin configuration in a human-readable, multi-line format.

By overloading the three-argument `show` method (with `MIME"text/plain"`),
the rich display is used by the REPL while collections such as
`Vector{SpinConfig}` and contexts like `print(config)` keep the default
compact representation, matching Julia's standard `show` conventions.

# Arguments
- `io::IO`: The output stream
- `::MIME"text/plain"`: Selects the rich multi-line representation
- `config::SpinConfig`: The spin configuration to display

# Output Format
```
energy (eV): <energy>
   atom  magmom     e_x     e_y     e_z     B_x         B_y         B_z         Bv_x        Bv_y        Bv_z
	  1  <value>  <value>    <value>    <value>    <value>     <value>     <value>     <value>     <value>     <value>
	  ...
```
"""
function show(io::IO, ::MIME"text/plain", config::SpinConfig)
	println(io, "energy (eV): $(config.energy)")
	println(
		io,
		"   atom  magmom     e_x     e_y     e_z     B_x         B_y         B_z         Bv_x        Bv_y        Bv_z",
	)

	num_atoms = length(config.magmom_size)
	@inbounds for i = 1:num_atoms
		println(io,
			lpad(i, 7),
			"  ",
			lpad(@sprintf("%.5f", config.magmom_size[i]), 7),
			"  ",
			lpad(@sprintf("%.6f", config.spin_directions[1, i]), 9),
			" ",
			lpad(@sprintf("%.6f", config.spin_directions[2, i]), 9),
			" ",
			lpad(@sprintf("%.6f", config.spin_directions[3, i]), 9),
			"  ",
			lpad(@sprintf("%.5e", config.local_magfield[1, i]), 12),
			" ",
			lpad(@sprintf("%.5e", config.local_magfield[2, i]), 12),
			" ",
			lpad(@sprintf("%.5e", config.local_magfield[3, i]), 12),
			" ",
			lpad(@sprintf("%.5e", config.local_magfield_vertical[1, i]), 12),
			" ",
			lpad(@sprintf("%.5e", config.local_magfield_vertical[2, i]), 12),
			" ",
			lpad(@sprintf("%.5e", config.local_magfield_vertical[3, i]), 12),
		)
	end
end

"""
	read_embset(filepath::AbstractString) -> Vector{SpinConfig}

Read spin configurations from an EMBSET file.

# Arguments
- `filepath::AbstractString`: Path to the EMBSET file

# Returns
- `Vector{SpinConfig}`: Array of spin configurations

# Throws
- `ErrorException` if the file format is invalid
- `ArgumentError` if the file does not exist
"""
function read_embset(filepath::AbstractString)::Vector{SpinConfig}
	if !isfile(filepath)
		throw(ArgumentError("File does not exist: $filepath"))
	end

	# Detect number of atoms from file
	num_atoms = detect_num_atoms(filepath)

	# Read and filter lines
	filtered_lines = String[]
	open(filepath, "r") do file
		for line in eachline(file)
			stripped_line = strip(line)
			if !isempty(stripped_line) && !startswith(stripped_line, "#")
				push!(filtered_lines, stripped_line)
			end
		end
	end

	# Validate file format
	if length(filtered_lines) % (num_atoms + 1) != 0
		throw(
			ErrorException(
				"Invalid EMBSET file format: Number of lines ($(length(filtered_lines))) " *
				"is not divisible by $(num_atoms + 1)",
			),
		)
	end

	# Process configurations
	num_configs = length(filtered_lines) ÷ (num_atoms + 1)
	configs = Vector{SpinConfig}(undef, num_configs)

	for i = 1:num_configs
		configs[i] = separate_embset(filtered_lines, num_atoms, i)
	end

	return configs
end

"""
	detect_num_atoms(filepath::AbstractString) -> Integer

Detect the number of atoms from an EMBSET file by counting the number of atom data lines in the first spin configuration.

# Arguments
- `filepath::AbstractString`: Path to the EMBSET file

# Returns
- `Integer`: Number of atoms detected from the file

# Throws
- `ErrorException` if the file format is invalid
"""
function detect_num_atoms(filepath::AbstractString)::Integer
	open(filepath, "r") do file
		# Find first comment line (start of first config)
		found_first_comment = false
		skipped_energy_line = false
		atom_count = 0
		
		for line in eachline(file)
			stripped_line = strip(line)
			
			# Skip empty lines
			if isempty(stripped_line)
				continue
			end
			
			# Find first comment line
			if startswith(stripped_line, "#")
				if found_first_comment && atom_count > 0
					# Reached next config, stop counting
					break
				end
				found_first_comment = true
				skipped_energy_line = false
				atom_count = 0
				continue
			end
			
			# After finding first comment, skip the energy line (first non-comment line)
			if found_first_comment && !skipped_energy_line
				skipped_energy_line = true
				continue
			end
			
			# Count atom data lines until next comment line
			if found_first_comment && skipped_energy_line
				atom_count += 1
			end
		end
		
		if atom_count == 0
			throw(ErrorException("Could not detect number of atoms from EMBSET file: no atom data found"))
		end
		
		return atom_count
	end
end

"""
	separate_embset(filtered_lines, num_atoms, data_index) -> SpinConfig

Extract a single spin configuration from filtered EMBSET lines.

# Arguments
- `filtered_lines::AbstractVector{<:AbstractString}`: Filtered lines from EMBSET file
- `num_atoms::Integer`: Number of atoms in the structure
- `data_index::Integer`: Index of the configuration to extract

# Returns
- `SpinConfig`: The extracted spin configuration

# Throws
- `ErrorException` if the file format is invalid
"""
function separate_embset(
	filtered_lines::AbstractVector{<:AbstractString},
	num_atoms::Integer,
	data_index::Integer,
)::SpinConfig
	start_line = (num_atoms + 1) * (data_index - 1) + 1
	end_line = (num_atoms + 1) * data_index

	# Parse energy
	energy = parse(Float64, filtered_lines[start_line])

	# Initialize arrays
	magmom_size = Vector{Float64}(undef, num_atoms)
	spin_directions = Matrix{Float64}(undef, 3, num_atoms)
	local_magfield = Matrix{Float64}(undef, 3, num_atoms)

	# Parse atom data
	for (i, line) in enumerate(filtered_lines[(start_line+1):end_line])
		line_split = split(line)
		if length(line_split) < 7
			throw(ErrorException("Invalid line format in EMBSET file"))
		end
		moment = parse.(Float64, line_split[2:4])
		magmom_size[i] = norm(moment)
		spin_directions[:, i] = moment ./ magmom_size[i]
		local_magfield[:, i] = parse.(Float64, line_split[5:7])
	end

	return SpinConfig(energy, magmom_size, spin_directions, local_magfield)
end

# Compact REPL display. Default Julia output dumps every per-atom matrix.
function Base.show(io::IO, sc::SpinConfig)
	print(io, "SpinConfig(num_atoms=", length(sc.magmom_size),
		", energy=", sc.energy, ")")
end

end