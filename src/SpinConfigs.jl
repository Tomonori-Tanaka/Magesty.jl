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
	for i in 1:size(local_magfield, 2)
		proj_B = dot(local_magfield[:, i], spin_directions[:, i]) * spin_directions[:, i]
		local_magfield_vertical[:, i] = local_magfield[:, i] - proj_B
	end
	return local_magfield_vertical
end

"""
	SpinConfig

A configuration of spins in a magnetic structure.

# Fields
- `energy::Float64`: The energy of the spin configuration [eV]
- `magmom_size::Vector{Float64}`: The magnitude of magnetic moments for each atom [μ_B]
- `spin_directions::Matrix{Float64}`: The direction cosines (unit vectors) of the spins [3 × num_atoms]
- `local_magfield::Matrix{Float64}`: The local magnetic field at each atom [T]
- `local_magfield_vertical::Matrix{Float64}`: The component of the local magnetic field perpendicular to the magnetic moments [T]

# Constructors
- `SpinConfig(energy, magmom_size, spin_directions, local_magfield)`: Create a spin configuration

# Throws
- `ArgumentError` if dimensions of input arrays do not match
- `ArgumentError` if any magnetic moment size is negative
"""
struct SpinConfig
	energy::Float64
	magmom_size::Vector{Float64}
	spin_directions::Matrix{Float64}
	local_magfield::Matrix{Float64}
	local_magfield_vertical::Matrix{Float64}

	function SpinConfig(
		energy::Real,
		magmom_size::AbstractVector{<:Real},
		spin_directions::AbstractMatrix{<:Real},
		local_magfield::AbstractMatrix{<:Real},
	)
		num_atoms = length(magmom_size)

		# Validate dimensions
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
		)
	end
end

"""
	show(io::IO, config::SpinConfig)

Display a spin configuration in a human-readable format.

# Arguments
- `io::IO`: The output stream
- `config::SpinConfig`: The spin configuration to display

# Output Format
```
energy (eV): <energy>
   atom  magmom     e_x     e_y     e_z     B_x         B_y         B_z         Bv_x        Bv_y        Bv_z
	  1  <value>  <value>    <value>    <value>    <value>     <value>     <value>     <value>     <value>     <value>
	  ...
```
"""
function show(io::IO, config::SpinConfig)
	println(io, "energy (eV): $(config.energy)")
	println(
		io,
		"   atom  magmom     e_x     e_y     e_z     B_x         B_y         B_z         Bv_x        Bv_y        Bv_z",
	)

	num_atoms = length(config.magmom_size)
	@inbounds for i in 1:num_atoms
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
	read_embset(filepath::AbstractString, num_atoms::Integer) -> Vector{SpinConfig}

Read spin configurations from an EMBSET file.

# Arguments
- `filepath::AbstractString`: Path to the EMBSET file
- `num_atoms::Integer`: Number of atoms in the structure

# Returns
- `Vector{SpinConfig}`: Array of spin configurations

# Throws
- `ErrorException` if the file format is invalid
- `ArgumentError` if the file does not exist
"""
function read_embset(filepath::AbstractString, num_atoms::Integer)::Vector{SpinConfig}
	if !isfile(filepath)
		throw(ArgumentError("File does not exist: $filepath"))
	end

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

	for i in 1:num_configs
		configs[i] = separate_embset(filtered_lines, num_atoms, i)
	end

	return configs
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
end