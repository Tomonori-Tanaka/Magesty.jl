"""
	SpinConfigs

A module for handling spin configurations in magnetic systems.

# Types
- `SpinConfig`: A structure representing a spin configuration with associated properties
- `DataSet`: A structure for managing training and validation data sets

# Functions
- `read_embset`: Read spin configurations from an EMBSET file
- `parse_embset`: Parse EMBSET file and create a DataSet with specified training ratio
"""
module SpinConfigs

using Printf
using LinearAlgebra
using StatsBase
import Base: show

export SpinConfig, read_embset, DataSet, parse_embset

"""
	calc_local_magfield_vertical(spin_directions, local_magfield) -> Matrix{Float64}

Calculate the component of the local magnetic field that is perpendicular to the magnetic moments.

# Arguments
- `spin_directions::Matrix{Float64}`: Spin direction vectors [3 × num_atoms]
- `local_magfield::Matrix{Float64}`: Local magnetic field vectors [3 × num_atoms]

# Returns
- `Matrix{Float64}`: Component of the local magnetic field perpendicular to the magnetic moments [3 × num_atoms]
"""
function calc_local_magfield_vertical(
	spin_directions::Matrix{Float64},
	local_magfield::Matrix{Float64},
)::Matrix{Float64}
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
   atom  magmom     spin_x     spin_y     spin_z     B_x         B_y         B_z         Bv_x        Bv_y        Bv_z
	  1  <value>  <value>    <value>    <value>    <value>     <value>     <value>     <value>     <value>     <value>
	  ...
```
"""
function show(io::IO, config::SpinConfig)
	println(io, "energy (eV): $(config.energy)")
	println(
		io,
		"   atom  magmom     spin_x     spin_y     spin_z     B_x         B_y         B_z         Bv_x        Bv_y        Bv_z",
	)

	num_atoms = length(config.magmom_size)
	for i in 1:num_atoms
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
	DataSet

A structure for managing training and validation data sets.

# Fields
- `spinconfigs::Vector{SpinConfig}`: All spin configurations
- `training_data_num::Int`: Number of training data points
- `validation_data_num::Int`: Number of validation data points
- `training_data_indices::Vector{Int}`: Indices of training data points
- `validation_data_indices::Vector{Int}`: Indices of validation data points
"""
struct DataSet
	spinconfigs::Vector{SpinConfig}
	training_data_num::Int
	validation_data_num::Int
	training_data_indices::Vector{Int}
	validation_data_indices::Vector{Int}

	function DataSet(
		spinconfigs::Vector{SpinConfig},
		training_data_indices::Vector{Int},
		validation_data_indices::Vector{Int},
	)
		return new(
			spinconfigs,
			length(training_data_indices),
			length(validation_data_indices),
			training_data_indices,
			validation_data_indices,
		)
	end
end

"""
	DataSet(spinconfigs::Vector{SpinConfig}, training_ratio::Real = 1.0)

Create a DataSet with specified training ratio.

# Arguments
- `spinconfigs::Vector{SpinConfig}`: Vector of spin configurations
- `training_ratio::Real`: Ratio of data to use for training (0.0 < ratio ≤ 1.0, default: 1.0)

# Returns
- `DataSet`: A DataSet containing the configurations split into training and validation sets

# Throws
- `ArgumentError` if training_ratio is not in the range (0, 1]
- `ArgumentError` if no configurations are provided
"""
function DataSet(spinconfigs::Vector{SpinConfig}, training_ratio::Real = 1.0)
	if training_ratio <= 0.0 || training_ratio > 1.0
		throw(ArgumentError("training_ratio must be in the range (0, 1]"))
	end
	num_configs = length(spinconfigs)
	if num_configs < 1
		throw(ArgumentError("At least one configuration is required"))
	end

	# When training_ratio is 1.0, use all data for training
	if training_ratio == 1.0
		training_data_indices = collect(1:num_configs)
		validation_data_indices = Int[]
	else
		training_data_num = max(1, Int(floor(num_configs * training_ratio)))
		training_data_indices = sample(1:num_configs, training_data_num, replace = false)
		validation_data_indices = setdiff(1:num_configs, training_data_indices)
	end

	return DataSet(spinconfigs, training_data_indices, validation_data_indices)
end

"""
	parse_embset(filename::AbstractString, num_atoms::Integer; training_ratio::Real = 1.0)

Parse EMBSET file and create a DataSet with specified training ratio.

# Arguments
- `filename::AbstractString`: Path to the EMBSET file
- `num_atoms::Integer`: Number of atoms in the system
- `training_ratio::Real`: Ratio of data to use for training (default: 1.0)

# Returns
- `DataSet`: A DataSet containing the parsed configurations

# Throws
- `ErrorException` if the file format is invalid
"""
function parse_embset(filename::AbstractString, num_atoms::Integer; training_ratio::Real = 1.0)
	spinconfigs = read_embset(filename, num_atoms)
	return DataSet(spinconfigs, training_ratio)
end

"""
	parse_embset(filename::AbstractString, num_atoms::Integer, use_data_indices::Vector{Int}; training_ratio::Real = 1.0)

Parse EMBSET file and create a DataSet using only specified configurations.

# Arguments
- `filename::AbstractString`: Path to the EMBSET file
- `num_atoms::Integer`: Number of atoms in the system
- `use_data_indices::Vector{Int}`: Indices of configurations to use
- `training_ratio::Real`: Ratio of data to use for training (default: 1.0)

# Returns
- `DataSet`: A DataSet containing the selected configurations

# Throws
- `ErrorException` if the file format is invalid
- `BoundsError` if any index in use_data_indices is out of range
"""
function parse_embset(
	filename::AbstractString,
	num_atoms::Integer,
	use_data_indices::Vector{Int};
	training_ratio::Real = 1.0,
)
	spinconfigs_orig = read_embset(filename, num_atoms)
	spinconfigs = spinconfigs_orig[use_data_indices]
	return DataSet(spinconfigs, training_ratio)
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
"""
function read_embset(filepath::AbstractString, num_atoms::Integer)::Vector{SpinConfig}
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
	for (i, line) in enumerate(filtered_lines[start_line+1:end_line])
		line_split = split(line)
		moment = parse.(Float64, line_split[2:4])
		magmom_size[i] = norm(moment)
		spin_directions[:, i] = moment ./ magmom_size[i]
		local_magfield[:, i] = parse.(Float64, line_split[5:7])
	end

	return SpinConfig(energy, magmom_size, spin_directions, local_magfield)
end
end