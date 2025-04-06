"""
	SpinConfigs

A module for spin configurations.
"""
module SpinConfigs

using Printf
using LinearAlgebra
import Base: show
export SpinConfig, read_embset
"""
	SpinConfig

A configuration of spins.

# Fields
- `energy::Float64`: The energy of the spin configuration. [num_atoms]
- `magmom_size::Vector{Float64}`: The size of the magnetic moments. [3, num_atoms]
- `spin_directions::Matrix{Float64}`: The direction cosines (unit vectors) of the spins. [3, num_atoms]
- `local_magfield::Matrix{Float64}`: The local magnetic field. [3, num_atoms]
- `torque::Matrix{Float64}`: The torque acting on each spin, calculated as cross product of spin direction and local magnetic field. [3, num_atoms]
"""
struct SpinConfig
	energy::Float64
	magmom_size::Vector{Float64}
	spin_directions::Matrix{Float64}
	local_magfield::Matrix{Float64}
	local_magfield_vertical::Matrix{Float64}
end

function SpinConfig(
	energy::Real,
	magmom_size::AbstractVector{<:Real},
	spin_directions::AbstractMatrix{<:Real},
	local_magfield::AbstractMatrix{<:Real},
)::SpinConfig
	num_atoms = length(magmom_size)

	if num_atoms != size(spin_directions, 2)
		throw(
			ArgumentError(
				"Sizes mismatch: spin_directions columns = $(size(spin_directions, 2)), magmom_size length = $(length(magmom_size))",
			),
		)
	elseif num_atoms != size(local_magfield, 2)
		throw(
			ArgumentError(
				"Sizes mismatch: local_magfield columns = $(size(local_magfield, 2)), magmom_size length = $(length(magmom_size))",
			),
		)
	end

	local_magfield_vertical = calc_local_magfield_vertical(spin_directions, local_magfield)

	return SpinConfig(energy, magmom_size, spin_directions, local_magfield, local_magfield_vertical)
end

function show(io::IO, config::SpinConfig)::Nothing
	println(io, "energy (eV): $(config.energy)")
	num_atoms = length(config.magmom_size)
	for i in 1:num_atoms
		println(
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

function read_embset(filepath::AbstractString, num_atoms::Integer)::Vector{SpinConfig}
	filtered_lines = String[]
	open(filepath, "r") do file
		for line in eachline(file)
			stripped_line = strip(line)
			if !isempty(stripped_line) && !startswith(stripped_line, "#")
				push!(filtered_lines, stripped_line)
			end
		end
	end

	if size(filtered_lines, 1) % (num_atoms + 1) != 0
		error(
			"EMBSET file '$filepath' is not compatible: Number of lines $(size(filtered_lines, 1)) is not divisible by $(num_atoms + 1).",
		)
	end
	num_configs::Int = size(filtered_lines, 1) ÷ (num_atoms + 1)

	configs = SpinConfig[]
	for i in 1:num_configs
		push!(configs, separate_embset(filtered_lines, num_atoms, i))
	end
	return configs
end

function separate_embset(
	filtered_lines::AbstractVector{<:AbstractString},
	num_atoms::Integer,
	data_index::Integer,
)::SpinConfig
	start_line::Int = (num_atoms + 1) * (data_index - 1) + 1
	end_line::Int = (num_atoms + 1) * data_index

	magmom_size = Vector{Float64}(undef, num_atoms)
	spin_directions = Matrix{Float64}(undef, 3, num_atoms)
	local_magfield = Matrix{Float64}(undef, 3, num_atoms)

	energy = parse(Float64, filtered_lines[start_line])
	for (i, line) in enumerate(filtered_lines[start_line+1:end_line])
		line_split = split(line)
		moment::Vector{Float64} = map(x -> parse(Float64, x), line_split[2:4])
		magmom_size[i] = norm(moment)
		spin_directions[:, i] = moment ./ magmom_size[i]
		magfield = map(x -> parse(Float64, x), line_split[5:7])
		local_magfield[:, i] = magfield
	end

	return SpinConfig(energy, magmom_size, spin_directions, local_magfield)
end

function calc_local_magfield_vertical(
	spin_directions::Matrix{Float64},
	local_magfield::Matrix{Float64}
)::Matrix{Float64}
	local_magfield_vertical = zeros(3, size(local_magfield, 2))
	for i in 1:size(local_magfield, 2)
		proj_B::Vector{Float64} = dot(local_magfield[:, i], spin_directions[:, i]) * spin_directions[:, i]
		local_magfield_vertical[:, i] = local_magfield[:, i] - proj_B
	end
	return local_magfield_vertical
end

end
