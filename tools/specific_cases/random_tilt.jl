include("../vasptools.jl")

using .VaspTools
using Random
using ArgParse
using LinearAlgebra

"""
	random_tilt(
		incar::String,
		num_atoms::Integer,
		theta_range::Tuple{Float64, Float64},
		num_patterns::Integer,
		start_number::Integer = 1,
		;
		file_prefix::String = "pattern",
	)

Generate random tilt MAGMOM and M_CONSTR patterns for a VASP calculation.

# Arguments
- `incar::String`: Path to the input INCAR file
- `num_atoms::Integer`: Number of atoms in the structure
- `theta_range::Tuple{Float64, Float64}`: Range of tilt angles in degrees from the original magnetic moment direction (min, max)
- `num_patterns::Integer`: Number of random patterns to generate
- `start_number::Integer`: Starting number for pattern files (default: 1)
- `file_prefix::String`: Prefix for output files (default: "pattern")

# Returns
- Generates `num_patterns` INCAR files with random magnetic moment orientations
- Each file is named as `{file_prefix}_{pattern_number}.incar`
- Pattern numbers are zero-padded based on the total number of patterns

# Notes
- The tilt angle θ is defined as the angle from the original magnetic moment direction
- The azimuthal angle φ is defined as the angle around the original magnetic moment direction
- The magnitude of each magnetic moment is preserved
- The rotation axis is randomly chosen perpendicular to the original magnetic moment
- Atoms with zero magnetic moment are kept unchanged
- The original magnetic moment direction is taken from M_CONSTR in the input INCAR file
- Both MAGMOM and M_CONSTR in the output files are set to the same values
"""
function random_tilt(
	incar::String,
	num_atoms::Integer,
	theta_range::Tuple{Float64, Float64},
	num_patterns::Integer,
	start_number::Integer = 1,
	;
	file_prefix::String = "pattern",
)
	# Parse INCAR file
	incar_dict = parse_incar(incar)

	# Get M_CONSTR values
	m_constr = incar_dict[:M_CONSTR]

	# Calculate required number of digits for zero-padding
	digits = length(string(start_number + num_patterns - 1))

	# Convert theta_range from degrees to radians
	theta_range_rad = (theta_range[1] * π / 180, theta_range[2] * π / 180)

	# Generate random patterns
	for i in start_number:(start_number + num_patterns - 1)
		new_magmom = Vector{Float64}()
		new_m_constr = Vector{Float64}()
		sizehint!(new_magmom, 3 * num_atoms)
		sizehint!(new_m_constr, 3 * num_atoms)

		# Generate random orientations for each atom
		for j in 1:num_atoms
			# Get original magnetic moment components
			m_orig = [m_constr[3*(j-1)+1], m_constr[3*(j-1)+2], m_constr[3*(j-1)+3]]

			# Calculate magnitude of original magnetic moment
			m = norm(m_orig)

			if m ≈ 0.0
				# Keep zero magnetic moment atoms unchanged
				append!(new_magmom, m_orig)
				append!(new_m_constr, m_orig)
			else
				# Normalize original magnetic moment
				m_hat = m_orig / m

				# Generate random tilt angle and azimuthal angle
				theta = rand() * (theta_range_rad[2] - theta_range_rad[1]) + theta_range_rad[1]
				phi = rand() * 2π

				# Find a vector perpendicular to m_hat
				if abs(m_hat[1]) < 0.5
					v = [1.0, 0.0, 0.0]
				else
					v = [0.0, 1.0, 0.0]
				end

				# Make it perpendicular to m_hat using Gram-Schmidt
				axis1 = v - dot(v, m_hat) * m_hat
				axis1 = axis1 / norm(axis1)

				# Create a second perpendicular axis
				axis2 = cross(m_hat, axis1)

				# Create rotation axis using azimuthal angle
				axis = cos(phi) * axis1 + sin(phi) * axis2

				# Rodrigues' rotation formula
				cos_theta = cos(theta)
				sin_theta = sin(theta)
				m_new =
					cos_theta * m_orig + sin_theta * cross(axis, m_orig) +
					(1 - cos_theta) * dot(axis, m_orig) * axis

				# Append new values
				append!(new_magmom, m_new)
				append!(new_m_constr, m_new)
			end
		end

		# Create new INCAR
		new_incar_dict = copy(incar_dict)
		new_incar_dict[:MAGMOM] = new_magmom
		new_incar_dict[:M_CONSTR] = new_m_constr

		# Write to file with zero-padded pattern number
		output_file = "$(file_prefix)_$(lpad(i, digits, "0")).incar"
		write_incar(output_file, new_incar_dict)
	end
end

"""
	write_conditions(
		output_file::String,
		incar::String,
		num_atoms::Integer,
		theta_range::Tuple{Float64, Float64},
		num_patterns::Integer,
		start_number::Integer,
		file_prefix::String,
	)

Write calculation conditions to a specified file.

# Arguments
- `output_file::String`: Path to the output file
- `incar::String`: Path to the input INCAR file
- `num_atoms::Integer`: Number of atoms in the structure
- `theta_range::Tuple{Float64, Float64}`: Range of tilt angles in degrees
- `num_patterns::Integer`: Number of random patterns to generate
- `start_number::Integer`: Starting number for pattern files
- `file_prefix::String`: Prefix for output files
"""
function write_conditions(
	output_file::String,
	incar::String,
	num_atoms::Integer,
	theta_range::Tuple{Float64, Float64},
	num_patterns::Integer,
	start_number::Integer,
	file_prefix::String,
)
	open(output_file, "w") do io
		write(io, "Calculation Conditions\n")
		write(io, "=====================\n\n")
		write(io, "Input INCAR file: $incar\n")
		write(io, "Number of atoms: $num_atoms\n")
		write(io, "Tilt angle range (degrees): $(theta_range[1]) to $(theta_range[2])\n")
		write(io, "Number of patterns: $num_patterns\n")
		write(io, "Start number: $start_number\n")
		write(io, "Output file prefix: $file_prefix\n")
		write(io, "\nGenerated files:\n")
		for i in start_number:(start_number + num_patterns - 1)
			digits = length(string(start_number + num_patterns - 1))
			write(io, "- $(file_prefix)_$(lpad(i, digits, "0")).incar\n")
		end
	end
end

function parse_commandline()
	s = ArgParseSettings(
		description = "Generate random tilt MAGMOM and M_CONSTR patterns for VASP calculations.",
		version = "1.0.0",
		add_version = true,
	)

	@add_arg_table! s begin
		"--incar", "-i"
			help = "Path to the input INCAR file"
			required = true
			arg_type = String

		"--num-atoms", "-n"
			help = "Number of atoms in the structure"
			required = true
			arg_type = Int

		"--theta-min", "-t"
			help = "Minimum tilt angle in degrees from the original magnetic moment direction"
			required = true
			arg_type = Float64

		"--theta-max", "-T"
			help = "Maximum tilt angle in degrees from the original magnetic moment direction"
			required = true
			arg_type = Float64

		"--num-patterns", "-p"
			help = "Number of random patterns to generate"
			required = true
			arg_type = Int

		"--start-number", "-s"
			help = "Starting number for pattern files (default: 1)"
			default = 1
			arg_type = Int

		"--prefix", "-f"
			help = "Prefix for output files"
			default = "pattern"
			arg_type = String

		"--conditions", "-c"
			help = "Path to write calculation conditions (optional)"
			arg_type = String
	end

	return parse_args(s)
end

function main()
	# Parse command line arguments
	args = parse_commandline()
	
	# Call random_tilt function with parsed arguments
	random_tilt(
		args["incar"],
		args["num-atoms"],
		(args["theta-min"], args["theta-max"]),
		args["num-patterns"],
		args["start-number"],
		file_prefix = args["prefix"],
	)

	# Write calculation conditions if specified
	if haskey(args, "conditions") && args["conditions"] !== nothing
		write_conditions(
			args["conditions"],
			args["incar"],
			args["num-atoms"],
			(args["theta-min"], args["theta-max"]),
			args["num-patterns"],
			args["start-number"],
			args["prefix"],
		)
	end
end

# Execute main function if script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
