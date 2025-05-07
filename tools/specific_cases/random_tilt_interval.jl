include("../vasptools.jl")

using .VaspTools
using Random
using ArgParse
using LinearAlgebra

"""
	random_tilt_interval(
		incar::String,
		num_atoms::Integer,
		theta_range::Tuple{Float64, Float64},
		interval::Float64,
		start_number::Integer = 1,
		;
		file_prefix::String = "pattern",
	)

Generate random tilt MAGMOM and M_CONSTR patterns for a VASP calculation.

# Arguments
- `incar::String`: Path to the input INCAR file
- `num_atoms::Integer`: Number of atoms in the structure
- `theta_range::Tuple{Float64, Float64}`: Range of tilt angles in degrees from the original magnetic moment direction (min, max)
- `interval::Float64`: Interval of tilt angles in degrees
- `start_number::Integer`: Starting number for pattern files (default: 1)
- `file_prefix::String`: Prefix for pattern files (default: "pattern")
- `random_axis::Bool`: Whether to randomize the axis of initial spins (default: false)
"""
function random_tilt_interval(
	incar::String,
	num_atoms::Integer,
	theta_range::Tuple{Float64, Float64},
	interval::Float64,
	start_number::Integer = 1,
	;
	file_prefix::String = "pattern_",
	random_axis::Bool = false,
)

	# Read INCAR file
	incar_dict = parse_incar(incar)

	# Get M_CONSTR values
	m_constr = incar_dict[:M_CONSTR]

	# Calculate the number of data points
	num_data = length(collect(theta_range[1]:interval:theta_range[2]))

	# convert theta list to radians
	theta_list_rad = deg2rad.(collect(theta_range[1]:interval:theta_range[2]))

	# Generate random patterns
	theta_min = theta_list_rad[begin]
	for (i, theta_max) in enumerate(theta_list_rad)
		new_magmom = Matrix{Float64}(undef, 3, num_atoms)

		# Generate random orientations for each atom
		for j in 1:num_atoms
			# Get original magnetic moment components
			m_orig = [m_constr[3*(j-1)+1], m_constr[3*(j-1)+2], m_constr[3*(j-1)+3]]

			# Calculate magnitude of original magnetic moment
			m = norm(m_orig)

			if m ≈ 0.0
				# Keep zero magnetic moment atoms unchanged
				new_magmom[:, j] = m_orig
			else
				# Normalize original magnetic moment
				m_hat = m_orig / m

				# Generate random tilt angle and azimuthal angle
				phi = rand() * 2π
				theta = rand() * (theta_max - theta_min) + theta_min

				# Generate random vector for Gram-Schmidt
				v = rand(3)
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
				new_magmom[:, j] = m_new
			end
		end

		if random_axis
			# Randomize the axis of initial spins
			desired_axis = rand(3)
			desired_axis = desired_axis / norm(desired_axis)
			# Calculate angle between desired_axis and [0,0,1]
			angle = acos(clamp(dot(desired_axis, [0.0, 0.0, 1.0]), -1.0, 1.0))
			# temporary axis for Rodrigues' rotation formula
			rotation_axis = cross([0.0, 0.0, 1.0], desired_axis)
			# Construct the rotation matrix using Rodrigues' rotation formula.
			K = [     0.0      -rotation_axis[3]   rotation_axis[2];
				rotation_axis[3]    0.0      -rotation_axis[1];
					-rotation_axis[2]   rotation_axis[1]    0.0]
			R = I + sin(angle)*K + (1-cos(angle))*(K*K)
			new_magmom = R * new_magmom
		end

		# Create new INCAR
		new_incar_dict = copy(incar_dict)
        # flatten the new_magmom matrix
        new_magmom_list = vec(new_magmom)
		new_incar_dict[:MAGMOM] = new_magmom_list
		new_incar_dict[:M_CONSTR] = new_magmom_list

		# Write to file with zero-padded pattern number
		output_file = "$(file_prefix)$(lpad(i+start_number-1, length(string(num_data+start_number-1)), "0")).incar"
		write_incar(output_file, new_incar_dict)
	end
end


function write_conditions(
	output_file::String,
	incar::String,
	num_atoms::Integer,
	theta_range::Tuple{Float64, Float64},
	interval::Float64,
	start_number::Integer,
	file_prefix::String,
	random_axis::Bool,
)
	# Input validation
	if !isfile(incar)
		error("Input INCAR file does not exist: $incar")
	end
	if num_atoms <= 0
		error("Number of atoms must be positive")
	end
	if theta_range[1] > theta_range[2]
		error("Minimum angle must be less than maximum angle")
	end
	if interval <= 0
		error("Interval must be positive")
	end

	try
		open(output_file, "w") do io
			write(io, "Calculation Conditions\n")
			write(io, "=====================\n\n")
			write(io, "Input INCAR file: $incar\n")
			write(io, "Number of atoms: $num_atoms\n")
			write(io, "Tilt angle range (degrees): $(theta_range[1]) to $(theta_range[2])\n")
			write(io, "Interval of tilt angles (degrees): $interval\n")
			write(io, "Start number: $start_number\n")
			write(io, "File prefix: $file_prefix\n")
			write(io, "Random axis: $random_axis\n")
		end
	catch
		error("Error writing to output file: $output_file")
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

		"--interval", "-I"
		help = "Interval of tilt angles in degrees"
		required = true
		arg_type = Float64

		"--start-number", "-s"
		help = "Starting number for pattern files"
		arg_type = Int
		default = 1

		"--file-prefix", "-f"
		help = "Prefix for pattern files"
		arg_type = String
		default = "pattern_"

		"--random-axis", "-r"
		help = "Randomize the axis of initial spins"
		action = :store_true

		"--conditions", "-c"
		help = "Write calculation conditions to a file"
		arg_type = String
	end
	return parse_args(s)
end

function main()
	# Parse command line arguments
	args = parse_commandline()

	# Call random_tilt_interval function with parsed arguments
	random_tilt_interval(
		args["incar"],
		args["num-atoms"],
		(args["theta-min"], args["theta-max"]),
		args["interval"],
		args["start-number"],
		file_prefix = args["file-prefix"],
		random_axis = args["random-axis"],
	)

	# Write calculation conditions if specified
	if haskey(args, "conditions") && args["conditions"] !== nothing
		write_conditions(
			args["conditions"],
			args["incar"],
			args["num-atoms"],
			(args["theta-min"], args["theta-max"]),
			args["interval"],
			args["start-number"],
			args["file-prefix"],
			args["random-axis"],
		)
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
