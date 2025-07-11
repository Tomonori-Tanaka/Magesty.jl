#!/usr/bin/env julia

using ArgParse
using Distributions
using LinearAlgebra
using Printf
using Random
using Roots
include("vasptools.jl")
using .VaspTools

const TOLERANCE = 1e-10
const MIN_TEMP = 0.00001
const MAX_TEMP = 0.99999

"""
	thermal_averaged_m(τ::Real)

Calculate the thermal averaged m for the given temperature by solving the equation:
m = coth(3m/τ) - τ/3m.

# Arguments
- `τ`: Scaled temperature (T/Tc)

# Returns
- Thermal averaged magnetization m

# Throws
- Error if temperature is negative
"""
function thermal_averaged_m(τ::Real)
	# Special cases
	if τ < MIN_TEMP
		return 0.0
	elseif τ > MAX_TEMP
		return 1.0
	end

	# Use the bisection method to solve the equation
	function f(m::Real, τ::Real)
		return m - coth(3m/τ) + τ/3m
	end

	m_brent = find_zero(m -> f(m, τ), (MIN_TEMP, MAX_TEMP), Roots.Brent())
	return m_brent
end

"""
	mfa_spin_sampler_with_τ(input_spin_matrix::AbstractMatrix{<:Real}, τ::Real)

Generate a single sample of spin configurations using von Mises-Fisher distribution.

# Arguments
- `input_spin_matrix`: Input spin configuration matrix (3 × num_atoms)
	- Each column represents a 3D spin vector
- `τ`: Scaled temperature (T/Tc) in range [0, 1]
	- τ = 0: Fully ordered state
	- τ = 1: Fully disordered state

# Returns
- Matrix of sampled spin configurations (3 × num_atoms)
	- Each column is a normalized 3D spin vector

# Notes
- Uses von Mises-Fisher distribution for sampling
- Preserves the magnitude of each spin vector
- Skips zero vectors in the input
"""
function mfa_spin_sampler_with_τ(input_spin_matrix::AbstractMatrix{<:Real}, τ::Real)
	output_spin_matrix = zeros(size(input_spin_matrix))
	for (i, spin) in enumerate(eachcol(input_spin_matrix))
		if isapprox(norm(spin), 0.0, atol = 1e-10)
			continue
		end
		# Prepare the mean direction and the concentration in the von Mises-Fisher distribution
		mean_direction = normalize(spin)  # more efficient than spin / norm(spin)
		if τ > MAX_TEMP
			concentration = 1e-6
		elseif τ < MIN_TEMP
			return input_spin_matrix
		else
			m = thermal_averaged_m(τ)
			concentration = 3m / τ
		end

		# Calculate the von Mises-Fisher distribution
		vmf = VonMisesFisher(mean_direction, concentration)
		output_spin_direction = rand(vmf)

		# Update the spin direction
		output_spin_matrix[:, i] = output_spin_direction * norm(spin)
	end
	return output_spin_matrix
end

"""
	calculate_τ_from_magnetization(m::Real)

Calculate the scaled temperature (τ = T/Tc) from the given magnetization.

# Arguments
- `m`: Magnetization in range [0, 1]

# Returns
- Scaled temperature (τ = T/Tc)

# Throws
- Error if magnetization is not in range [0, 1]
"""
function calculate_τ_from_magnetization(m::Real)
	# Early return if m is 0 or 1
	if m ≈ 0.0
		return 1.0
	elseif m ≈ 1.0
		return 0.0
	end

	# Use the bisection method to solve the equation
	function f(τ::Real, m::Real)
		return m - coth(3m/τ) + τ/3m
	end

	τ_brent = find_zero(τ -> f(τ, m), (MIN_TEMP, MAX_TEMP), Roots.Brent())
	return τ_brent
end

"""
	mfa_spin_sampler_with_magnetization(input_spin_matrix::AbstractMatrix{<:Real}, magnetization::Real)

Generate a single sample of spin configurations using the given magnetization.

# Arguments
- `input_spin_matrix`: Input spin configuration matrix (3 × num_atoms)
	- Each column represents a 3D spin vector
"""
function mfa_spin_sampler_with_magnetization(
	input_spin_matrix::AbstractMatrix{<:Real},
	magnetization::Real,
)
	# Calculate τ from the given magnetization
	τ = calculate_τ_from_magnetization(magnetization)

	# Generate a single sample of spin configurations using the given magnetization
	output_spin_matrix = mfa_spin_sampler_with_τ(input_spin_matrix, τ)

	return output_spin_matrix
end

"""
	randomize_quantization_axis(spin_matrix::AbstractMatrix{<:Real})

Randomize the quantization axis of the spin configuration.

# Arguments
- `spin_matrix`: Input spin configuration matrix (3 × num_atoms)

# Returns
- Rotated spin configuration matrix (3 × num_atoms)
"""
function randomize_quantization_axis(spin_matrix::AbstractMatrix{<:Real})
	# Randomly choose a quantization axis
	axis = randn(3)
	axis = axis / norm(axis)  # Normalize the axis

	# Rotate the spin matrix to align with the new quantization axis
	rotation_matrix = rotation_matrix_from_vectors([0.0, 0.0, 1.0], axis)
	return rotation_matrix * spin_matrix
end

"""
	rotation_matrix_from_vectors(v1::Vector{Float64}, v2::Vector{Float64})

Calculate the rotation matrix that rotates vector v1 to align with vector v2.

# Arguments
- `v1`: Source vector (vector to be rotated)
- `v2`: Target vector (vector to align with)

# Returns
- 3×3 rotation matrix R such that R * v1 = v2

# Notes
- Both vectors are automatically normalized
- If vectors are parallel:
  - Returns identity matrix if vectors point in the same direction
  - Returns negative identity matrix if vectors point in opposite directions
- Uses Rodrigues' rotation formula for the calculation
"""
function rotation_matrix_from_vectors(v1::Vector{Float64}, v2::Vector{Float64})
	# Normalize the unit vectors
	v1_normalized = v1 / norm(v1)
	v2_normalized = v2 / norm(v2)

	# Calculate the cross product of v1 and v2
	cross_product = cross(v1_normalized, v2_normalized)
	# If the two vectors are parallel
	if norm(cross_product) ≈ 0
		if dot(v1_normalized, v2_normalized) > 0
			return Matrix{Float64}(I, 3, 3)  # If they are in the same direction, return the identity matrix
		else
			return -Matrix{Float64}(I, 3, 3)  # If they are in the opposite direction, return the identity matrix with a minus sign
		end
	end
	# Normalize the rotation axis
	v = cross_product / norm(cross_product)

	# Calculate cos(θ)
	c = dot(v1_normalized, v2_normalized)

	# Skew-symmetric matrix
	s = [          0 -v[3] v[2];
		v[3] 0 -v[1];
		-v[2] v[1] 0]

	# Calculate sin(θ) and cos(θ)
	cos_θ = c
	sin_θ = sqrt(1 - c^2)  # Since we're dealing with unit vectors

	# Rodrigues' rotation formula: R = I + sin(θ) * K + (1 - cos(θ)) * K²
	R = Matrix{Float64}(I, 3, 3) + sin_θ * s + (1 - cos_θ) * (s^2)

	return R
end


"""
	sampling_mfa(input_path::AbstractString, value::Real, num_sample::Integer, randomize::Bool, output_index::Integer, digits::Integer, variable::AbstractString)

Sample spin configurations using Mean Field Approximation (MFA) with von Mises-Fisher distribution.

# Arguments
- `input_path`: Path to the input INCAR file containing initial spin configuration
- `value`: Value of the variable to be sampled (τ or m)
	- For τ: Scaled temperature (T/Tc) in range [0, 1]
	- For m: Magnetization in range [0, 1]
- `num_sample`: Number of samples to generate for each value
- `randomize`: Whether to randomize the quantization axis
	- If true, applies random rotation to the spin configuration
- `output_index`: Index of the output file for naming
- `digits`: Number of digits in the output file name
- `variable`: Variable name to be sampled
	- "tau": Sample using scaled temperature
	- "m": Sample using magnetization

# Returns
- Nothing, but writes sampled configurations to INCAR files
	- Output files are named as "sample-{index}.INCAR"
	- Each file contains the sampled spin configuration in MAGMOM or M_CONSTR

# Throws
- `ArgumentError`: If input file is not found
- `ArgumentError`: If MAGMOM or M_CONSTR is not found in INCAR
- `ArgumentError`: If MAGMOM/M_CONSTR length is not a multiple of 3
- `ArgumentError`: If variable is not "tau" or "m"

# Examples
```julia
# Sample using temperature
sampling_mfa("input.INCAR", 0.5, 10, true, 0, 3, "tau")

# Sample using magnetization
sampling_mfa("input.INCAR", 0.8, 10, true, 0, 3, "m")
```
"""
function sampling_mfa(
	input_path::AbstractString,
	value::Real,
	num_sample::Integer,
	randomize::Bool,
	output_index::Integer,
	digits::Integer,
	variable::AbstractString,
)
	# Parse the INCAR file
	incar = parse_incar(input_path)

	if :MAGMOM in keys(incar)
		magmom = incar[:MAGMOM]
	elseif :M_CONSTR in keys(incar)
		magmom = incar[:M_CONSTR]
	else
		error("MAGMOM or M_CONSTR not found in INCAR file")
	end

	if length(magmom) % 3 != 0
		error("""
			Invalid MAGMOM/M_CONSTR length: $(length(magmom))
			The length must be a multiple of 3 for 3D spin vectors
			""")
	end

	println("input magmom: ", magmom)

	num_atoms = length(magmom) ÷ 3
	# Read the spin configuration
	spin_matrix = zeros(3, num_atoms)
	for i in 1:num_atoms
		spin_matrix[:, i] = magmom[(3i-2):3i]
	end

	# Sample the spin configurations
	for i in 1:num_sample
		if variable == "tau"
			output_spin_matrix = mfa_spin_sampler_with_τ(spin_matrix, value)
		elseif variable == "m"
			output_spin_matrix = mfa_spin_sampler_with_magnetization(spin_matrix, value)
		else
			error("""
				Invalid variable: $variable
				Expected one of: ["tau", "m"]
				""")
		end

		if randomize# Randomize the quantization axis if specified
			output_spin_matrix = randomize_quantization_axis(output_spin_matrix)
		end

		# Write the spin configuration to the output file
		output_spin_flattened = reshape(output_spin_matrix, :)  # more efficient than output_spin_matrix[:]
		output_incar = copy(incar)
		output_incar[:MAGMOM] = output_spin_flattened
		if :M_CONSTR in keys(incar)
			output_incar[:M_CONSTR] = output_spin_flattened
		end

		# Write the INCAR file with sample number
		output_path = @sprintf("sample-%0*d.INCAR", digits, output_index*num_sample + i)
		try
			write_incar(output_path, output_incar)
		catch e
			@error "Failed to write INCAR file" path=output_path exception=(e, catch_backtrace())
			rethrow(e)
		end
	end
end


"""
	Main function to parse arguments and execute sampling.
"""
function main()
	s = ArgParseSettings(description = """
	A program for sampling spin configurations using Mean Field Approximation (MFA).
	The distribution function becomes the von Mises-Fisher distribution.

	Main features:
	- Sampling spin configurations at specified temperature range (T/Tc)
	- Generate specified number of samples for each temperature
	""")

	@add_arg_table s begin
		"input"
		help = "The input file path"
		required = true

		"variable"
		help = "The variable name to be sampled (tau or m)"
		range_tester = in(["tau", "m"])
		required = true

		"--start", "-s"
		help = "The starting temperature (T/Tc)"
		required = true
		arg_type = Float64

		"--end", "-e"
		help = "The ending temperature (T/Tc)"
		required = true
		arg_type = Float64

		"--step", "-w"
		help = "The temperature step size"
		required = true
		arg_type = Float64

		"--num_samples", "-n"
		help = "The number of samples in each step"
		arg_type = Int64
		default = 1

		"--randomize", "-r"
		help = "Randomize the quantization axis"
		action = "store_true"

	end

	args = parse_args(s)

	# Generate temperature range
	sampling_list = args["start"]:args["step"]:args["end"]

	digits = length(string(length(collect(sampling_list)) * args["num_samples"]))

	# check the temperature range
	if args["start"] > args["end"]
		error("start must be less than or equal to end")
	end

	if args["step"] <= 0
		error("step must be positive")
	end

	for (i, value) in enumerate(sampling_list)
		sampling_mfa(
			args["input"],
			value,
			args["num_samples"],
			args["randomize"],
			i-1,
			digits,
			args["variable"],
		)
	end

	print_info(args)
end

function print_info(args)
	@printf("Input file: %s\n", args["input"])
	@printf("Variable: %s\n", args["variable"])
	@printf(
		"Sampling list: %.2f to %.2f with step %.2f\n",
		args["start"],
		args["end"],
		args["step"]
	)
	@printf("Number of samples per step: %d\n", args["num_samples"])
	@printf("Randomize quantization axis: %s\n", args["randomize"] ? "Yes" : "No")
end

# Execute the main function
if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
