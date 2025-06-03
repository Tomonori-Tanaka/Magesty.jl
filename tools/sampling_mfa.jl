#!/usr/bin/env julia

using ArgParse
using Distributions
using LinearAlgebra
using Printf
using Random
using Roots
include("vasptools.jl")
using .VaspTools

"""
	sampling_mfa(input_path::AbstractString, τ::Real, num_sample::Int64)

Sample spin configurations using Mean Field Approximation.

# Arguments
- `input_path`: Path to the input INCAR file
- `τ`: Scaled temperature (T/Tc)
- `num_sample`: Number of samples to generate

# Returns
- Nothing, but writes sampled configurations to INCAR files
"""
function sampling_mfa(input_path::AbstractString, τ::Real, num_sample::Int64, randomize::Bool, output_index::Int64, digits::Int64)
	# Parse the INCAR file
	incar = parse_incar(input_path)
	m_constr = get(incar, :M_CONSTR, nothing)
	if isnothing(m_constr)
		error("M_CONSTR not found in INCAR file")
	end
	
	if length(m_constr) % 3 != 0
		error("The length of M_CONSTR must be a multiple of 3")
	end
	
	num_atoms = length(m_constr) ÷ 3
	# Read the spin configuration
	spin_matrix = zeros(3, num_atoms)
	for i in 1:num_atoms
		spin_matrix[:, i] = m_constr[(3i-2):3i]
	end
	
	# Sample the spin configurations
	for i in 1:num_sample
		output_spin_matrix = sampling_mfa_for_each_sample(spin_matrix, τ)
		if randomize	# Randomize the quantization axis if specified
			output_spin_matrix = randomize_quantization_axis(output_spin_matrix)
		end

		# Write the spin configuration to the output file
		output_spin_flattened = vec(output_spin_matrix)  # more efficient than output_spin_matrix[:]
		output_incar = copy(incar)
		output_incar[:M_CONSTR] = output_spin_flattened
		output_incar[:MAGMOM] = output_spin_flattened

		# Write the INCAR file with sample number
		output_path = @sprintf("sample-%0*d.INCAR", digits, output_index*num_sample + i)
		write_incar(output_path, output_incar)
	end
end

"""
	sampling_mfa_for_each_sample(input_spin_matrix::AbstractMatrix{<:Real}, τ::Real)

Generate a single sample of spin configurations using von Mises-Fisher distribution.

# Arguments
- `input_spin_matrix`: Input spin configuration matrix (3 × num_atoms)
- `τ`: Scaled temperature (T/Tc)

# Returns
- Matrix of sampled spin configurations
"""
function sampling_mfa_for_each_sample(input_spin_matrix::AbstractMatrix{<:Real}, τ::Real)
	output_spin_matrix = zeros(size(input_spin_matrix))
	for (i, spin) in enumerate(eachcol(input_spin_matrix))
		if isapprox(norm(spin), 0.0, atol=1e-10)
			continue
		end
		# Prepare the mean direction and the concentration in the von Mises-Fisher distribution
		mean_direction = normalize(spin)  # more efficient than spin / norm(spin)
        if τ > 0.99999
            concentration = 0.000001
        elseif τ < 0.00001
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

	# Use the bisection method to solve the equation
	m_min = 0.00001
	m_max = 0.99999

	function f(m::Real, τ::Real)
		return m - coth(3m/τ) + τ/3m
	end

	m_brent = find_zero(m -> f(m, τ), (m_min, m_max), Roots.Brent())
	return m_brent
end

function randomize_quantization_axis(spin_matrix::AbstractMatrix{<:Real})
	# Randomly choose a quantization axis
	axis = randn(3)
	axis = axis / norm(axis)  # Normalize the axis

	# Rotate the spin matrix to align with the new quantization axis
	rotation_matrix = rotation_matrix_from_vectors([0.0, 0.0, 1.0], axis)
	return rotation_matrix * spin_matrix
end

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
	s = [        0 -v[3] v[2];
		v[3] 0 -v[1];
		-v[2] v[1] 0]

	# Rodrigues' rotation formula
	R = Matrix{Float64}(I, 3, 3) + s + s^2 * (1 - c) / (norm(v)^2)

	return R
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

		"--temp_start", "-s"
		help = "The starting temperature (T/Tc)"
		required = true
		arg_type = Float64

		"--temp_end", "-e"
		help = "The ending temperature (T/Tc)"
		required = true
		arg_type = Float64

		"--temp_step", "-w"
		help = "The temperature step size"
		required = true
		arg_type = Float64

		"--num_samples", "-n"
		help = "The number of samples for each temperature"
		arg_type = Int64
        default = 1

		"--randomize", "-r"
		help = "Randomize the quantization axis"
		action = "store_true"	

	end

	args = parse_args(s)
	
	# Generate temperature range
	temps = args["temp_start"]:args["temp_step"]:args["temp_end"]

	digits = length(string(length(collect(temps)) * args["num_samples"]))
	
	for (i, τ) in enumerate(temps)
		sampling_mfa(args["input"], τ, args["num_samples"], args["randomize"], i-1, digits)
	end

	print_info(args)
end

function print_info(args)
	@printf("Input file: %s\n", args["input"])
	@printf("Temperature range: %.2f to %.2f with step %.2f\n", args["temp_start"], args["temp_end"], args["temp_step"])
	@printf("Number of samples per temperature: %d\n", args["num_samples"])
	@printf("Randomize quantization axis: %s\n", args["randomize"] ? "Yes" : "No")
end

# Execute the main function
if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
