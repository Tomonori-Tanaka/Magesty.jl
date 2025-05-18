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
function sampling_mfa(input_path::AbstractString, τ::Real, num_sample::Int64, output_index::Int64, digits::Int64)
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

	end

	args = parse_args(s)
	
	# Generate temperature range
	temps = args["temp_start"]:args["temp_step"]:args["temp_end"]

	digits = length(string(length(collect(temps)) * args["num_samples"]))
	
	for (i, τ) in enumerate(temps)
		sampling_mfa(args["input"], τ, args["num_samples"], i-1, digits)
	end
end

# Execute the main function
if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
