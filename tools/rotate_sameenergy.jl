"""
	rotate_sameenergy.jl

This script rotates the spin directions and magnetic fields in EMBSET with keeping the same energy by using the Fibonacci sphere.

This is intended to generate multiple data sets with different spin orientations at the same energy level,
specifically for systems where spin-orbit coupling is not considered.
"""
module RotateSameEnergy

using LinearAlgebra
using Statistics
using Printf
include("../src/SpinConfigs.jl")
using .SpinConfigs

export rotate_spin_and_field
"""
	fibonacci_sphere(n_points::Int) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Generate Fibonacci lattice points on the sphere.

Arguments:
	n_points: Number of points to generate

Returns:
	(x, y, z): 3D coordinates of each point on the unit sphere
"""
function fibonacci_sphere(n_points::Int)
	# Golden ratio
	φ = (1 + √5) / 2

	# Array to store coordinates
	x = Vector{Float64}(undef, n_points)
	y = Vector{Float64}(undef, n_points)
	z = Vector{Float64}(undef, n_points)

	for i in 0:n_points-1
		z[i+1] = 1 - (2i + 1) / n_points
		r = √(1 - z[i+1]^2)
		θ = 2π * i / φ
		x[i+1] = r * cos(θ)
		y[i+1] = r * sin(θ)
	end

	return x, y, z
end

"""
	rotation_matrix_from_vectors(v1::Vector{Float64}, v2::Vector{Float64}) -> Matrix{Float64}

Calculate the rotation matrix from two unit vectors.

Returns the 3x3 rotation matrix that represents the rotation from v1 to v2.

Arguments:
	v1: The starting unit vector
	v2: The ending unit vector

Returns:
	3x3 rotation matrix
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
	s = [        0 -v[3] v[2];
		v[3] 0 -v[1];
		-v[2] v[1] 0]

	# Rodrigues' rotation formula
	R = Matrix{Float64}(I, 3, 3) + s + s^2 * (1 - c) / (norm(v)^2)

	return R
end

"""
	rotate_spin_and_field(embset::AbstractString, saxis::Vector{Float64}, n_points::Int)

Rotate the spin directions and magnetic fields in EMBSET with keeping the same energy by using the Fibonacci sphere.

Arguments:
	embset: The path to the EMBSET file
	num_atoms: The number of atoms in the structure
	saxis: The reference axis of the spin quantization
	n_points: The number of points to generate on the Fibonacci sphere
"""
function rotate_spin_and_field(
	embset::AbstractString,
	num_atoms::Int,
	saxis::Vector{Float64},
	n_points::Int,
	check::Bool,
)
	# Read the EMBSET file
	spin_configs::Vector{SpinConfig} = read_embset(embset, num_atoms)

	# Generate the Fibonacci sphere points
	x_list::Vector{Float64}, y_list::Vector{Float64}, z_list::Vector{Float64} =
		fibonacci_sphere(n_points)

	if check
		println("the number of the Fibonacci sphere points: $(n_points)")
		println("mean values of the Fibonacci sphere points")
		println("x: $(mean(x_list)), y: $(mean(y_list)), z: $(mean(z_list))")
		return
	end

	# Rotate the spin directions and magnetic fields
	for i in 1:n_points
		x = x_list[i]
		y = y_list[i]
		z = z_list[i]

		R::Matrix{Float64} = rotation_matrix_from_vectors(saxis, [x, y, z])

		# Rotate the spin and magnetic field directions
		for (num_config, config) in enumerate(spin_configs)
			energy = config.energy
			spin_direction = R * config.spin_directions # [3, num_atoms]
			field_direction = R * config.local_magfield # [3, num_atoms]

			concated_matrix = hcat(spin_direction', field_direction')   # [num_atoms, 6]

			println(
				"# the number of original config:$num_config, the number of the Fibonacci sphere point:$i, energy unit = eV, magmom unit = Bohr magneton, magnetic field unit = T",
			)
			@printf("%.5f\n", energy)
			for (row_idx, row) in enumerate(eachrow(concated_matrix))
				println(
					lpad(row_idx, 7),
					"  ",
					join([lpad(@sprintf("%.7f", x), 12) for x in row[1:3]], " "),
					"  ",
					join([lpad(@sprintf("%.5e", x), 14) for x in row[4:6]], " "),
				)
			end
		end
	end
end
end

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
	"embset"
	help = "The path to the EMBSET file"
	required = true

	"--num_atoms", "-n"
	help = "The number of atoms in the structure"
	required = true
	arg_type = Int

	"--saxis"
	help = "The reference axis of the spin quantization"
	nargs = 3
	arg_type = Float64
	default = [0.0, 0.0, 1.0]

	"--n_points"
	help = "The number of points to generate on the Fibonacci sphere"
	default = 1000
	arg_type = Int

	"--check", "-c"
	help = "Switch check mode"
	action = :store_true
end

parsed_args = parse_args(s)

RotateSameEnergy.rotate_spin_and_field(
	parsed_args["embset"],
	parsed_args["num_atoms"],
	parsed_args["saxis"],
	parsed_args["n_points"],
	parsed_args["check"],
)

