"""
	module Structures

A module for managing crystal structures and periodic systems.
"""
module Structures
using LinearAlgebra
using Printf
using StaticArrays

using ..ConfigParser

import Base: show

# Constants
const NUM_CELLS = 27  # Total number of cells: center cell and its neighboring virtual cells
const DIMENSIONS = 3  # Number of spatial dimensions

export Structure

"""
	Cell

Represents the unit cell of a crystal structure with lattice vectors, atomic positions, and related properties.

# Fields

- `lattice_vectors::SMatrix{3,3,Float64}`  
  The 3×3 matrix of lattice vectors in Cartesian coordinates.
[a1_x a2_x a3_x
 a1_y a2_y a3_y
 a1_z a2_z a3_z]

- `reciprocal_vectors::SMatrix{3,3,Float64}`  
  The 3×3 matrix of reciprocal lattice vectors, derived from `lattice_vectors`.
[b1_x b1_y b1_z
 b2_x b2_y b2_z
 b3_x b3_y b3_z]

- `num_atoms::Int`  
  The total number of atoms in the cell, matching the length of `kd_int_list`.

- `num_elements::Int`  
  The number of distinct atomic species, inferred from the unique values in `kd_int_list`.

- `kd_int_list::Vector{Int}`  
  An integer list representing the atomic species or "kind" identifiers for each atom.

- `x_frac::Matrix{Float64}`  
  A matrix holding the fractional coordinates of the atoms within the cell (`[3 × num_atoms]` or `[num_atoms × 3]` depending on usage).


# Notes
The matrix representing reciprocal vectors are defined to be the inverse matrix of the lattice vectors.
"""
struct Cell
	lattice_vectors::SMatrix{3, 3, Float64}
	reciprocal_vectors::SMatrix{3, 3, Float64}
	num_atoms::Int
	num_elements::Int
	kd_int_list::Vector{Int}
	x_frac::Matrix{Float64}

end

function Cell(
	lattice_vectors::AbstractMatrix{<:Real},
	kd_int_list::AbstractVector{<:Integer},
	x_frac::AbstractMatrix{<:Real},
)
	validate_lattice_vectors(lattice_vectors)

	return Cell(lattice_vectors,
		calc_reciprocal_vectors(lattice_vectors),
		length(kd_int_list),
		length(Set(kd_int_list)),
		kd_int_list,
		x_frac,
	)
end

function show(io::IO, cell::Cell)
	println(io, "\tlattice_vectors: ", cell.lattice_vectors)
	println(io, "\treciprocal_vectors: ", cell.reciprocal_vectors)
	println(io, "\tnum_atoms: ", cell.num_atoms)
	println(io, "\tnum_elements: ", cell.num_elements)
	println(io, "\tkd_ind_list: ", cell.kd_int_list)
	println(io, "\tx_frac: ", cell.x_frac)
end

function validate_lattice_vectors(lattice_vectors::AbstractMatrix{<:Real})
	# Check linear independence and right-handed coordinate system
	det_value = det(lattice_vectors)
	if det_value ≈ 0
		error("Lattice vectors are linearly dependent. det(lattice_vectors) = $det_value")
	elseif det_value < 0
		error("Lattice vectors do not form a right-handed coordinate system. det(lattice_vectors) = $det_value")
	end
end

function calc_reciprocal_vectors(lattice_vectors::AbstractMatrix{<:Real})
	return inv(lattice_vectors)
end

"""
	Structure

Represents a periodic structure built from a Cell, with information about periodicity, atom types, and neighboring images.

# Fields

- `supercell::Cell`  
  The unit cell of the structure, containing lattice vectors, reciprocal vectors, atomic positions, and magnetic moments.

- `is_periodic::Vector{Bool}`  
  A vector indicating periodicity along each of the three principal axes. Each element corresponds to whether the structure is periodic (`true`) or non-periodic (`false`) along that axis.

- `kd_name::Vector{String}`  
  A list of element names present in the structure (e.g., `["Fe", "Co", "Ni"]`).

- `x_image_frac::Array{Float64, 3}`  
  Fractional coordinates of atoms in neighboring (imaginary) cells. The dimensions typically represent the number of images, number of atoms, and the three fractional coordinates.

- `x_image_cart::Array{Float64, 3}`  
  Cartesian coordinates of atoms in neighboring (imaginary) cells. Similar to `x_image_frac`, the dimensions represent the number of images, number of atoms, and the three Cartesian coordinates.

- `exist_image::Vector{Bool}`  
  Indicates the existence of neighboring cells based on the structure's periodicity. Each element corresponds to whether a particular neighboring cell exists (`true`) or not (`false`).

- `atomtype_group::Vector{Vector{Int}}`  
  Groups of atom indices categorized by their types. Each sub-vector contains the indices of atoms belonging to a specific element type.

- `elapsed_time::Float64`  
  Time taken to create the structure in seconds.

"""
struct Structure
	supercell::Cell
	is_periodic::SVector{3, Bool}
	kd_name::Vector{String}
	x_image_frac::Array{Float64, 3}
	x_image_cart::Array{Float64, 3}
	exist_image::Vector{Bool}
	atomtype_group::Vector{Vector{Int}}
	elapsed_time::Float64

	function Structure(
		lattice_vectors::AbstractMatrix{<:Real},
		is_periodic::AbstractVector{Bool},
		kd_name::AbstractVector{<:AbstractString},
		kd_int_list::AbstractVector{<:Integer},
		x_frac::AbstractMatrix{<:Real},
	)
		start_time::UInt64 = time_ns()
		supercell::Cell = Cell(lattice_vectors, kd_int_list, x_frac)
		x_image_frac::Array{Float64, 3} = zeros(Float64, DIMENSIONS, size(x_frac, 2), NUM_CELLS)
		x_image_cart::Array{Float64, 3} = zeros(Float64, DIMENSIONS, size(x_frac, 2), NUM_CELLS)
		x_image_frac, x_image_cart = calc_x_images(lattice_vectors, x_frac)
		exist_image::Vector{Bool} = calc_exist_image(is_periodic)
		atomtype_group::Vector{Vector{Int}} = calc_atomtype_group(kd_int_list)
		elapsed_time::Float64 = (time_ns() - start_time) / 1e9

		return new(
			supercell,
			is_periodic,
			kd_name,
			x_image_frac,
			x_image_cart,
			exist_image,
			atomtype_group,
			elapsed_time,
		)
	end
end

"""
	Structure(config::Config4System) -> Structure

Create a Structure from a Config4System object.
"""
function Structure(config::Config4System)::Structure
	lattice_vectors::SMatrix{3, 3, Float64} = config.lattice_vectors
	is_periodic::SVector{3, Bool} = config.is_periodic
	kd_name::Vector{String} = config.kd_name
	kd_int_list::Vector{Int} = config.kd_int_list
	x_frac::Matrix{Float64} = config.x_fractional

	return Structure(
		lattice_vectors,
		is_periodic,
		kd_name,
		kd_int_list,
		x_frac,
	)
end

"""
	calc_atomtype_group(kd_int_list) -> Vector{Vector{Int}}

Group atom indices by their types.
"""
function calc_atomtype_group(kd_int_list::AbstractVector{<:Integer})::Vector{Vector{Int}}
	unique_vals::Vector{Int} = unique(kd_int_list)
	return [findall(x -> x == val, kd_int_list) for val in unique_vals]
end

"""
	calc_x_images(lattice_vectors, x_frac) -> Tuple{Array{Float64, 3}, Array{Float64, 3}}

Calculate fractional and Cartesian coordinates of atoms in neighboring cells.
"""
function calc_x_images(
	lattice_vectors::AbstractMatrix{<:Real},
	x_frac::AbstractMatrix{<:Real},
)::Tuple{Array{Float64, 3}, Array{Float64, 3}}
	num_atoms::Int = size(x_frac, 2)
	x_image_frac::Array{Float64, 3} = zeros(Float64, DIMENSIONS, num_atoms, NUM_CELLS)
	x_image_cart::Array{Float64, 3} = zeros(Float64, DIMENSIONS, num_atoms, NUM_CELLS)
	x_image_check::Array{Bool, 3} = fill(false, DIMENSIONS, num_atoms, NUM_CELLS)

	# Set up the center cell
	x_image_frac[:, :, 1] = x_frac
	x_image_cart[:, :, 1] = frac2cart(lattice_vectors, x_frac)
	x_image_check[:, :, 1] .= true

	# Pre-allocate static vectors for better performance
	x_image_tmp = MVector{3, Float64}(undef)
	
	# Calculate virtual cells
	cell::Int = 1
	for k in -1:1
		for j in -1:1
			for i in -1:1
				if i == j == k == 0
					continue  # Skip the center cell
				end
				cell += 1

				for iat in 1:num_atoms
					# Use static vector for temporary calculations
					x_image_tmp[1] = x_frac[1, iat] + convert(Float64, i)
					x_image_tmp[2] = x_frac[2, iat] + convert(Float64, j)
					x_image_tmp[3] = x_frac[3, iat] + convert(Float64, k)
					
					x_image_frac[:, iat, cell] = x_image_tmp
					x_image_check[:, iat, cell] .= true
				end
				x_image_cart[:, :, cell] = frac2cart(lattice_vectors, @view(x_image_frac[:, :, cell]))
			end
		end
	end

	# Check completeness of calculations
	if false in x_image_check
		indices::Vector{CartesianIndex{3}} = findall(x -> x == false, x_image_check)
		error("""
			Error in `calc_x_images`: Incomplete calculation detected.
			Missing coordinates at indices: $indices
			Please check the calculation of neighboring cell coordinates.
		""")
	end

	return x_image_frac, x_image_cart
end

"""
	frac2cart(lattice_vectors, x_frac) -> Matrix{Float64}

Convert fractional coordinates to Cartesian coordinates.
"""
function frac2cart(
	lattice_vectors::AbstractMatrix{<:Real},
	x_frac::AbstractMatrix{<:Real},
)::Matrix{Float64}
	# Convert to SMatrix for better performance with matrix multiplication
	static_lattice = SMatrix{3,3,Float64}(lattice_vectors)
	return Matrix(static_lattice * x_frac)
end

"""
	calc_exist_image(is_periodic) -> Vector{Bool}

Determine which neighboring cells exist based on periodicity.
"""
function calc_exist_image(is_periodic::SVector{3, Bool})::Vector{Bool}
	exist_image::Vector{Bool} = fill(true, NUM_CELLS)
	# Cell index 1 represents the central supercell
	# Other indices (2-27) represent neighboring virtual cells
	cell::Int = 1
	
	# Use static vector for offset calculations
	offset = MVector{3, Int}(undef)
	
	for i in -1:1
		offset[1] = i
		for j in -1:1
			offset[2] = j
			for k in -1:1
				offset[3] = k
				if i == j == k == 0
					continue  # Skip the center cell
				end
				cell += 1

				if (
					(abs(offset[1]) == 1 && !is_periodic[1])
					|| (abs(offset[2]) == 1 && !is_periodic[2])
					|| (abs(offset[3]) == 1 && !is_periodic[3])
				)
					exist_image[cell] = false
				end
			end
		end
	end
	return exist_image
end

function show(io::IO, structure::Structure)
	println(io, "supercell:\n ", structure.supercell)
	println(io, "is_periodic: ", structure.is_periodic)
	println(io, "kd_name: ", structure.kd_name)
	println(io, "x_image_frac: \n", structure.x_image_frac)
	println(io, "x_image_cart: \n", structure.x_image_cart)
	println(io, "exist_image: ", structure.exist_image)
	println(io, "atomtype_group: ", structure.atomtype_group)
end

function print_info(structure::Structure)
	supercell = structure.supercell
	println("""
	======
	SYSTEM
	======
	""")
	println("Total Number of atoms: ", supercell.num_atoms)
	println("Number of atomic species: ", supercell.num_elements)

	println("Lattice vector (in Angstrom):")
	for (i, label) in enumerate(["a1", "a2", "a3"])
		println(
			@sprintf("  %12.8e   %12.8e   %12.8e : %s",
				supercell.lattice_vectors[:, i]..., label)
		)
	end
	println("")
	println("Periodicity: ", Int.(structure.is_periodic), "\n")

	digits = length(string(abs(supercell.num_atoms))) + 1
	println("Atomic species:")
	for (i, species) in enumerate(structure.kd_name)
		println(@sprintf("  %*d: %s", digits, i, species))
	end
	println("")

	println("Atomic positions in fractional coordinates and atomic species:")
	for i in 1:supercell.num_atoms
		println(
			@sprintf("  %*d: %12.8e   %12.8e   %12.8e   %*d",
				digits, i,
				supercell.x_frac[1, i],
				supercell.x_frac[2, i],
				supercell.x_frac[3, i],
				digits, supercell.kd_int_list[i])
		)
	end
	println("")
	println(@sprintf("Elapsed time: %.6f seconds", structure.elapsed_time))
	println("-------------------------------------------------------------------")
end
end
