"""
	module Systems

A module for managing crystal structures and periodic systems. It defines `Cell` and `System` types to handle lattice vectors, atomic positions, and related properties.

# Types
- **`Cell`**: Represents the unit cell of a crystal structure.
- **`System`**: Represents a periodic system built from a `Cell`, with additional information about periodicity, atom types, and neighboring images.

"""

module Systems
using LinearAlgebra
using Printf
using StaticArrays

import Base: show

export System

"""
	Cell

A structure that stores crystallographic information for a simulation cell.

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

function Cell(lattice_vectors::AbstractMatrix{<:Real},
	kd_int_list::AbstractVector{<:Integer},
	x_frac::AbstractMatrix{<:Real},
)
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

function calc_reciprocal_vectors(lattice_vectors::AbstractMatrix{<:Real})
	return inv(lattice_vectors)
end

"""
   System(lattice_vectors, is_periodic, kd_name, kd_int_list, x_frac)

A structure that represents a periodic system constructed from a `Cell`.

# Fields

- `supercell::Cell`  
  The unit cell of the system, containing lattice vectors, reciprocal vectors, atomic positions, and magnetic moments.

- `is_periodic::Vector{Bool}`  
  A vector indicating periodicity along each of the three principal axes. Each element corresponds to whether the system is periodic (`true`) or non-periodic (`false`) along that axis.

- `kd_name::Vector{String}`  
  A list of element names present in the system (e.g., `["Fe", "Co", "Ni"]`).

- `x_image_frac::Array{Float64, 3}`  
  Fractional coordinates of atoms in neighboring (imaginary) cells. The dimensions typically represent the number of images, number of atoms, and the three fractional coordinates.

- `x_image_cart::Array{Float64, 3}`  
  Cartesian coordinates of atoms in neighboring (imaginary) cells. Similar to `x_image_frac`, the dimensions represent the number of images, number of atoms, and the three Cartesian coordinates.

- `exist_image::Vector{Bool}`  
  Indicates the existence of neighboring cells based on the system's periodicity. Each element corresponds to whether a particular neighboring cell exists (`true`) or not (`false`).

- `atomtype_group::Vector{Vector{Int}}`  
  Groups of atom indices categorized by their types. Each sub-vector contains the indices of atoms belonging to a specific element type.

"""
struct System
	supercell::Cell
	is_periodic::Vector{Bool}  # Periodicity flags for x, y, z directions
	kd_name::Vector{String}    # Element names (e.g., ["Fe", "Co", "Ni"])
	x_image_frac::Array{Float64, 3}  # Fractional coordinates of atoms in neighboring cells [3, num_atoms, 27]
	x_image_cart::Array{Float64, 3}  # Cartesian coordinates of atoms in neighboring cells [3, num_atoms, 27]
	exist_image::Vector{Bool}  # Flags indicating existence of neighboring cells based on periodicity
	atomtype_group::Vector{Vector{Int}}  # Groups of atom indices by element type

end

function System(
	lattice_vectors::AbstractMatrix{<:Real},
	is_periodic::AbstractVector{Bool},
	kd_name::AbstractVector{<:AbstractString},
	kd_int_list::AbstractVector{<:Integer},
	x_frac::AbstractMatrix{<:Real},
)
	supercell = Cell(lattice_vectors, kd_int_list, x_frac)
	x_image_frac, x_image_cart = calc_x_images(lattice_vectors, x_frac)
	exist_image::Vector{Bool} = calc_exist_image(is_periodic)
	atomtype_group::Vector{Vector{Int}} = calc_atomtype_group(kd_int_list)
	return System(
		supercell,
		is_periodic,
		kd_name,
		x_image_frac,
		x_image_cart,
		exist_image,
		atomtype_group,
	)
end

function calc_atomtype_group(kd_int_list::AbstractVector{<:Integer})::Vector{Vector{Int}}
	unique_vals = unique(kd_int_list)
	return [findall(x -> x == val, kd_int_list) for val in unique_vals]
end

function calc_x_images(
	lattice_vectors::AbstractMatrix{<:Real},
	x_frac::AbstractMatrix{<:Real},
)::Tuple{Array{Float64, 3}, Array{Float64, 3}}
	num_atoms::Int = size(x_frac, 2)
	num_cell::Int = 27  # Total number of cells: center cell and its neighboring virtual cells
	x_image_frac::Array{Float64, 3} = zeros(Float64, 3, num_atoms, num_cell)
	x_image_cart::Array{Float64, 3} = zeros(Float64, 3, num_atoms, num_cell)
	x_image_check::Array{Bool, 3} = fill(false, 3, num_atoms, num_cell)

	# Set up the center cell
	x_image_frac[:, :, 1] = x_frac
	x_image_cart[:, :, 1] = frac2cart(lattice_vectors, x_frac)
	x_image_check[:, :, 1] .= true

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
					x_image_frac[1, iat, cell] = x_frac[1, iat] + convert(Float64, i)
					x_image_frac[2, iat, cell] = x_frac[2, iat] + convert(Float64, j)
					x_image_frac[3, iat, cell] = x_frac[3, iat] + convert(Float64, k)

					x_image_check[:, iat, cell] .= true
				end
				x_image_cart[:, :, cell] = frac2cart(lattice_vectors, x_image_frac[:, :, cell])
			end
		end
	end

	# Check completeness of calculations
	if false in x_image_check
		indices = findall(x -> x == false, x_image_check)
		error("""
			Error in `calc_x_images`: Incomplete calculation detected.
			Missing coordinates at indices: $indices
			Please check the calculation of neighboring cell coordinates.
		""")
	end

	return x_image_frac, x_image_cart
end

function frac2cart(
	lattice_vectors::AbstractMatrix{<:Real},
	x_frac::AbstractMatrix{<:Real})::Matrix{Float64}

	return lattice_vectors * x_frac
end

function calc_exist_image(is_periodic::AbstractVector{Bool})::Vector{Bool}
	num_cell = 27
	exist_image = fill(true, num_cell)
	# Cell index 1 represents the central supercell
	# Other indices (2-27) represent neighboring virtual cells
	cell = 1
	for i in -1:1
		for j in -1:1
			for k in -1:1
				if i == j == k == 0
					continue  # Skip the center cell
				end
				cell += 1

				if (
					(abs(i) == 1 && is_periodic[1] == false)
					|| (abs(j) == 1 && is_periodic[2] == false)
					|| (abs(k) == 1 && is_periodic[3] == false)
				)
					exist_image[cell] = false
				end
			end
		end
	end
	return exist_image
end

function show(io::IO, system::System)
	println(io, "supercell:\n ", system.supercell)
	println(io, "is_periodic: ", system.is_periodic)
	println(io, "kd_name: ", system.kd_name)
	println(io, "x_image_frac: \n", system.x_image_frac)
	println(io, "x_image_cart: \n", system.x_image_cart)
	println(io, "exist_image: ", system.exist_image)
	println(io, "atomtype_group: ", system.atomtype_group)
end

function print_info(system::System)
	supercell = system.supercell
	println("""
	======
	SYSTEM
	======
	""")
	println("Total Number of atoms: ", supercell.num_atoms)
	println("Number of atomic species: \n", supercell.num_elements)

	println("Lattice vector (in Angstrom):")
	for (i, label) in enumerate(["a1", "a2", "a3"])
		println(@sprintf("  %12.8e   %12.8e   %12.8e : %s", 
			supercell.lattice_vectors[:, i]..., label))
	end
	println("")
	println("Periodicity: ", Int.(system.is_periodic), "\n")
	
	digits = length(string(abs(supercell.num_atoms))) + 1
	println("Atomic species:")
	for (i, species) in enumerate(system.kd_name)
		println(@sprintf("  %*d: %s", digits, i, species))
	end
	println("")
	
	println("Atomic positions in fractional coordinates and atomic species:")
	for i in 1:supercell.num_atoms
		println(@sprintf("  %*d: %12.8e   %12.8e   %12.8e   %*d",
			digits, i,
			supercell.x_frac[1, i],
			supercell.x_frac[2, i],
			supercell.x_frac[3, i],
			digits, supercell.kd_int_list[i]))
	end
	println("")
	println("-------------------------------------------------------------------")
end
end
