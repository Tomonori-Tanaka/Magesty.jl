"""
	module Systems

A module for managing crystal structures and periodic systems. It defines `Cell` and `System` types to handle lattice vectors, atomic positions, and related properties.

# Types
- **`Cell`**: Represents the unit cell of a crystal structure.
- **`System`**: Represents a periodic system built from a `Cell`, with additional information about periodicity, atom types, and neighboring images.

# Functions
- `calc_volume(lattice_vectors)`: Calculate the volume of the unit cell.
- `calc_reciprocal_vectors(lattice_vectors)`: Calculate the reciprocal lattice vectors.
- `calc_x_image_cart(lattice_vectors, x_frac)`: Compute atomic positions in neighboring cells.
- `calc_exist_image(is_periodic)`: Determine the existence of neighboring images based on periodicity.
- `calc_atomtype_group(kd_int_list)`: Group atoms by type.
- `frac2cart(lattice_vectors, xf)`: Convert fractional coordinates to Cartesian coordinates.
"""

module Systems

using LinearAlgebra
using StaticArrays

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

- `volume::Real`  
  The volume of the simulation cell, computed from the lattice vectors.

- `num_atoms::Int`  
  The total number of atoms in the cell, matching the length of `kd_int_list`.

- `num_elements::Int`  
  The number of distinct atomic species, inferred from the unique values in `kd_int_list`.

- `kd_int_list::Vector{Int}`  
  An integer list representing the atomic species or "kind" identifiers for each atom.

- `x_frac::Matrix{Float64}`  
  A matrix holding the fractional coordinates of the atoms within the cell (`[3 × num_atoms]` or `[num_atoms × 3]` depending on usage).

- `magmom::Matrix{Float64}`  
  A matrix storing the magnetic moments for each atom (`[3 × num_atoms]` or `[num_atoms × 3]` depending on usage).

# Notes
The matrix representing reciprocal vectors are defined to be the inverse matrix of the lattice vectors.
"""
struct Cell
	lattice_vectors::SMatrix{3, 3, Float64}
	reciprocal_vectors::SMatrix{3, 3, Float64}
	volume::Real
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
		calc_volume(lattice_vectors),
		length(kd_int_list),
		length(Set(kd_int_list)),
		kd_int_list,
		x_frac,
	)
end

function calc_reciprocal_vectors(lattice_vectors::AbstractMatrix{<:Real})
	return inv(lattice_vectors)
end

function calc_volume(lattice_vectors::AbstractMatrix{<:Real})
	a1 = lattice_vectors[:, 1]
	a2 = lattice_vectors[:, 2]
	a3 = lattice_vectors[:, 3]
	return abs(LinearAlgebra.dot(a1, LinearAlgebra.cross(a2, a3)))
end


"""
   System(lattice_vectors, is_periodic, kd_name, kd_int_list, x_frac, magmom)

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
	is_periodic::Vector{Bool}
	kd_name::Vector{String}# [number of elements]
	x_image_frac::Array{Float64, 3} # [≤ 3 (means x, y, z), ≤ nat, cell]
	x_image_cart::Array{Float64, 3} # the same with right above
	exist_image::Vector{Bool}
	atomtype_group::Vector{Vector{Int}}

end

function System(
	lattice_vectors::AbstractMatrix{<:Real},
	is_periodic::AbstractVector{Bool},
	kd_name::AbstractVector{<:AbstractString},
	kd_int_list::AbstractVector{<:Integer},
	x_frac::AbstractMatrix{<:Real},
)

	supercell = Cell(lattice_vectors, kd_int_list, x_frac, magmom)
	x_image_frac, x_image_cart = calc_x_images(lattice_vectors, x_frac)# [≤ 3, ≤ nat, cell]
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
	atomtype_group = Dict{Int, Vector{Int}}()
	for (iat, kd) in enumerate(kd_int_list)
		push!(get!(atomtype_group, kd, []), iat)
	end
	return collect(values(atomtype_group))
end

function calc_x_images(
	lattice_vectors::AbstractMatrix{<:Real},
	x_frac::AbstractMatrix{<:Real},
)::Tuple{Array{Float64, 3}, Array{Float64, 3}}

	num_atoms::Int = size(x_frac, 2)
	num_cell::Int = 27# the sum of centering and its neighboring imaginary cells.
	x_image_frac::Array{Float64, 3} = zeros{Float64}(3, num_atoms, num_cell)
	x_image_cart::Array{Float64, 3} = zeros{Float64}(3, num_atoms, num_cell)
	x_image_check = fill(false, 3, num_atoms, num_cell)

	x_image_frac[:, :, 1] = x_frac
	x_image_cart[:, :, 1] =
		frac2cart(lattice_vectors, x_frac)
	x_image_check[:, :, 1] .= true
	# cell means index of virtual neighboring cell.
	cell = 1
	for i in -1:1
		for j in -1:1
			for k in -1:1
				if i == j == k == 0
					# centering cell
					continue
				end
				cell += 1

				for iat in 1:num_atoms
					x_image_frac[1, iat, cell] = x_frac[1, iat] + convert(Float64, k)
					x_image_frac[2, iat, cell] = x_frac[2, iat] + convert(Float64, j)
					x_image_frac[3, iat, cell] = x_frac[3, iat] + convert(Float64, i)

					x_image_check[:, iat, cell] .= true
				end
				x_image_cart[:, :, cell] =
					frac2cart(lattice_vectors, x_image_frac[:, :, cell])
			end
		end
	end

	if false in x_image_check
		indices = findall(x -> x == false, x_image_check)
		error(
			"Error in `calc_x_image_cart`: Incomplete x_image_cart calculation detected: \n$indices",
		)
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
	#= cell means index of imaginary (virtual) neighboring cell.
	1 means centrally located supercell. =#
	cell = 1
	for i in -1:1
		for j in -1:1
			for k in -1:1
				if i == j == k == 0
					continue
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

end