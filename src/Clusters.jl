"""
	module Clusters

This module provides data structures and functions for managing and analyzing clusters of interacting atoms within a crystal structure.

# Types
- `DistInfo`: Represents distance information between two atoms
- `Cluster`: Manages multiple interaction clusters based on structure parameters
- `DistList`: Internal structure for sorting distance information

# Functions
- `Cluster(structure, symmetry, nbody, cutoff_radii)`: Constructs a `Cluster` instance
- `set_mindist_pairs(num_atoms, cartesian_coords, cell_exists)`: Computes minimum distance pairs
- `is_within_cutoff(atomcell_list, kd_int_list, cutoff_radii, body, x_image_cart, min_distance_pairs)`: Checks cutoff conditions
- `distance_atomcells(atomcell1, atomcell2, x_image_cart)`: Calculates distance between two atoms in different cells
- `generate_clusters(structure, symmetry, cutoff_radii, nbody)`: Generates interaction clusters
- `print_cluster_stdout(min_distance_pairs, atoms_in_prim, kd_name, kd_int_list)`: Prints cluster information to stdout
"""

module Clusters

using Combinat
using DataStructures
using LinearAlgebra
using OffsetArrays
using Printf

using ..CountingContainer
using ..SortedContainer
using ..AtomCells
using ..ConfigParser
using ..Structures
using ..Symmetries

import Base: isless, ==

export Cluster

# Constants
const NUM_VIRTUAL_CELLS = 27  # Number of virtual cells in 3x3x3 supercell
const DEFAULT_TOLERANCE = 1e-5  # Default tolerance for floating-point comparisons

"""
	struct DistInfo

Represents distance information between two atoms.

# Fields
- `cell_index::Int`: Index of the virtual cell containing the second atom (1-27)
- `distance::Float64`: Cartesian distance between atoms in Angstroms
- `relative_vector::Vector{Float64}`: Relative Cartesian vector from first to second atom

# Constructor
	DistInfo(cell::Integer, distance::Real, relvec::AbstractVector{<:Real})

Creates a new `DistInfo` instance. Ensures that `relvec` has length 3.

# Examples
```julia
# Create a DistInfo instance for atoms in the same cell
dist = DistInfo(1, 2.5, [1.0, 0.0, 0.0])

# Create a DistInfo instance for atoms in different cells
dist = DistInfo(2, 3.0, [2.0, 1.0, 0.0])
```
"""
struct DistInfo
	cell_index::Int
	distance::Float64
	relative_vector::Vector{Float64}

	function DistInfo(cell::Integer, distance::Real, relvec::AbstractVector{<:Real})
		@assert length(relvec) == 3 "The length of \"relvec\" must be 3."
		return new(Int(cell), Float64(distance), Vector{Float64}(relvec))
	end
end

function Base.isless(distinfo1::DistInfo, distinfo2::DistInfo)::Bool
	@inbounds begin
		if distinfo1.distance < distinfo2.distance
			return true
		elseif distinfo1.distance ≈ distinfo2.distance
			return distinfo1.cell_index < distinfo2.cell_index
		end
		return false
	end
end


"""
	struct Cluster

Represents a collection of interaction clusters based on the specified number of bodies and cutoff radii.

# Fields
- `num_bodies::Int`: Number of interacting bodies
- `cutoff_radii::OffsetArray{Float64, 3}`: Cutoff radii for each atomic element pair and interaction body
- `min_distance_pairs::Matrix{Vector{DistInfo}}`: Matrix of minimum distance pairs between atoms
- `cluster_dict::Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}}`: Dictionary of interaction clusters organized by body and primitive atom index

# Constructor
	Cluster(structure, symmetry, nbody, cutoff_radii)

Creates a new `Cluster` instance based on the provided structure, symmetry information, number of bodies, and cutoff radii.

# Example
```julia
cluster = Cluster(structure, symmetry, 3, cutoff_radii)
```
"""
struct Cluster
	num_bodies::Int
	cutoff_radii::OffsetArray{Float64, 3}
	min_distance_pairs::Matrix{Vector{DistInfo}}
	cluster_dict::Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}}

	function Cluster(
		structure::Structure,
		symmetry::Symmetry,
		nbody::Integer,
		cutoff_radii::AbstractArray{<:Real, 3},
		;
		verbosity::Bool = true,
	)

		start_time = time_ns()

		cluster_dict::Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}} =
			generate_clusters(structure, symmetry, cutoff_radii, nbody)

		min_distance_pairs = set_mindist_pairs(
			structure.supercell.num_atoms,
			structure.x_image_cart,
			structure.exist_image,
			tol = symmetry.tol,
		)


		if verbosity
			print_cluster_stdout(
				min_distance_pairs,
				symmetry.atoms_in_prim,
				structure.kd_name,
				structure.supercell.kd_int_list,
			)
			elapsed_time = (time_ns() - start_time) / 1e9
			println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
			println("-------------------------------------------------------------------")
		end


		return new(
			nbody,
			cutoff_radii,
			min_distance_pairs,
			cluster_dict,
		)
	end
end

function Cluster(
	structure::Structure,
	symmetry::Symmetry,
	config::Config4System;
	verbosity::Bool = true,
)
	return Cluster(structure, symmetry, config.nbody, config.bodyn_cutoff; verbosity = verbosity)
end

"""
	set_mindist_pairs(num_atoms, cartesian_coords, cell_exists; tol=DEFAULT_TOLERANCE)

Computes minimum distance pairs between atoms in a crystal structure.

# Arguments
- `num_atoms::Integer`: Number of atoms in the structure
- `cartesian_coords::AbstractArray{<:AbstractFloat, 3}`: Cartesian coordinates of atoms
- `cell_exists::AbstractVector{Bool}`: Indicates which virtual cells exist
- `tol::Real`: Tolerance for distance comparison (default: DEFAULT_TOLERANCE)

# Returns
- `Matrix{Vector{DistInfo}}`: Matrix containing minimum distance pairs between atoms

# Throws
- `AssertionError`: If input parameters are invalid
"""
function set_mindist_pairs(
	num_atoms::Integer,
	cartesian_coords::AbstractArray{<:AbstractFloat, 3},
	cell_exists::AbstractVector{Bool};
	tol::Real = DEFAULT_TOLERANCE,
)::Matrix{Vector{DistInfo}}
	@assert num_atoms > 0 "Number of atoms must be positive."
	@assert size(cartesian_coords, 2) == num_atoms "Number of atoms in cartesian_coords must match num_atoms."
	@assert length(cell_exists) == NUM_VIRTUAL_CELLS "Length of cell_exists must match NUM_VIRTUAL_CELLS."

	distance_all = Matrix{Vector{DistInfo}}(undef, num_atoms, num_atoms)
	initialized = falses(size(distance_all))

	@inbounds for i in 1:num_atoms, j in 1:num_atoms
		distinfo_list = Vector{DistInfo}()
		for cell_index in 1:NUM_VIRTUAL_CELLS
			cell_exists[cell_index] || continue

			distance = norm(cartesian_coords[:, i, 1] - cartesian_coords[:, j, cell_index])
			relative_vector = cartesian_coords[:, j, cell_index] - cartesian_coords[:, i, 1]
			push!(distinfo_list, DistInfo(cell_index, distance, relative_vector))
		end
		sort!(distinfo_list)
		distance_all[i, j] = distinfo_list
		initialized[i, j] = true
	end

	@assert all(initialized) "Unassigned indices in distance_all variable"

	min_distance_pairs = Matrix{Vector{DistInfo}}(undef, num_atoms, num_atoms)
	initialized = falses(size(min_distance_pairs))

	@inbounds for i in 1:num_atoms, j in 1:num_atoms
		dist_vec_tmp = Vector{DistInfo}()
		min_dist = distance_all[i, j][1].distance
		for distinfo in distance_all[i, j]
			isapprox(distinfo.distance, min_dist, atol = tol) || break
			push!(dist_vec_tmp, distinfo)
		end
		min_distance_pairs[i, j] = dist_vec_tmp
		initialized[i, j] = true
	end

	@assert all(initialized) "Unassigned indices in min_distance_pairs variable"
	return min_distance_pairs
end


"""
	is_within_cutoff(atomcell_list, kd_int_list, cutoff_radii, body, x_image_cart, min_distance_pairs) -> Bool

Checks if all pairs of atoms in the given atom cell list are within the cutoff radius.

# Arguments
- `atomcell_list::Vector{AtomCell}`: List of atom-cell pairs to check
- `kd_int_list::Vector{Int}`: List of atom types
- `cutoff_radii::AbstractArray{<:Real, 3}`: Cutoff radii for interactions
- `body::Integer`: Number of interacting bodies
- `x_image_cart::AbstractArray{<:Real, 3}`: Cartesian coordinates of atoms in different image cells
- `min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}}`: Minimum distance pairs between atoms

# Returns
- `Bool`: `true` if all atom pairs are within cutoff, `false` otherwise
"""
function is_within_cutoff(
	atomcell_list::Vector{AtomCell},
	kd_int_list::Vector{Int},
	cutoff_radii::AbstractArray{<:Real, 3},
	body::Integer,
	x_image_cart::AbstractArray{<:Real, 3},
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Bool
	for comb::Vector{AtomCell} in collect(combinations(atomcell_list, 2))
		rc = cutoff_radii[body, kd_int_list[comb[1].atom], kd_int_list[comb[2].atom]]
		distance = distance_atomcells(comb[1], comb[2], x_image_cart)
		if (rc < 0.0 || distance ≤ rc)
			if min_distance_pairs[comb[1].atom, comb[2].atom][1].distance + 0.00001 > distance
				continue
			else
				return false
			end
		else
			return false
		end
	end
	return true
end


"""
	distance_atomcells(atomcell1, atomcell2, x_image_cart) -> Float64

Calculates the Cartesian distance between two atoms in different cells.

# Arguments
- `atomcell1::AtomCell`: First atom-cell pair
- `atomcell2::AtomCell`: Second atom-cell pair
- `x_image_cart::AbstractArray{<:Real, 3}`: Cartesian coordinates of atoms in different image cells

# Returns
- `Float64`: Distance between the two atoms in Angstroms
"""
function distance_atomcells(
	atomcell1::AtomCell,
	atomcell2::AtomCell,
	x_image_cart::AbstractArray{<:Real, 3},
)::Float64
	return norm(
		x_image_cart[:, atomcell1.atom, atomcell1.cell] -
		x_image_cart[:, atomcell2.atom, atomcell2.cell],
	)
end



"""
	generate_clusters(structure, symmetry, cutoff_radii, nbody) -> Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}}

Generates interaction clusters based on cutoff radii and structure information.

# Arguments
- `structure::Structure`: Crystal structure information
- `symmetry::Symmetry`: Symmetry operations
- `cutoff_radii::AbstractArray{<:Real, 3}`: Cutoff radii for interactions
- `nbody::Integer`: Number of interacting bodies

# Returns
- `Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}}`: Dictionary of clusters organized by body and primitive atom index
"""
function generate_clusters(
	structure::Structure,
	symmetry::Symmetry,
	cutoff_radii::AbstractArray{<:Real, 3},
	nbody::Integer,
)::Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}}

	min_distance_pairs = set_mindist_pairs(
		structure.supercell.num_atoms,
		structure.x_image_cart,
		structure.exist_image,
		tol = symmetry.tol,
	)
	interaction_cutoff_dict = Dict{Int, Dict{Int, SortedVector{AtomCell}}}()
	for body in 2:nbody
		interaction_cutoff_dict[body] = Dict{Int, SortedVector{AtomCell}}()
		for prim_atom_sc in symmetry.atoms_in_prim
			prim_atom_type = structure.supercell.kd_int_list[prim_atom_sc]
			interaction_cutoff_dict[body][prim_atom_sc] = SortedVector{AtomCell}()
			for other_atom_sc in 1:structure.supercell.num_atoms
				if prim_atom_sc == other_atom_sc
					;
					continue;
				end
				other_atom_type = structure.supercell.kd_int_list[other_atom_sc]
				rc = cutoff_radii[body, prim_atom_type, other_atom_type]

				if rc < 0.0 || min_distance_pairs[prim_atom_sc, other_atom_sc][1].distance ≤ rc
					for distinfo in min_distance_pairs[prim_atom_sc, other_atom_sc]
						push!(
							interaction_cutoff_dict[body][prim_atom_sc],
							AtomCell(other_atom_sc, distinfo.cell_index),
						)
					end
				end
			end
		end
	end


	# interaction_clusters[body][prim_atom_sc] = SortedVector{SortedVector{AtomCell}}()
	interaction_clusters = Dict{Int, Dict{Int, Vector{SortedVector{AtomCell}}}}()
	for body in 2:nbody
		interaction_clusters[body] = Dict{Int, SortedVector{AtomCell}}()
		for prim_atom_sc in symmetry.atoms_in_prim
			interaction_clusters[body][prim_atom_sc] = SortedVector{AtomCell}()
		end
	end


	for body in 2:nbody, prim_atom_sc in symmetry.atoms_in_prim
		prim_atom_ac::AtomCell = AtomCell(prim_atom_sc, 1)
		interactiong_atoms::SortedVector = interaction_cutoff_dict[body][prim_atom_sc]
		if body == 2
			for other_atom_ac::AtomCell in interactiong_atoms
				distance = distance_atomcells(prim_atom_ac, other_atom_ac, structure.x_image_cart)
				rc = cutoff_radii[
					body,
					structure.supercell.kd_int_list[prim_atom_sc],
					structure.supercell.kd_int_list[other_atom_ac.atom],
				]
				if rc < 0.0 || distance ≤ rc
					push!(interaction_clusters[body][prim_atom_sc], SortedVector([other_atom_ac]))
				end
			end
		else
			for atom_combination::Vector{AtomCell} in
				collect(combinations(interactiong_atoms, body - 1))
				atom_cell_list_all = vcat([prim_atom_ac], atom_combination)
				atom_list_all = [atom_cell.atom for atom_cell in atom_cell_list_all]
				if (
					is_within_cutoff(
						atom_cell_list_all,
						structure.supercell.kd_int_list,
						cutoff_radii,
						body,
						structure.x_image_cart,
						min_distance_pairs,
					) && (length(atom_list_all) == length(unique(atom_list_all)))
				)
					sorted_vector = SortedVector(atom_combination)
					push!(interaction_clusters[body][prim_atom_sc], sorted_vector)
				end
			end
		end
	end

	result = Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}}()
	for body in 2:nbody
		result[body] = Dict{Int, CountingUniqueVector{Vector{Int}}}()
		for prim_atom_sc in symmetry.atoms_in_prim
			result[body][prim_atom_sc] = CountingUniqueVector{Vector{Int}}()
		end
	end

	for body in 2:nbody
		for prim_atom_sc in symmetry.atoms_in_prim
			counting_unique_vector = CountingUniqueVector{Vector{Int}}()
			for cluster::SortedVector{AtomCell} in interaction_clusters[body][prim_atom_sc]
				atom_list = [atom_cell.atom for atom_cell in cluster]
				atom_list = vcat([prim_atom_sc], atom_list)
				push!(counting_unique_vector, atom_list)
			end
			result[body][prim_atom_sc] = counting_unique_vector
		end
	end

	return result
end

"""
	struct DistList

Internal structure used to sort distance information in `print_cluster_stdout`.

# Fields
- `atom::Int`: Atom index
- `dist::Float64`: Distance value
"""
struct DistList
	atom::Int
	dist::Float64
end

function Base.isless(a::DistList, b::DistList)
	return a.dist < b.dist
end

"""
	print_cluster_stdout(min_distance_pairs, atoms_in_prim, kd_name, kd_int_list)

Prints the list of neighboring atoms and distances for each atom in the primitive cell.

# Arguments
- `min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}}`: Matrix of minimum distance pairs
- `atoms_in_prim::AbstractVector{<:Integer}`: Indices of atoms in the primitive cell
- `kd_name::AbstractVector{<:AbstractString}`: Names of atomic elements
- `kd_int_list::AbstractVector{<:Integer}`: List of atom types

# Throws
- `AssertionError`: If input parameters are invalid
"""
function print_cluster_stdout(
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
	atoms_in_prim::AbstractVector{<:Integer},
	kd_name::AbstractVector{<:AbstractString},
	kd_int_list::AbstractVector{<:Integer},
)
	@assert length(atoms_in_prim) > 0 "atoms_in_prim must not be empty."
	@assert length(kd_name) > 0 "kd_name must not be empty."
	@assert length(kd_int_list) > 0 "kd_int_list must not be empty."

	println("""

	INTERACTION
	===========
	""")

	num_atoms::Int = size(min_distance_pairs, 1)
	neighbor_list = Vector{Vector{DistList}}(undef, length(atoms_in_prim))

	# Create neighbor list for each primitive atom
	for (i_prim, i_prim_atom) in enumerate(atoms_in_prim)
		neighbor_list[i_prim] = DistList[]
		for j in 1:num_atoms
			push!(
				neighbor_list[i_prim],
				DistList(j, min_distance_pairs[i_prim_atom, j][1].distance),
			)
		end
		sort!(neighbor_list[i_prim])
	end

	println(" List of neighboring atoms below.")
	println(" Format [N th-nearest shell, distance (Number of atoms on the shell)]")
	println()

	# Print neighbor information for each primitive atom
	for (i_prim, i_prim_atom) in enumerate(atoms_in_prim)
		nth_nearest = 0
		atom_list = Int[]

		@printf("%5d (%3s): ", i_prim_atom, kd_name[kd_int_list[i_prim_atom]])

		dist_tmp = 0.0

		for j in 1:num_atoms
			# Skip if distance is zero (same atom)
			if neighbor_list[i_prim][j].dist < 1e-8
				continue
			end

			# Check if this is a new distance shell
			if abs(neighbor_list[i_prim][j].dist - dist_tmp) > 1e-6
				# Print previous shell if not empty
				if !isempty(atom_list)
					nth_nearest += 1

					if nth_nearest > 1
						print(" " ^ 13)
					end

					@printf("%3d%10.6f (%3d) -", nth_nearest, dist_tmp, length(atom_list))

					# Print atoms in this shell
					for (k, atom_idx) in enumerate(atom_list)
						if k > 1 && k % 4 == 1
							println()
							print(" " ^ 34)
						end
						@printf("%4d(%3s)", atom_idx, kd_name[kd_int_list[atom_idx]])
					end
					println()
				end

				# Start new shell
				dist_tmp = neighbor_list[i_prim][j].dist
				atom_list = [neighbor_list[i_prim][j].atom]
			else
				# Add to current shell
				push!(atom_list, neighbor_list[i_prim][j].atom)
			end
		end

		# Print the last shell if not empty
		if !isempty(atom_list)
			nth_nearest += 1

			if nth_nearest > 1
				print(" " ^ 13)
			end

			@printf("%3d%10.6f (%3d) -", nth_nearest, dist_tmp, length(atom_list))

			# Print atoms in this shell
			for (k, atom_idx) in enumerate(atom_list)
				if k > 1 && k % 4 == 1
					println()
					print(" " ^ 34)
				end
				@printf("%4d(%3s)", atom_idx, kd_name[kd_int_list[atom_idx]])
			end
			println()
		end
		println("\n")
	end
end


end
