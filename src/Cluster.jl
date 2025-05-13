"""
	module Clusters

This module provides data structures and functions for managing and analyzing clusters of interacting atoms within a crystal structure.

# Types
- `DistInfo`: Represents distance information between two atoms
- `InteractionCluster`: Represents a cluster of interacting atoms
- `Cluster`: Manages multiple interaction clusters based on structure parameters

# Functions
- `Cluster(structure, symmetry, nbody, cutoff_radii)`: Constructs a `Cluster` instance
- `set_mindist_pairs(num_atoms, cartesian_coords, cell_exists)`: Computes minimum distance pairs
- `set_interaction_by_cutoff(num_atoms, num_primitive_atoms, atom_type_list, ...)`: Determines interacting pairs
- `set_interaction_clusters(num_primitive_atoms, atom_type_list, ...)`: Organizes interaction clusters
- `is_within_cutoff(atom_indices, atom_types, cutoff_radii, min_distance_pairs)`: Checks cutoff conditions
- `calc_distmax(atom_indices, min_distance_pairs)`: Calculates maximum distance
- `generate_pairs(primitive_atom_indices, interaction_clusters, num_bodies)`: Generates atom pairs
- `classify_equivalent_atoms(primitive_atom_indices, symmetry_map)`: Classifies equivalent atoms
- `generate_combinations(atom_indices, cell_indices)`: Generates combinations
- `all_atomlist_by_symop(atom_indices, symmetry_map)`: Generates atom lists by symmetry
- `print_info(cluster)`: Prints cluster information
"""

module Clusters

using Combinatorics: combinations
using DataStructures
using LinearAlgebra
using Printf

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
	struct InteractionCluster

Represents a cluster of interacting atoms within a crystal structure.

# Fields
- `atom_indices::Vector{Int}`: List of partner atom indices (length ≤ body-1)
- `cell_indices::Vector{Vector{Int}}`: List of cell indices for each partner atom
- `max_distance::Float64`: Maximum distance between any pair of atoms in the cluster

# Examples
```julia
# Create a 2-body interaction cluster
cluster = InteractionCluster([2], [[1]], 2.5)

# Create a 3-body interaction cluster
cluster = InteractionCluster([2, 3], [[1], [2]], 3.0)
```
"""
struct InteractionCluster
	atom_indices::Vector{Int}
	cell_indices::Vector{Vector{Int}}
	max_distance::Float64

	function InteractionCluster(atoms::Vector{Int}, cells::Vector{Vector{Int}}, distmax::Real)
		@assert length(atoms) > 0 "Atoms vector must not be empty."
		@assert length(atoms) == length(cells) "Length of atoms and cells must match."
		@assert distmax ≥ 0.0 "Maximum distance must be non-negative."
		return new(atoms, cells, Float64(distmax))
	end
end

Base.isless(intclus1::InteractionCluster, intclus2::InteractionCluster) = intclus1.atom_indices < intclus2.atom_indices

function ==(intclus1::InteractionCluster, intclus2::InteractionCluster)
	return intclus1.atom_indices == intclus2.atom_indices &&
	       intclus1.cell_indices == intclus2.cell_indices &&
	       intclus1.max_distance == intclus2.max_distance
end

"""
	struct Cluster

Represents a collection of interaction clusters based on the specified number of bodies and cutoff radii.

# Fields
- `num_bodies::Int`: Number of interacting bodies
- `cutoff_radii::Array{Float64, 3}`: Cutoff radii for each atomic element pair and interaction body
- `min_distance_pairs::Matrix{Vector{DistInfo}}`: Matrix of minimum distance pairs between atoms
- `interaction_clusters::Matrix{OrderedSet{InteractionCluster}}`: Matrix of interaction clusters
- `cluster_list::Vector{SortedVector{Vector{AtomCell}}}`: List of interacting atom clusters
- `equivalent_atom_list::Vector{Vector{Int}}`: List of equivalent atom groups
- `elapsed_time::Float64`: Time taken to create the cluster in seconds

# Constructor
	Cluster(structure, symmetry, nbody, cutoff_radii)

Creates a new `Cluster` instance based on the provided structure, symmetry information, number of bodies, and cutoff radii.

# Example
```julia
cluster = Cluster(structure, symmetry, 3, cutoff_radii)
"""
struct Cluster
	num_bodies::Int
	cutoff_radii::Array{Float64, 3}
	min_distance_pairs::Matrix{Vector{DistInfo}}
	# interaction_clusters::Matrix{OrderedSet{InteractionCluster}}
	cluster_list::Vector{SortedVector{Vector{AtomCell}}}
	equivalent_atom_list::Vector{Vector{Int}}
	elapsed_time::Float64

	function Cluster(
		structure::Structure,
		symmetry::Symmetry,
		nbody::Integer,
		cutoff_radii::AbstractArray{<:Real, 3},
	)
		@assert size(cutoff_radii, 3) == nbody "Cutoff radii dimensions must match nbody."

		start_time = time_ns()
		
		min_distance_pairs = set_mindist_pairs(
			structure.supercell.num_atoms,
			structure.x_image_cart,
			structure.exist_image,
			tol = symmetry.tol,
		)

		interaction_pairs = set_interaction_by_cutoff(
			structure.supercell.num_atoms,
			symmetry.nat_prim,
			structure.supercell.kd_int_list,
			symmetry.map_p2s,
			cutoff_radii,
			nbody,
			min_distance_pairs,
		)

		interaction_clusters = set_interaction_clusters(
			symmetry.nat_prim,
			structure.supercell.kd_int_list,
			symmetry.map_p2s,
			structure.x_image_cart,
			structure.exist_image,
			interaction_pairs,
			min_distance_pairs,
			nbody,
		)

		cluster_list = generate_pairs(symmetry.atoms_in_prim, interaction_clusters, nbody)
		equivalent_atom_list = classify_equivalent_atoms(symmetry.atoms_in_prim, symmetry.map_sym)
		elapsed_time = (time_ns() - start_time) / 1e9

		return new(
			nbody,
			cutoff_radii,
			min_distance_pairs,
			# interaction_clusters,
			cluster_list,
			equivalent_atom_list,
			elapsed_time,
		)
	end
end

function Cluster(structure::Structure, symmetry::Symmetry, config::Config4System)
	return Cluster(structure, symmetry, config.nbody, config.cutoff_radii)
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
	set_interaction_by_cutoff(num_atoms, num_primitive_atoms, atom_type_list, ...)

Determines interacting pairs of atoms based on cutoff radii.

# Arguments
- `num_atoms::Integer`: Number of atoms in the structure
- `num_primitive_atoms::Integer`: Number of atoms in the primitive cell
- `atom_type_list::AbstractVector{<:Integer}`: List of atom types
- `primitive_to_supercell_map::AbstractMatrix{<:Integer}`: Mapping from primitive to supercell
- `cutoff_radii::AbstractArray{<:Real, 3}`: Cutoff radii for interactions
- `num_bodies::Integer`: Number of interacting bodies
- `min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}}`: Minimum distance pairs

# Returns
- `Vector{Vector{Vector{Int}}}`: List of interacting pairs for each body

# Throws
- `AssertionError`: If input parameters are invalid
"""
function set_interaction_by_cutoff(
	num_atoms::Integer,
	num_primitive_atoms::Integer,
	atom_type_list::AbstractVector{<:Integer},
	primitive_to_supercell_map::AbstractMatrix{<:Integer},
	cutoff_radii::AbstractArray{<:Real, 3},
	num_bodies::Integer,
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Vector{Vector{Vector{Int}}}
	@assert num_atoms > 0 "Number of atoms must be positive."
	@assert num_primitive_atoms > 0 "Number of primitive atoms must be positive."
	@assert length(atom_type_list) == num_atoms "Length of atom_type_list must match num_atoms."
	@assert size(primitive_to_supercell_map, 1) == num_primitive_atoms "Number of rows in primitive_to_supercell_map must match num_primitive_atoms."

	interaction_pairs = Vector{Vector{Vector{Int}}}()

	@inbounds for body in 2:num_bodies
		pairs_per_body = Vector{Vector{Int}}()
		for prim_atom_index in 1:num_primitive_atoms
			pair_list = Int[]
			atom_index = primitive_to_supercell_map[prim_atom_index, 1]
			atom_type = atom_type_list[atom_index]

			for other_atom_index in 1:num_atoms
				atom_index == other_atom_index && continue
				other_atom_type = atom_type_list[other_atom_index]
				cutoff = cutoff_radii[atom_type, other_atom_type, body]
				if cutoff < 0.0 || min_distance_pairs[atom_index, other_atom_index][1].distance ≤ cutoff
					push!(pair_list, other_atom_index)
				end
			end
			push!(pairs_per_body, pair_list)
		end
		push!(interaction_pairs, pairs_per_body)
	end
	return interaction_pairs
end

"""
	set_interaction_clusters(num_primitive_atoms, atom_type_list, ...)

Organizes interaction clusters based on cutoff radii and minimum distances.

# Arguments
- `num_primitive_atoms::Integer`: Number of atoms in the primitive cell
- `atom_type_list::AbstractVector{<:Integer}`: List of atom types
- `primitive_to_supercell_map::AbstractMatrix{<:Integer}`: Mapping from primitive to supercell
- `cartesian_coords::AbstractArray{<:Real, 3}`: Cartesian coordinates of atoms
- `cell_exists::AbstractVector{Bool}`: Indicates which virtual cells exist
- `interaction_pairs::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}}`: Interacting pairs
- `cutoff_radii::AbstractArray{<:Real, 3}`: Cutoff radii for interactions
- `min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}}`: Minimum distance pairs
- `num_bodies::Integer`: Number of interacting bodies

# Returns
- `Matrix{OrderedSet{InteractionCluster}}`: Matrix of interaction clusters

# Throws
- `AssertionError`: If input parameters are invalid
"""
function set_interaction_clusters(
	num_primitive_atoms::Integer,
	atom_type_list::AbstractVector{<:Integer},
	primitive_to_supercell_map::AbstractMatrix{<:Integer},
	cartesian_coords::AbstractArray{<:Real, 3},
	cell_exists::AbstractVector{Bool},
	interaction_pairs::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}},
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
	num_bodies::Integer,
)::Matrix{OrderedSet{InteractionCluster}}
	@assert num_primitive_atoms > 0 "Number of primitive atoms must be positive."
	@assert length(atom_type_list) > 0 "atom_type_list must not be empty."
	@assert size(primitive_to_supercell_map, 1) == num_primitive_atoms "Number of rows in primitive_to_supercell_map must match num_primitive_atoms."

	interaction_clusters = Matrix{OrderedSet{InteractionCluster}}(undef, num_primitive_atoms, num_bodies - 1)
	initialized = falses(size(interaction_clusters))

	@inbounds for body in 2:num_bodies, prim_atom_index in 1:num_primitive_atoms
		cluster_set = OrderedSet{InteractionCluster}()
		atom_index = primitive_to_supercell_map[prim_atom_index, 1]
		interacting_atoms = interaction_pairs[body-1][prim_atom_index]
		sort!(interacting_atoms)

		if body == 2
			for other_atom_index in interacting_atoms
				atom_indices = [other_atom_index]
				cell_indices = Int[]
				for distinfo in min_distance_pairs[atom_index, other_atom_index]
					push!(cell_indices, distinfo.cell_index)
				end
				distance = min_distance_pairs[atom_index, other_atom_index][1].distance
				push!(cluster_set, InteractionCluster(atom_indices, [cell_indices], distance))
			end
		else
			for atom_combination in collect(combinations(interacting_atoms, body - 1))
				cell_indices = Vector{Vector{Int}}()
				distance = calc_distmax(vcat([atom_index], atom_combination), min_distance_pairs)
				for other_atom_index in atom_combination
					cell_list = Int[]
					for distinfo in min_distance_pairs[atom_index, other_atom_index]
						push!(cell_list, distinfo.cell_index)
					end
					push!(cell_indices, cell_list)
				end
				push!(cluster_set, InteractionCluster(atom_combination, cell_indices, distance))
			end
		end
		interaction_clusters[prim_atom_index, body-1] = cluster_set
		initialized[prim_atom_index, body-1] = true
	end

	@assert all(initialized) "Unassigned indices in interaction_clusters variable"
	return interaction_clusters
end

"""
	calc_distmax(atom_indices, min_distance_pairs)

Calculates the maximum distance between any pair of atoms in a cluster.

# Arguments
- `atom_indices::AbstractVector{<:Integer}`: Indices of atoms in the cluster
- `min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}}`: Minimum distance pairs

# Returns
- `Float64`: Maximum distance between any pair of atoms

# Throws
- `AssertionError`: If input parameters are invalid
"""
function calc_distmax(
	atom_indices::AbstractVector{<:Integer},
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Float64
	@assert length(atom_indices) > 1 "atom_indices must contain at least 2 elements."

	max_distance = -1.0
	@inbounds for pair in collect(combinations(atom_indices, 2))
		max_distance = max(max_distance, min_distance_pairs[pair[1], pair[2]][1].distance)
	end
	return max_distance
end

"""
	generate_pairs(primitive_atom_indices, interaction_clusters, num_bodies)

Generates pairs of interacting atoms for each interaction body.

# Arguments
- `primitive_atom_indices::AbstractVector{<:Integer}`: Indices of atoms in the primitive cell
- `interaction_clusters::AbstractMatrix{<:AbstractSet{InteractionCluster}}`: Interaction clusters
- `num_bodies::Integer`: Number of interacting bodies

# Returns
- `Vector{SortedVector{Vector{AtomCell}}}`: List of interacting atom pairs

# Throws
- `AssertionError`: If input parameters are invalid
"""
function generate_pairs(
	primitive_atom_indices::AbstractVector{<:Integer},
	interaction_clusters::AbstractMatrix{<:AbstractSet{InteractionCluster}},
	num_bodies::Integer,
)::Vector{SortedVector{Vector{AtomCell}}}
	@assert length(primitive_atom_indices) > 0 "primitive_atom_indices must not be empty."

	cluster_list = Vector{SortedVector{Vector{AtomCell}}}()

	@inbounds for body in 2:num_bodies
		sorted_clusters = SortedVector{Vector{AtomCell}}()
		for (i, prim_atom_index) in enumerate(primitive_atom_indices)
			for cluster in interaction_clusters[i, body-1]	# target is i-th primitive atom
				pairs::Vector{Vector{Tuple{Int, Int}}} = generate_combinations(cluster.atom_indices, cluster.cell_indices)
				for partners::Vector{Tuple{Int, Int}} in pairs
					atom_cells = Vector{AtomCell}([AtomCell(prim_atom_index, 1)])
					for (atom_index, cell_index) in partners
						push!(atom_cells, AtomCell(atom_index, cell_index))
					end
					push!(sorted_clusters, atom_cells)
				end
			end
		end
		push!(cluster_list, sorted_clusters)
	end
	return cluster_list
end

"""
	classify_equivalent_atoms(primitive_atom_indices, symmetry_map)

Classifies atoms in the primitive cell into equivalent groups based on symmetry operations.

# Arguments
- `primitive_atom_indices::AbstractVector{<:Integer}`: Indices of atoms in the primitive cell
- `symmetry_map::AbstractMatrix{<:Integer}`: Symmetry mapping matrix

# Returns
- `Vector{Vector{Int}}`: List of equivalent atom groups

# Throws
- `AssertionError`: If input parameters are invalid
"""
function classify_equivalent_atoms(
	primitive_atom_indices::AbstractVector{<:Integer},
	symmetry_map::AbstractMatrix{<:Integer},
)::Vector{Vector{Int}}
	@assert length(primitive_atom_indices) > 0 "primitive_atom_indices must not be empty."
	@assert size(symmetry_map, 1) ≥ maximum(primitive_atom_indices) "symmetry_map must cover all atoms in primitive_atom_indices."

	equivalent_groups = Vector{Vector{Int}}()
	checked = falses(size(primitive_atom_indices))
	num_symmetry_ops = size(symmetry_map, 2)

	@inbounds for (idx, atom_index) in enumerate(primitive_atom_indices)
		checked[idx] && continue

		equivalent_atoms = Int[atom_index]
		for sym_op in 1:num_symmetry_ops
			equivalent_atom = symmetry_map[atom_index, sym_op]
			atom_index == equivalent_atom && continue
			if equivalent_atom in primitive_atom_indices && !(equivalent_atom in equivalent_atoms)
				eq_idx = findfirst(x -> x == equivalent_atom, primitive_atom_indices)
				checked[eq_idx] = true
				push!(equivalent_atoms, equivalent_atom)
			end
		end
		push!(equivalent_groups, equivalent_atoms)
	end

	return equivalent_groups
end

"""
	generate_combinations(atom_indices, cell_indices)

Generates all possible combinations of atom indices and cell indices.

# Arguments
- `atom_indices::AbstractVector`: List of atom indices
- `cell_indices::AbstractVector{<:AbstractVector}`: List of cell indices for each atom

# Returns
- `Vector{Vector{Tuple{Int, Int}}}`: List of tuples containing atom and cell index pairs

# Examples
```julia
julia> generate_combinations([2, 3], [[1, 4], [5, 6]])
4-element Vector{Vector{Tuple{Int64, Int64}}}:
 [(2, 1), (3, 5)]
 [(2, 1), (3, 6)]
 [(2, 4), (3, 5)]
 [(2, 4), (3, 6)]
```

# Throws
- `AssertionError`: If input parameters are invalid
"""
function generate_combinations(
	atom_indices::AbstractVector,
	cell_indices::AbstractVector{<:AbstractVector},
)::Vector{Vector{Tuple{Int, Int}}}
	@assert length(atom_indices) == length(cell_indices) "Vectors must have the same length."
	@assert length(atom_indices) > 0 "Vectors must not be empty."
	
	result = Vector{Vector{Tuple{Int, Int}}}()
	for comb in Iterators.product(cell_indices...)
		pairs = [(atom_idx, cell_idx) for (atom_idx, cell_idx) in zip(atom_indices, comb)]
		push!(result, pairs)
	end
	return result
end

"""
	all_atomlist_by_symop(atom_indices, symmetry_map)

Generates all possible atom lists by applying symmetry operations.

# Arguments
- `atom_indices::AbstractVector{<:Integer}`: List of atom indices
- `symmetry_map::AbstractMatrix{<:Integer}`: Symmetry mapping matrix

# Returns
- `Vector{Vector{Int}}`: List of unique atom lists

# Throws
- `AssertionError`: If input parameters are invalid
"""
function all_atomlist_by_symop(
	atom_indices::AbstractVector{<:Integer},
	symmetry_map::AbstractMatrix{<:Integer},
)::Vector{Vector{Int}}
	@assert length(atom_indices) > 0 "atom_indices must not be empty."
	@assert size(symmetry_map, 1) ≥ maximum(atom_indices) "symmetry_map must cover all atoms in atom_indices."

	num_symmetry_ops = size(symmetry_map, 2)
	atom_lists = Vector{Vector{Int}}()

	@inbounds for sym_op in 1:num_symmetry_ops
		transformed_atoms = Int[]
		for atom_index in atom_indices
			push!(transformed_atoms, symmetry_map[atom_index, sym_op])
		end
		push!(atom_lists, transformed_atoms)
	end

	return unique(atom_lists)
end

"""
	print_info(cluster)

Prints information about a cluster.

# Arguments
- `cluster::Cluster`: The cluster to print information about
"""
function print_info(cluster::Cluster)
	println("""
	========
	CLUSTER
	========
	""")
	println("Number of interacting bodies (nbody): ", cluster.num_bodies)
	println("Number of interaction clusters:")
	for i in 2:cluster.num_bodies
		println("\t$i-body: ", length(cluster.cluster_list[i-1]))
		for j in eachindex(cluster.cluster_list[i-1])
			cluster_str = join(["(atom: $(ac.atom), cell: $(ac.cell))" for ac in cluster.cluster_list[i-1][j]], ", ")
			println("\t\t[$cluster_str]")
		end
	end
	println("\nElapsed time: ", cluster.elapsed_time, " seconds")
	println("-------------------------------------------------------------------")
end

end
