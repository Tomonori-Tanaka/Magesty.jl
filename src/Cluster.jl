"""
	module Clusters

This module provides data structures and functions for managing and analyzing clusters of interacting atoms within a crystal structure. It defines `DistInfo` for storing distance information between atoms, `InteractionCluster` for representing clusters of interacting atoms, and `Cluster` for handling multiple interaction clusters based on specified cutoff radii.

# Types
- **`DistInfo`**: Represents distance information between two atoms.
- **`InteractionCluster`**: Represents a cluster of interacting atoms.
- **`Cluster`**: Manages multiple interaction clusters based on structure parameters.

# Functions
- `Cluster(structure, symmetry, nbody::Int, cutoff_radii)`: Constructs a `Cluster` instance based on the given structure and symmetry information.
- `set_mindist_pairs(nat::Int, xc_in::AbstractArray{<:AbstractFloat, 3}, exist_cell::AbstractVector{Bool})`: Computes minimum distance pairs between atoms.
- `set_interaction_by_cutoff(nat::Int, nat_prim::Int, kd_int_list::AbstractVector{Int}, map_p2s::AbstractMatrix{Int}, cutoff_radii::AbstractArray{Float64, 3}, nbody::Int, mindist_pairs::Matrix{<:Vector{<:DistInfo}})`: Determines interacting pairs based on cutoff radii.
- `set_interaction_clusters(nat_prim::Int, kd_int_list::AbstractVector{Int}, map_p2s::AbstractMatrix{Int}, x_image_cart::AbstractArray{Float64, 3}, exist_image::AbstractVector{Bool}, interaction_pairs::Vector{<:Vector{Int}}, cutoff_radii::AbstractArray{Float64, 3}, mindist_pairs::Matrix{<:Vector{<:DistInfo}}, nbody::Int)`: Organizes interaction clusters.
- `is_within_cutoff(atoms::AbstractVector{<:Int}, kd_int::AbstractVector{<:Int}, rc::AbstractArray{Float64, 3}, mindist_pairs::Matrix{<:AbstractVector{<:DistInfo}})`: Checks if a cluster is within cutoff radii.
- `calc_distmax(atoms::AbstractVector{<:Int}, mindist_pairs::Matrix{<:AbstractVector{<:DistInfo}})`: Calculates the maximum distance within a cluster.

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
const NUM_VIRTUAL_CELLS = 27
const DEFAULT_TOLERANCE = 1e-5
const MIN_NBODY = 2

"""
	struct DistInfo <: MyAbstractDistInfo

Represents distance information between two atoms.

# Fields
- `cell_index::Int`: The virtual cell index of the second atom (1-27 for 3x3x3 supercell).
- `distance::Float64`: The Cartesian distance between the two atoms in Angstroms.
- `relative_vector::Vector{Float64}`: The relative Cartesian vector from the first atom to the second atom in Angstroms.

# Constructor
	DistInfo(cell::Int, distance::Real, relvec::Vector{<:Real})

Creates a new `DistInfo` instance. Ensures that `relvec` has a length of 3.

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

isless(distinfo1::DistInfo, distinfo2::DistInfo) = distinfo1.distance < distinfo2.distance

"""
	struct InteractionCluster

Represents a cluster of interacting atoms within a crystal structure.

# Fields
- `atom_indices::Vector{Int}`: List of partner atom indices (length ≤ body-1)
- `cell_indices::Vector{Vector{Int}}`: List of cell indices for each partner atom (length ≤ body-1)
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

isless(intclus1::InteractionCluster, intclus2::InteractionCluster) =
	intclus1.atom_indices < intclus2.atom_indices

function ==(intclus1::InteractionCluster, intclus2::InteractionCluster)
	return intclus1.atom_indices == intclus2.atom_indices &&
	       intclus1.cell_indices == intclus2.cell_indices &&
	       intclus1.max_distance == intclus2.max_distance
end


"""
	struct Cluster

Represents a collection of interaction clusters based on the specified number of bodies and cutoff radii.

# Fields
- `num_bodies::Int`: The number of interacting bodies.
- `cutoff_radii::Array{Float64, 3}`: Cutoff radii for each atomic element pair and interaction body. Dimensions: [nkd, nkd, nbody].
- `min_distance_pairs::Matrix{Vector{DistInfo}}`: Matrix containing the minimum distance pairs between atoms. Dimensions: [num_atoms, num_atoms].
- `interaction_clusters::Matrix{OrderedSet{InteractionCluster}}`: Matrix of interaction clusters for each primitive atom and interaction body. Dimensions: [nat_prim, nbody-1].
- `cluster_list::Vector{SortedVector{Vector{AtomCell}}}`: List of interacting atom clusters for each interaction body.
- `equivalent_atom_list::Vector{Vector{Int}}`: List of equivalent atom groups based on symmetry operations.
- `elapsed_time::Float64`: Time taken to create the cluster in seconds.

# Constructor
	Cluster(structure, symmetry, nbody::Int, cutoff_radii)

Creates a new `Cluster` instance based on the provided structure, symmetry information, number of bodies, and cutoff radii.

# Example
```julia
cluster = Cluster(structure, symmetry, 3, cutoff_radii)
"""
struct Cluster
	num_bodies::Int
	cutoff_radii::Array{Float64, 3}
	min_distance_pairs::Matrix{Vector{DistInfo}} # [≤ num_atoms, ≤ num_atoms]
	interaction_clusters::Matrix{OrderedSet{InteractionCluster}}# [≤ nat_prim, ≤ nbody-1]
	cluster_list::Vector{SortedVector{Vector{AtomCell}}} # [≤ nbody-1][≤ number of clusters][≤ nbody]
	equivalent_atom_list::Vector{Vector{Int}}
	elapsed_time::Float64  # Time taken to create the cluster in seconds

	function Cluster(
		structure::Structure,
		symmetry::Symmetry,
		nbody::Integer,
		cutoff_radii::AbstractArray{<:Real, 3},
	)
		@assert nbody ≥ MIN_NBODY "Number of bodies must be at least $MIN_NBODY."
		@assert size(cutoff_radii, 3) == nbody "Cutoff radii dimensions must match nbody."

		# Start timing
		start_time = time_ns()
		
		min_distance_pairs =
			set_mindist_pairs(
				structure.supercell.num_atoms,
				structure.x_image_cart,
				structure.exist_image,
				tol = symmetry.tol,
			)
		# interaction_pairs[≤ nat_prim][nbody-1]
		interaction_pairs::Vector{Vector{Vector{Int}}} = set_interaction_by_cutoff(
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
			cutoff_radii,
			min_distance_pairs,
			nbody,
		)
		cluster_list =
			generate_pairs(symmetry.atoms_in_prim, interaction_clusters, nbody)

		equivalent_atom_list =
			classify_equivalent_atoms(symmetry.atoms_in_prim, symmetry.map_sym)

		# End timing
		elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds

		return new(
			nbody,
			cutoff_radii,
			min_distance_pairs,
			interaction_clusters,
			cluster_list,
			equivalent_atom_list,
			elapsed_time
		)
	end
end

function Cluster(structure::Structure, symmetry::Symmetry, config::Config4System)
	return Cluster(structure, symmetry, config.nbody, config.cutoff_radii)
end

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
	for iat in 1:num_atoms
		for jat in 1:num_atoms
			distinfo_list = Vector{DistInfo}()
			for icell in 1:NUM_VIRTUAL_CELLS
				if !(cell_exists[icell])
					continue
				end

				dist_tmp::Float64 = norm(cartesian_coords[:, iat, 1] - cartesian_coords[:, jat, icell])
				relvec::Vector{Float64} = cartesian_coords[:, jat, icell] - cartesian_coords[:, iat, 1]
				push!(distinfo_list, DistInfo(icell, dist_tmp, relvec))
			end
			sort!(distinfo_list)
			distance_all[iat, jat] = distinfo_list
			initialized[iat, jat] = true
		end
	end
	# check assigned or not
	undef_indices = findall(x -> x == false, initialized)
	if length(undef_indices) >= 1
		error("unassigned indices in distance_all variable: $undef_indices")
	end

	min_distance_pairs = Matrix{Vector{DistInfo}}(undef, num_atoms, num_atoms)
	initialized = falses(size(min_distance_pairs))

	for iat in 1:num_atoms
		for jat in 1:num_atoms
			dist_vec_tmp = Vector{DistInfo}()
			dist_min = distance_all[iat, jat][1].distance
			for distinfo in distance_all[iat, jat]
				if isapprox(distinfo.distance, dist_min, atol=tol)
					push!(dist_vec_tmp, distinfo)
				end
			end
			min_distance_pairs[iat, jat] = dist_vec_tmp
			initialized[iat, jat] = true
		end
	end
	# check assigned or not
	undef_indices = findall(x -> x == false, initialized)
	if length(undef_indices) >= 1
		error("unassigned indices in min_distance_pairs variable: $undef_indices")
	end

	return min_distance_pairs
end

function set_interaction_by_cutoff(
	num_atoms::Integer,
	num_primitive_atoms::Integer,
	atom_type_list::AbstractVector{<:Integer},
	primitive_to_supercell_map::AbstractMatrix{<:Integer},
	cutoff_radii::AbstractArray{<:Real, 3},
	nbody::Integer,
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Vector{Vector{Vector{Int}}}

	interaction_pairs = Vector{Vector{Vector{Int}}}()
	for body in 2:nbody # 1 is skipped cuz it's 1-body term.
		pairs_tmp_a_body = Vector{Vector{Int}}()
		for iat_prim in 1:num_primitive_atoms
			pair_tmp = Int[]
			iat = primitive_to_supercell_map[iat_prim, 1]  # index in the central cell
			ikd = atom_type_list[iat]

			for jat in 1:num_atoms
				jkd = atom_type_list[jat]
				cutoff_tmp::Float64 = cutoff_radii[ikd, jkd, body]
				if iat == jat
					continue
				elseif cutoff_tmp < 0.0
					append!(pair_tmp, jat)
				elseif min_distance_pairs[iat, jat][1].distance ≤ cutoff_tmp
					append!(pair_tmp, jat)
				end
			end
			push!(pairs_tmp_a_body, pair_tmp)
		end
		push!(interaction_pairs, pairs_tmp_a_body)
	end
	return interaction_pairs
end

function set_interaction_clusters(
	num_primitive_atoms::Integer,
	atom_type_list::AbstractVector{<:Integer},
	primitive_to_supercell_map::AbstractMatrix{<:Integer},
	cartesian_coords::AbstractArray{<:Real, 3},
	cell_exists::AbstractVector{Bool},
	interaction_pairs::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}},
	cutoff_radii::AbstractArray{<:Real, 3},
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
	nbody::Integer,
)::Matrix{OrderedSet{InteractionCluster}}
	interaction_clusters =
		Matrix{OrderedSet{InteractionCluster}}(undef, num_primitive_atoms, nbody - 1)
	initialized = falses(size(interaction_clusters))
	for body in 2:nbody
		for iat_prim in 1:num_primitive_atoms
			inter_cluster_set = OrderedSet{InteractionCluster}()
			iat = primitive_to_supercell_map[iat_prim, 1]
			# intlist = interaction_pairs[iat_prim, body-1]
			intlist::Vector{Int} = interaction_pairs[body-1][iat_prim]
			sort!(intlist)

			if body == 2# 2-body terms
				for jat in intlist
					atoms_tmp = [jat]
					cell_tmp = Int[]
					for distinfo::DistInfo in min_distance_pairs[iat, jat]
						append!(cell_tmp, distinfo.cell_index)
					end
					dist_max = min_distance_pairs[iat, jat][1].distance
					push!(
						inter_cluster_set,
						InteractionCluster(atoms_tmp, [cell_tmp], dist_max),
					)
				end
			else# 3- or more body terms
				for combi in collect(combinations(intlist, body - 1))
					atoms_tmp = combi
					cell_tmp = Vector{Vector{Int}}()
					if !is_within_cutoff(
						vcat([iat], combi),
						atom_type_list,
						cutoff_radii,
						min_distance_pairs,
					)
						continue
					end
					dist_max = calc_distmax(vcat([iat], combi), min_distance_pairs)
					for jat in combi
						cell_j_tmp = []
						for distinfo in min_distance_pairs[iat, jat]
							push!(cell_j_tmp, distinfo.cell_index)
						end
						push!(cell_tmp, cell_j_tmp)
					end
					push!(
						inter_cluster_set,
						InteractionCluster(atoms_tmp, cell_tmp, dist_max),
					)
				end
			end
			interaction_clusters[iat_prim, body-1] = deepcopy(inter_cluster_set)
			initialized[iat_prim, body-1] = true
		end
	end

	unassigned_indices = findall(x -> x == false, initialized)
	if length(unassigned_indices) >= 1
		error("unassigned indices in `interaction_clusters` variable: $unassigned_indices")
	end
	return interaction_clusters
end

"""
For combinations of three or more atoms, determine whether all of them are within the cutoff radius.
"""
function is_within_cutoff(
	atom_indices::AbstractVector{<:Integer},
	atom_types::AbstractVector{<:Integer},
	cutoff_radii::AbstractArray{<:Real, 3},
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Bool
	for pair in collect(combinations(atom_indices, 2))
		if (min_distance_pairs[pair[1], pair[2]])[1].distance ≥ cutoff_radii[atom_types[pair[1]], atom_types[pair[2]], 1]
			return false
		end
	end
	return true
end

function calc_distmax(
	atom_indices::AbstractVector{<:Integer},
	min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Float64
	dist_tmp = -1.0
	for pair in collect(combinations(atom_indices, 2))
		if dist_tmp < min_distance_pairs[pair[1], pair[2]].distance
			dist_tmp = min_distance_pairs[pair[1], pair[2]].distance
		end
	end
	return dist_tmp
end

function generate_pairs(
	primitive_atom_indices::AbstractVector{<:Integer},
	interactoin_clusters::AbstractMatrix{<:AbstractSet{InteractionCluster}},# [≤ nat_prim, ≤ nbody-1]
	nbody::Integer,
)::Vector{SortedVector{Vector{AtomCell}}}
	cluster_list = Vector{SortedVector{Vector{AtomCell}}}()
	vsv = Vector{SortedVector{Vector{AtomCell}}}()
	for body in 2:nbody
		sv = SortedVector{Vector{AtomCell}}()
		for (i, iat) in enumerate(primitive_atom_indices)
			for intclus_i::InteractionCluster in interactoin_clusters[i, body-1]
				pairs = generate_combinations(intclus_i.atom_indices, intclus_i.cell_indices)
				# examples of partners variable
				# ((2, 1),) for 2-body case
				# ((2, 1), (3, 1)) for 3-body case
				for partners in pairs
					vec_tmp = Vector{AtomCell}([AtomCell(iat, 1)])
					for tuple in partners
						push!(vec_tmp, AtomCell(tuple[1], tuple[2]))
					end
					push!(sv, vec_tmp)
				end
			end
		end
		push!(vsv, sv)
	end
	return vsv
end

"""
	classify_equivalent_atoms(atoms_in_prim::AbstractVector{<:Integer}, map_sym::AbstractMatrix{<:Integer}) -> Vector{Vector{Int}}

Classifies atoms in the primitive cell into equivalent groups based on symmetry operations.

# Arguments
- `atoms_in_prim::AbstractVector{<:Integer}`: List of atom indices in the primitive cell
- `map_sym::AbstractMatrix{<:Integer}`: Symmetry mapping matrix where map_sym[i,j] gives the image of atom i under symmetry operation j

# Returns
- `Vector{Vector{Int}}`: List of equivalent atom groups, where each group contains atoms that are related by symmetry operations

# Examples
```julia
# For a structure with 4 atoms and 2 symmetry operations
atoms = [1, 2, 3, 4]
map_sym = [1 2; 2 1; 3 4; 4 3]  # Example symmetry mapping

# Classify equivalent atoms
groups = classify_equivalent_atoms(atoms, map_sym)

# Result will be:
# [[1, 2], [3, 4]]  # Atoms 1,2 are equivalent and 3,4 are equivalent
```
"""
function classify_equivalent_atoms(
	primitive_atom_indices::AbstractVector{<:Integer},
	map_sym::AbstractMatrix{<:Integer},
)::Vector{Vector{Int}}
	group = Vector{Vector{Int}}()
	checked = falses(size(primitive_atom_indices))

	nsym = size(map_sym, 2)
	for (idx, atom) in enumerate(primitive_atom_indices)
		if checked[idx]
			continue
		end

		eqlist = Int[atom]

		for isym in 1:nsym
			eqatom = map_sym[atom, isym]
			if atom == eqatom
				continue
			elseif eqatom in primitive_atom_indices && !(eqatom in eqlist)
				eq_idx = findfirst(x -> x == eqatom, primitive_atom_indices)
				checked[eq_idx] = true
				push!(eqlist, eqatom)
			end
		end
		push!(group, eqlist)
	end

	return group
end

"""
	generate_combinations(vec::AbstractVector{<:Any}, vecofvec::AbstractVector{<:AbstractVector{<:Any}}) -> Vector{Tuple}

Generates all possible combinations by pairing each element of `vec` with each element of the corresponding sublist in `vecofvec`.

# Arguments
- `vec::AbstractVector{<:Any}`: A vector containing elements to be paired
- `vecofvec::AbstractVector{<:AbstractVector{<:Any}}`: A vector of vectors, where each subvector contains elements to be paired with the corresponding element in `vec`

# Returns
- `Vector{Tuple}`: A vector of tuples, where each tuple contains paired elements from `vec` and `vecofvec`

# Examples
```julia
# Define the vectors
atoms = [1, 2]
cells = [[1, 2], [3]]

# Generate combinations
result = generate_combinations(atoms, cells)

# Result will be:
# [((1, 1), (2, 3)),
#  ((1, 2), (2, 3))]
```
"""
function generate_combinations(
	vec::AbstractVector,
	vecofvec::AbstractVector{<:AbstractVector},
)::Vector{Tuple}
	if length(vec) != length(vecofvec)
		error("Vectors A and B must have the same length.")
	end
	combination = Iterators.product(vecofvec...)
	return [Tuple(zip(vec, comb)) for comb in combination]
end

"""
	all_atomlist_by_symop(atomlist::AbstractVector{<:Integer}, map_sym::AbstractMatrix{<:Integer}) -> Vector{Vector{Int}}

Generates all possible atom lists by applying symmetry operations to the input atom list.

# Arguments
- `atomlist::AbstractVector{<:Integer}`: List of atom indices to transform
- `map_sym::AbstractMatrix{<:Integer}`: Symmetry mapping matrix where map_sym[i,j] gives the image of atom i under symmetry operation j

# Returns
- `Vector{Vector{Int}}`: List of unique atom lists obtained by applying all symmetry operations

# Examples
```julia
# For a structure with 3 atoms and 2 symmetry operations
atoms = [1, 2, 3]
map_sym = [1 2; 2 1; 3 3]  # Example symmetry mapping

# Generate all possible atom lists
lists = all_atomlist_by_symop(atoms, map_sym)

# Result will be:
# [[1, 2, 3], [2, 1, 3]]  # Two unique configurations
```
"""
function all_atomlist_by_symop(
	atom_indices::AbstractVector{<:Integer},
	map_sym::AbstractMatrix{<:Integer},
)::Vector{Vector{Int}}
	# total number of symmetry operations
	nsym = size(map_sym, 2)
	atomlist_list = Vector{Vector{Int}}()

	for i in 1:nsym
		atomlist_tmp = Int[]
		for atom in atom_indices
			push!(atomlist_tmp, map_sym[atom, i])
		end
		push!(atomlist_list, atomlist_tmp)
	end

	return unique(atomlist_list)

end

function print_info(cluster::Cluster)
	println("""
	========
	CLUSTER
	========
	""")
	println("Number of bodies: ", cluster.num_bodies)
	println("Number of interaction clusters: ", length(cluster.cluster_list))
	println("")
	println("Elapsed time: ", cluster.elapsed_time, " seconds")
	println("-------------------------------------------------------------------")
end

end
