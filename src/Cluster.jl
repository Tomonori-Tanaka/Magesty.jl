"""
	module Clusters

This module provides data structures and functions for managing and analyzing clusters of interacting atoms within a crystal structure. It defines `DistInfo` for storing distance information between atoms, `InteractionCluster` for representing clusters of interacting atoms, and `Cluster` for handling multiple interaction clusters based on specified cutoff radii.

# Types
- **`DistInfo`**: Represents distance information between two atoms.
- **`InteractionCluster`**: Represents a cluster of interacting atoms.
- **`Cluster`**: Manages multiple interaction clusters based on system parameters.

# Functions
- `Cluster(system, symmetry, nbody::Int, cutoff_radii)`: Constructs a `Cluster` instance based on the given system and symmetry information.
- `set_mindist_pairs(nat::Int, xc_in::AbstractArray{<:AbstractFloat, 3}, exist_cell::AbstractVector{Bool})`: Computes minimum distance pairs between atoms.
- `set_interaction_by_cutoff(nat::Int, nat_prim::Int, kd_int_list::AbstractVector{Int}, map_p2s::AbstractMatrix{Int}, cutoff_radii::AbstractArray{Float64, 3}, nbody::Int, mindist_pairs::Matrix{<:Vector{<:DistInfo}})`: Determines interacting pairs based on cutoff radii.
- `set_interaction_clusters(nat_prim::Int, kd_int_list::AbstractVector{Int}, map_p2s::AbstractMatrix{Int}, x_image_cart::AbstractArray{Float64, 3}, exist_image::AbstractVector{Bool}, interaction_pairs::Vector{<:Vector{Int}}, cutoff_radii::AbstractArray{Float64, 3}, mindist_pairs::Matrix{<:Vector{<:DistInfo}}, nbody::Int)`: Organizes interaction clusters.
- `is_within_cutoff(atoms::AbstractVector{<:Int}, kd_int::AbstractVector{<:Int}, rc::AbstractArray{Float64, 3}, mindist_pairs::Matrix{<:AbstractVector{<:DistInfo}})`: Checks if a cluster is within cutoff radii.
- `calc_distmax(atoms::AbstractVector{<:Int}, mindist_pairs::Matrix{<:AbstractVector{<:DistInfo}})`: Calculates the maximum distance within a cluster.
- `generate_pairs(nat_prim::Int, map_p2s::Matrix{<:Int}, interaction_clusters::Matrix{OrderedSet{InteractionCluster}}, nbody::Int)`: Generates a list of interaction clusters.

"""

module Clusters

using Combinatorics
using DataStructures
using LinearAlgebra

using ..SortedContainer
using ..AtomCells
using ..Systems
using ..Symmetries

import Base: isless, ==

export Cluster

"""
	struct DistInfo <: MyAbstractDistInfo

Represents distance information between an atom in the central cell and an atom in a central or virtual neighboring cell.

# Fields
- `cell::Int`: The virtual cell index of the second atom.
- `distance::Float64`: The Cartesian distance between the two atoms.
- `relvec::Vector{Float64}`: The relative Cartesian vector from the first atom to the second atom.

# Constructor
	DistInfo(cell::Int, distance::Real, relvec::Vector{<:Real})

Creates a new `DistInfo` instance. Ensures that `relvec` has a length of 3.
"""
struct DistInfo
	cell::Int   # virtual cell index of atom j
	distance::Float64   # in Cartesian coordinates
	relvec::Vector{Float64} # relative vector in Cartesian coordinates.

	function DistInfo(cell::Integer, distance::Real, relvec::AbstractVector{<:Real})
		if length(relvec) != 3
			error("The length of \"relvec\" is not 3.")
		end
		return new(Int(cell), Float64(distance), Vector{Float64}(relvec))
	end
end

isless(distinfo1::DistInfo, distinfo2::DistInfo) = distinfo1.distance < distinfo2.distance

"""
Cluster of interactiong atoms. Note that this structure contains only partner atoms.
- `atom_list::Vector{Int}`: [≤ body-1]
- `cell_list::`: [≤ body-1][icell ≤ the number of equivalent atoms]
- `distmax::Float64`:
"""
struct InteractionCluster
	atoms::Vector{Int}
	cells::Vector{Vector{Int}}
	distmax::Float64
end

isless(intclus1::InteractionCluster, intclus2::InteractionCluster) =
	intclus1.atoms < intclus2.atoms

function ==(intclus1::InteractionCluster, intclus2::InteractionCluster)
	if intclus1.atoms != intclus2.atoms
		return false
	elseif intclus1.cells != intclus2.cells
		return false
	elseif intclus1.distmax != intclus2.distmax
		return false
	else
		return true
	end
end


"""
	struct Cluster

Represents a collection of interaction clusters based on the specified number of bodies and cutoff radii.

# Fields
- `nbody::Int`: The number of interacting bodies.
- `cutoff_radii::Array{Float64, 3}`: Cutoff radii for each atomic element pair and interaction body. Dimensions: [nkd, nkd, nbody].
- `mindist_pairs::Matrix{Vector{DistInfo}}`: Matrix containing the minimum distance pairs between atoms. Dimensions: [num_atoms, num_atoms].
- `interaction_clusters::Matrix{OrderedSet{InteractionCluster}}`: Matrix of interaction clusters for each primitive atom and interaction body. Dimensions: [nat_prim, nbody-1].
- `cluster_list::Vector{SortedVector{Vector{Int}}}`: List of interacting atom clusters for each interaction body.

# Constructor
	Cluster(system, symmetry, nbody::Int, cutoff_radii)

Creates a new `Cluster` instance based on the provided system, symmetry information, number of bodies, and cutoff radii.

# Example
```julia
cluster = Cluster(system, symmetry, 3, cutoff_radii)
"""
struct Cluster
	nbody::Int
	cutoff_radii::Array{Float64, 3}
	mindist_pairs::Matrix{Vector{DistInfo}} # [≤ num_atoms, ≤ num_atoms]
	interaction_clusters::Matrix{OrderedSet{InteractionCluster}}# [≤ nat_prim, ≤ nbody-1]
	cluster_list::Vector{SortedVector{Vector{Int}}}# [≤ nbody-1][≤ number of clusters][≤ nbody]
	cluster_list_with_cell::Vector{SortedVector{Vector{AtomCell}}} # [≤ nbody-1][≤ number of clusters][≤ nbody]
	equivalent_atom_list::Vector{Vector{Int}}

	function Cluster(
		system::System,
		symmetry::Symmetry,
		nbody::Integer,
		cutoff_radii::AbstractArray{<:Real},
	)
		mindist_pairs =
			set_mindist_pairs(
				system.supercell.num_atoms,
				system.x_image_cart,
				system.exist_image,
			)
		# interaction_pairs[≤ nat_prim][nbody-1]
		interaction_pairs::Vector{Vector{Int}} = set_interaction_by_cutoff(
			system.supercell.num_atoms,
			symmetry.nat_prim,
			system.supercell.kd_int_list,
			symmetry.map_p2s,
			cutoff_radii,
			nbody,
			mindist_pairs,
		)
		interaction_clusters = set_interaction_clusters(
			symmetry.nat_prim,
			system.supercell.kd_int_list,
			symmetry.map_p2s,
			system.x_image_cart,
			system.exist_image,
			interaction_pairs,
			cutoff_radii,
			mindist_pairs,
			nbody,
		)
		cluster_list = generate_pairs(
			symmetry.nat_prim,
			symmetry.map_p2s,
			interaction_clusters,
			nbody,
		)
		cluster_list_with_cell =
			generate_pairs_with_icells(symmetry.atoms_in_prim, interaction_clusters, nbody)

		equivalent_atom_list =
			classify_equivalent_atoms(symmetry.atoms_in_prim, symmetry.map_sym)

		return new(
			nbody,
			cutoff_radii,
			mindist_pairs,
			interaction_clusters,
			cluster_list,
			cluster_list_with_cell,
			equivalent_atom_list,
		)
	end
end

function set_mindist_pairs(
	nat::Integer,
	xc_in::AbstractArray{<:AbstractFloat, 3}, # x_image_cart in `System`
	exist_cell::AbstractVector{Bool},
)::Matrix{Vector{DistInfo}}

	distance_all = Matrix{Vector{DistInfo}}(undef, nat, nat)
	initialized = falses(size(distance_all))
	for iat in 1:nat
		for jat in 1:nat
			distinfo_list = Vector{DistInfo}()
			for icell in 1:27   # 27 is the number of virtual neighboring cells
				if !(exist_cell[icell])
					continue
				end

				dist_tmp::Float64 = norm(xc_in[:, iat, 1] - xc_in[:, jat, icell])
				relvec::Vector{Float64} = xc_in[:, jat, icell] - xc_in[:, iat, 1]
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

	mindist_pairs = Matrix{Vector{DistInfo}}(undef, nat, nat)
	initialized = falses(size(mindist_pairs))

	for iat in 1:nat
		for jat in 1:nat
			dist_vec_tmp = Vector{DistInfo}()
			dist_min = distance_all[iat, jat][1].distance
			for distinfo in distance_all[iat, jat]
				if distinfo.distance ≈ dist_min
					push!(dist_vec_tmp, distinfo)
				end
			end
			mindist_pairs[iat, jat] = dist_vec_tmp
			initialized[iat, jat] = true
		end
	end
	# check assigned or not
	undef_indices = findall(x -> x == false, initialized)
	if length(undef_indices) >= 1
		error("unassigned indices in mindist_pairs variable: $undef_indices")
	end

	return mindist_pairs
end

function set_interaction_by_cutoff(
	nat::Integer,
	nat_prim::Integer,
	kd_int_list::AbstractVector{<:Integer},
	map_p2s::AbstractMatrix{<:Integer},
	cutoff_radii::AbstractArray{<:Real, 3},
	nbody::Integer,
	mindist_pairs::Matrix{<:Vector{<:DistInfo}},
)::Vector{Vector{Int}}

	interaction_pairs = Vector{Vector{Int}}()
	for body in 2:nbody # 1 is skipped cuz it's 1-body term.
		for iat_prim in 1:nat_prim
			pair_tmp = Int[]
			iat = map_p2s[iat_prim, 1]  # index in the central cell
			ikd = kd_int_list[iat]

			for jat in 1:nat
				jkd = kd_int_list[jat]
				cutoff_tmp::Float64 = cutoff_radii[ikd, jkd, body]
				if iat == jat
					continue
				elseif cutoff_tmp < 0.0
					append!(pair_tmp, jat)
				elseif mindist_pairs[iat, jat][1].distance ≤ cutoff_tmp
					append!(pair_tmp, jat)
				end
			end
			push!(interaction_pairs, pair_tmp)
		end
	end
	return interaction_pairs
end

function set_interaction_clusters(
	nat_prim::Integer,
	kd_int_list::AbstractVector{<:Integer},
	map_p2s::AbstractMatrix{<:Integer},
	x_image_cart::AbstractArray{<:Real, 3},
	exist_image::AbstractVector{Bool},
	interaction_pairs::AbstractVector{<:AbstractVector{<:Integer}},
	cutoff_radii::AbstractArray{<:Real, 3},
	mindist_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
	nbody::Integer,
)::Matrix{OrderedSet{InteractionCluster}}
	interaction_clusters =
		Matrix{OrderedSet{InteractionCluster}}(undef, nat_prim, nbody - 1)
	initialized = falses(size(interaction_clusters))
	for body in 2:nbody
		for iat_prim in 1:nat_prim
			inter_cluster_set = OrderedSet{InteractionCluster}()
			iat = map_p2s[iat_prim, 1]
			intlist = interaction_pairs[iat_prim, body-1]
			sort!(intlist)

			if body == 2# 2-body terms
				for jat in intlist
					atoms_tmp = [jat]
					cell_tmp = Int[]
					for distinfo::DistInfo in mindist_pairs[iat, jat]
						append!(cell_tmp, distinfo.cell)
					end
					dist_max = mindist_pairs[iat, jat][1].distance
					push!(
						inter_cluster_set,
						InteractionCluster(atoms_tmp, [cell_tmp], dist_max),
					)
				end
			else# 3- or more body terms
				combinations = collect(combinations(intlist), body - 1)
				for combi in combinations
					atoms_tmp = combi
					cell_tmp = Vector{Vector{Int}}()
					if !is_within_cutoff(
						vcat([iat], combi),
						kd_int_list,
						cutoff_radii,
						mindist_pairs,
					)
						continue
					end
					dist_max = calc_distmax(vcat([iat], combi), mindist_pairs)
					for jat in combi
						cell_j_tmp = []
						for distinfo in mindist_pairs[iat, jat]
							push!(cell_j_tmp, distinfo.cell)
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

	undef_indices = findall(x -> x == false, initialized)
	if length(undef_indices) >= 1
		error("undef is detected in `interaction_clusters` variable: $undef_indices")
	end
	return interaction_clusters
end

"""
For combinations of three or more atoms, determine whether all of them are within the cutoff radius.
"""
function is_within_cutoff(
	atoms::AbstractVector{<:Integer},
	kd_int::AbstractVector{<:Integer},
	rc::AbstractArray{<:Real, 3},
	mindist_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Bool
	combinations = collect(combinations(atoms, 2))
	for pair in combinations
		if mindist_pairs[pair[1], pair[2]].distance ≥ rc[kd_int[pair[1]], kd_int[pair[2]]]
			return false
		end
	end
	return true
end

function calc_distmax(
	atoms::AbstractVector{<:Integer},
	mindist_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Float64
	combinations = collect(combinations(atoms, 2))
	dist_tmp = -1.0
	for pair in combinations
		if dist_tmp < mindist_pairs[pair[1], pair[2]].distance
			dist_tmp = mindist_pairs[pair[1], pair[2]].distance
		end
	end
	return dist_tmp
end

function generate_pairs(
	nat_prim::Integer,
	map_p2s::AbstractMatrix{<:Integer},
	interaction_clusters::AbstractMatrix{<:AbstractSet{InteractionCluster}},# [≤ nat_prim, ≤ nbody-1]
	nbody::Integer,
)::Vector{SortedVector{Vector{Int}}}
	cluster_list = Vector{SortedVector{Vector{Int}}}()
	for body in 2:nbody
		ordered_set = SortedVector{Vector{Int}}()
		for iat_prim in 1:nat_prim
			iat = map_p2s[iat_prim][1]
			for inter_clus::InteractionCluster in interaction_clusters[iat_prim, body-1]
				partners_tmp = sort(inter_clus.atoms)
				vcated = vcat(iat, partners_tmp)

				push!(ordered_set, vcated)
			end
		end
		push!(cluster_list, ordered_set)
	end
	return cluster_list
end

function generate_pairs_with_icells(
	atoms_in_prim::AbstractVector{<:Integer},
	interactoin_clusters::AbstractMatrix{<:AbstractSet{InteractionCluster}},# [≤ nat_prim, ≤ nbody-1]
	nbody::Integer,
)::Vector{SortedVector{Vector{AtomCell}}}
	cluster_list_with_cell = Vector{SortedVector{Vector{AtomCell}}}()
	vsv = Vector{SortedVector{Vector{AtomCell}}}()
	for body in 2:nbody
		sv = SortedVector{Vector{AtomCell}}()
		for (i, iat) in enumerate(atoms_in_prim)
			for intclus_i::InteractionCluster in interactoin_clusters[i, body-1]
				pairs = generate_combinations(intclus_i.atoms, intclus_i.cells)
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
classify the atoms in the primitive cell into equivalent group.
"""
function classify_equivalent_atoms(
	atoms_in_prim::AbstractVector{<:Integer},
	map_sym::AbstractMatrix{<:Integer},
)::Vector{Vector{Int}}
	group = Vector{Vector{Int}}()
	checked = falses(size(atoms_in_prim))

	nsym = size(map_sym, 2)
	for (idx, atom) in enumerate(atoms_in_prim)
		if checked[idx]
			continue
		end

		eqlist = Int[atom]

		for isym in 1:nsym
			eqatom = map_sym[atom, isym]
			if atom == eqatom
				continue
			elseif eqatom in atoms_in_prim && !(eqatom in eqlist)
				eq_idx = findfirst(x -> x == eqatom, atoms_in_prim)
				checked[eq_idx] = true
				push!(eqlist, eqatom)
			end
		end
		push!(group, eqlist)
	end

	return group
end

function classify_equivalent_clusters(
	interactoin_clusters::AbstractMatrix{<:AbstractSet{InteractionCluster}},
	kd_int_list::AbstractVector{<:Integer},
	atoms_in_prim::AbstractVector{<:Integer},
	map_sym::AbstractMatrix{<:Integer},
	lmax::AbstractMatrix{<:Integer},
)

end

"""
	generate_combinations(vec::AbstractVector{<:Any}, vecofvec::AbstractVector{<:AbstractVector{<:Any}}) -> Vector{Tuple}

Generates all possible combinations by pairing each element of `vec` with each element of the corresponding sublist in `vecofvec`.

# Arguments

- `vec::AbstractVector{<:Any}`: A vector containing elements to be paired.
- `vecofvec::AbstractVector{<:AbstractVector{<:Any}}`: A vector of vectors, where each subvector contains elements to be paired with the corresponding element in `vec`.

# Returns

- `Vector{Tuple}`: A vector of tuples, where each tuple contains paired elements from `vec` and `vecofvec`. Each tuple consists of pairs corresponding to the elements in `vec` and each subvector in `vecofvec`.

# Errors

- Throws an `ArgumentError` if the lengths of `vec` and `vecofvec` are not equal.

# Examples

```julia
# Define the vectors
A = [1, 2]
B = [['a', 'b', 'c'], ['d', 'e']]

# Generate combinations
result = generate_combinations(A, B)

# Display the result
println(result)
# Output:
# [((1, 'a'), (2, 'd')),
#  ((1, 'a'), (2, 'e')),
#  ((1, 'b'), (2, 'd')),
#  ((1, 'b'), (2, 'e')),
#  ((1, 'c'), (2, 'd')),
#  ((1, 'c'), (2, 'e'))]
"""
function generate_combinations(
	vec::AbstractVector,
	vecofvec::AbstractVector{<:AbstractVector},
)::Vector{Tuple}
	if length(vec) != length(vecofvec)
		error("Vectors A and B must have the same length.")
	end
	combinations = Iterators.product(vecofvec...)
	return [Tuple(zip(vec, comb)) for comb in combinations]
end

function all_atomlist_by_symop(
	atomlist::AbstractVector{<:Integer},
	map_sym::AbstractMatrix{<:Integer},
)::Vector{Vector{Int}}
	# total number of symmetry operations
	nsym = size(map_sym, 2)
	atomlist_list = Vector{Vector{Int}}()

	for i in 1:nsym
		atomlist_tmp = Int[]
		for atom in atomlist
			push!(atomlist_tmp, map_sym[atom, i])
		end
		push!(atomlist_list, atomlist_tmp)
	end

	return unique(atomlist_list)

end

end
