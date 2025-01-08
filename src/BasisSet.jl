"""
This module provides functionality for storing and handling of basis set.
"""
module BasisSets

using Arpack
using Combinatorics
using LinearAlgebra

using ..SortedContainer
using ..AtomCells
using ..AtomicIndices
using ..Systems
using ..Symmetries
using ..Clusters

include("./utils/Projection.jl")

export BasisSet

"""
 - 'salc_coeffs:', coefficient vectors to represent symmetry-adapted linear combinations.
"""
struct BasisSet
	basislist::SortedCountingUniqueVector{IndicesUniqueList}
	salc_coeffs::Vector{Vector{Float64}}
end

function BasisSet(
	system::System,
	symmetry::Symmetry,
	cluster::Cluster,
	lmax::AbstractMatrix{<:Integer},    # [≤ nkd, ≤ nbody]
	bodymax::Integer)

	basislist = construct_basislist(
		system.supercell.kd_int_list,
		symmetry.atoms_in_prim,
		cluster.cluster_list_with_cell,
		lmax,
		bodymax,
	)

	projection_matrix, each_matrix_list =
		construct_projectionmatrix(
			basislist,
			symmetry.symdata,
			symmetry.map_sym,
		)
	projection_matrix = Matrix(projection_matrix)
	eigenval, eigenvec = eigen(projection_matrix)
	# eigenval, eigenvec = eigen(Matrix(projection_matrix))
	eigenvec = round.(eigenvec, digits = 6)
	eigenval = round.(eigenval, digits = 6)
	eigenval = real.(eigenval)
	eigenvec = real.(eigenvec)
	@show eigenval
	# for (val, basis) in zip(eigenvec[:, end], basislist)
	# 	println(val, "\t", basis)
	# end
	idx = 2
	display(symmetry.symdata[idx].rotation_frac)
	display(symmetry.symdata[idx].translation_frac)
	display(each_matrix_list[idx])

	tmp = [[]]

	return BasisSet(basislist, tmp)

end

function construct_basislist(
	kd_int_list::AbstractVector{<:Integer},
	atoms_in_prim::AbstractVector{<:Integer},
	cluster_list::AbstractVector{<:AbstractVector{<:AbstractVector{AtomCell}}}, # Vector{SortedVector{Vector{AtomCell}}}  2025-01-06
	lmax_mat::AbstractMatrix{<:Integer},
	bodymax::Integer,
)::SortedCountingUniqueVector{IndicesUniqueList}

	basislist = SortedCountingUniqueVector{IndicesUniqueList}()

	# firstly treat 1-body case which needs special treatments.
	for iat in atoms_in_prim
		lmax = lmax_mat[kd_int_list[iat], 1]
		if lmax == 0
			continue
		end
		iul::IndicesUniqueList = AtomicIndices.indices_singleatom(iat, lmax)
		for indices::Indices in iul
			push!(basislist, IndicesUniqueList(indices))
		end
	end

	# multi-body cases
	for body in 2:bodymax
		for cluster in cluster_list[body-1]
			atomlist, llist =
				get_atomsls_from_cluster(cluster, lmax_mat, kd_int_list)
			for iul in AtomicIndices.product_indices(atomlist, llist)
				push!(basislist, iul)
			end
		end
	end

	return basislist
end

function get_atomsls_from_cluster(
	cluster::AbstractVector{AtomCell}, # [≤ nbody]
	lmax::AbstractMatrix{<:Integer},
	kd_int_list::AbstractVector{<:Integer},
)::Tuple{Vector{Int}, Vector{Int}}

	atomlist = [atomcell.atom for atomcell in cluster]

	body = length(cluster)
	llist = [lmax[kd_int_list[atomcell.atom], body] for atomcell in cluster]

	return atomlist, llist
end

end
