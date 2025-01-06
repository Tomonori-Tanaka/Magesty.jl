"""
This module provides functionality for storing and handling of basis set.
"""
module BasisSets

using Combinatorics

using ..SortedContainer
using ..AtomCells
using ..AtomicIndices
using ..Systems
using ..Symmetries
using ..Clusters

export BasisSet

"""
 - 'salc_coeffs:', coefficient vectors to represent symmetry-adapted linear combinations.
"""
struct BasisSet
	basislist::SortedCountingVector{IndicesUniqueList}
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

	# to pass compile
	tmp = [[]]

	return BasisSet(basislist, tmp)

end

function construct_basislist(
	kd_int_list::AbstractVector{<:Integer},
	atoms_in_prim::AbstractVector{<:Integer},
	cluster_list::AbstractVector{<:AbstractVector{<:AbstractVector{AtomCell}}}, # Vector{SortedVector{Vector{AtomCell}}}  2025-01-06
	lmax_mat::AbstractMatrix{<:Integer},
	bodymax::Integer,
)::SortedCountingVector{IndicesUniqueList}
	basislist = SortedCountingVector{IndicesUniqueList}()

	# firstly treat 1-body case which needs special treatments.
	for iat in atoms_in_prim
		lmax = lmax_mat[kd_int_list[iat], 1]
		if lmax == 0
			continue
		end
		iul::IndicesUniqueList = AtomicIndices.indices_singleatom(iat, 1, lmax) # where 1 means virtual cell index (i.e. centering cell)
		append!(basislist, iul)
	end

	for body in 2:bodymax

	end

	return basislist
end


end
