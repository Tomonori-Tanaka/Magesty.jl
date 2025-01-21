"""
This module provides functionality for storing and handling of basis set.
"""
module BasisSets

using Combinatorics
using DataStructures
using LinearAlgebra
using DelimitedFiles# for test

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
	projection_dict::Dict{Int, Matrix{Float64}}
	each_projection_dict::Any
	salc_coeffs::Vector{Vector{Float64}}
end

function BasisSet(
	system::System,
	symmetry::Symmetry,
	cluster::Cluster,
	lmax::AbstractMatrix{<:Integer},    # [≤ nkd, ≤ nbody]
	bodymax::Integer)

	basislist = construct_basislist(
		system,
		symmetry,
		cluster,
		lmax,
		bodymax,
	)

	classified_basisdict = classify_basislist(basislist, symmetry.map_sym)
	println(classified_basisdict)

	projection_dict::Dict{Int, Matrix{Float64}},
	each_projection_dict =
		construct_projectionmatrix(
			classified_basisdict,
			system,
			symmetry,
		)

	display(SparseMatrixCSC(projection_dict[1]))
	#= 	for (idx, mat) in enumerate(each_projection_dict[1])
			if idx > 0
				println(idx)
				display(mat)
			end
		end =#
	for idx in 1:maximum(keys(projection_dict))
		eigenval, eigenvec = eigen(projection_dict[idx])
		eigenval = round.(eigenval, digits = 6)
		eigenvec = round.(eigenvec, digits = 6)
		println(idx, "\t", eigenval)
		display(eigenvec[:, end-1])
		display(eigenvec[:, end])
	end

	# projection_matrix = Matrix(projection_matrix)
	# eigenval, eigenvec = eigen(projection_matrix)
	# eigenval, eigenvec = eigen(Matrix(projection_matrix))
	# eigenvec = round.(eigenvec, digits = 6)
	# eigenval = round.(eigenval, digits = 6)
	# eigenval = real.(eigenval)
	# eigenvec = real.(eigenvec)
	# @show eigenval
	# for (val, basis) in zip(eigenvec[:, end-7], basislist)
	# 	println(val, "\t", basis)
	# end

	tmp = [[]]

	return BasisSet(basislist, projection_dict, each_projection_dict, tmp)

end

function construct_basislist(
	system::System,
	symmetry::Symmetry,
	cluster::Cluster,
	lmax_mat::AbstractMatrix{<:Integer},
	bodymax::Integer,
)::SortedCountingUniqueVector{IndicesUniqueList}

	basislist = SortedCountingUniqueVector{IndicesUniqueList}()

	# aliases
	kd_int_list = system.supercell.kd_int_list
	cluster_list = cluster.cluster_list_with_cell

	# firstly treat 1-body case which needs special treatments.
	for iat in symmetry.atoms_in_prim
		lmax = lmax_mat[kd_int_list[iat], 1]

		for l in 1:lmax
			iul::Vector{Indices} = indices_singleatom(iat, l, 1)
			for indices::Indices in iul
				push!(basislist, IndicesUniqueList(indices))
			end
		end
	end

	# multi-body cases
	for body in 2:bodymax
		for cluster in cluster_list[body-1]
			atomlist, llist, celllist =
				get_atomsls_from_cluster(cluster, lmax_mat, kd_int_list)
			for iul in product_indices(atomlist, llist, celllist)
				for basis in basislist
					if equivalent(basis, iul)
						basislist.counts[basis] += 1
						@goto skip
					elseif is_translationally_equiv_basis(iul, basis, symmetry, system)
						basislist.counts[basis] += 1
						@goto skip
					end
				end
				push!(basislist, iul)
				@label skip
			end
		end
	end

	# basislist = merge_duplicated_elements(basislist, symmetry)

	return basislist
end

function get_atomsls_from_cluster(
	cluster::AbstractVector{AtomCell}, # [≤ nbody]
	lmax::AbstractMatrix{<:Integer},
	kd_int_list::AbstractVector{<:Integer},
)::Tuple{Vector{Int}, Vector{Int}, Vector{Int}}

	atomlist = [atomcell.atom for atomcell in cluster]
	celllist = [atomcell.cell for atomcell in cluster]

	body = length(cluster)
	llist = [lmax[kd_int_list[atomcell.atom], body] for atomcell in cluster]

	return atomlist, llist, celllist
end

function classify_basislist(
	basislist::AbstractVector{IndicesUniqueList},
	map_sym::AbstractMatrix{<:Integer},
)

	count = 1
	label_list = zeros(Int, size(basislist))
	for (idx, basis) in enumerate(basislist)
		if label_list[idx] != 0
			continue
		end

		atom_l_list_base = get_atom_l_list(basis)

		for isym in 1:size(map_sym, 2)
			mapped_list = map_atom_l_list(atom_l_list_base, map_sym, isym)
			for (idx2, basis2) in enumerate(basislist)
				if sort(mapped_list) == sort(get_atom_l_list(basis2))
					label_list[idx2] = count
				end
			end
		end
		count += 1
	end


	dict = OrderedDict{Int, SortedCountingUniqueVector}()

	for idx in 1:maximum(label_list)
		if !haskey(dict, idx)
			dict[idx] = SortedCountingUniqueVector{IndicesUniqueList}()
		end
	end

	for (basis, label) in zip(basislist, label_list)
		if !(haskey(dict, label))
			dict[label] = SortedCountingUniqueVector{IndicesUniqueList}()
		end
		push!(dict[label], basis, getcount(basislist, basis))
	end

	return dict
end

function map_atom_l_list(
	atom_l_list::AbstractVector{<:AbstractVector{<:Integer}},
	map_sym::AbstractMatrix{<:Integer},
	isym::Integer,
)::Vector{Vector{Int}}
	mapped_atom_l_list = Vector{Vector{Int}}()
	for atom_l_vec in atom_l_list
		mapped_atom = map_sym[atom_l_vec[1], isym]
		push!(mapped_atom_l_list, [mapped_atom, atom_l_vec[2]])
	end

	return mapped_atom_l_list
end

function is_translationally_equiv_basis(
	basis_target::IndicesUniqueList,
	basis_ref::IndicesUniqueList,
	symmetry::Symmetry,
	system::System,
)::Bool

	atomlist = get_atomlist(basis_target)
	celllist = get_celllist(basis_target)

	if atomlist == get_atomlist(basis_ref)
		return false
	end

	for i in 2:length(basis_target)
		iatom = get_atomlist(basis_target)[i]
		icell = get_celllist(basis_target)[i]
		iatom_in_prim = symmetry.atoms_in_prim[symmetry.map_s2p[iatom].atom]

		# cartesian relative vector b/w iatom and iatom_in_prim
		relvec::Vector{Float64} =
			calc_relvec_in_cart((iatom, icell), (iatom_in_prim, 1), system.x_image_cart)

		moved_atomlist = Int[]
		moved_celllist = Int[]
		for indices::Indices in basis_target
			# corresponding atom and cell obtained by adding relvec
			crrsp_atom, crrsp_cell = find_corresponding_atom(
				(indices.atom, indices.cell),
				relvec,
				system.x_image_cart,
			)
			push!(moved_atomlist, crrsp_atom)
			push!(moved_celllist, crrsp_cell)
		end

		iul = IndicesUniqueList()
		for (idx, (atom, cell)) in enumerate(zip(moved_atomlist, moved_celllist))
			push!(iul, Indices(atom, basis_target[idx].l, basis_target[idx].m, cell))
		end

		if equivalent(iul, basis_ref)
			return true
		end
	end
	return false
end

"""
Finds the cartesian relative vector b/w 2 atoms specified by (atom, cell) tuples, where cell means virtual cell index (1 <= cell <= 27).
The equation is
r(atom2) - r(atom1)
"""
function calc_relvec_in_cart(
	atom1::NTuple{2, Integer},# (atom, cell)
	atom2::NTuple{2, Integer},
	x_image_cart::AbstractArray{<:Real, 3},
)::Vector{Float64}
	relvec::Vector{Float64} =
		x_image_cart[:, atom2[1], atom2[2]] - x_image_cart[:, atom1[1], atom1[2]]
	return relvec
end

"""
Move an atom by a cartesian relative vector (calculated by `calc_relvec_in_cart` function) and find corresponding atom and cell indices as a tuple.
"""
function find_corresponding_atom(
	atom::NTuple{2, Int},# (atom, cell)
	relvec::AbstractVector{<:Real},
	x_image_cart::AbstractArray{<:Real, 3},
)::NTuple{2, Int}

	moved_coords::Vector{Float64} = x_image_cart[:, atom[1], atom[2]] + relvec

	num_atoms = size(x_image_cart, 2)
	num_cells = size(x_image_cart, 3)
	for iatom in 1:num_atoms
		for icell in 1:num_cells
			if isapprox(x_image_cart[:, iatom, icell], moved_coords, atol = 1e-8)
				return (iatom, icell)
			end
		end
	end
	error("No matching (atom, cell) indices found.")
end


function merge_duplicated_elements(
	basislist::SortedCountingUniqueVector,
	symmetry::Symmetry,
)
	basislist_copy = copy(basislist)
	duplication_list = Vector{IndicesUniqueList}()

	for (i, iul_outer) in enumerate(basislist)
		if iul_outer in duplication_list
			continue
		end

		for (j, iul_inner) in enumerate(basislist)
			if j ≤ i
				continue
			end

			for itrans in symmetry.symnum_translation[2:end]
				iul_mapped = IndicesUniqueList()
				for indices in iul_inner
					push!(
						iul_mapped,
						Indices(
							symmetry.map_sym[indices.atom, itrans],
							indices.l,
							indices.m,
						),
					)
				end
				if equivalent(iul_outer, iul_mapped)
					addcount!(basislist_copy, iul_outer, getcount(basislist, iul_inner))
					push!(duplication_list, iul_inner)
					break
				end
			end
		end
	end

	@show duplication_list

	for iul_delete in duplication_list
		delete!(basislist_copy, iul_delete)
	end

	return basislist_copy
end

function __write_martix(
	matrix::AbstractMatrix,
)::Nothing
	open("matrix.csv", "w") do io
		writedlm(io, Matrix(matrix), ',')
	end
end

end
