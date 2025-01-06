using SparseArrays

function construct_projectionmatrix(basislist::AbstractVector{IndicesUniqueList},
	symdata::AbstractVector{SymmetryOperation},
	map_sym::AbstractArray{AtomCell},
)
	dimension = length(basislist)

	projection_mat = spzeros(Float64, dimension, dimension)
	projection_mat_list = Vector{SparseMatrixCSC}()

	nneq = 0

	projection_mat_per_symop = spzeros(Float64, dimension, dimension)

	for (n, symop) in enumerate(symdata)
		projection_mat_per_symop =
			calc_projection(basislist, symop, n, map_sym; threshold = 1e-12)
		if nnz(projection_mat_per_symop) != 0
			nneq += 1
			projection_mat += projection_mat_per_symop
		end
		push!(projection_mat_list, projection_mat_per_symop)

	end

	time_reversal_sym = true
	if time_reversal_sym
		projection_mat_per_symop .= 0
		nneq += 1
		for (i, basis) in enumerate(basislist)
			totall = gettotall(basis)
			factor = (-1)^totall
			projection_mat_per_symop[i, i] = factor
		end
		push!(projection_mat_list, projection_mat_per_symop)
		projection_mat += projection_mat_per_symop
	end

end

function calc_projection(
	basislist::AbstractVector{IndicesUniqueList},
	symop::SymmetryOperation,
	isym::Integer,
	map_sym::AbstractArray{AtomCell};
	threshold = 1e-12,
)::SparseMatrixCSC{Float64, Int}

	projection_matrix = spzeros(Float64, length(basislist), length(basislist))
	for (ir, rbasis::IndicesUniqueList) in enumerate(basislist)  # right-hand basis
		moved_atomlist, moved_celllist, moved_llist =
			move_atomscellsls(rbasis, isym, map_sym)
		for (il, lbasis::IndicesUniqueList) in enumerate(basislist)  # left-hand basis

			if length(rbasis) != length(lbasis)
				continue
			end

		end
	end
	return projection_matrix
end

function move_atomscellsls(
	iul::IndicesUniqueList,
	isym::Integer,
	map_sym::AbstractArray{AtomCell},   # [≤ num_atoms, ≤ 27, ≤ nsym]
)::Tuple{Vector{Int}, Vector{Int}, Vector{Int}}
	moved_atomlist = Int[]
	moved_celllist = Int[]
	moved_llist = Int[]
	for indices::Indices in iul
		atomcell_tmp = map_sym[indices.atom, indices.cell, isym]
		push!(moved_atomlist, atomcell_tmp.atom)
		push!(moved_celllist, atomcell_tmp.cell)
		push!(moved_llist, indices.l)
	end
	return moved_atomlist, moved_celllist, moved_llist
end
