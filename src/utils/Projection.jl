using SparseArrays

using ..RotationMatrices


function construct_projectionmatrix(basislist::AbstractVector{IndicesSortedUniqueList},
	symdata::AbstractVector{SymmetryOperation},
	map_sym::AbstractMatrix{<:Integer},
)::Tuple{SparseMatrixCSC, Vector{SparseMatrixCSC}}
	dimension = length(basislist)

	projection_mat = spzeros(Float64, dimension, dimension)
	projection_mat_list = Vector{SparseMatrixCSC}()

	nneq = 0

	projection_mat_per_symop = spzeros(Float64, dimension, dimension)

	for (n, symop) in enumerate(symdata)
		projection_mat_per_symop =
			calc_projection(
				basislist,
				symop,
				n,
				map_sym;
				threshold_digits = 10,
				time_reversal_sym = false,
			)
		if nnz(projection_mat_per_symop) != 0
			nneq += 1
			projection_mat += projection_mat_per_symop
		end
		push!(projection_mat_list, projection_mat_per_symop)
	end

	# time_reversal_sym will be optional keyword
	time_reversal_sym = true
	if time_reversal_sym
		for (n, symop) in enumerate(symdata)
			projection_mat_per_symop =
				calc_projection(
					basislist,
					symop,
					n,
					map_sym;
					threshold_digits = 10,
					time_reversal_sym = time_reversal_sym,
				)
			if nnz(projection_mat_per_symop) != 0
				nneq += 1
				projection_mat += projection_mat_per_symop
			end
			push!(projection_mat_list, projection_mat_per_symop)
		end
	end


	projection_mat = projection_mat / nneq
	# projection_mat_list = projection_mat_list ./ nneq

	return projection_mat, projection_mat_list
end

function calc_projection(
	basislist::AbstractVector{IndicesSortedUniqueList},
	symop::SymmetryOperation,
	isym::Integer,
	map_sym::AbstractMatrix{<:Integer};
	threshold_digits::Integer = 10,
	time_reversal_sym::Bool = false,
)::SparseMatrixCSC{Float64, Int}

	projection_matrix = spzeros(Float64, length(basislist), length(basislist))
	for (ir, rbasis::IndicesSortedUniqueList) in enumerate(basislist)  # right-hand basis
		moved_atomlist, llist =
			move_atoms(rbasis, isym, map_sym)
		# moved_rbasis will be used later to determine a matrix element
		moved_rbasis = IndicesSortedUniqueList()
		for (idx, atom) in enumerate(moved_atomlist)
			indices = Indices(atom, rbasis[idx].l, rbasis[idx].m)
			push!(moved_rbasis, indices)
		end

		partial_moved_basis::Vector{IndicesSortedUniqueList} =
			AtomicIndices.product_indices(moved_atomlist, llist)
		partial_r_idx = findfirst(x -> x == moved_rbasis, partial_moved_basis)
		if isnothing(partial_r_idx)
			error("Something is wrong at partial_r_idx variable.")
		end

		# calculate rotation matrix
		is_proper::Bool = symop.is_proper
		multiplier::Float64 = 1.0
		if is_proper
			rotmat = symop.rotation_cart
		else
			rotmat = -1 * symop.rotation_cart
			multiplier = (-1)^(gettotall(rbasis))
		end

		if time_reversal_sym
			multiplier *= (-1)^(gettotall(rbasis))
		end

		euler_angles::Tuple{Float64, Float64, Float64} = rotmat2euler(rotmat)
		rotation_list = Vector{Matrix{Float64}}()
		for indices in moved_rbasis
			push!(rotation_list, Δl(indices.l, euler_angles...))
		end

		if length(moved_rbasis) == 1    # 1-body term
			rotmat_kron = multiplier * rotation_list[begin]
		else
			rotmat_kron = multiplier * kron(rotation_list...)
		end


		for (il, lbasis::IndicesSortedUniqueList) in enumerate(basislist)  # left-hand basis
			partial_l_idx = findfirst(x -> x == lbasis, partial_moved_basis)
			if isnothing(partial_l_idx)
				continue
			end

			projection_matrix[il, ir] = rotmat_kron[partial_l_idx, partial_r_idx]
		end
	end

	projection_matrix = round.(projection_matrix, digits = threshold_digits)

	return projection_matrix
end

function move_atoms(
	iul::IndicesSortedUniqueList,
	isym::Integer,
	map_sym::AbstractMatrix{<:Integer},   # [≤ num_atoms, ≤ 27, ≤ nsym]
)::Tuple{Vector{Int}, Vector{Int}}
	moved_atomlist = Int[]
	moved_llist = Int[]
	for indices::Indices in iul
		moved_atom = map_sym[indices.atom, isym]
		push!(moved_atomlist, moved_atom)
		push!(moved_llist, indices.l)
	end
	return moved_atomlist, moved_llist
end
