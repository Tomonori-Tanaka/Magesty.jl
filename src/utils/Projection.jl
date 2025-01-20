using SparseArrays

using ..RotationMatrices


function construct_projectionmatrix(
	basisdict::AbstractDict{<:Integer, <:AbstractVector},
	system::System,
	symmetry::Symmetry,
)::Tuple{Dict{Int, Matrix{Float64}}, Dict{Int, Any}}

	# aliases
	num_atoms::Int = system.supercell.num_atoms# total number of atoms in the supercell
	symdata::Vector{SymmetryOperation} = symmetry.symdata
	map_sym::Matrix{Int} = symmetry.map_sym
	map_s2p::Vector{Maps} = symmetry.map_s2p
	atoms_in_prim::Vector{Int} = symmetry.atoms_in_prim# atoms in the primitive cell
	symnum_translation::Vector{Int} = symmetry.symnum_translation

	dict = Dict{Int, Matrix{Float64}}()
	dict_each_matrix = Dict{Int, Vector{SparseMatrixCSC{Float64, Int}}}()
	for (idx::Int, basislist::SortedCountingUniqueVector) in basisdict
		println("calculating projection matrix of $idx-th basis list.")
		dim = length(basislist)
		projection_mat = spzeros(Float64, dim, dim)
		dict_each_matrix[idx] = []

		nneq = 0

		for (n, symop) in enumerate(symdata)
			projection_mat_per_symop =
				calc_projection(
					basislist,
					num_atoms,
					symop,
					n,
					map_sym,
					map_s2p,
					atoms_in_prim,
					symnum_translation,
					threshold_digits = 10,
					time_reversal_sym = false,
				)
			push!(dict_each_matrix[idx], projection_mat_per_symop)
			if nnz(projection_mat_per_symop) != 0
				nneq += 1
				projection_mat += projection_mat_per_symop
			end
		end

		# time_reversal_sym will be optional keyword
		time_reversal_sym = true
		if time_reversal_sym
			for (n, symop) in enumerate(symdata)
				projection_mat_per_symop =
					calc_projection(
						basislist,
						num_atoms,
						symop,
						n,
						map_sym,
						map_s2p,
						atoms_in_prim,
						symnum_translation,
						threshold_digits = 10,
						time_reversal_sym = time_reversal_sym,
					)
				if nnz(projection_mat_per_symop) != 0
					nneq += 1
					projection_mat += projection_mat_per_symop
				end
			end
		end

		projection_mat = projection_mat / nneq
		dict[idx] = projection_mat
	end

	return dict, dict_each_matrix
end

function calc_projection(
	basislist::AbstractVector{IndicesUniqueList},
	num_atoms::Integer,
	symop::SymmetryOperation,
	isym::Integer,
	map_sym::AbstractMatrix{<:Integer},
	map_s2p,
	atoms_in_prim,
	symnum_translation,
	;
	threshold_digits::Integer = 10,
	time_reversal_sym::Bool = false,
)::SparseMatrixCSC{Float64, Int}

	projection_matrix = spzeros(Float64, length(basislist), length(basislist))
	if symop.is_translation_included == true
		return projection_matrix
	end
	for (ir, rbasis::IndicesUniqueList) in enumerate(basislist)  # right-hand basis
		moved_atomlist, llist =
			move_atoms(rbasis, isym, map_sym)
		moved_atomlist = translate_atomlist2primitive(
			moved_atomlist,
			num_atoms,
			map_sym,
			map_s2p,
			atoms_in_prim,
			symnum_translation,
		)
		# moved_rbasis will be used later to determine a matrix element
		moved_rbasis = IndicesUniqueList()
		for (idx, atom) in enumerate(moved_atomlist)
			indices = Indices(atom, rbasis[idx].l, rbasis[idx].m)
			push!(moved_rbasis, indices)
		end

		partial_moved_basis::Vector{IndicesUniqueList} =
			AtomicIndices.product_indices_fixed_l(
				get_atomlist(moved_rbasis),
				get_llist(moved_rbasis),
			)
		partial_r_idx = findfirst(x -> equivalent(x, moved_rbasis), partial_moved_basis)
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
			multiplier = (-1)^(get_totalL(moved_rbasis))
		end

		if time_reversal_sym
			multiplier *= (-1)^(get_totalL(moved_rbasis))
		end

		euler_angles::Tuple{Float64, Float64, Float64} = rotmat2euler(rotmat)
		rotation_list = [Δl(indices.l, euler_angles...) for indices in moved_rbasis]

		if length(moved_rbasis) == 1    # 1-body term
			rotmat_kron = multiplier * rotation_list[begin]
		else
			rotmat_kron = multiplier * kron(rotation_list...)
		end


		for (il, lbasis::IndicesUniqueList) in enumerate(basislist)  # left-hand basis
			partial_l_idx = findfirst(x -> equivalent(x, lbasis), partial_moved_basis)
			if isnothing(partial_l_idx)
				continue
			end
			projection_matrix[il, ir] = rotmat_kron[partial_l_idx, partial_r_idx]
		end
	end

	if !(RotationMatrices.is_orthogonal(projection_matrix, tol = 1e-6))
		println("symmetry operation index: $isym")
		println(symop)
		display(projection_matrix)
		error("not orthogonal")
	end
	projection_matrix = round.(projection_matrix, digits = threshold_digits)


	return projection_matrix
end

function move_atoms(
	iul::IndicesUniqueList,
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

"""
	translate_atomlist2primitive

This function is designed to transform the given atoms_list using translational operations
	so that the first atom in the list is positioned within the primitive cell.
"""
function translate_atomlist2primitive(
	atom_list::AbstractVector{<:Integer},
	num_atoms::Integer,
	map_sym::AbstractMatrix{<:Integer},
	map_s2p::AbstractVector{Symmetries.Maps},
	atoms_in_prim::AbstractVector{<:Integer},
	symnum_translation::AbstractVector{<:Integer},
)::Vector{Int}
	header_atom = first(atom_list)
	header_atom_in_prim = atoms_in_prim[map_s2p[header_atom].atom]

	# identify corresponding translational operation
	trans_op_idx = 0
	for idx_trans in symnum_translation
		if map_sym[header_atom_in_prim, idx_trans] == header_atom
			trans_op_idx = idx_trans
			break
		end
	end
	if trans_op_idx == 0
		error("Something is wrong.")
	end

	moved_atomlist = Int[]
	for atom in atom_list
		for iat in 1:num_atoms
			if map_sym[iat, trans_op_idx] == atom
				push!(moved_atomlist, iat)
			end
		end
	end

	if length(moved_atomlist) != length(atom_list)
		error("Something is wrong.")
	end

	return moved_atomlist
end
