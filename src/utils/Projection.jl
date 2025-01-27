using SparseArrays

using ..RotationMatrices


function construct_projectionmatrix(
	basisdict::AbstractDict{<:Integer, <:AbstractVector},
	system::System,
	symmetry::Symmetry,
)::Tuple{Dict{Int, Matrix{Float64}}, Dict{Int, Any}}

	# aliases
	num_atoms::Int = system.supercell.num_atoms# total number of atoms in the supercell
	x_image_frac::Array{Float64, 3} = system.x_image_frac
	symdata::Vector{SymmetryOperation} = symmetry.symdata
	map_sym::Matrix{Int} = symmetry.map_sym
	map_sym_cell::Array{AtomCell} = symmetry.map_sym_cell
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

		for time_rev_sym in [false, true]
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
						map_sym_cell,
						x_image_frac,
						threshold_digits = 10,
						time_reversal_sym = time_rev_sym,
					)
				push!(dict_each_matrix[idx], projection_mat_per_symop)
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
	map_sym_cell::AbstractArray{AtomCell},
	x_image_frac
	;
	threshold_digits::Integer = 10,
	time_reversal_sym::Bool = false,
)::SparseMatrixCSC{Float64, Int}

	projection_matrix = spzeros(Float64, length(basislist), length(basislist))
	if symop.is_translation_included == true
		return projection_matrix
	end
	for (ir, rbasis::IndicesUniqueList) in enumerate(basislist)  # right-hand basis

		moved_atomlist, moved_celllist =
			apply_symop_to_basis(
				rbasis,
				isym,
				symop,
				map_sym_cell,
				map_s2p,
				atoms_in_prim,
				x_image_frac,
			)
		# moved_rbasis will be used later to determine a matrix element
		moved_rbasis = IndicesUniqueList()
		for (idx, (atom, cell)) in enumerate(zip(moved_atomlist, moved_celllist))
			indices = Indices(atom, rbasis[idx].l, rbasis[idx].m, cell)
			push!(moved_rbasis, indices)
		end

		partial_moved_basis::Vector{IndicesUniqueList} =
			AtomicIndices.product_indices_fixed_l(
				get_atomlist(moved_rbasis),
				get_llist(moved_rbasis),
				get_celllist(moved_rbasis),
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

"""
Apply the symmetry operation and return the list of correspoinding atom and cell list.
Note that the first atom is keeped to be located at the primitive cell. 
"""
function apply_symop_to_basis(
	basis::IndicesUniqueList,
	isym::Integer,
	symop::SymmetryOperation,
	map_sym_cell::AbstractArray{AtomCell},
	map_s2p::AbstractVector,
	atoms_in_prim::AbstractVector{<:Integer},
	x_image_frac::AbstractArray,
)::NTuple{2, Vector{Int}}

	atom_list = get_atomlist(basis)
	cell_list = get_celllist(basis)

	header_atom = atom_list[begin]
	header_cell = cell_list[begin]

	# firstly, apply the symmetry operation to atoms
	moved_coords = Vector{Vector{Float64}}()
	for (atom, cell) in zip(atom_list, cell_list)
		coords_tmp =
			symop.rotation_frac * x_image_frac[:, atom, cell] + symop.translation_frac
		push!(moved_coords, coords_tmp)
	end

	# secondly, calculate translation vector to shift first atom to the primitive cell
	moved_header_indices::NTuple{2, Int} =
		(
			map_sym_cell[header_atom, header_cell, isym].atom,
			map_sym_cell[header_atom, header_cell, isym].cell,
		)
	moved_header_prim::NTuple{2, Int} =
		(atoms_in_prim[map_s2p[moved_header_indices[1]].atom], 1)

	if map_sym_cell[header_atom, header_cell, isym].cell < 0
		@show header_atom
		@show header_cell
		@show isym
		@show moved_header_indices
	end
	translation_vec =
		calc_relvec_in_frac(moved_header_indices, moved_header_prim, x_image_frac)

	# thirdly, shift all atoms by using the translation vector
	shifted_coords = Vector{Vector{Float64}}()
	for coords in moved_coords
		push!(shifted_coords, coords + translation_vec)
	end

	# finally, find corresponding atom indices
	result_atom_list = Int[]
	result_cell_list = Int[]
	for coords in shifted_coords
		for iatom in axes(x_image_frac, 2)
			for icell in axes(x_image_frac, 3)
				if isapprox(coords, x_image_frac[:, iatom, icell], atol = 1e-6)
					push!(result_atom_list, iatom)
					push!(result_cell_list, icell)
					@goto found
				end
			end
		end

		# hints of error
		@show atom_list
		@show moved_coords
		@show translation_vec
		@show result_atom_list
		error("Something is wrong.")
		@label found
	end

	return result_atom_list, result_cell_list
end

function calc_relvec_in_frac(atom1::NTuple{2, Integer},# (atom, cell)
	atom2::NTuple{2, Integer},
	x_image_frac::AbstractArray{<:Real, 3},
)::Vector{Float64}
	result::Vector{Float64} =
		x_image_frac[:, atom2[1], atom2[2]] - x_image_frac[:, atom1[1], atom1[2]]

	return result
end

