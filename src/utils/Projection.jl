using SparseArrays

using ..RotationMatrices

"""
Construct projection matrices for basis functions

Parameters
----------
basisdict : Dictionary of basis functions
system : System information
symmetry : Symmetry information

Returns
-------
Tuple{Dict{Int, Matrix{Float64}}, Dict{Int, Any}}
	Dictionary of projection matrices and dictionary of matrices for each symmetry operation
"""
function construct_projectionmatrix(
	basisdict::AbstractDict{<:Integer, <:AbstractVector},
	system::System,
	symmetry::Symmetry,
)::Tuple{Dict{Int, Matrix{Float64}}, Dict{Int, Any}}

	# aliases
	num_atoms::Int = system.supercell.num_atoms# total number of atoms in the supercell
	x_image_cart::Array{Float64, 3} = system.x_image_cart
	lattice_vectors::Array{Float64, 2} = system.supercell.lattice_vectors
	symdata::Vector{SymmetryOperation} = symmetry.symdata
	map_sym::Matrix{Int} = symmetry.map_sym
	map_sym_cell::Array{AtomCell} = symmetry.map_sym_cell
	map_s2p::Vector{Maps} = symmetry.map_s2p
	atoms_in_prim::Vector{Int} = symmetry.atoms_in_prim# atoms in the primitive cell
	symnum_translation::Vector{Int} = symmetry.symnum_translation

	result_projection = Dict{Int, Matrix{Float64}}()
	result_each_matrix = Dict{Int, Vector{SparseMatrixCSC{Float64, Int}}}()
	for (idx::Int, basislist::SortedCountingUniqueVector) in basisdict
		println("calculating projection matrix of $idx-th basis list.")
		dim = length(basislist)
		projection_mat = spzeros(Float64, dim, dim)
		result_each_matrix[idx] = []

		num_nonzero_projections = 0

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
						x_image_cart,
						lattice_vectors,
						threshold_digits = 10,
						time_reversal_sym = time_rev_sym,
					)
				push!(result_each_matrix[idx], projection_mat_per_symop)
				if nnz(projection_mat_per_symop) != 0
					num_nonzero_projections += 1
					projection_mat += projection_mat_per_symop
				end
			end
		end

		projection_mat = projection_mat / num_nonzero_projections
		result_projection[idx] = projection_mat
	end

	return result_projection, result_each_matrix
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
	x_image_cart,
	lattice_vectors,
	;
	threshold_digits::Integer = 6,
	time_reversal_sym::Bool = false,
)::SparseMatrixCSC{Float64, Int}

	projection_matrix = spzeros(Float64, length(basislist), length(basislist))
	if symop.is_translation_included == true
		return projection_matrix
	end
	for (ir, rbasis::IndicesUniqueList) in enumerate(basislist)  # right-hand basis
	
		moved_atomlist, moved_celllist =
			apply_symop_to_basis_with_shift(
				rbasis,
				isym,
				symop,
				map_sym_cell,
				map_s2p,
				atoms_in_prim,
				x_image_cart,
				lattice_vectors,
			)
		# moved_rbasis will be used later to determine a matrix element
		moved_rbasis = IndicesUniqueList()
		for (idx, (atom, cell)) in enumerate(zip(moved_atomlist, moved_celllist))
			indices = Indices(atom, rbasis[idx].l, rbasis[idx].m, cell)
			push!(moved_rbasis, indices)
		end

		found = false
		for basis in basislist
			if equivalent(basis, moved_rbasis)
				moved_rbasis = basis
				found = true
				break
			elseif is_translationally_equiv_basis(moved_rbasis, basis, atoms_in_prim, map_s2p, x_image_cart)
				moved_rbasis = basis
				found = true
				break
			end
		end
		if !found
			@show isym
			@show symop
			@show rbasis
			@show moved_rbasis
			error("Failed to find moved basis in basis projection.")
		end

		partial_moved_basis::Vector{IndicesUniqueList} =
			AtomicIndices.product_indices_fixed_l(
				get_atomlist(moved_rbasis),
				get_llist(moved_rbasis),
				get_celllist(moved_rbasis),
			)
		partial_r_idx = findfirst(x -> equivalent(x, moved_rbasis), partial_moved_basis)
		if isnothing(partial_r_idx)
			error("Failed to find partial basis index in basis projection.")
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
			# @show lbasis
			# @show moved_rbasis
			partial_l_idx = findfirst(x -> equivalent(x, lbasis), partial_moved_basis)
			if isnothing(partial_l_idx)
				continue
			end
			projection_matrix[il, ir] = rotmat_kron[partial_l_idx, partial_r_idx]
		end
	end

	projection_matrix = round.(projection_matrix, digits = threshold_digits)
	if !(RotationMatrices.is_orthogonal(projection_matrix, tol = 1e-7))
		println("symmetry operation index: $isym")
		println(symop)
		display(projection_matrix)
		# display(Matrix(projection_matrix))
		error("not orthogonal")
	end

	return projection_matrix
end

function apply_symop_to_basis(
	basis::IndicesUniqueList,
	isym::Integer,
	map_sym_cell::AbstractArray{AtomCell},
)::NTuple{2, Vector{Int}}
	atom_list = get_atomlist(basis)
	cell_list = get_celllist(basis)

	moved_atom_list = Vector{Int}()
	moved_cell_list = Vector{Int}()

	for (atom, cell) in zip(atom_list, cell_list)
		push!(moved_atom_list, map_sym_cell[atom, cell, isym].atom)
		push!(moved_cell_list, map_sym_cell[atom, cell, isym].cell)
	end

	return moved_atom_list, moved_cell_list

end

"""
Apply the symmetry operation and return the list of correspoinding atom and cell list.
Note that the first atom is keeped to be located at the primitive cell. 
"""
function apply_symop_to_basis_with_shift(
	basis::IndicesUniqueList,
	isym::Integer,
	symop::SymmetryOperation,
	map_sym_cell::AbstractArray{AtomCell},
	map_s2p::AbstractVector,
	atoms_in_prim::AbstractVector{<:Integer},
	x_image_cart::AbstractArray,
	lattice_vectors::AbstractArray,
)::NTuple{2, Vector{Int}}

	atom_list = get_atomlist(basis)
	cell_list = get_celllist(basis)

	header_atom = atom_list[begin]
	header_cell = cell_list[begin]

	# firstly, apply the symmetry operation to atoms
	moved_coords = Vector{Vector{Float64}}()
	for (atom, cell) in zip(atom_list, cell_list)
		translation_cart = lattice_vectors * symop.translation_frac
		coords_tmp =
			symop.rotation_cart * x_image_cart[:, atom, cell] + translation_cart
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
		error("Something is wrong.")
	end
	translation_vec =
		calc_relvec_in_cart(moved_header_prim, moved_header_indices, x_image_cart)

	# thirdly, shift all atoms by using the translation vector
	shifted_coords = Vector{Vector{Float64}}()
	for coords in moved_coords
		push!(shifted_coords, coords + translation_vec)
	end

	# finally, find corresponding atom indices
	result_atom_list = Int[]
	result_cell_list = Int[]
	for coords in shifted_coords
		found = false
		for iatom in axes(x_image_cart, 2)
			for icell in axes(x_image_cart, 3)
				if isapprox(coords, x_image_cart[:, iatom, icell], atol = 1e-5)
					push!(result_atom_list, iatom)
					push!(result_cell_list, icell)
					found = true
					break
				end
			end
			found && break
		end
		if !found
			error("Could not find atomic position after symmetry operation.")
		end
	end

	return result_atom_list, result_cell_list
end

function calc_relvec_in_frac(atom1::NTuple{2, Integer},# (atom, cell)
	atom2::NTuple{2, Integer},
	x_image_frac::AbstractArray{<:Real, 3},
)::Vector{Float64}
	result::Vector{Float64} =
		x_image_frac[:, atom1[1], atom1[2]] - x_image_frac[:, atom2[1], atom2[2]]

	return result
end

