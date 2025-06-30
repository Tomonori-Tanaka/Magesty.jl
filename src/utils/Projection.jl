using SparseArrays
using Base.Threads
using StaticArrays

using ..RotationMatrix

"""
Construct projection matrices for basis functions

Parameters
----------
basisdict : Dictionary of basis functions
structure : Structure information
symmetry : Symmetry information

Returns
-------
Tuple{Dict{Int, Matrix{Float64}}, Dict{Int, Any}}
	Dictionary of projection matrices and dictionary of matrices for each symmetry operation
"""
function construct_projectionmatrix(
	basisdict::AbstractDict{<:Integer, <:AbstractVector},
	structure::Structure,
	symmetry::Symmetry,
)::Tuple{Vector{Matrix{Float64}}, Vector{Int}}

	# aliases
	x_image_cart::Array{Float64, 3} = structure.x_image_cart
	# Convert to SMatrix for better performance
	lattice_vectors::SMatrix{3, 3, Float64} =
		SMatrix{3, 3, Float64}(structure.supercell.lattice_vectors)
	symdata::Vector{SymmetryOperation} = symmetry.symdata
	map_sym = symmetry.map_sym
	map_s2p::Vector{Maps} = symmetry.map_s2p
	atoms_in_prim::Vector{Int} = symmetry.atoms_in_prim
	tol::Real = symmetry.tol
	result_projection = Vector{Matrix{Float64}}(undef, length(basisdict))
	num_nonzero_projection_list = Vector{Int}(undef, length(basisdict))

	idx_list = collect(keys(basisdict))
	# Process each basis function in parallel
	@threads for idx in idx_list
		# Create thread-local storage for calculations
		basislist = basisdict[idx]
		dim = length(basislist)
		local_projection_mat = zeros(Float64, dim, dim)
		local_num_nonzero_projections = 0

		# Calculate projection matrix for each symmetry operation
		for (n, symop) in enumerate(symdata), time_rev_sym in [false, true]
			projection_mat_per_symop = proj_matrix_a_symop(
				basislist,
				symop,
				n,
				map_s2p,
				atoms_in_prim,
				map_sym,
				x_image_cart,
				lattice_vectors,
				threshold_digits = 10,
				tol = tol,
				time_reversal_sym = time_rev_sym,
			)
			if nnz(projection_mat_per_symop) != 0
				local_num_nonzero_projections += 1
				local_projection_mat += projection_mat_per_symop
			end
		end

		# Write results to shared arrays
		@inbounds begin
			num_nonzero_projection_list[idx] = local_num_nonzero_projections
			result_projection[idx] = local_projection_mat
		end
	end

	return result_projection, num_nonzero_projection_list
end

function proj_matrix_a_symop(
	basislist::AbstractVector{IndicesUniqueList},
	symop::SymmetryOperation,
	isym::Integer,
	map_s2p,
	atoms_in_prim,
	map_sym,
	x_image_cart,
	lattice_vectors,
	;
	tol::Real = 1e-5,
	threshold_digits::Integer = 6,
	time_reversal_sym::Bool = false,
)::SparseMatrixCSC{Float64, Int}

	result_matrix = spzeros(Float64, length(basislist), length(basislist))

	# Early return if the symmetry operation includes translation part beyond the primitive cell.
	# if symop.is_translation_included == true
	# 	return result_matrix
	# end

	# aliases
	rotation_cart = SMatrix{3, 3, Float64}(symop.rotation_cart)

	for (ir, rbasis::IndicesUniqueList) in enumerate(basislist)  # right-hand basis

		moved_atomlist, moved_celllist =
			apply_symop_to_basis_with_shift(
				rbasis,
				isym,
				symop,
				map_sym,
				map_s2p,
				atoms_in_prim,
				x_image_cart,
				lattice_vectors,
				tol = tol,
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
			elseif is_translationally_equiv_basis(
				moved_rbasis,
				basis,
				atoms_in_prim,
				map_s2p,
				x_image_cart,
				tol = tol,
			)
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

		atom_list = [indices.atom for indices in moved_rbasis]
		l_list = [indices.l for indices in moved_rbasis]
		cell_list = [indices.cell for indices in moved_rbasis]
		partial_moved_basis::Vector{IndicesUniqueList} =
			AtomicIndices.product_indices(
				atom_list,
				l_list,
				cell_list,
			)
		partial_r_idx = findfirst(x -> equivalent(x, moved_rbasis), partial_moved_basis)
		if isnothing(partial_r_idx)
			error("Failed to find partial basis index in basis projection.")
		end

		# calculate rotation matrix
		is_proper::Bool = symop.is_proper
		multiplier::Float64 = 1.0
		if is_proper
			rotmat = rotation_cart
		else
			rotmat = -1 * rotation_cart
			# multiplier = (-1)^(get_total_L(moved_rbasis))
		end

		if time_reversal_sym
			multiplier *= (-1)^(get_total_L(moved_rbasis))
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
			result_matrix[il, ir] = rotmat_kron[partial_l_idx, partial_r_idx]
		end
	end

	result_matrix = round.(result_matrix, digits = threshold_digits)
	if !(RotationMatrix.is_orthogonal(result_matrix, tol = 1e-7))
		println("symmetry operation index: $isym")
		println(symop)
		display(result_matrix)
		error("not orthogonal")
	end

	return result_matrix
end

"""
Apply the symmetry operation and return the list of correspoinding atom and cell list.
Note that the first atom is keeped to be located at the primitive cell. 
"""
function apply_symop_to_basis_with_shift(
	basis::IndicesUniqueList,
	isym::Integer,
	symop::SymmetryOperation,
	map_sym::AbstractArray{<:Integer},
	map_s2p::AbstractVector,
	atoms_in_prim::AbstractVector{<:Integer},
	x_image_cart::AbstractArray,
	lattice_vectors::AbstractMatrix,
	;
	tol::Real = 1e-5,
)::NTuple{2, Vector{Int}}

	atom_list = [indices.atom for indices in basis]
	cell_list = [indices.cell for indices in basis]

	header_atom = atom_list[begin]
	header_cell = cell_list[begin]

	# firstly, apply the symmetry operation to atoms
	translation_cart = SVector{3, Float64}(lattice_vectors * symop.translation_frac)
	moved_coords = Matrix{Float64}(undef, 3, length(atom_list))
	for (i, (atom, cell)) in enumerate(zip(atom_list, cell_list))
		moved_coords[:, i] =
			symop.rotation_cart * x_image_cart[:, atom, cell] + translation_cart
	end

	# secondly, calculate translation vector to shift first atom to the primitive cell
	moved_header_atom::Int = map_sym[header_atom, isym]
	moved_header_atom_in_prim::Int = atoms_in_prim[map_s2p[moved_header_atom].atom]

	relative_vec_cart::SVector{3, Float64} =
		x_image_cart[:, moved_header_atom_in_prim, 1] - moved_coords[:, 1]

	# shift the moved_coords to the primitive cell
	shifted_coords = moved_coords .+ relative_vec_cart


	# finally, find corresponding atom indices
	result_atom_list = Int[]
	result_cell_list = Int[]

	for coords in eachcol(shifted_coords)
		found = false
		for iatom in axes(x_image_cart, 2)
			for icell in axes(x_image_cart, 3)
				if isapprox(
					coords,
					SVector{3, Float64}(x_image_cart[:, iatom, icell]),
					atol = tol,
				)
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

function calc_relvec_in_frac(
	atom1::NTuple{2, Integer},# (atom, cell)
	atom2::NTuple{2, Integer},
	x_image_frac::AbstractArray{<:Real, 3},
)::SVector{3, Float64}
	# Use SVector for better performance
	result = SVector{3, Float64}(
		x_image_frac[:, atom1[1], atom1[2]] - x_image_frac[:, atom2[1], atom2[2]],
	)
	return result
end


