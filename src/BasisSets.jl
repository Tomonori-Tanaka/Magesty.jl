"""
This module provides functionality for handling basis sets in crystal structure calculations.
It includes tools for constructing, classifying, and manipulating basis functions that are
adapted to the symmetry of the crystal structure.
"""
module BasisSets

using Combinatorics
using DataStructures
using LinearAlgebra
using Printf
using StaticArrays
using ..SortedContainer
using ..AtomCells
using ..AtomicIndices
using ..ConfigParser
using ..Structures
using ..Symmetries
using ..Clusters
using ..SALCs

include("./utils/Projection.jl")

export BasisSet

"""
	struct BasisSet

Represents a set of basis functions for atomic interactions in a crystal structure.
This structure is used to store and manage basis functions that are adapted to the symmetry of the crystal.

# Fields
- `basislist::SortedCountingUniqueVector{IndicesUniqueList}`: List of unique basis indices with their counts
- `classified_basisdict::Dict{Int, SortedCountingUniqueVector}`: Dictionary mapping symmetry labels to basis sets
- `projection_dict::Dict{Int, Matrix{Float64}}`: Dictionary of projection matrices for each symmetry label
- `each_projection_dict::Any`: Dictionary containing individual projection information
- `salc_list::Vector{SALC}`: List of symmetry-adapted linear combinations

# Constructors
	BasisSet(structure::Structure, symmetry::Symmetry, cluster::Cluster, lmax::AbstractMatrix{<:Integer}, bodymax::Integer)

Constructs a new `BasisSet` instance for atomic interactions in a crystal structure.

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information for atomic interactions
- `lmax::AbstractMatrix{<:Integer}`: Matrix of maximum angular momentum values for each atomic species and body [nkd × nbody]
- `bodymax::Integer`: Maximum number of bodies in interactions

# Returns
- `BasisSet`: A new basis set instance containing:
  - List of unique basis functions
  - Symmetry-classified basis dictionary
  - Projection matrices
  - Symmetry-adapted linear combinations

# Examples
```julia
# Create a basis set for a structure with 2 atomic species and 3-body interactions
lmax_matrix = [2 3; 3 2]  # lmax for each species and body
basis = BasisSet(structure, symmetry, cluster, lmax_matrix, 3)
```

# Note
The constructor performs the following steps:
1. Constructs the basis list by considering all possible combinations of atoms and angular momenta
2. Classifies basis functions by symmetry operations
3. Constructs projection matrices for each symmetry label
4. Generates symmetry-adapted linear combinations (SALCs)
"""
struct BasisSet
	basislist::SortedCountingUniqueVector{IndicesUniqueList}
	classified_basisdict::Dict{Int, SortedCountingUniqueVector}
	# projection_list::Vector{Matrix{Float64}}
	salc_list::Vector{SALC}

	function BasisSet(
		structure::Structure,
		symmetry::Symmetry,
		cluster::Cluster,
		lmax::AbstractMatrix{<:Integer},    # [≤ nkd, ≤ nbody]
		bodymax::Integer,
		;
		verbosity::Bool = true,
	)
		# Start timing
		start_time = time_ns()

		if verbosity
			println(
				"""

				BASIS SET
				=========
				""",
			)
		end

		# Validate input parameters
		nkd, nbody = size(lmax)
		bodymax > 0 || throw(ArgumentError("bodymax must be positive"))
		nbody ≥ bodymax ||
			throw(ArgumentError("lmax matrix must have at least bodymax columns"))

		# Construct basis list
		# basislist consists of all possible basis functions which is the product of spherical harmonics.
		if verbosity
			println("Constructing basis list...")
		end
		basislist = construct_basislist(
			structure,
			symmetry,
			cluster,
			lmax,
			bodymax,
		)

		# Classify basis functions by symmetry
		classified_basisdict = classify_basislist(basislist, symmetry.map_sym)

		# display(classified_basisdict)

		# Construct projection matrices
		if verbosity
			println("Constructing projection matrices...")
		end
		projection_list, num_nonzero_projection_list = construct_projectionmatrix(
			classified_basisdict,
			structure,
			symmetry,
		)

		# for (idx, proj) in enumerate(projection_list)
		# 	if proj != proj'
		# 		throw(
		# 			DomainError(
		# 				"Projection matrix at index $idx is not Hermitian)",
		# 			),
		# 		)
		# 	end
		# end

		# Generate symmetry-adapted linear combinations
		if verbosity
			println("Deriving SALCs...\n")
		end
		salc_list = Vector{SALC}()
		for idx in eachindex(projection_list)
			eigenval, eigenvec = eigen(Symmetric(projection_list[idx]))
			# eigenval, eigenvec = eigs(projection_list[idx], nev = )
			eigenval = eigenval ./ num_nonzero_projection_list[idx]

			# Set tolerance constants
			tol_zero = 1e-5
			tol_eigen = 1e-8

			# Process eigenvalues and eigenvectors
			eigenval = real.(round.(eigenval, digits = 6))
			eigenvec = round.(eigenvec .* (abs.(eigenvec) .≥ tol_zero), digits = 10)

			!check_eigenval(eigenval, tol = tol_eigen) &&
				throw(
					DomainError(
						"Critical error: Eigenvalues must be either 0 or 1. index: $idx",
					),
				)

			# Process vectors corresponding to eigenvalue 1
			for idx_eigenval in findall(x -> isapprox(x, 1.0, atol = tol_eigen), eigenval)
				eigenvec_real = to_real_vector(eigenvec[:, idx_eigenval])
				# Filter small values in the eigenvector
				# eigenvec_real = filter_eigenvec(eigenvec_real, rtol = 1e-3)
				eigenvec_real = flip_vector_if_negative_sum(eigenvec_real)
				push!(
					salc_list,
					SALC(classified_basisdict[idx], eigenvec_real / norm(eigenvec_real)),
				)
			end
		end

		salc_list = filter(salc -> !is_identically_zero(salc), salc_list)

		if verbosity
			print_basisset_stdout(salc_list)
			elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds
			println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
			println("-------------------------------------------------------------------")
		end


		return new(
			basislist,
			classified_basisdict,
			salc_list,
		)
	end
end

function BasisSet(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	config::Config4System,
	;
	verbosity::Bool = true,
)
	return BasisSet(structure, symmetry, cluster, config.lmax, config.nbody, verbosity = verbosity)
end

function construct_basislist(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	lmax_mat::AbstractMatrix{<:Integer},
	bodymax::Integer,
)::SortedCountingUniqueVector{IndicesUniqueList}

	result_basislist = SortedCountingUniqueVector{IndicesUniqueList}()

	# Get aliases for better readability
	kd_int_list = structure.supercell.kd_int_list
	cluster_list = cluster.cluster_list

	# Handle 1-body case
	for iat in symmetry.atoms_in_prim
		# Use @view for better performance when accessing matrix row
		lmax = @view(lmax_mat[kd_int_list[iat], :])

		for l in 1:lmax[1]
			iul::Vector{Indices} = indices_singleatom(iat, l, 1)
			for indices::Indices in iul
				push!(result_basislist, IndicesUniqueList(indices))
			end
		end
	end

	# Process multi-body cases
	for body in 2:bodymax
		for cluster in cluster_list[body-1]
			# Convert cluster into atomlist, llist, and celllist
			atomlist, llist, celllist =
				get_atomsls_from_cluster(cluster, lmax_mat, kd_int_list)

			for iul in product_indices_of_all_comb(atomlist, llist, celllist)
				for basis in result_basislist
					# Check for equivalent clusters in primitive cell
					if equivalent(basis, iul)
						result_basislist.counts[basis] += 1
						@goto skip
						# Check for translationally equivalent clusters
					elseif is_translationally_equiv_basis(
						iul,
						basis,
						symmetry.atoms_in_prim,
						symmetry.map_s2p,
						structure.x_image_cart,
						tol = symmetry.tol,
					)
						result_basislist.counts[basis] += 1
						@goto skip
					end
				end
				push!(result_basislist, iul)
				@label skip
			end
		end
	end

	return result_basislist
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
)::AbstractDict{Int, SortedCountingUniqueVector}

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

"""
Checks if two basis sets are translationally equivalent.

# Arguments
- `basis_target::IndicesUniqueList`: The target basis set to check
- `basis_ref::IndicesUniqueList`: The reference basis set
- `atoms_in_prim::AbstractVector{<:Integer}`: List of atoms in the primitive cell
- `map_s2p::AbstractVector`: Mapping from supercell to primitive cell
- `x_image_cart::AbstractArray{<:Real, 3}`: Cartesian coordinates of atoms in the supercell
- `tol::Real = 1e-5`: Tolerance for floating-point comparisons


"""
function is_translationally_equiv_basis(
	basis_target::IndicesUniqueList,
	basis_ref::IndicesUniqueList,
	atoms_in_prim::AbstractVector{<:Integer},
	map_s2p::AbstractVector,
	x_image_cart::AbstractArray{<:Real, 3},
	;
	tol::Real = 1e-5,
)::Bool

	# Early return if the atomlist is the same including the order
	if [indices.atom for indices in basis_target] == [indices.atom for indices in basis_ref]
		return false
		# Early return if the first atom is the same
		# because this function is intended to be used for different first atoms but translationally equivalent clusters
	elseif basis_target[1].atom == basis_ref[1].atom
		return false
	end

	for i in eachindex(basis_target)
		iatom = basis_target[i].atom
		icell = basis_target[i].cell
		iatom_in_prim = atoms_in_prim[map_s2p[iatom].atom]

		# cartesian relative vector b/w iatom and iatom_in_prim
		relvec::SVector{3, Float64} =
			calc_relvec_in_cart((iatom_in_prim, 1), (iatom, icell), x_image_cart)

		moved_atomlist = Int[]
		moved_celllist = Int[]
		for indices::Indices in basis_target
			# corresponding atom and cell obtained by adding relvec
			crrsp_atom, crrsp_cell = find_corresponding_atom(
				(indices.atom, indices.cell),
				relvec,
				x_image_cart,
				tol = tol,
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
Finds the cartesian relative vector between 2 atoms specified by (atom, cell) tuples, where cell means virtual cell index (1 <= cell <= 27).
The equation is
r(atom2) - r(atom1)
"""
function calc_relvec_in_cart(
	atom1::NTuple{2, Integer},# (atom, cell)
	atom2::NTuple{2, Integer},
	x_image_cart::AbstractArray{<:Real, 3},
)::SVector{3, Float64}
	# Use SVector for better performance with vector operations
	relvec = SVector{3, Float64}(
		x_image_cart[:, atom1[1], atom1[2]] - x_image_cart[:, atom2[1], atom2[2]],
	)
	return relvec
end

"""
Move an atom by a cartesian relative vector (calculated by `calc_relvec_in_cart` function) and find corresponding atom and cell indices as a tuple.
"""
function find_corresponding_atom(
	atom::NTuple{2, Int},# (atom, cell)
	relvec::AbstractVector{<:Real},
	x_image_cart::AbstractArray{<:Real, 3},
	;
	tol::Real = 1e-5,
)::NTuple{2, Int}
	# Use SVector for better performance with vector operations
	moved_coords =
		SVector{3, Float64}(x_image_cart[:, atom[1], atom[2]]) + SVector{3, Float64}(relvec)

	num_atoms = size(x_image_cart, 2)
	num_cells = size(x_image_cart, 3)
	for iatom in 1:num_atoms
		for icell in 1:num_cells
			if isapprox(
				SVector{3, Float64}(x_image_cart[:, iatom, icell]),
				moved_coords,
				atol = tol,
			)
				return (iatom, icell)
			end
		end
	end
	error("atom: $(atom[1]), cell: $(atom[2]) \n
	relvec: $relvec \n
	No matching (atom, cell) indices found.")
end

"""
	check_eigenval(eigenval::AbstractVector; tol = 1e-8)::Bool

Checks whether all eigenvalues in the given vector are approximately 0 or 1.

# Arguments
- `eigenval::AbstractVector`: Vector of eigenvalues to check
- `tol::Real = 1e-8`: Tolerance for floating-point comparisons

# Returns
- `Bool`: `true` if all eigenvalues are approximately 0 or 1, `false` otherwise

# Examples
```julia
# Check eigenvalues from a projection matrix
eigenvals = [0.0, 1.0, 0.0, 1.0]
is_valid = check_eigenval(eigenvals)  # true

# With some tolerance
eigenvals = [0.0, 1.0 + 1e-9, 0.0, 1.0]
is_valid = check_eigenval(eigenvals)  # true

# Invalid eigenvalues
eigenvals = [0.0, 1.0, 0.5, 1.0]
is_valid = check_eigenval(eigenvals)  # false
```
"""
function check_eigenval(eigenval::AbstractVector; tol = 1e-8)::Bool
	invalid_values =
		filter(x -> !isapprox(x, 0, atol = tol) && !isapprox(x, 1, atol = tol), eigenval)
	if !isempty(invalid_values)
		throw(BasisSetError("Invalid eigenvalues found: $invalid_values"))
	end
	return true
end

function to_real_vector(v::AbstractVector, atol::Real = 1e-12)
	if eltype(v) <: Complex && any(zi -> !isapprox(imag(zi), 0, atol = atol), v)
		idx = findfirst(zi -> !isapprox(imag(zi), 0, atol = atol), v)
		throw(
			DomainError(
				v[idx],
				"Vector contains complex numbers with significant imaginary parts",
			),
		)
	end
	return real.(v)
end

"""
	is_identically_zero(salc::SALC, atol::Real = 1e-8)::Bool

Check if a SALC (Symmetry-Adapted Linear Combination) always returns zero for arbitrary spin configuration.

# Arguments
- `salc::SALC`: The SALC to check
- `atol::Real = 1e-8`: Absolute tolerance for floating-point comparisons

# Returns
- `Bool`: `true` if the SALC is identically zero, `false` otherwise

"""
function is_identically_zero(salc::SALC, atol::Real = 1e-6)::Bool
	group_lists::Vector{Vector{Int}} = group_same_basis(salc)
	# println("group_lists:", group_lists)
	for index_list in group_lists
		partial_sum::Float64 = 0.0
		for index in index_list
			partial_sum += salc.multiplicity[index] * salc.coeffs[index]
		end
		if !isapprox(partial_sum, 0.0, atol = atol)
			return false
		end
	end
	return true
end

"""
	group_same_basis(salc::SALC)::Vector{Vector{Int}}

Group basis functions in a SALC that have the same atom, l, and m values (cell is not considered).

# Arguments
- `salc::SALC`: The SALC to group

# Returns
- `Vector{Vector{Int}}`: List of groups, where each group contains indices of basis functions with the same atom, l, and m values

# Examples
# For a SALC with basis functions:
# [(1, 1, 1, 1), (2, 1, -1, 1)]
# [(1, 1, 1, 1), (2, 1, -1, 6)]
# [(1, 1, 1, 1), (2, 1,  0, 1)]
# Returns: [[1, 2], [3]]
"""
function group_same_basis(salc::SALC)::Vector{Vector{Int}}
	group_dict = OrderedDict{Tuple{Vararg{NTuple{3, Int}}}, Vector{Int}}()

	for (i, basis::IndicesUniqueList) in enumerate(salc.basisset)
		atom_l_m_lists::Vector{Vector{Int}} = AtomicIndices.get_atom_l_m_list(basis)
		basisset_tuple::Tuple{Vararg{NTuple{3, Int}}} =
			Tuple(Tuple(atom_l_m) for atom_l_m in atom_l_m_lists)

		if haskey(group_dict, basisset_tuple)
			push!(group_dict[basisset_tuple], i)
		else
			group_dict[basisset_tuple] = [i]
		end
	end

	result_list = Vector{Vector{Int}}()
	for (key, value) in group_dict
		push!(result_list, value)
	end
	return result_list
end

"""
	filter_eigenvec(eigenvec::AbstractVector{<:Real}, rtol::Real = 1e-3)::Vector{Float64}
Filter the eigenvector to remove small values based on a relative tolerance and re-normalize the filtered vector.

"""
function filter_eigenvec(eigenvec::AbstractVector{<:Real}; rtol::Real = 1e-3)
	# Get the maximum absolute value in the eigenvector
	max_abs = maximum(abs.(eigenvec))
	# Filter the eigenvector to keep only values that are significant relative to the maximum
	filtered_vec = deepcopy(eigenvec)
	for i in eachindex(eigenvec)
		if abs(eigenvec[i]) < rtol * max_abs
			filtered_vec[i] = 0.0
		end
	end
	# Re-normalize the filtered vector
	return filtered_vec ./ norm(filtered_vec)  # Ensure the filtered vector is normalized

end

"""
	flip_vector_if_negative_sum(v::AbstractVector{<:Real}; tol::Real=1e-10)::AbstractVector{<:Real}

Flip the sign of the vector if the sum of the elements is negative.
If the sum is approximately zero (within tolerance), return the original vector.

# Arguments
- `v::AbstractVector{<:Real}`: Input vector
- `tol::Real=1e-10`: Tolerance for considering sum as zero

# Returns
- `AbstractVector{<:Real}`: Vector with sign flipped if sum is negative
"""
function flip_vector_if_negative_sum(
	v::AbstractVector{<:Real};
	tol::Real = 1e-10,
)::AbstractVector{<:Real}
	sum_v = sum(v)
	if abs(sum_v) < tol
		return v
	end
	return sum_v < 0 ? -v : v
end


function print_basisset_stdout(salc_list::AbstractVector{<:SALC})
	println(" Number of symmetry-adapted basis functions: $(length(salc_list))\n")
	println(" List of symmetry-adapted basis functions:")
	println(" # multiplicity  coefficient  basis")
	for (i, salc) in enumerate(salc_list)
		println(" $i-th salc")
		display(salc)
	end
	println("")
end

# Custom exception type for basis set errors
struct BasisSetError <: Exception
	msg::String
end

end