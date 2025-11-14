"""
This module provides functionality for handling basis sets in crystal structure calculations.
It includes tools for constructing, classifying, and manipulating basis functions that are
adapted to the symmetry of the crystal structure.
"""
module BasisSets

using Combinat
using Combinatorics
using DataStructures
using LinearAlgebra
using OffsetArrays
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
using ..RotationMatrix

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
		body1_lmax::Vector{Int},
		bodyn_lsum::OffsetArray{Int, 1},
		nbody::Integer,
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
		# Construct basis list
		# basislist consists of all possible basis functions which is the product of spherical harmonics.
		if verbosity
			println("Constructing basis list...")
		end
		basislist = construct_basislist(
			structure,
			symmetry,
			cluster,
			body1_lmax,
			bodyn_lsum,
			nbody,
		)

		basislist_simple = construct_basislist_simple(
			structure,
			symmetry,
			cluster,
			body1_lmax,
			bodyn_lsum,
			nbody,
		)
		# @show basislist_simple

		# Classify basis functions by symmetry
		classified_basisdict = classify_basislist(basislist, symmetry.map_sym)
		classified_basisdict_simple = classify_basislist_simple(basislist_simple, symmetry.map_sym)
		# display(classified_basisdict_simple)

		# display(classified_basisdict)

		# Construct projection matrices
		if verbosity
			println("Constructing projection matrices...")
		end
		projection_list_simple = projection_matrix_simple(classified_basisdict_simple, symmetry)

		for idx in eachindex(projection_list_simple)
			eigenvals, eigenvecs = eigen(projection_list_simple[idx])
			eigenvals = real.(round.(eigenvals, digits = 6))
			eigenvecs = round.(eigenvecs .* (abs.(eigenvecs) .≥ 1e-8), digits = 10)
			if !is_proper_eigenvals(eigenvals)
				display(classified_basisdict_simple[idx-1])
				display(classified_basisdict_simple[idx])
				error("Critical error: Eigenvalues must be either 0 or 1. index: $idx")
			end
			for idx_eigenval in findall(x -> isapprox(x, 1.0, atol = 1e-8), eigenvals)
				eigenvec = eigenvecs[:, idx_eigenval]
				eigenvec = flip_vector_if_negative_sum(eigenvec)
			end
		end
		print("Done")


		# for (idx, proj) in enumerate(projection_list)
		# 	if proj != proj'
		# 		throw(
		# 			DomainError(
		# 				"Projection matrix at index $idx is not Hermitian)",
		# 			),
		# 		)
		# 	end
		# end

		projection_list, num_nonzero_projection_list = construct_projectionmatrix(
			classified_basisdict,
			structure,
			symmetry,
		)

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

			!is_proper_eigenvals(eigenval, tol = tol_eigen) &&
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
	return BasisSet(
		structure,
		symmetry,
		cluster,
		config.body1_lmax,
		config.bodyn_lsum,
		config.nbody,
		verbosity = verbosity,
	)
end

function construct_basislist(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	body1_lmax::Vector{Int},
	bodyn_lsum::OffsetArray{Int, 1},
	nbody::Integer,
)::SortedCountingUniqueVector{IndicesUniqueList}

	result_basislist = SortedCountingUniqueVector{IndicesUniqueList}()

	# Get aliases for better readability
	kd_int_list = structure.supercell.kd_int_list
	cluster_list = cluster.cluster_list

	# Handle 1-body case (single-thread for safety; cost is small)
	for iat in symmetry.atoms_in_prim
		# Use @view for better performance when accessing matrix row
		lmax = body1_lmax[kd_int_list[iat]]

		for l in 1:lmax[1]
			iul::Vector{Indices} = indices_singleatom(iat, l, 1)
			for indices::Indices in iul
				push!(result_basislist, IndicesUniqueList(indices))
			end
		end
	end

	# Process multi-body cases
	for body in 2:nbody
		for cluster in cluster_list[body-1]
			for lsum in 1:bodyn_lsum[body]

				# skip odd lsum cases due to the time-reversal symmetry
				if mod(lsum, 2) == 1
					continue
				end

				iul_list = listup_basislist(cluster, lsum)

				for iul in iul_list
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
	end

	return result_basislist
end

function construct_basislist_simple(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	body1_lmax::Vector{Int},
	bodyn_lsum::OffsetArray{Int, 1},
	nbody::Integer,
)::SortedCountingUniqueVector{SHProduct}

	result_basislist = SortedCountingUniqueVector{SHProduct}()
	cluster_dict::Dict{Int, SortedCountingUniqueVector{SortedVector{Int}}} = cluster.cluster_dict

	# Handle 1-body case
	for iat in symmetry.atoms_in_prim
		lmax = body1_lmax[structure.supercell.kd_int_list[iat]]
		for l in 2:lmax[1] # skip l = 1 because it is prohibited by the time-reversal symmetry
			if l % 2 == 1 # skip odd l cases due to the time-reversal symmetry
				continue
			end
			shsi_list::Vector{SHSiteIndex} = shsiteindex_singleatom(iat, l)
			for shsi::SHSiteIndex in shsi_list
				push!(result_basislist, SHProduct([shsi]))
			end
		end
	end

	# Process multi-body cases
	for body in 2:nbody
		for atom_list::SortedVector{Int} in cluster_dict[body]
			count = cluster_dict[body].counts[atom_list]

			shp_list::Vector{SHProduct} = listup_basislist_simple(atom_list, bodyn_lsum[body])
			for shp::SHProduct in shp_list
				push!(result_basislist, shp, count)
			end
		end
	end
	return result_basislist
end

function listup_basislist_simple(
	atom_list::SortedVector{<:Integer},
	lsum::Integer,
)::Vector{SHProduct}

	result_basislist = Vector{SHProduct}()
	l_list = Combinat.compositions(lsum, length(atom_list); min = 1)

	for l_vec::Vector{Int} in l_list
		shp_list = product_shsiteindex(atom_list, l_vec)
		for shp::SHProduct in shp_list
			push!(result_basislist, shp)
		end
	end

	return result_basislist
end


function listup_basislist(
	cluster::AbstractVector{AtomCell},
	lsum::Integer,
)::Vector{IndicesUniqueList}

	result_basislist = Vector{IndicesUniqueList}()
	nbody = length(cluster)

	atomlist = [atomcell.atom for atomcell in cluster]
	celllist = [atomcell.cell for atomcell in cluster]

	l_list = Combinat.compositions(lsum, nbody; min = 1)
	for l_vec::Vector{Int} in l_list
		iul_list = product_indices(atomlist, l_vec, celllist)
		for iul::IndicesUniqueList in iul_list
			push!(result_basislist, iul)
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

function classify_basislist_simple(
	basislist::AbstractVector{SHProduct},
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
			dict[idx] = SortedCountingUniqueVector{SHProduct}()
		end
	end

	for (basis, label) in zip(basislist, label_list)
		if !(haskey(dict, label))
			dict[label] = SortedCountingUniqueVector{SHProduct}()
		end
		push!(dict[label], basis, getcount(basislist, basis))
	end

	return dict
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
function is_proper_eigenvals(eigenval::AbstractVector; tol = 1e-8)::Bool
	invalid_values =
		filter(x -> !isapprox(x, 0, atol = tol) && !isapprox(x, 1, atol = tol), eigenval)
	if !isempty(invalid_values)
		return false
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

function projection_matrix_simple(
	basisdict::AbstractDict,
	symmetry::Symmetry,
)::Vector{Matrix{Float64}}

	result_projections = Vector{Matrix{Float64}}(undef, length(basisdict))

	idx_list = sort(collect(keys(basisdict)))
	for idx in idx_list
		basislist::SortedCountingUniqueVector{SHProduct} = basisdict[idx]
		dim = length(basislist)
		local_projection_mat = zeros(Float64, dim, dim)
		for (n, symop) in enumerate(symmetry.symdata), time_rev_sym in [false, true]
			projection_mat_per_symop = proj_matrix_a_symop_simple(
				basislist,
				symop,
				@view(symmetry.map_sym[:, n]),
				symmetry.map_s2p,
				symmetry.map_sym,
				symmetry.symnum_translation,
				time_rev_sym,
			)
			local_projection_mat += projection_mat_per_symop
		end
		local_projection_mat = local_projection_mat ./ (2 * symmetry.nsym)
		local_projection_mat = hermitianpart(local_projection_mat)
		# local_projection_mat = round.(local_projection_mat, digits = 6)

		if !is_symmetric(local_projection_mat, tol = 1e-10)
			error("Projection matrix is not symmetric. index: $idx")
		end

		result_projections[idx] = local_projection_mat
	end
	return result_projections
end

function proj_matrix_a_symop_simple(
	basislist::SortedCountingUniqueVector{SHProduct},
	symop::SymmetryOperation,
	map_sym_per_symop::AbstractVector{<:Integer},
	map_s2p::AbstractVector{Maps},
	map_sym::AbstractMatrix{<:Integer},
	symnum_translation::AbstractVector{<:Integer},
	time_rev_sym::Bool,
)::Matrix{Float64}

	# collect atom list in basislist used for symmetry operation
	atom_list = [[shsi.i for shsi in basis] for basis in basislist]

	projection_mat = zeros(Float64, length(basislist), length(basislist))

	for (j, basis_j::SHProduct) in enumerate(basislist)
		lco_j = operate_symop(
			basis_j,
			atom_list,
			symop,
			map_sym_per_symop,
			map_sym,
			map_s2p,
			symnum_translation,
			time_rev_sym,
		)
		for (i, basis_i::SHProduct) in enumerate(basislist)
			projection_mat[i, j] = inner_product(basis_i, lco_j)
		end
	end
	if all(
		isapprox(projection_mat[i, j], 0.0, atol = 1e-8) for i in eachindex(basislist) for
		j in eachindex(basislist)
	)
		@assert false "Projection matrix is zero matrix"
	end
	return projection_mat
end

function operate_symop(
	basis::SHProduct,
	atom_list::AbstractVector{<:AbstractVector{<:Integer}},
	symop::SymmetryOperation,
	map_sym_per_symop::AbstractVector{<:Integer},
	map_sym::AbstractMatrix{<:Integer},
	map_s2p::AbstractVector{Maps},
	symnum_translation::AbstractVector{<:Integer},
	time_rev_sym::Bool,
)::LinearCombo
	# Apply the symmetry operation to atoms
	new_atom_list = [map_sym_per_symop[shsi.i] for shsi in basis]

	# Shift the new atom list to the primitive cell
	# Find the source atoms in supercell that map to new_atom_list under translation
	prim_atom_list = similar(new_atom_list)
	found = false
	for i in eachindex(new_atom_list)
		# Get translation operation for atom i
		map_s2p_i = map_s2p[new_atom_list[i]]
		itran = map_s2p_i.translation  # 1 <= itran <= ntran
		symnum_translation_i = symnum_translation[itran]

		# Find source atom that maps to new_atom_list[i] under this translation
		# i.e., find k such that map_sym[k, symnum_translation_i] == new_atom_list[i]
		source_atom_i = findfirst(k -> map_sym[k, symnum_translation_i] == new_atom_list[i],
			1:size(map_sym, 1))
		if source_atom_i === nothing
			continue  # Try next reference atom
		end

		# Set source atom for position i
		prim_atom_list[i] = source_atom_i

		# Find source atoms for all other positions
		all_found = true
		for j in eachindex(new_atom_list)
			if j == i
				continue
			end
			# Find source atom that maps to new_atom_list[j] under the same translation
			source_atom_j = findfirst(k -> map_sym[k, symnum_translation_i] == new_atom_list[j],
				1:size(map_sym, 1))
			if source_atom_j === nothing
				all_found = false
				break
			end
			prim_atom_list[j] = source_atom_j
		end

		# Check if all atoms were found and this atom list exists in the basis list
		if all_found && sort(prim_atom_list) in atom_list
			new_atom_list = sort(prim_atom_list)
			found = true
			break
		end
	end
	if !found
		@show [shsi.i for shsi in basis]
		@show prim_atom_list
		@show new_atom_list
		@show atom_list
		error("Failed to find corresponding atom in the primitive cell.")
	end

	new_basis = replace_atom(basis, new_atom_list)# replace atom only (l and m are kept)
	idx = corresponding_idx(new_basis)

	is_proper = symop.is_proper
	multiplier = 1.0
	rotmat = similar(symop.rotation_cart)
	if is_proper
		rotmat = symop.rotation_cart
	else
		rotmat = -1 * symop.rotation_cart
	end
	if time_rev_sym
		multiplier *= (-1)^(sum([shsi.l for shsi in new_basis]))
	end
	rotation_list = [Δl(shsi.l, rotmat2euler(rotmat)...) for shsi in new_basis]
	if length(new_basis) == 1
		rotmat_kron = multiplier * rotation_list[begin]
	else
		rotmat_kron = multiplier * kron(rotation_list...)
	end

	l_list = [shsi.l for shsi in new_basis]

	shp_list = product_shsiteindex(new_atom_list, l_list)
	coeffs = rotmat_kron[:, idx]
	return LinearCombo(shp_list, coeffs)
end

"""
	corresponding_idx(shp::SHProduct)::Int

Find the index of the rotation matrix element corresponding to a given `SHProduct`.

# Example
((l=1, m=-1), (l=1, m=-1)) -> 1
((l=1, m=-1), (l=1, m= 1)) -> 3
((l=1, m=1),  (l=1, m= 1)) -> 9
((l=2, m=2),  (l=2, m=2))  -> 25

# Arguments
- `shp::SHProduct`: The basis function (product of spherical harmonics)

# Returns
- `Int`: The 1-based linear index corresponding to `shp`
"""
function corresponding_idx(shp::SHProduct)::Int
	# Extract l and m from each factor
	l_list = [shsi.l for shsi in shp]
	m_list = [shsi.m for shsi in shp]

	@assert !isempty(l_list) "SHProduct must have at least one factor"

	# Local dimension of each site: d_i = 2l_i + 1
	dims = [2*l + 1 for l in l_list]

	# Map m_i = -l_i,…,l_i → p_i = 1,…,2l_i+1
	pos = [m + l + 1 for (m, l) in zip(m_list, l_list)]

	# Optional range check
	@assert all(1 .<= pos .<= dims) "m is out of range for given l"

	# Compute 1-based linear index with the last factor varying fastest
	idx    = 1
	stride = 1
	@inbounds for i in length(pos):-1:1
		idx += (pos[i] - 1) * stride
		stride *= dims[i]
	end

	return idx
end

"""
	is_symmetric(mat::AbstractMatrix{<:Real}; tol::Float64 = 1e-10) -> Bool

Check if a matrix is symmetric (i.e., `A ≈ A'`).

# Arguments
- `mat::AbstractMatrix{<:Real}`: The matrix to check
- `tol::Float64`: Tolerance for floating-point comparison (default: 1e-10)

# Returns
- `Bool`: `true` if the matrix is symmetric within the tolerance, `false` otherwise
"""
function is_symmetric(mat::AbstractMatrix{<:Real}; tol::Float64 = 1e-10)::Bool
	# Check if matrix is square
	if size(mat, 1) != size(mat, 2)
		return false
	end
	# Check if matrix is symmetric: A ≈ A'
	return isapprox(mat, mat', atol = tol)
end

end
