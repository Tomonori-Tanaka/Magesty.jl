"""
This module provides functionality for handling basis sets in crystal structure calculations.
It includes tools for constructing, classifying, and manipulating basis functions that are
adapted to the symmetry of the crystal structure.
"""
module BasisSets

using Base.Threads
using Combinat
using DataStructures
using LinearAlgebra
using OffsetArrays
using Printf

using ..CountingContainer
using ..SortedContainer
using ..AtomicIndices
using ..ConfigParser
using ..Structures
using ..Symmetries
using ..Clusters
using ..RotationMatrix
using ..Basis

export BasisSet


"""
	struct BasisSet

Represents a set of basis functions for atomic interactions in a crystal structure.
This structure is used to store and manage basis functions that are adapted to the symmetry of the crystal.

# Fields
- `coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis}`: List of coupled angular momentum basis functions with their multiplicities
- `salc_list::Vector{Vector{Basis.CoupledBasis_with_coefficient}}`: List of symmetry-adapted linear combinations (SALCs), where each element is a vector of coupled basis functions belonging to the same key group

# Constructors
	BasisSet(structure, symmetry, cluster, body1_lmax, bodyn_lsum, nbody; isotropy=false, verbosity=true)
	BasisSet(structure, symmetry, cluster, config; verbosity=true)

Constructs a new `BasisSet` instance for atomic interactions in a crystal structure.

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information for atomic interactions
- `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions for each atomic species
- `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
- `nbody::Integer`: Maximum number of bodies in interactions
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`
- `verbosity::Bool`: Whether to print progress information (default: `true`)

# Returns
- `BasisSet`: A new basis set instance containing coupled basis functions and symmetry-adapted linear combinations

# Examples
```julia
# Create a basis set using explicit parameters
body1_lmax = [2, 3]  # lmax for each atomic species
bodyn_lsum = OffsetArray([0, 0, 4, 6], 0:3)  # lsum for each body
basis = BasisSet(structure, symmetry, cluster, body1_lmax, bodyn_lsum, 3)

# Create a basis set using configuration
basis = BasisSet(structure, symmetry, cluster, config)
```

# Note
The constructor performs the following steps:
1. Constructs and classifies coupled basis functions using orbit information
2. Constructs projection matrices for each symmetry label
3. Generates symmetry-adapted linear combinations (SALCs) of `CoupledBasis_with_coefficient` objects
"""
struct BasisSet
	coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis}
	salc_list::Vector{Vector{Basis.CoupledBasis_with_coefficient}}
	angular_momentum_couplings::Vector{Basis.AngularMomentumCouplingResult}
end

function BasisSet(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	body1_lmax::Vector{Int},
	bodyn_lsum::OffsetArray{Int, 1},
	nbody::Integer,
	isotropy::Bool = false,
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
	if verbosity
		print("Constructing and classifying coupled basis list...")
	end
	classified_coupled_basisdict::Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}} =
		construct_and_classify_coupled_basislist(
			structure,
			symmetry,
			cluster,
			body1_lmax,
			bodyn_lsum,
			nbody,
			isotropy = isotropy,
		)
	classified_coupled_basisdict = filter_basisdict(classified_coupled_basisdict, symmetry)


	if verbosity
		println(" Done.")
	end

	if verbosity
		print("Constructing projection matrix...")
	end
	keys_list = sort(collect(keys(classified_coupled_basisdict)))
	num_keys = length(keys_list)
	# Pre-allocate array to store results for each key (preserving order)
	# Each key can have multiple SALC groups (one per eigenvector)
	salc_list_per_key = Vector{Vector{Vector{Basis.CoupledBasis_with_coefficient}}}(undef, num_keys)

	@threads for idx in 1:num_keys
		key = keys_list[idx]
		coupled_basislist = classified_coupled_basisdict[key]
		# display(coupled_basislist)
		projection_mat =
			projection_matrix_coupled_basis(coupled_basislist, symmetry)
		# projection_mat is real symmetric → use Hermitian + eigen! to save memory
		h_projection = Hermitian(projection_mat)
		# Free projection_mat memory before eigen decomposition
		projection_mat = nothing
		eigenvals, eigenvecs = eigen!(h_projection)
		eigenvals = real.(round.(eigenvals, digits = 6))
		eigenvecs = round.(eigenvecs .* (abs.(eigenvecs) .≥ 1e-8), digits = 10)
		if !is_proper_eigenvals(eigenvals)
			@show eigenvals
			display(coupled_basislist)
			# error("Critical error: Eigenvalues must be either 0 or 1. index: $key")
			@warn "Critical error: Eigenvalues must be either 0 or 1. index: $key"
		end
		Lf = coupled_basislist[1].Lf
		submatrix_dim = 2 * Lf + 1
		nbasis = length(coupled_basislist)

		# Create a list of SALC groups for this key (one per eigenvector)
		key_salc_groups = Vector{Vector{Basis.CoupledBasis_with_coefficient}}()

		for idx_eigenval in findall(x -> isapprox(x, 1.0, atol = 1e-8), eigenvals)
			eigenvec = eigenvecs[:, idx_eigenval]
			eigenvec = real.(eigenvec)
			eigenvec = round.(eigenvec .* (abs.(eigenvec) .≥ 1e-8), digits = 10)
			eigenvec = eigenvec / norm(eigenvec)
			# Flip sign to make sum(eigenvec) non-negative for determinism
			if sum(eigenvec) < 0
				eigenvec .= -eigenvec
			end

			# Create a new SALC group for this eigenvector
			salc_group = Vector{Basis.CoupledBasis_with_coefficient}()

			# Create CoupledBasis_with_coefficient for each basis in coupled_basislist
			for (idx_basis, cb) in enumerate(coupled_basislist)
				# Extract coefficients for this basis (Mf-dependent)
				coeff_start = (idx_basis - 1) * submatrix_dim + 1
				coeff_end = idx_basis * submatrix_dim
				coefficient = eigenvec[coeff_start:coeff_end]

				# Skip if coefficient is all zeros (or nearly zero)
				if isapprox(norm(coefficient), 0.0, atol = 1e-10)
					continue
				end

				# Get multiplicity from counts
				multiplicity = coupled_basislist.counts[cb]

				# Create CoupledBasis_with_coefficient and add to salc_group
				cbc = Basis.CoupledBasis_with_coefficient(cb, coefficient, multiplicity)
				push!(salc_group, cbc)
			end

			# Add this SALC group if it's not empty
			if !isempty(salc_group)
				push!(key_salc_groups, salc_group)
			end
		end

		# Store the list of SALC groups for this key (preserving order by index)
		salc_list_per_key[idx] = key_salc_groups
		# Free large matrices after processing each key
		h_projection = nothing
		eigenvecs = nothing
	end

	# Collect all SALC groups in order
	salc_list = Vector{Vector{Basis.CoupledBasis_with_coefficient}}()
	for key_salc_groups in salc_list_per_key
		for salc_group in key_salc_groups
			push!(salc_list, salc_group)
		end
	end
	if verbosity
		println(" Done.")
	end
	if verbosity
		elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds
		println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
		println("-------------------------------------------------------------------")
	end




	# Reconstruct basislist from classified dictionary for BasisSet storage
	result_basislist = SortedCountingUniqueVector{Basis.CoupledBasis}()
	for (key, classified_basislist) in classified_coupled_basisdict
		for cb in classified_basislist
			count = classified_basislist.counts[cb]
			push!(result_basislist, cb, count)
		end
	end

	# Collect angular momentum coupling results from cache
	angular_momentum_couplings = Vector{Basis.AngularMomentumCouplingResult}()
	# Collect unique ls combinations from coupled_basislist (preserve original order)
	ls_combinations_set = Set{Vector{Int}}()
	for cb in result_basislist
		push!(ls_combinations_set, cb.ls)
	end

	# Extract angular momentum coupling results from cache
	for ls_vec in ls_combinations_set
		# Check cache for this ls combination (try both isotropy options)
		# Cache key uses original ls order (not sorted)
		for isotropy_flag in [false, true]
			cache_key = (ls_vec, :none, isotropy_flag)
			if haskey(Basis._angular_momentum_cache, cache_key)
				bases_by_L, paths_by_L = Basis._angular_momentum_cache[cache_key]
				for Lf in sort(collect(keys(bases_by_L)))
					tensors = bases_by_L[Lf]
					Lseqs = paths_by_L[Lf]
					for (tensor, Lseq) in zip(tensors, Lseqs)
						push!(
							angular_momentum_couplings,
							Basis.AngularMomentumCouplingResult(
								ls_vec,  # Use original ls order
								Lseq,
								Lf,
								copy(tensor),  # Make a copy to avoid reference issues
							),
						)
					end
				end
				break  # Found in cache, no need to check other isotropy flag
			end
		end
	end

	return BasisSet(result_basislist, salc_list, angular_momentum_couplings)
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
		config.isotropy,
		verbosity = verbosity,
	)
end

"""
	construct_and_classify_coupled_basislist(
		structure::Structure,
		symmetry::Symmetry,
		cluster::Cluster,
		body1_lmax::Vector{Int},
		bodyn_lsum::OffsetArray{Int, 1},
		nbody::Integer;
		isotropy::Bool = false,
	) -> Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}

Construct coupled basis functions and classify them simultaneously using orbit information.

This function combines the functionality of `construct_coupled_basislist` and `classify_coupled_basislist_test`
by generating basis functions orbit-by-orbit and classifying them on-the-fly. This approach is more
memory-efficient and allows for better organization of basis functions by symmetry.

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information containing orbit classification
- `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions
- `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
- `nbody::Integer`: Maximum number of bodies in interactions
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`

# Returns
- `Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}`: Dictionary keyed by classification labels
  (based on `(nbody, Lf, sum(ls), Tuple(sort(ls)...))`), containing classified basis functions
"""
function construct_and_classify_coupled_basislist(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	body1_lmax::Vector{Int},
	bodyn_lsum::OffsetArray{Int, 1},
	nbody::Integer;
	isotropy::Bool = false,
)::Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}
	# Result dictionary for classified basis functions
	classified_dict = OrderedDict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}()
	label_map = Dict{Any, Int}()
	next_label = 0

	irreducible_cluster_dict::Dict{Int, SortedCountingUniqueVector{Vector{Int}}} =
		cluster.irreducible_cluster_dict
	cluster_orbits_dict::Dict{Int, Dict{Int, Vector{Vector{Int}}}} =
		cluster.cluster_orbits_dict

	# Handle 1-body case
	for iat in symmetry.atoms_in_prim
		lmax = body1_lmax[structure.supercell.kd_int_list[iat]]
		for l in 2:lmax # skip l = 1 because it is prohibited by the time-reversal symmetry
			if l % 2 == 1 # skip odd l cases due to the time-reversal symmetry
				continue
			end
			# For 1-body case, create CoupledBasis with single atom
			cb_list = tesseral_coupled_bases_from_tesseral_bases(
				[l],
				[iat];
				isotropy = isotropy,
			)
			for cb::Basis.CoupledBasis in cb_list
				# Classify on-the-fly: key is (nbody, Lf, sum(ls), Tuple(sort(ls)...))
				nbody_val = length(cb.ls)
				ls_sorted = Tuple(sort(cb.ls))
				key = (nbody_val, cb.Lf, sum(cb.ls), ls_sorted)

				label = get(label_map, key, 0)
				if label == 0
					next_label += 1
					label = next_label
					label_map[key] = label
					classified_dict[label] = SortedCountingUniqueVector{Basis.CoupledBasis}()
				end
				push!(classified_dict[label], cb, 1)
			end
		end
	end

	# Process multi-body cases using cluster_orbits_dict for efficient processing
	# Group by orbit first, then classify within each orbit
	for body in 2:nbody
		# Use cluster_orbits_dict to process clusters grouped by symmetry orbits
		if haskey(cluster_orbits_dict, body)
			for (orbit_index, orbit_clusters) in cluster_orbits_dict[body]
				# Process all clusters in this orbit
				# Clusters in the same orbit have the same projection matrix structure
				# Collect all basis functions from this orbit first
				orbit_basis_list = Vector{Basis.CoupledBasis}()
				# Use IdDict instead of Dict since CoupledBasis doesn't implement hash
				orbit_basis_counts = IdDict{Basis.CoupledBasis, Int}()

				for atom_list::Vector{Int} in orbit_clusters
					# Get multiplicity from irreducible_cluster_dict
					count = irreducible_cluster_dict[body].counts[atom_list]
					sorted_atom_list = sort(atom_list)
					cb_list::Vector{Basis.CoupledBasis} = listup_coupled_basislist(
						sorted_atom_list,
						bodyn_lsum[body];
						isotropy = isotropy,
					)

					# Collect basis functions from this cluster
					for cb::Basis.CoupledBasis in cb_list
						# Check for translationally equivalent within orbit
						found_equivalent = false
						for existing_cb in orbit_basis_list
							if is_translationally_equivalent_coupled_basis(
								cb,
								existing_cb,
								symmetry,
							)
								found_equivalent = true
								# Safe access: get existing count or 0, then add
								orbit_basis_counts[existing_cb] =
									get(orbit_basis_counts, existing_cb, 0) + count
								break
							end
						end
						if !found_equivalent
							push!(orbit_basis_list, cb)
							orbit_basis_counts[cb] = count
						end
					end
				end

				# Classify basis functions from this orbit
				# Use orbit index in classification key to group by orbit
				for cb::Basis.CoupledBasis in orbit_basis_list
					ls_sorted = Tuple(sort(cb.ls))
					# Include orbit_index in classification key to utilize orbit information
					key = (body, orbit_index, cb.Lf, sum(cb.ls), ls_sorted)

					label = get(label_map, key, 0)
					if label == 0
						next_label += 1
						label = next_label
						label_map[key] = label
						classified_dict[label] = SortedCountingUniqueVector{Basis.CoupledBasis}()
					end

					count = orbit_basis_counts[cb]
					push!(classified_dict[label], cb, count)
				end
			end
		end
	end

	return classified_dict
end

"""
	construct_coupled_basislist(
		structure::Structure,
		symmetry::Symmetry,
		cluster::Cluster,
		body1_lmax::Vector{Int},
		bodyn_lsum::OffsetArray{Int, 1},
		nbody::Integer;
		isotropy::Bool = false,
	) -> SortedCountingUniqueVector{Basis.CoupledBasis}

Construct a list of coupled angular momentum basis functions.

This function generates all coupled basis functions for atomic interactions, handling both 1-body
and multi-body cases. It uses cluster orbit information for efficient processing.

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information containing orbit classification
- `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions
- `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
- `nbody::Integer`: Maximum number of bodies in interactions
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`

# Returns
- `SortedCountingUniqueVector{Basis.CoupledBasis}`: List of coupled basis functions with multiplicities

# Note
- Skips l=1 and odd l values due to time-reversal symmetry constraints
- Uses cluster orbits for efficient processing of multi-body cases
- **Unused**: This function is currently not used. Its functionality has been integrated into `construct_and_classify_coupled_basislist`.
  Kept for reference purposes.
"""
function construct_coupled_basislist(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	body1_lmax::Vector{Int},
	bodyn_lsum::OffsetArray{Int, 1},
	nbody::Integer;
	isotropy::Bool = false,
)
	result_coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis} =
		SortedCountingUniqueVector{Basis.CoupledBasis}()
	irreducible_cluster_dict::Dict{Int, SortedCountingUniqueVector{Vector{Int}}} =
		cluster.irreducible_cluster_dict
	cluster_orbits_dict::Dict{Int, Dict{Int, Vector{Vector{Int}}}} =
		cluster.cluster_orbits_dict

	# Handle 1-body case
	for iat in symmetry.atoms_in_prim
		lmax = body1_lmax[structure.supercell.kd_int_list[iat]]
		for l in 2:lmax # skip l = 1 because it is prohibited by the time-reversal symmetry
			if l % 2 == 1 # skip odd l cases due to the time-reversal symmetry
				continue
			end
			# For 1-body case, create LinearCombo with single atom
			cb_list = tesseral_coupled_bases_from_tesseral_bases(
				[l],
				[iat];
				isotropy = isotropy,
			)
			for cb::Basis.CoupledBasis in cb_list
				push!(result_coupled_basislist, cb, 1)
			end
		end
	end


	# Process multi-body cases using cluster_orbits_dict for efficient processing
	for body in 2:nbody
		body_coupled_basislist = SortedCountingUniqueVector{Basis.CoupledBasis}()
		# Use cluster_orbits_dict to process clusters grouped by symmetry orbits
		if haskey(cluster_orbits_dict, body)
			for (orbit_index, orbit_clusters) in cluster_orbits_dict[body]
				# Process all clusters in this orbit
				# Clusters in the same orbit have the same projection matrix structure
				for atom_list::Vector{Int} in orbit_clusters
					# Get multiplicity from irreducible_cluster_dict
					count = irreducible_cluster_dict[body].counts[atom_list]
					sorted_atom_list = sort(atom_list)
					cb_list::Vector{Basis.CoupledBasis} = listup_coupled_basislist(
						sorted_atom_list,
						bodyn_lsum[body];
						isotropy = isotropy,
					)
					for cb::Basis.CoupledBasis in cb_list
						push_unique_coupled_basis!(body_coupled_basislist, cb, count, symmetry)
					end
				end
			end
		end
		for cb in body_coupled_basislist
			push!(result_coupled_basislist, cb, body_coupled_basislist.counts[cb])
		end
	end
	return result_coupled_basislist
end

"""
	listup_coupled_basislist(atom_list, lsum; isotropy=false) -> Vector{Basis.CoupledBasis}

List up all coupled angular momentum basis functions for a given atom list and maximum angular momentum sum.

# Arguments
- `atom_list::Vector{<:Integer}`: List of atom indices
- `lsum::Integer`: Maximum sum of angular momentum values
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`

# Returns
- `Vector{Basis.CoupledBasis}`: List of coupled basis functions

# Note
- Skips l values less than the number of atoms or odd l values
- Generates all possible combinations of angular momenta that sum to `lsum` or less
"""
function listup_coupled_basislist(
	atom_list::Vector{<:Integer},
	lsum::Integer;
	isotropy::Bool = false,
)::Vector{Basis.CoupledBasis}
	result_basislist = Vector{Basis.CoupledBasis}()
	for l in 2:lsum
		# Skip if l is less than number of atoms (each atom needs at least l=1)
		# or if l is odd (prohibited by time-reversal symmetry)
		if l < length(atom_list) || isodd(l)
			continue
		end
		l_list = Combinat.compositions(l, length(atom_list); min = 1)
		for l_vec::Vector{Int} in l_list
			cb_list =
				tesseral_coupled_bases_from_tesseral_bases(l_vec, atom_list; isotropy = isotropy)
			append!(result_basislist, cb_list)
		end
	end
	return result_basislist
end

"""
Check if the product of two `CoupledBasis` objects is obviously zero.
This is a quick check for the step constructing projection matrix.
"""
function is_obviously_zero_coupled_basis_product(
	cb1::Basis.CoupledBasis,
	cb2::Basis.CoupledBasis,
)::Bool
	if cb1.Lf != cb2.Lf
		return true
	end
	if cb1.ls != cb2.ls
		return true
	end
	if cb1.atoms != cb2.atoms
		return true
	end
	return false
end

"""
	tensor_inner_product(tensor1::AbstractArray{T, N}, tensor2::AbstractArray{T, N}) where {T, N} -> Float64

Compute the inner product of two tensors.

Both tensors must have the same element type `T` and the same number of dimensions `N`.

# Arguments
- `tensor1::AbstractArray{T, N}`: First tensor
- `tensor2::AbstractArray{T, N}`: Second tensor (must have same element type and dimensions as tensor1)

# Returns
- `Float64`: Inner product of the two tensors

# Examples
```julia
t1 = [1.0 2.0; 3.0 4.0]
t2 = [5.0 6.0; 7.0 8.0]
tensor_inner_product(t1, t2)  # 70.0
```
"""
function tensor_inner_product(
	tensor1::AbstractArray{T, N},
	tensor2::AbstractArray{T, N},
) where {T, N}
	return sum(conj.(tensor1) .* tensor2)
end


"""
	collect_cluster_atoms(coupled_basislist) -> Set{Vector{Int}}

Collect all unique atom lists from a list of coupled basis functions.

# Arguments
- `coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis}`: List of coupled basis functions

# Returns
- `Set{Vector{Int}}`: Set of unique atom lists (as vectors)
"""
function collect_cluster_atoms(
	coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis},
)::Set{Vector{Int}}
	result_set = Set{Vector{Int}}()
	for cb in coupled_basislist
		push!(result_set, cb.atoms)
	end
	return result_set
end

"""
	find_translation_atoms(atom_list, cluster_atoms, symmetry) -> Vector{Int}

Find the translationally equivalent atom list in the primitive cell for a given atom list.

# Arguments
- `atom_list::Vector{<:Integer}`: Input atom list
- `cluster_atoms::Set{Vector{Int}}`: Set of atom lists from cluster basis functions
- `symmetry::Symmetry`: Symmetry information containing translation mappings

# Returns
- `Vector{Int}`: Translationally equivalent atom list in the primitive cell

# Throws
- `ErrorException`: If no matching atom list is found or multiple matches are found
"""
function find_translation_atoms(
	atom_list::Vector{<:Integer},
	cluster_atoms::Set{Vector{Int}},
	symmetry::Symmetry,
)::Vector{Int}
	matched_set = Set{Vector{Int}}()
	for sym_tran in symmetry.symnum_translation
		atom_list_shifted = Vector{Int}([symmetry.map_sym[atom, sym_tran] for atom in atom_list])
		for itran in symmetry.symnum_translation
			atom_list_shifted_shifted =
				Vector{Int}([symmetry.map_sym[atom, itran] for atom in atom_list_shifted])
			sorted_atom_list_shifted_shifted = sort(atom_list_shifted_shifted)
			if sorted_atom_list_shifted_shifted in cluster_atoms
				push!(matched_set, atom_list_shifted_shifted)
			end
		end
	end
	if isempty(matched_set)
		error("Failed to find translation atoms")
	elseif length(matched_set) > 1
		sorted_matched_set = Set([sort(vec) for vec in matched_set])
		if length(sorted_matched_set) > 1
			@show atom_list
			@show cluster_atoms
			@show matched_set
			error("Multiple translation atoms found")
		end
	end
	return first(matched_set)
end


"""
	projection_matrix_coupled_basis(coupled_basislist, symmetry) -> Matrix{Float64}

Construct the projection matrix for coupled basis functions under symmetry operations.

The projection matrix is computed by averaging over all symmetry operations (including time-reversal)
and projects onto the subspace of basis functions that are invariant under the symmetry group.

# Arguments
- `coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis}`: List of coupled basis functions
- `symmetry::Symmetry`: Symmetry information containing all symmetry operations

# Returns
- `Matrix{Float64}`: Projection matrix of size (nbasis * submatrix_dim) × (nbasis * submatrix_dim),
  where submatrix_dim = 2*Lf + 1 for the total angular momentum Lf

# Note
- The matrix is averaged over all symmetry operations and time-reversal symmetry
- Each basis function contributes a submatrix of size (2*Lf + 1) × (2*Lf + 1)
"""
function projection_matrix_coupled_basis(
	coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis},
	symmetry::Symmetry,
)::Matrix{Float64}
	Lf = coupled_basislist[1].Lf
	submatrix_dim = 2 * Lf + 1
	cluster_atoms::Set{Vector{Int}} = collect_cluster_atoms(coupled_basislist)
	nbasis = length(coupled_basislist)
	full_matrix_dim = nbasis * submatrix_dim

	projection_mat = zeros(Float64, full_matrix_dim, full_matrix_dim)
	for (n, symop) in enumerate(symmetry.symdata), time_rev_sym in [false, true]
		# Calculate rotation matrix
		is_proper = symop.is_proper
		rotmat = is_proper ? symop.rotation_cart : -1 * symop.rotation_cart

		submat_in_mat::Matrix{Union{Matrix{Float64}, Nothing}} =
			Matrix{Union{Matrix{Float64}, Nothing}}(undef, nbasis, nbasis)
		for (i, cb1) in enumerate(coupled_basislist)
			atoms_shifted_list = [symmetry.map_sym[atom, n] for atom in cb1.atoms]
			primitive_atoms = find_translation_atoms(atoms_shifted_list, cluster_atoms, symmetry)
			reordered_cb = reorder_atoms(cb1, primitive_atoms)
			for (j, cb2) in enumerate(coupled_basislist)
				if is_obviously_zero_coupled_basis_product(reordered_cb, cb2)
					submat_in_mat[j, i] = nothing
					continue
				else
					phase = tensor_inner_product(cb2.coeff_tensor, reordered_cb.coeff_tensor) / (2*Lf+1)
					# Calculate rotation matrix for this basis
					rot_mat = Δl(reordered_cb.Lf, rotmat2euler(rotmat)...)
					# Apply time reversal symmetry multiplier if needed
					if time_rev_sym
						total_l = sum(reordered_cb.ls)
						multiplier = (-1)^total_l
						rot_mat = multiplier * rot_mat
					end
					submat_in_mat[j, i] = rot_mat * phase
				end
			end
		end

		# Expand submat_in_mat into a temporary Matrix{Float64}
		temp_projection_mat = zeros(Float64, full_matrix_dim, full_matrix_dim)
		for j in 1:nbasis, i in 1:nbasis
			submat = submat_in_mat[j, i]
			if submat !== nothing
				row_range = ((j-1)*submatrix_dim+1):(j*submatrix_dim)
				col_range = ((i-1)*submatrix_dim+1):(i*submatrix_dim)
				temp_projection_mat[row_range, col_range] = submat
			end
			# If submat is nothing, the corresponding block remains zero
		end

		# Accumulate the projection matrix
		projection_mat += temp_projection_mat
	end

	# Average over all symmetry operations (2 for time reversal)
	return projection_mat ./ (2 * symmetry.nsym)
end


"""
	is_translationally_equivalent_coupled_basis(
		cb1::Basis.CoupledBasis,
		cb2::Basis.CoupledBasis,
		symmetry::Symmetry,
	) -> Bool

Check if two `CoupledBasis` objects are translationally equivalent.

Two `CoupledBasis` objects are translationally equivalent if:
- They are physically equivalent (same `Lf`, `Lseq`, and `(atom, l)` pairs)
- Their atom lists are related by a translation operation in the supercell

This function checks if the atom lists can be mapped to each other via translation operations
defined in `symmetry.symnum_translation`.

# Arguments
- `cb1::Basis.CoupledBasis`: First `CoupledBasis` to compare
- `cb2::Basis.CoupledBasis`: Second `CoupledBasis` to compare
- `symmetry::Symmetry`: Symmetry information containing translation mappings

# Returns
- `Bool`: `true` if the `CoupledBasis` objects are translationally equivalent, `false` otherwise
"""
function is_translationally_equivalent_coupled_basis(
	cb1::Basis.CoupledBasis,
	cb2::Basis.CoupledBasis,
	symmetry::Symmetry,
)::Bool
	# Different number of sites
	length(cb1.ls) != length(cb2.ls) && return false

	# Different Lf
	cb1.Lf != cb2.Lf && return false

	# Different Lseq
	cb1.Lseq != cb2.Lseq && return false

	# Different ls values (as multisets) means different basis functions
	ls1_sorted = sort(cb1.ls)
	ls2_sorted = sort(cb2.ls)
	if ls1_sorted != ls2_sorted
		return false
	end

	# Check if (atom, l) pairs match as multisets (required for physical equivalence)
	atom_l_pairs1 = collect(zip(cb1.atoms, cb1.ls))
	atom_l_pairs2 = collect(zip(cb2.atoms, cb2.ls))
	if sort(atom_l_pairs1) != sort(atom_l_pairs2)
		return false
	end

	# Check if coeff_tensor has the same size and values
	if size(cb1.coeff_tensor) != size(cb2.coeff_tensor)
		return false
	end
	if !isapprox(cb1.coeff_tensor, cb2.coeff_tensor, atol = 1e-10)
		return false
	end

	# Check if atom lists are translationally equivalent
	atom_list1 = cb1.atoms
	atom_list2 = cb2.atoms

	# Early return if atom lists are the same
	if atom_list1 == atom_list2
		return false
	end

	# Early return if first atom is the same
	# because this function is intended to be used for different first atoms but translationally equivalent clusters
	if atom_list1[1] == atom_list2[1]
		return false
	end

	# Check translation operations
	for itran in symmetry.symnum_translation
		# Method 1: Apply forward translation (map_sym) to atom_list1
		atom_list1_shifted = [symmetry.map_sym[atom, itran] for atom in atom_list1]
		# Sort both lists to compare as multisets (order doesn't matter)
		if sort(atom_list1_shifted) == sort(atom_list2)
			return true
		end

		# Method 2: Apply inverse translation (map_sym_inv) to atom_list1
		atom_list1_shifted = [symmetry.map_sym_inv[atom, itran] for atom in atom_list1]
		# Sort both lists to compare as multisets (order doesn't matter)
		if sort(atom_list1_shifted) == sort(atom_list2)
			return true
		end
	end

	return false
end

"""
	push_unique_coupled_basis!(target, cb, count, symmetry)

Add a `CoupledBasis` to the target list only if it is not translationally equivalent
to any existing `CoupledBasis` in the list.

# Arguments
- `target::SortedCountingUniqueVector{Basis.CoupledBasis}`: Target list to add to
- `cb::Basis.CoupledBasis`: Coupled basis function to add
- `count::Integer`: Multiplicity count for the basis function
- `symmetry::Symmetry`: Symmetry information for checking translational equivalence

# Note
- If an equivalent basis function is found, the function returns without adding
- Otherwise, adds the basis function with the given count
"""
function push_unique_coupled_basis!(
	target::SortedCountingUniqueVector{Basis.CoupledBasis},
	cb::Basis.CoupledBasis,
	count::Integer,
	symmetry::Symmetry,
)
	# Quick check: sum of l values
	lsum_cb = sum(collect(cb.ls))
	for existing_cb in target
		lsum_existing = sum(collect(existing_cb.ls))
		if lsum_cb != lsum_existing
			continue
		end
		# Check if physically equivalent (same Lf, Lseq, (atom, l) pairs, and coeff_list)
		# This checks for exact matches (same atoms, same ls order)
		if is_translationally_equivalent_coupled_basis(cb, existing_cb, symmetry)
			return
		end
	end
	# No equivalent LinearCombo found: add with the given count
	push!(target, cb, count)
end


"""
	construct_basislist(structure, symmetry, cluster, body1_lmax, bodyn_lsum, nbody) -> SortedCountingUniqueVector{SHProduct}

Construct a list of basis functions using `SHProduct` objects.

This function generates basis functions for atomic interactions, handling both 1-body and multi-body cases.
It is a legacy function that uses `SHProduct` instead of `CoupledBasis`.

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information for atomic interactions
- `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions
- `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
- `nbody::Integer`: Maximum number of bodies in interactions

# Returns
- `SortedCountingUniqueVector{SHProduct}`: List of basis functions with multiplicities

# Note
- This function is kept for backward compatibility
- New code should use `construct_coupled_basislist` instead
"""
function construct_basislist(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	body1_lmax::Vector{Int},
	bodyn_lsum::OffsetArray{Int, 1},
	nbody::Integer,
)::SortedCountingUniqueVector{SHProduct}

	result_basislist = SortedCountingUniqueVector{SHProduct}()
	cluster_dict::Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}} =
		cluster.cluster_dict

	# Handle 1-body case
	for iat in symmetry.atoms_in_prim
		lmax = body1_lmax[structure.supercell.kd_int_list[iat]]
		for l in 2:lmax # skip l = 1 because it is prohibited by the time-reversal symmetry
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
		body_basislist = SortedCountingUniqueVector{SHProduct}()
		for prim_atom_sc in symmetry.atoms_in_prim
			cuv::CountingUniqueVector{Vector{Int}} = cluster_dict[body][prim_atom_sc]
			for atom_list::Vector{Int} in cuv
				count = cuv.counts[atom_list]
				shp_list::Vector{SHProduct} =
					listup_basislist(atom_list, bodyn_lsum[body])
				for shp::SHProduct in shp_list
					push_unique_body!(body_basislist, shp, count, symmetry)
				end
			end
		end
		for basis in body_basislist
			push!(result_basislist, basis, body_basislist.counts[basis])
		end
	end

	return result_basislist
end

"""
	push_unique_body!(target, shp, count, symmetry)

Add a `SHProduct` to the target list only if it is not translationally equivalent
to any existing `SHProduct` in the list.

# Arguments
- `target::SortedCountingUniqueVector{SHProduct}`: Target list to add to
- `shp::SHProduct`: Basis function to add
- `count::Integer`: Multiplicity count for the basis function
- `symmetry::Symmetry`: Symmetry information for checking translational equivalence

# Note
- If an equivalent basis function is found, the function returns without adding
- Otherwise, adds the basis function with the given count
"""
function push_unique_body!(
	target::SortedCountingUniqueVector{SHProduct},
	shp::SHProduct,
	count::Integer,
	symmetry::Symmetry,
)
	shp_sorted = sort(shp)
	for basis in target
		if sum([shsi.l for shsi in shp]) != sum([shsi.l for shsi in basis])
			continue
		end
		basis_sorted = sort(basis)
		if shp_sorted == basis_sorted || is_translationally_equiv_basis(basis, shp, symmetry)
			return
		end
	end
	push!(target, shp, count)
end


"""
	listup_basislist(atom_list, lsum) -> Vector{SHProduct}

List up all basis functions as `SHProduct` objects for a given atom list and maximum angular momentum sum.

# Arguments
- `atom_list::Vector{<:Integer}`: List of atom indices
- `lsum::Integer`: Maximum sum of angular momentum values

# Returns
- `Vector{SHProduct}`: List of basis functions

# Note
- Skips l values less than the number of atoms or odd l values
- Generates all possible combinations of angular momenta that sum to `lsum` or less
"""
function listup_basislist(
	atom_list::Vector{<:Integer},
	lsum::Integer,
)::Vector{SHProduct}
	result_basislist = Vector{SHProduct}()
	for l in 2:lsum
		if l < length(atom_list) || isodd(l)
			continue
		end
		l_list = Combinat.compositions(l, length(atom_list); min = 1)
		for l_vec::Vector{Int} in l_list
			shp_list = product_shsiteindex(atom_list, l_vec)
			append!(result_basislist, shp_list)
		end
	end
	return result_basislist
end




function map_atom_l_list(
	atom_l_list::AbstractVector{<:AbstractVector{<:Integer}},
	map_sym::AbstractMatrix{<:Integer},
	isym::Integer,
)::Vector{Vector{Int}}
	return [[map_sym[atom_l_vec[1], isym], atom_l_vec[2]] for atom_l_vec in atom_l_list]
end


function classify_basislist(
	basislist::AbstractVector{SHProduct},
	map_sym::AbstractMatrix{<:Integer},
)::AbstractDict{Int, SortedCountingUniqueVector}

	count = 1
	label_list = zeros(Int, size(basislist))
	for (idx, basis) in enumerate(basislist)
		if label_list[idx] != 0
			continue
		end
		atom_l_list_base = AtomicIndices.get_atom_l_list(basis)
		for isym in 1:size(map_sym, 2)
			mapped_list = map_atom_l_list(atom_l_list_base, map_sym, isym)
			sorted_mapped = sort(mapped_list)
			for (idx2, basis2) in enumerate(basislist)
				if label_list[idx2] == 0 &&
				   sort(AtomicIndices.get_atom_l_list(basis2)) == sorted_mapped
					label_list[idx2] = count
				end
			end
		end
		count += 1
	end

	dict = OrderedDict{Int, SortedCountingUniqueVector}()
	max_label = maximum(label_list)
	if max_label > 0
		for idx in 1:max_label
			dict[idx] = SortedCountingUniqueVector{SHProduct}()
		end
	end

	for (basis, label) in zip(basislist, label_list)
		push!(dict[label], basis, getcount(basislist, basis))
	end

	return dict
end

"""
	classify_coupled_basislist_test(
		coupled_basislist::AbstractVector{Basis.CoupledBasis},
	) -> Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}

Simplified classifier for `CoupledBasis` objects used in tests.

This version ignores spatial symmetry and groups basis functions solely by
interaction order (number of sites), final angular momentum `Lf`, sum of `ls`, and sorted `ls`.
It trades efficiency for robustness so that test fixtures can rely on deterministic
grouping without depending on symmetry metadata.

# Arguments
- `coupled_basislist::AbstractVector{Basis.CoupledBasis}`: List of CoupledBasis objects

# Returns
- `Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}`: Dictionary keyed by
  labels assigned per `(nbody, Lf, sum(ls), Tuple(sort(ls)...))` tuple
"""
function classify_coupled_basislist_test(
	coupled_basislist::AbstractVector{<:Basis.CoupledBasis},
)::Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}
	if isempty(coupled_basislist)
		return OrderedDict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}()
	end

	label_map = Dict{Any, Int}()
	label_list = Vector{Int}(undef, length(coupled_basislist))
	next_label = 0

	for (idx, cb) in enumerate(coupled_basislist)
		ls_sorted = Tuple(sort(cb.ls))
		key = (length(cb.atoms), cb.Lf, sum(cb.ls), ls_sorted)
		label = get(label_map, key, 0)
		if label == 0
			next_label += 1
			label = next_label
			label_map[key] = label
		end
		label_list[idx] = label
	end

	dict = OrderedDict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}()
	for label in 1:next_label
		dict[label] = SortedCountingUniqueVector{Basis.CoupledBasis}()
	end

	for (cb, label) in zip(coupled_basislist, label_list)
		if coupled_basislist isa SortedCountingUniqueVector
			count_val = get(coupled_basislist.counts, cb, 1)
			push!(dict[label], cb, count_val)
		else
			push!(dict[label], cb, 1)
		end
	end

	return dict
end


"""
Checks if two basis sets (SHProduct) are translationally equivalent.

# Arguments
- `basis_target::SHProduct`: The target basis set to check
- `basis_ref::SHProduct`: The reference basis set
- `symmetry::Symmetry`: Symmetry information containing translation mappings

# Returns
- `Bool`: `true` if the basis sets are translationally equivalent, `false` otherwise
"""
function is_translationally_equiv_basis(
	basis_target::SHProduct,
	basis_ref::SHProduct,
	symmetry::Symmetry,
)::Bool
	if length(basis_target) != length(basis_ref)
		return false
	end
	# Early return if the atomlist is the same including the order
	atom_list_target = [shsi.i for shsi in basis_target]
	atom_list_ref = [shsi.i for shsi in basis_ref]
	if atom_list_target == atom_list_ref
		return false
	end
	# Early return if the first atom is the same
	# because this function is intended to be used for different first atoms but translationally equivalent clusters
	if basis_target[1].i == basis_ref[1].i
		return false
	end

	for itran in symmetry.symnum_translation
		# Method 1: Apply forward translation (map_sym) to new_atom_list
		atom_list_target_shifted = [symmetry.map_sym[shsi.i, itran] for shsi in basis_target]
		forward_candidate = replace_atom(
			SHProduct([shsi for shsi in basis_target]),
			atom_list_target_shifted,
		)
		if sort(forward_candidate) == sort(basis_ref)
			return true
		end

		# Method 2: Apply inverse translation (map_sym_inv) to new_atom_list
		atom_list_target_shifted = [symmetry.map_sym_inv[shsi.i, itran] for shsi in basis_target]
		inverse_candidate = replace_atom(
			SHProduct([shsi for shsi in basis_target]),
			atom_list_target_shifted,
		)
		if sort(inverse_candidate) == sort(basis_ref)
			return true
		end

	end

	return false
end

"""
	is_proper_eigenvals(eigenval::AbstractVector; tol = 1e-8)::Bool

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
is_valid = is_proper_eigenvals(eigenvals)  # true

# With some tolerance
eigenvals = [0.0, 1.0 + 1e-9, 0.0, 1.0]
is_valid = is_proper_eigenvals(eigenvals)  # true

# Invalid eigenvalues
eigenvals = [0.0, 1.0, 0.5, 1.0]
is_valid = is_proper_eigenvals(eigenvals)  # false
```
"""
function is_proper_eigenvals(eigenval::AbstractVector; tol = 1e-8)::Bool
	return all(x -> isapprox(x, 0, atol = tol) || isapprox(x, 1, atol = tol), eigenval)
end




"""
	projection_matrix(basisdict, symmetry) -> Vector{Matrix{Float64}}

Construct projection matrices for each classification label in the `basisdict`.

This function computes the average projection matrix over all symmetry operations
(including time-reversal symmetry) for each group of basis functions.

# Arguments
- `basisdict::AbstractDict`: Dictionary mapping classification labels to groups of basis functions
- `symmetry::Symmetry`: Symmetry information containing all symmetry operations

# Returns
- `Vector{Matrix{Float64}}`: List of projection matrices, one for each classification label
"""
function projection_matrix(
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
			projection_mat_per_symop = proj_matrix_a_symop(
				basislist,
				symop,
				@view(symmetry.map_sym[:, n]),
				symmetry.map_sym_inv,
				symmetry.map_s2p,
				symmetry.atoms_in_prim,
				symmetry.symnum_translation,
				time_rev_sym,
			)
			local_projection_mat += projection_mat_per_symop
		end
		local_projection_mat = local_projection_mat ./ (2 * symmetry.nsym)
		local_projection_mat = hermitianpart(local_projection_mat)

		if !is_symmetric(local_projection_mat, tol = 1e-10)
			error("Projection matrix is not symmetric. index: $idx")
		end

		result_projections[idx] = local_projection_mat
	end
	return result_projections
end

"""
	proj_matrix_a_symop(basislist, symop, map_sym_per_symop, map_sym_inv, map_s2p, atoms_in_prim, symnum_translation, time_rev_sym) -> Matrix{Float64}

Compute the projection matrix for a single symmetry operation applied to a list of basis functions.

# Arguments
- `basislist::SortedCountingUniqueVector{SHProduct}`: List of basis functions
- `symop::SymmetryOperation`: The symmetry operation to apply
- `map_sym_per_symop::AbstractVector{<:Integer}`: Atom mapping for this symmetry operation
- `map_sym_inv::AbstractMatrix{<:Integer}`: Inverse atom mapping matrix
- `map_s2p::AbstractVector{<:Maps}`: Mapping from supercell to primitive cell
- `atoms_in_prim::AbstractVector{<:Integer}`: List of atoms in primitive cell
- `symnum_translation::AbstractVector{<:Integer}`: Translation symmetry numbers
- `time_rev_sym::Bool`: Whether to apply time-reversal symmetry

# Returns
- `Matrix{Float64}`: The projection matrix for this symmetry operation
"""
function proj_matrix_a_symop(
	basislist::SortedCountingUniqueVector{SHProduct},
	symop::SymmetryOperation,
	map_sym_per_symop::AbstractVector{<:Integer},
	map_sym_inv::AbstractMatrix{<:Integer},
	map_s2p::AbstractVector{<:Maps},
	atoms_in_prim::AbstractVector{<:Integer},
	symnum_translation::AbstractVector{<:Integer},
	time_rev_sym::Bool,
)::Matrix{Float64}
	projection_mat = zeros(Float64, length(basislist), length(basislist))

	for (j, basis_j::SHProduct) in enumerate(basislist)
		lco_j = operate_symop(
			basislist,
			basis_j,
			symop,
			map_sym_per_symop,
			map_sym_inv,
			map_s2p,
			atoms_in_prim,
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
		display(basislist)
		display(symop)
		@assert false "Projection matrix is zero matrix"
	end
	# if !is_unitary(projection_mat, tol = 1e-10)
	# 	error("Projection matrix is not unitary")
	# end
	return projection_mat
end

"""
	operate_symop(basislist, basis, symop, map_sym_per_symop, map_sym_inv, map_s2p, atoms_in_prim, symnum_translation, time_rev_sym) -> AtomicIndices.LinearCombo

Apply a symmetry operation to a basis function and return it as a linear combination.

# Arguments
- `basislist::SortedCountingUniqueVector{SHProduct}`: List of basis functions
- `basis::SHProduct`: The basis function to transform
- `symop::SymmetryOperation`: The symmetry operation to apply
- `map_sym_per_symop::AbstractVector{<:Integer}`: Atom mapping for this symmetry operation
- `map_sym_inv::AbstractMatrix{<:Integer}`: Inverse atom mapping matrix
- `map_s2p::AbstractVector{<:Maps}`: Mapping from supercell to primitive cell
- `atoms_in_prim::AbstractVector{<:Integer}`: List of atoms in primitive cell
- `symnum_translation::AbstractVector{<:Integer}`: Translation symmetry numbers
- `time_rev_sym::Bool`: Whether to apply time-reversal symmetry

# Returns
- `AtomicIndices.LinearCombo`: The transformed basis function as a linear combination

# Throws
- `ErrorException`: If no corresponding basis function is found in the primitive cell
"""
function operate_symop(
	basislist::SortedCountingUniqueVector{SHProduct},
	basis::SHProduct,
	symop::SymmetryOperation,
	map_sym_per_symop::AbstractVector{<:Integer},
	map_sym_inv::AbstractMatrix{<:Integer},
	map_s2p::AbstractVector{<:Maps},
	atoms_in_prim::AbstractVector{<:Integer},
	symnum_translation::AbstractVector{<:Integer},
	time_rev_sym::Bool,
)::AtomicIndices.LinearCombo
	# Apply the symmetry operation to atoms
	new_atom_list = [map_sym_per_symop[shsi.i] for shsi in basis]

	# Shift the new atom list to the primitive cell
	new_basis_found = SHProduct()
	translated_atom_list = similar(new_atom_list)

	found = false
	for i in eachindex(new_atom_list)
		ref_atom = new_atom_list[i]
		ref_atom_in_prim = atoms_in_prim[map_s2p[ref_atom].atom]
		corresponding_translation = symnum_translation[map_s2p[ref_atom].translation]
		atom_translated = similar(new_atom_list)
		atom_translated[i] = ref_atom_in_prim
		for (n, new_atom) in enumerate(new_atom_list)
			if n == i
				continue
			end
			atom_translated[n] = map_sym_inv[new_atom, corresponding_translation]
		end
		new_basis_candidate = replace_atom(basis, atom_translated)
		for basis_i in basislist
			if sort(new_basis_candidate) == sort(basis_i)
				found = true
				translated_atom_list = [shsi.i for shsi in basis_i]
				new_basis_found = basis_i
				break
			end
		end
		if found
			break
		end
	end
	if !found
		error("Failed to find corresponding basis in the primitive cell.")
	end

	idx = corresponding_idx(new_basis_found)

	is_proper = symop.is_proper
	multiplier = 1.0
	if is_proper
		rotmat = symop.rotation_cart
	else
		rotmat = -1 * symop.rotation_cart
	end
	if time_rev_sym
		multiplier *= (-1)^(sum([shsi.l for shsi in new_basis_found]))
	end
	rotation_list = [Δl(shsi.l, rotmat2euler(rotmat)...) for shsi in new_basis_found]
	if length(new_basis_found) == 1
		rotmat_kron = multiplier * rotation_list[begin]
	else
		rotmat_kron = multiplier * kron(rotation_list...)
	end

	l_list = [shsi.l for shsi in new_basis_found]

	shp_list::Vector{SHProduct} = product_shsiteindex(translated_atom_list, l_list)
	coeffs = rotmat_kron[:, idx]
	return AtomicIndices.LinearCombo(shp_list, coeffs)
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
	# This matches the order of kron(rotation_list...) used in operate_symop
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

"""
	is_unitary(mat::AbstractMatrix{<:Real}; tol::Float64 = 1e-10) -> Bool

Check if a matrix is unitary (i.e., UᵀU ≈ I and UUᵀ ≈ I within tolerance).

# Arguments
- `mat::AbstractMatrix{<:Real}`: Matrix to check
- `tol::Float64`: Tolerance for floating-point comparison (default: 1e-10)

# Returns
- `Bool`: `true` if the matrix is unitary within the tolerance, `false` otherwise
"""
function is_unitary(mat::AbstractMatrix{<:Real}; tol::Float64 = 1e-10)::Bool
	size(mat, 1) == size(mat, 2) || return false
	return isapprox(mat' * mat, I, atol = tol) && isapprox(mat * mat', I, atol = tol)
end


"""
	filter_basisdict(
		basisdict::Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}},
		symmetry::Symmetry,
	) -> Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}

Filter out basis entries that correspond to clusters where all atoms map to the same
primitive atom, have odd total angular momentum Lf, and have multiplicity greater than 1.

# Description
For clusters consisting of atom sites located at cell boundaries that are connected by
pure translations, the symmetry-adapted linear combinations (SALCs) vanish when the total
angular momentum Lf is odd. This function filters out such basis entries beforehand to
avoid unnecessary computation in the SALC construction process.

A cluster is filtered if all of the following conditions are satisfied:
1. All atoms in the cluster map to the same atom in the primitive cell (i.e., they are
   translationally equivalent)
2. The total angular momentum Lf of the coupled basis is odd
3. The multiplicity (count) of the basis is greater than 1

The third condition ensures that only clusters at cell boundaries are filtered, as
translationally equivalent sites at cell boundaries always have multiplicity > 1 due to
periodic boundary conditions.

# Arguments
- `basisdict::Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}`: Dictionary of
  classified coupled basis sets, where keys are labels and values are lists of basis
  functions belonging to the same symmetry class. The `counts` field of each
  `SortedCountingUniqueVector` stores the multiplicity of each basis.
- `symmetry::Symmetry`: Symmetry information containing the mapping from supercell atoms
  to primitive cell atoms (`symmetry.map_s2p`)

# Returns
- `Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}`: Filtered dictionary with
  renumbered labels. Entries that satisfy all filtering conditions are removed, and
  remaining entries are assigned new consecutive labels starting from 1.

# Notes
- The function checks only the first basis in each `SortedCountingUniqueVector` to determine
  the filtering condition, as all bases in the same list share the same Lf value.
- The multiplicity is accessed via `basislist.counts[first_basis]`, which counts how many
  times the same basis appears due to translational symmetry.
- The input dictionary keys are sorted before processing to ensure deterministic output.
"""
function filter_basisdict(
	basisdict::Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}},
	symmetry::Symmetry,
)::Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}
	result_basisdict = Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}}()
	new_label = 1
	for label in sort!(collect(keys(basisdict)))
		basislist::SortedCountingUniqueVector{Basis.CoupledBasis} = basisdict[label]
		first_basis = basislist[begin]
		atom_list = first_basis.atoms
		atom_list_in_prim = [symmetry.map_s2p[atom].atom for atom in atom_list]
		unique_atom_list_in_prim = unique(atom_list_in_prim)
		if length(unique_atom_list_in_prim) == 1 && isodd(first_basis.Lf) && basislist.counts[first_basis] > 1
			continue
		end

		result_basisdict[new_label] = deepcopy(basislist)
		new_label += 1
	end
	return result_basisdict
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

"""
	print_basisset_stdout(salc_list)

Print symmetry-adapted basis functions to stdout.

# Arguments
- `salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}}`: List of SALC groups, where each group is a vector of `CoupledBasis_with_coefficient` objects
"""
function print_basisset_stdout(
	salc_list::AbstractVector{Vector{Basis.CoupledBasis_with_coefficient}},
)
	println(" Number of symmetry-adapted basis functions: $(length(salc_list))\n")
	println(" List of symmetry-adapted basis functions:")
	for (i, key_group) in enumerate(salc_list)
		println(" $i-th salc")
		println(" number of terms: $(length(key_group))")
		for cbc in key_group
			# Format coefficient vector as space-separated string
			coeff_str = join([@sprintf("% 15.10f", x) for x in cbc.coefficient], " ")
			line_str = @sprintf("%2d  [%s]  %s", cbc.multiplicity, coeff_str, cbc)
			println(line_str)
		end
		println("")
	end
end


end # module BasisSets
