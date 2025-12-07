"""
This module provides functionality for handling basis sets in crystal structure calculations.
It includes tools for constructing, classifying, and manipulating basis functions that are
adapted to the symmetry of the crystal structure.
"""
module BasisSets

using Combinat
using DataStructures
using LinearAlgebra
using OffsetArrays
using Printf
using StaticArrays

using ..CountingContainer
using ..SortedContainer
using ..AtomCells
using ..AtomicIndices
using ..ConfigParser
using ..Structures
using ..Symmetries
using ..Clusters
using ..SALCs
using ..RotationMatrix
using ..Basis

export BasisSet

"""
	SALC_LinearCombo

A structure representing a symmetry-adapted linear combination of `LinearCombo` objects.

# Fields
- `basisset::Vector{Basis.LinearCombo}`: The basis set of `LinearCombo` objects
- `coeffs::Vector{Float64}`: The coefficients of the linear combination
- `multiplicity::Vector{Int}`: The multiplicity of each basis element

# Constructors
- `SALC_LinearCombo(basisset::Vector{Basis.LinearCombo}, coeffs::Vector{Float64}, multiplicity::Vector{Int})`:
  Create a SALC_LinearCombo from a basis set, coefficients, and multiplicity
- `SALC_LinearCombo(basislist::SortedCountingUniqueVector{Basis.LinearCombo}, coeffs::Vector{<:Real})`:
  Create a SALC_LinearCombo from a basis list and coefficients
"""
# struct SALC_LinearCombo
# 	basisset::Vector{Basis.LinearCombo}
# 	coeffs::Vector{Float64}
# 	multiplicity::Vector{Int}

# 	function SALC_LinearCombo(
# 		basisset::Vector{Basis.LinearCombo},
# 		coeffs::Vector{Float64},
# 		multiplicity::Vector{Int},
# 	)
# 		length(basisset) == length(coeffs) == length(multiplicity) ||
# 			throw(ArgumentError("All vectors must have the same length"))
# 		all(x -> x > 0, multiplicity) ||
# 			throw(ArgumentError("Multiplicity must be positive"))
# 		new(basisset, coeffs, multiplicity)
# 	end
# end

"""
	SALC_LinearCombo(
		basislist::SortedCountingUniqueVector{Basis.LinearCombo},
		coeffs::Vector{<:Real},
	) -> SALC_LinearCombo

Create a SALC_LinearCombo from a basis list and coefficients.

# Arguments
- `basislist::SortedCountingUniqueVector{Basis.LinearCombo}`: The basis list with counts
- `coeffs::Vector{<:Real}`: The coefficients for each basis element

# Returns
- `SALC_LinearCombo`: A new symmetry-adapted linear combination

# Throws
- `ArgumentError` if the lengths of `basislist` and `coeffs` differ
- `ArgumentError` if the resulting coefficient vector has zero norm
"""
# function SALC_LinearCombo(
# 	basislist::SortedCountingUniqueVector{Basis.LinearCombo},
# 	coeffs::Vector{<:Real},
# )
# 	length(basislist) == length(coeffs) ||
# 		throw(ArgumentError("The length of basislist and coeffs must be the same"))

# 	result_basisset = Vector{Basis.LinearCombo}()
# 	result_coeffs = Vector{Float64}()
# 	result_multiplicity = Vector{Int}()

# 	for (idx, basis) in enumerate(basislist)
# 		count = basislist.counts[basis]
# 		if !isapprox(coeffs[idx], 0.0, atol = 1e-8)
# 			push!(result_basisset, basis)
# 			push!(result_coeffs, coeffs[idx])
# 			push!(result_multiplicity, count)
# 		end
# 	end

# 	# normalize coefficient vector
# 	norm_coeffs = norm(result_coeffs)
# 	if isapprox(norm_coeffs, 0.0, atol = 1e-8)
# 		throw(ArgumentError("The norm of the coefficient vector is zero."))
# 	end
# 	result_coeffs ./= norm_coeffs

# 	return SALC_LinearCombo(result_basisset, result_coeffs, result_multiplicity)
# end

"""
	show(io::IO, salc::SALC_LinearCombo)

Display a SALC_LinearCombo in a human-readable format.

# Arguments
- `io::IO`: The output stream
- `salc::SALC_LinearCombo`: The SALC_LinearCombo to display

# Output Format
```
number of terms: N
M  COEFFICIENT  LinearCombo(...)
```
where:
- N is the number of terms
- M is the multiplicity
- COEFFICIENT is the normalized coefficient
- LinearCombo(...) is the basis element
"""
# function Base.show(io::IO, salc::SALC_LinearCombo)
# 	println(io, "number of terms: ", length(salc.basisset))
# 	for (basis, coeff, multiplicity) in zip(salc.basisset, salc.coeffs, salc.multiplicity)
# 		println(io, @sprintf("%2d  % 15.10f  %s", multiplicity, coeff, basis))
# 	end
# end

"""
	print_basisset_stdout_linearcombo(salc_list::AbstractVector{<:SALC_LinearCombo})

Print symmetry-adapted basis functions constructed from `LinearCombo` objects.

# Arguments
- `salc_list::AbstractVector{<:SALC_LinearCombo}`: List of SALC_LinearCombo objects

# Output Format
```
 Number of symmetry-adapted basis functions: N

 List of symmetry-adapted basis functions:
 # multiplicity  coefficient  LinearCombo
 1-th salc
 ...
```
"""
# function print_basisset_stdout_linearcombo(salc_list::AbstractVector{<:SALC_LinearCombo})
# 	println(" Number of symmetry-adapted basis functions: $(length(salc_list))\n")
# 	println(" List of symmetry-adapted basis functions:")
# 	println(" # multiplicity  coefficient  LinearCombo")
# 	for (i, salc) in enumerate(salc_list)
# 		println(" $i-th salc")
# 		display(salc)
# 	end
# 	println("")
# end

"""
	struct BasisSet

Represents a set of basis functions for atomic interactions in a crystal structure.
This structure is used to store and manage basis functions that are adapted to the symmetry of the crystal.

# Fields
- `salc_linearcombo_list::Vector{SALC_LinearCombo}`: List of symmetry-adapted linear combinations of `LinearCombo` objects

# Constructors
	BasisSet(structure, symmetry, cluster, body1_lmax, bodyn_lsum, nbody; verbosity=true)
	BasisSet(structure, symmetry, cluster, config; verbosity=true)

Constructs a new `BasisSet` instance for atomic interactions in a crystal structure.

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information for atomic interactions
- `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions for each atomic species
- `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
- `nbody::Integer`: Maximum number of bodies in interactions
- `verbosity::Bool`: Whether to print progress information (default: `true`)

# Returns
- `BasisSet`: A new basis set instance containing symmetry-adapted linear combinations (SALCs) of `LinearCombo` objects

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
1. Constructs the tesseral basis list using `LinearCombo` objects
2. Classifies basis functions by symmetry operations
3. Constructs projection matrices for each symmetry label
4. Generates symmetry-adapted linear combinations (SALCs) of `LinearCombo` objects
"""
struct BasisSet
	# salc_linearcombo_list::Vector{SALC_LinearCombo}
	coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis}
	salc_list::Vector{Vector{Basis.CoupledBasis_with_coefficient}}

	function BasisSet(
		coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis},
		salc_list::Vector{Vector{Basis.CoupledBasis_with_coefficient}},
	)
		return new(coupled_basislist, salc_list)
	end
end

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

	print("Constructing and classifying coupled basis list...")
	classified_coupled_basisdict::Dict{Int, SortedCountingUniqueVector{Basis.CoupledBasis}} =
		construct_and_classify_coupled_basislist(
			structure,
			symmetry,
			cluster,
			body1_lmax,
			bodyn_lsum,
			nbody,
		)
	println(" Done.")

	if verbosity
		elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds
		println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
		println("-------------------------------------------------------------------")
	end

	salc_list = Vector{Vector{Basis.CoupledBasis_with_coefficient}}()

	print("Constructing projection matrix...")
	keys_list = sort(collect(keys(classified_coupled_basisdict)))
	for key in keys_list
		coupled_basislist = classified_coupled_basisdict[key]
		projection_mat =
			projection_matrix_coupled_basis(coupled_basislist, symmetry)
		println("key : $key")
		# projection_mat is real symmetric → use Hermitian + eigen! to save memory
		h_projection = Hermitian(projection_mat)
		# Free projection_mat memory before eigen decomposition
		projection_mat = nothing
		eigenvals, eigenvecs = eigen!(h_projection)
		eigenvals = real.(round.(eigenvals, digits = 6))
		eigenvecs = round.(eigenvecs .* (abs.(eigenvecs) .≥ 1e-8), digits = 10)
		if !is_proper_eigenvals(eigenvals)
			@show eigenvals
			error("Critical error: Eigenvalues must be either 0 or 1. index: $key")
		end
		Lf = coupled_basislist[1].Lf
		submatrix_dim = 2 * Lf + 1
		nbasis = length(coupled_basislist)
		
		# Create a list for this key
		key_salc_list = Vector{Basis.CoupledBasis_with_coefficient}()
		
		for idx_eigenval in findall(x -> isapprox(x, 1.0, atol = 1e-8), eigenvals)
			eigenvec = eigenvecs[:, idx_eigenval]
			eigenvec = real.(eigenvec)
			eigenvec = round.(eigenvec .* (abs.(eigenvec) .≥ 1e-8), digits = 10)
			eigenvec = eigenvec / norm(eigenvec)

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

				# Create CoupledBasis_with_coefficient and add to key_salc_list
				cbc = Basis.CoupledBasis_with_coefficient(cb, coefficient, multiplicity)
				push!(key_salc_list, cbc)
			end
		end
		
		# Add the list for this key to salc_list only if it's not empty
		if !isempty(key_salc_list)
			push!(salc_list, key_salc_list)
		end
		# Free large matrices after processing each key
		h_projection = nothing
		eigenvecs = nothing
	end
	println(" Done.")
	if verbosity
		elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds
		println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
		println("-------------------------------------------------------------------")
	end

	for salc in salc_list
		println("salc : $salc")
	end



	# Validate input parameters
	# Construct basis list
	# basislist consists of all possible basis functions which is the product of spherical harmonics.
	# if verbosity
	# 	print("Constructing basis list...")
	# end
	# tesseral_basislist::SortedCountingUniqueVector{Basis.LinearCombo} =
	# 	construct_tesseral_basislist(
	# 		structure,
	# 		symmetry,
	# 		cluster,
	# 		body1_lmax,
	# 		bodyn_lsum,
	# 		nbody,
	# 	)

	# Classify tesseral_basislist by symmetry operations
	# classified_tesseral_basisdict::Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}} =
	# 	classify_tesseral_basislist(tesseral_basislist, symmetry.map_sym)

	# classified_tesseral_basisdict = classify_tesseral_basislist_test(tesseral_basislist)
	# keys_list = sort(collect(keys(classified_tesseral_basisdict)))
	# for key in keys_list
	# 	println("key : $key")
	# 	lc_list = classified_tesseral_basisdict[key]
	# 	for lc in lc_list
	# 		println("  lc : $lc")
	# 	end
	# end

	# projection_list = projection_matrix_linearcombo(classified_tesseral_basisdict, symmetry)
	# salc_linearcombo_list = Vector{SALC_LinearCombo}()
	# for (idx, projection_mat) in enumerate(projection_list)
	# 	basislist = classified_tesseral_basisdict[idx]
	# 	eigenvals, eigenvecs = eigen(projection_mat)
	# 	eigenvals = real.(round.(eigenvals, digits = 6))
	# 	eigenvecs = round.(eigenvecs .* (abs.(eigenvecs) .≥ 1e-8), digits = 10)
	# 	if !is_proper_eigenvals(eigenvals)
	# 		println(eigenvals)
	# 		@warn "Critical error: Eigenvalues must be either 0 or 1. index: $idx"
	# 	end
	# 	for idx_eigenval in findall(x -> isapprox(x, 1.0, atol = 1e-8), eigenvals)
	# 		eigenvec = eigenvecs[:, idx_eigenval]
	# 		eigenvec = flip_vector_if_negative_sum(eigenvec)
	# 		eigenvec = round.(eigenvec .* (abs.(eigenvec) .≥ 1e-8), digits = 10)
	# 		push!(salc_linearcombo_list, SALC_LinearCombo(basislist, eigenvec / norm(eigenvec)))
	# 	end
	# end

	# if verbosity
	# 	print_basisset_stdout_linearcombo(salc_linearcombo_list)
	# 	elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds
	# 	println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
	# 	println("-------------------------------------------------------------------")
	# end

	# return BasisSet(salc_linearcombo_list)
	# Reconstruct basislist from classified dictionary for BasisSet storage
	result_basislist = SortedCountingUniqueVector{Basis.CoupledBasis}()
	for (key, classified_basislist) in classified_coupled_basisdict
		for cb in classified_basislist
			count = classified_basislist.counts[cb]
			push!(result_basislist, cb, count)
		end
	end
	return BasisSet(result_basislist, salc_list)
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
							if is_translationally_equivalent_coupled_basis(cb, existing_cb, symmetry)
								found_equivalent = true
								# Safe access: get existing count or 0, then add
								orbit_basis_counts[existing_cb] = get(orbit_basis_counts, existing_cb, 0) + count
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
	if cb1.Lseq != cb2.Lseq
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




function collect_cluster_atoms(
	coupled_basislist::SortedCountingUniqueVector{Basis.CoupledBasis},
)::Set{Vector{Int}}
	result_set = Set{Vector{Int}}()
	for cb in coupled_basislist
		push!(result_set, cb.atoms)
	end
	return result_set
end

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
				push!(matched_set, sorted_atom_list_shifted_shifted)
			end
		end
	end
	if isempty(matched_set)
		error("Failed to find translation atoms")
	elseif length(matched_set) > 1
		@show atom_list
		@show cluster_atoms
		@show matched_set
		error("Multiple translation atoms found")
	end
	return first(matched_set)
end


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
					# Calculate rotation matrix for this basis
					rot_mat = Δl(reordered_cb.Lf, rotmat2euler(rotmat)...)
					# Apply time reversal symmetry multiplier if needed
					if time_rev_sym
						total_l = sum(reordered_cb.ls)
						multiplier = (-1)^total_l
						rot_mat = multiplier * rot_mat
					end
					submat_in_mat[j, i] = rot_mat
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

# """
# 	listup_tesseral_basislist(atom_list, lsum; normalize=:none, isotropy::Bool=false)

# List up all tesseral (real) basis functions as `LinearCombo` objects for a given atom list and maximum angular momentum sum.

# This function generates all possible combinations of angular momenta that sum to `lsum` or less,
# and constructs tesseral basis functions using `tesseral_linear_combos_from_tesseral_bases`.

# # Arguments
# - `atom_list::Vector{<:Integer}`: List of atom indices
# - `lsum::Integer`: Maximum sum of angular momentum values
# - `normalize::Symbol`: Normalization option (`:none` or `:fro`, default: `:none`)
# - `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`

# # Returns
# - `Vector{LinearCombo}`: List of `LinearCombo` objects representing tesseral basis functions

# # Examples
# ```julia
# atom_list = [1, 2]
# lsum = 4
# basis_list = listup_tesseral_basislist(atom_list, lsum)
# ```
# """
# function listup_tesseral_basislist(
# 	atom_list::Vector{<:Integer},
# 	lsum::Integer;
# 	normalize::Symbol = :none,
# 	isotropy::Bool = false,
# )::Vector{Basis.LinearCombo}
# 	result_basislist = Vector{Basis.LinearCombo}()
# 	for l in 2:lsum
# 		if l < length(atom_list) || isodd(l)
# 			continue
# 		end
# 		l_list = Combinat.compositions(l, length(atom_list); min = 1)
# 		for l_vec::Vector{Int} in l_list
# 			lc_list = tesseral_linear_combos_from_tesseral_bases(
# 				l_vec,
# 				atom_list;
# 				normalize = normalize,
# 				isotropy = isotropy,
# 			)
# 			append!(result_basislist, lc_list)
# 		end
# 	end
# 	# Remove physically equivalent LinearCombos
# 	return Basis.remove_duplicate_linear_combos(result_basislist)
# end

# """
# 	construct_tesseral_basislist(
# 		structure::Structure,
# 		symmetry::Symmetry,
# 		cluster::Cluster,
# 		body1_lmax::Vector{Int},
# 		bodyn_lsum::OffsetArray{Int, 1},
# 		nbody::Integer;
# 		normalize::Symbol = :none,
# 		isotropy::Bool = false,
# 	)

# Construct a list of tesseral (real) basis functions as `LinearCombo` objects, similar to `construct_basislist`
# but using `tesseral_linear_combos_from_tesseral_bases` instead of `SHProduct`.

# This function follows the same structure as `construct_basislist`:
# 1. Handles 1-body case (skipping l=1 and odd l due to time-reversal symmetry)
# 2. Processes multi-body cases (body=2 to nbody)
# 3. Uses `listup_tesseral_basislist` to generate basis functions for each atom list

# # Arguments
# - `structure::Structure`: Structure information containing atomic positions and species
# - `symmetry::Symmetry`: Symmetry information for the crystal structure
# - `cluster::Cluster`: Cluster information for atomic interactions
# - `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions for each atomic species
# - `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
# - `nbody::Integer`: Maximum number of bodies in interactions
# - `normalize::Symbol`: Normalization option (`:none` or `:fro`, default: `:none`)
# - `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`

# # Returns
# - `Vector{LinearCombo}`: List of `LinearCombo` objects representing tesseral basis functions

# # Examples
# ```julia
# body1_lmax = [2, 3]
# bodyn_lsum = OffsetArray([0, 0, 4, 6], 0:3)
# basis_list = construct_tesseral_basislist(structure, symmetry, cluster, body1_lmax, bodyn_lsum, 3)
# ```
# """
# function construct_tesseral_basislist(
# 	structure::Structure,
# 	symmetry::Symmetry,
# 	cluster::Cluster,
# 	body1_lmax::Vector{Int},
# 	bodyn_lsum::OffsetArray{Int, 1},
# 	nbody::Integer;
# 	normalize::Symbol = :none,
# 	isotropy::Bool = false,
# )::SortedCountingUniqueVector{Basis.LinearCombo}
# 	result_basislist = SortedCountingUniqueVector{Basis.LinearCombo}()
# 	cluster_dict::Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}} =
# 		cluster.cluster_dict

# 	# Handle 1-body case
# 	for iat in symmetry.atoms_in_prim
# 		lmax = body1_lmax[structure.supercell.kd_int_list[iat]]
# 		for l in 2:lmax[1] # skip l = 1 because it is prohibited by the time-reversal symmetry
# 			if l % 2 == 1 # skip odd l cases due to the time-reversal symmetry
# 				continue
# 			end
# 			# For 1-body case, create LinearCombo with single atom
# 			lc_list = tesseral_linear_combos_from_tesseral_bases(
# 				[l],
# 				[iat];
# 				normalize = normalize,
# 				isotropy = isotropy,
# 			)
# 			for lc::Basis.LinearCombo in lc_list
# 				push!(result_basislist, lc, 1)  # multiplicity = 1 for 1-body case
# 			end
# 		end
# 	end

# 	# Process multi-body cases
# 	for body in 2:nbody
# 		body_basislist = SortedCountingUniqueVector{Basis.LinearCombo}()
# 		for prim_atom_sc in symmetry.atoms_in_prim
# 			cuv::CountingUniqueVector{Vector{Int}} = cluster_dict[body][prim_atom_sc]
# 			for atom_list::Vector{Int} in cuv
# 				count = cuv.counts[atom_list]  # Get multiplicity from cluster
# 				lc_list::Vector{Basis.LinearCombo} =
# 					listup_tesseral_basislist(atom_list, bodyn_lsum[body];
# 						normalize = normalize,
# 						isotropy = isotropy,
# 					)
# 				for lc::Basis.LinearCombo in lc_list
# 					push_unique_tesseral_body!(body_basislist, lc, count, symmetry)
# 				end
# 			end
# 		end
# 		for lc in body_basislist
# 			push!(result_basislist, lc, body_basislist.counts[lc])
# 		end
# 	end

# 	return result_basislist
# end

# """
# 	is_translationally_equivalent_linear_combo(lc1::Basis.LinearCombo, lc2::Basis.LinearCombo, symmetry::Symmetry) -> Bool

# Check if two `LinearCombo` objects are translationally equivalent.

# Two `LinearCombo` objects are translationally equivalent if:
# - They are physically equivalent (same `Lf`, `Lseq`, and `(atom, l)` pairs)
# - Their atom lists are related by a translation operation in the supercell

# This function checks if the atom lists can be mapped to each other via translation operations
# defined in `symmetry.symnum_translation`.

# # Arguments
# - `lc1::Basis.LinearCombo`: First `LinearCombo` to compare
# - `lc2::Basis.LinearCombo`: Second `LinearCombo` to compare
# - `symmetry::Symmetry`: Symmetry information containing translation mappings

# # Returns
# - `Bool`: `true` if the `LinearCombo` objects are translationally equivalent, `false` otherwise
# """
# function is_translationally_equivalent_linear_combo(
# 	lc1::Basis.LinearCombo{T1, N1},
# 	lc2::Basis.LinearCombo{T2, N2},
# 	symmetry::Symmetry,
# ) where {T1, T2, N1, N2}
# 	# Different number of sites
# 	N1 != N2 && return false

# 	# Different Lf
# 	lc1.Lf != lc2.Lf && return false

# 	# Different Lseq
# 	lc1.Lseq != lc2.Lseq && return false

# 	# Different coeff_list means different basis functions
# 	if lc1.coeff_list != lc2.coeff_list
# 		return false
# 	end

# 	# Different ls values (as multisets) means different basis functions
# 	ls1_sorted = sort(collect(lc1.ls))
# 	ls2_sorted = sort(collect(lc2.ls))
# 	if ls1_sorted != ls2_sorted
# 		return false
# 	end

# 	# Check if atom lists are translationally equivalent
# 	atom_list1 = lc1.atoms
# 	atom_list2 = lc2.atoms

# 	# Early return if atom lists are the same
# 	if atom_list1 == atom_list2
# 		return false
# 	end

# 	# Early return if first atom is the same
# 	# because this function is intended to be used for different first atoms but translationally equivalent clusters
# 	if atom_list1[1] == atom_list2[1]
# 		return false
# 	end

# 	# Check translation operations
# 	for itran in symmetry.symnum_translation
# 		# Method 1: Apply forward translation (map_sym) to atom_list1
# 		atom_list1_shifted = [symmetry.map_sym[atom, itran] for atom in atom_list1]
# 		# Sort both lists to compare as multisets (order doesn't matter)
# 		if sort(atom_list1_shifted) == sort(atom_list2)
# 			return true
# 		end

# 		# Method 2: Apply inverse translation (map_sym_inv) to atom_list1
# 		atom_list1_shifted = [symmetry.map_sym_inv[atom, itran] for atom in atom_list1]
# 		# Sort both lists to compare as multisets (order doesn't matter)
# 		if sort(atom_list1_shifted) == sort(atom_list2)
# 			return true
# 		end
# 	end

# 	return false
# end

# """
# 	push_unique_tesseral_body!(target::SortedCountingUniqueVector{Basis.LinearCombo}, lc::Basis.LinearCombo, count::Integer, symmetry::Symmetry)

# Add a `LinearCombo` to the target list only if it is not physically or translationally equivalent
# to any existing `LinearCombo` in the list. If equivalent, add the count to the existing entry.

# This function checks:
# 1. If the sum of `l` values matches
# 2. If the `LinearCombo` is physically equivalent (same `Lf`, `Lseq`, `(atom, l)` pairs, and `coeff_list`)
# 3. If the `LinearCombo` is translationally equivalent (same cluster shifted by translation)

# Similar to `push_unique_body!` but for `LinearCombo` objects.
# Note: `coeff_list` must also match for two `LinearCombo` objects to be considered equivalent,
# since different `coeff_list` values represent different basis functions.
# """
# function push_unique_tesseral_body!(
# 	target::SortedCountingUniqueVector{Basis.LinearCombo},
# 	lc::Basis.LinearCombo,
# 	count::Integer,
# 	symmetry::Symmetry,
# )
# 	# Quick check: sum of l values
# 	lsum_lc = sum(collect(lc.ls))
# 	for existing_lc in target
# 		lsum_existing = sum(collect(existing_lc.ls))
# 		if lsum_lc != lsum_existing
# 			continue
# 		end
# 		# Check if physically equivalent (same Lf, Lseq, (atom, l) pairs, and coeff_list)
# 		# This checks for exact matches (same atoms, same ls order)
# 		if Basis.is_physically_equivalent(lc, existing_lc)
# 			# If physically equivalent and atoms are the same, they are duplicates
# 			# Don't add count, just return (already counted)
# 			if lc.atoms == existing_lc.atoms
# 				return
# 			end
# 			# If physically equivalent but atoms differ, check if translationally equivalent
# 			if is_translationally_equivalent_linear_combo(lc, existing_lc, symmetry)
# 				# Translationally equivalent: don't add count, just return (already counted)
# 				return
# 			end
# 		else
# 			# If not physically equivalent, still check if translationally equivalent
# 			# (e.g., [1, 10] vs [9, 2] with same ls and coeff_list)
# 			if is_translationally_equivalent_linear_combo(lc, existing_lc, symmetry)
# 				# Translationally equivalent: don't add count, just return (already counted)
# 				return
# 			end
# 		end
# 	end
# 	# No equivalent LinearCombo found: add with the given count
# 	push!(target, lc, count)
# end


function map_atom_l_list(
	atom_l_list::AbstractVector{<:AbstractVector{<:Integer}},
	map_sym::AbstractMatrix{<:Integer},
	isym::Integer,
)::Vector{Vector{Int}}
	return [[map_sym[atom_l_vec[1], isym], atom_l_vec[2]] for atom_l_vec in atom_l_list]
end

# """
# 	get_atom_l_list_linearcombo(lc::Basis.LinearCombo) -> Vector{Vector{Int}}

# Get the atom-l list from a `LinearCombo` object.

# # Arguments
# - `lc::Basis.LinearCombo`: The `LinearCombo` object

# # Returns
# - `Vector{Vector{Int}}`: List of `[atom, l]` pairs
# """
# function get_atom_l_list_linearcombo(lc::Basis.LinearCombo)::Vector{Vector{Int}}
# 	atom_list = lc.atoms
# 	ls_list = collect(lc.ls)
# 	vec = Vector{Vector{Int}}()
# 	for (atom, l) in zip(atom_list, ls_list)
# 		push!(vec, Int[atom, l])
# 	end
# 	return vec
# end

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

# """
# 	classify_tesseral_basislist(
# 		tesseral_basislist::AbstractVector{Basis.LinearCombo},
# 		map_sym::AbstractMatrix{<:Integer},
# 	) -> Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}

# Classify `tesseral_basislist` (a list of `LinearCombo` objects) by symmetry operations.
# Only combinations that can have finite matrix elements under symmetry operations are grouped together.

# This function groups `LinearCombo` objects that are related by symmetry operations,
# meaning they can have non-zero matrix elements between them.

# Classification criteria:
# - `atom_l_list` (atom and l pairs) must be related by symmetry operations
# - `Lf` (total angular momentum) must be the same (invariant under rotations)
# - `Lseq` (intermediate L values) must be the same (invariant under rotations)

# # Arguments
# - `tesseral_basislist::AbstractVector{Basis.LinearCombo}`: List of `LinearCombo` objects to classify
# - `map_sym::AbstractMatrix{<:Integer}`: Symmetry mapping matrix

# # Returns
# - `Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}`: Dictionary mapping classification labels to groups of `LinearCombo` objects

# # Examples
# ```julia
# classified_dict = classify_tesseral_basislist(tesseral_basislist, symmetry.map_sym)
# ```
# """
# function classify_tesseral_basislist(
# 	tesseral_basislist::AbstractVector{<:Basis.LinearCombo},
# 	map_sym::AbstractMatrix{<:Integer},
# )::Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}
# 	count = 1
# 	label_list = zeros(Int, length(tesseral_basislist))
# 	for (idx, lc) in enumerate(tesseral_basislist)
# 		if label_list[idx] != 0
# 			continue
# 		end
# 		atom_l_list_base = get_atom_l_list_linearcombo(lc)
# 		Lf_base = lc.Lf
# 		Lseq_base = lc.Lseq
# 		for isym in 1:size(map_sym, 2)
# 			mapped_list = map_atom_l_list(atom_l_list_base, map_sym, isym)
# 			sorted_mapped = sort(mapped_list)
# 			for (idx2, lc2) in enumerate(tesseral_basislist)
# 				if label_list[idx2] == 0 &&
# 				   sort(get_atom_l_list_linearcombo(lc2)) == sorted_mapped &&
# 				   lc2.Lf == Lf_base &&
# 				   lc2.Lseq == Lseq_base
# 					label_list[idx2] = count
# 				end
# 			end
# 		end
# 		count += 1
# 	end

# 	dict = OrderedDict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}()
# 	max_label = maximum(label_list)
# 	if max_label > 0
# 		for idx in 1:max_label
# 			dict[idx] = SortedCountingUniqueVector{Basis.LinearCombo}()
# 		end
# 	end

# 	for (lc, label) in zip(tesseral_basislist, label_list)
# 		# Get count from the original tesseral_basislist if it's a SortedCountingUniqueVector
# 		if tesseral_basislist isa SortedCountingUniqueVector
# 			count_val = tesseral_basislist.counts[lc]
# 			push!(dict[label], lc, count_val)
# 		else
# 			push!(dict[label], lc, 1)
# 		end
# 	end

# 	return dict
# end

# """
# 	classify_tesseral_basislist_test(
# 		tesseral_basislist::AbstractVector{Basis.LinearCombo},
# 	) -> Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}

# Simplified classifier for `LinearCombo` objects used in tests.

# This version ignores spatial symmetry and groups basis functions solely by
# interaction order (number of sites) and final angular momentum `Lf`. It trades
# efficiency for robustness so that test fixtures can rely on deterministic
# grouping without depending on symmetry metadata.

# # Arguments
# - `tesseral_basislist::AbstractVector{Basis.LinearCombo}`: List of LinearCombo objects

# # Returns
# - `Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}`: Dictionary keyed by
#   labels assigned per `(nbody, Lf)` pair
# """
# function classify_tesseral_basislist_test(
# 	tesseral_basislist::AbstractVector{<:Basis.LinearCombo},
# )::Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}
# 	if isempty(tesseral_basislist)
# 		return OrderedDict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}()
# 	end

# 	label_map = Dict{Tuple{Int, Int}, Int}()
# 	label_list = Vector{Int}(undef, length(tesseral_basislist))
# 	next_label = 0

# 	for (idx, lc) in enumerate(tesseral_basislist)
# 		key = (length(lc.atoms), lc.Lf)
# 		label = get(label_map, key, 0)
# 		if label == 0
# 			next_label += 1
# 			label = next_label
# 			label_map[key] = label
# 		end
# 		label_list[idx] = label
# 	end

# 	dict = OrderedDict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}()
# 	for label in 1:next_label
# 		dict[label] = SortedCountingUniqueVector{Basis.LinearCombo}()
# 	end

# 	for (lc, label) in zip(tesseral_basislist, label_list)
# 		if tesseral_basislist isa SortedCountingUniqueVector
# 			count_val = get(tesseral_basislist.counts, lc, 1)
# 			push!(dict[label], lc, count_val)
# 		else
# 			push!(dict[label], lc, 1)
# 		end
# 	end

# 	return dict
# end

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
	atom_list_target = [indices.atom for indices in basis_target]
	atom_list_ref = [indices.atom for indices in basis_ref]
	if atom_list_target == atom_list_ref
		return false
	end
	# Early return if the first atom is the same
	# because this function is intended to be used for different first atoms but translationally equivalent clusters
	if basis_target[1].atom == basis_ref[1].atom
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
Finds the cartesian relative vector between 2 atoms specified by (atom, cell) tuples, where cell means virtual cell index (1 <= cell <= 27).
The equation is
r(atom2) - r(atom1)
"""
function calc_relvec_in_cart(
	atom1::NTuple{2, Integer},# (atom, cell)
	atom2::NTuple{2, Integer},
	x_image_cart::AbstractArray{<:Real, 3},
)::SVector{3, Float64}
	return SVector{3, Float64}(
		x_image_cart[:, atom1[1], atom1[2]] - x_image_cart[:, atom2[1], atom2[2]],
	)
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
		basisdict::Dict{Int, SortedCountingUniqueVector{SHProduct}},
	) -> Dict{Int, SortedCountingUniqueVector{SHProduct}}

Filter out basis entries whose `SortedCountingUniqueVector` contains at least
one `SHProduct` with (a) multiplicity greater than 1 and (b) the sum of `l`
quantum numbers exceeding 2. The surviving lists are deep-copied to avoid
aliasing and relabeled sequentially (following ascending original keys) so the
resulting dictionary has deterministic, compact indices.
"""
function filter_basisdict(
	basisdict::Dict{Int, SortedCountingUniqueVector{SHProduct}},
)::Dict{Int, SortedCountingUniqueVector{SHProduct}}
	result_basisdict = Dict{Int, SortedCountingUniqueVector{SHProduct}}()
	new_label = 1
	for label in sort!(collect(keys(basisdict)))
		basislist::SortedCountingUniqueVector{SHProduct} = basisdict[label]
		has_forbidden_basis = any(basislist) do basis::SHProduct
			basislist.counts[basis] > 1 && sum(shsi.l for shsi in basis) > 2
		end
		if has_forbidden_basis
			continue
		end
		result_basisdict[new_label] = deepcopy(basislist)
		new_label += 1
	end
	return result_basisdict
end

"""
	operate_symop_linearcombo(
		basislist::AbstractVector{Basis.LinearCombo},
		lc::Basis.LinearCombo,
		symop::SymmetryOperation,
		map_sym_per_symop::AbstractVector{<:Integer},
		map_sym_inv::AbstractMatrix{<:Integer},
		map_s2p::AbstractVector{<:Maps},
		atoms_in_prim::AbstractVector{<:Integer},
		symnum_translation::AbstractVector{<:Integer},
		time_rev_sym::Bool,
	) -> Basis.LinearCombo

Apply a symmetry operation to a `LinearCombo` object.

This is a more refined implementation that:
1. Applies the symmetry operation to atom indices
2. Tries translation operations directly (more efficient than trying each atom as reference)
3. Finds matching LinearCombo by comparing (atom, l) pairs as multisets
4. Applies rotation matrix to coeff_list (for Lf)
5. Applies time-reversal symmetry if needed

# Arguments
- `basislist::AbstractVector{Basis.LinearCombo}`: List of `LinearCombo` objects
- `lc::Basis.LinearCombo`: The `LinearCombo` to transform
- `symop::SymmetryOperation`: The symmetry operation to apply
- `map_sym_per_symop::AbstractVector{<:Integer}`: Atom mapping for this symmetry operation
- `map_sym_inv::AbstractMatrix{<:Integer}`: Inverse atom mapping matrix
- `map_s2p::AbstractVector{<:Maps}`: Mapping from supercell to primitive cell
- `atoms_in_prim::AbstractVector{<:Integer}`: List of atoms in primitive cell
- `symnum_translation::AbstractVector{<:Integer}`: Translation symmetry numbers
- `time_rev_sym::Bool`: Whether to apply time-reversal symmetry

# Returns
- `Basis.LinearCombo`: The transformed `LinearCombo`
"""
# function operate_symop_linearcombo(
# 	basislist::AbstractVector{Basis.LinearCombo},
# 	lc::Basis.LinearCombo,
# 	symop::SymmetryOperation,
# 	map_sym_per_symop::AbstractVector{<:Integer},
# 	map_sym_inv::AbstractMatrix{<:Integer},
# 	map_s2p::AbstractVector{<:Maps},
# 	atoms_in_prim::AbstractVector{<:Integer},
# 	symnum_translation::AbstractVector{<:Integer},
# 	time_rev_sym::Bool,
# )::Basis.LinearCombo
# 	# Apply symmetry operation to atoms
# 	new_atom_list = [map_sym_per_symop[atom] for atom in lc.atoms]

# 	# Prepare (atom, l) pairs for matching (as multisets, order doesn't matter)
# 	ls_vec = collect(lc.ls)

# 	# Try translation operations to find matching LinearCombo in primitive cell
# 	new_lc_found = nothing
# 	found = false

# 	for itran in symnum_translation
# 		# Apply translation to shift atoms to primitive cell
# 		atom_list_shifted = [map_sym_inv[atom, itran] for atom in new_atom_list]
# 		atom_l_pairs_shifted = collect(zip(atom_list_shifted, ls_vec))

# 		# Find matching LinearCombo in basislist
# 		for lc_candidate in basislist
# 			# Quick checks: Lf and Lseq must match (they are rotation-invariant)
# 			if lc_candidate.Lf != lc.Lf || lc_candidate.Lseq != lc.Lseq
# 				continue
# 			end

# 			# Check if (atom, l) pairs match as multisets
# 			atom_l_pairs_candidate = collect(zip(lc_candidate.atoms, collect(lc_candidate.ls)))
# 			if sort(atom_l_pairs_candidate) == sort(atom_l_pairs_shifted)
# 				found = true
# 				new_lc_found = lc_candidate
# 				break
# 			end
# 		end
# 		if found
# 			break
# 		end
# 	end

# 	if !found
# 		error("Failed to find corresponding LinearCombo in the primitive cell.")
# 	end

# 	# Calculate rotation matrix for Lf and apply to coeff_list
# 	is_proper = symop.is_proper
# 	rotmat = is_proper ? symop.rotation_cart : -1 * symop.rotation_cart

# 	multiplier = time_rev_sym ? (-1)^(sum(ls_vec)) : 1.0

# 	euler_angles = rotmat2euler(rotmat)
# 	rotmat_Lf = multiplier * Δl(lc.Lf, euler_angles...)
# 	new_coeff_list = rotmat_Lf * lc.coeff_list

# 	# Return transformed LinearCombo
# 	return Basis.LinearCombo(
# 		new_lc_found.ls,
# 		new_lc_found.Lf,
# 		new_lc_found.Lseq,
# 		new_lc_found.atoms,
# 		new_coeff_list,
# 		new_lc_found.coeff_tensor,
# 	)
# end

# """
# 	proj_matrix_a_symop_linearcombo(
# 		basislist::AbstractVector{Basis.LinearCombo},
# 		symop::SymmetryOperation,
# 		map_sym_per_symop::AbstractVector{<:Integer},
# 		map_sym_inv::AbstractMatrix{<:Integer},
# 		map_s2p::AbstractVector{<:Maps},
# 		atoms_in_prim::AbstractVector{<:Integer},
# 		symnum_translation::AbstractVector{<:Integer},
# 		time_rev_sym::Bool,
# 	) -> Matrix{Float64}

# Compute the projection matrix for a single symmetry operation applied to a list of `LinearCombo` objects.

# # Arguments
# - `basislist::AbstractVector{Basis.LinearCombo}`: List of `LinearCombo` objects
# - `symop::SymmetryOperation`: The symmetry operation
# - `map_sym_per_symop::AbstractVector{<:Integer}`: Atom mapping for this symmetry operation
# - `map_sym_inv::AbstractMatrix{<:Integer}`: Inverse atom mapping matrix
# - `map_s2p::AbstractVector{<:Maps}`: Mapping from supercell to primitive cell
# - `atoms_in_prim::AbstractVector{<:Integer}`: List of atoms in primitive cell
# - `symnum_translation::AbstractVector{<:Integer}`: Translation symmetry numbers
# - `time_rev_sym::Bool`: Whether to apply time-reversal symmetry

# # Returns
# - `Matrix{Float64}`: The projection matrix for this symmetry operation
# """
# function proj_matrix_a_symop_linearcombo(
# 	basislist::AbstractVector{Basis.LinearCombo},
# 	symop::SymmetryOperation,
# 	map_sym_per_symop::AbstractVector{<:Integer},
# 	map_sym_inv::AbstractMatrix{<:Integer},
# 	map_s2p::AbstractVector{<:Maps},
# 	atoms_in_prim::AbstractVector{<:Integer},
# 	symnum_translation::AbstractVector{<:Integer},
# 	time_rev_sym::Bool,
# )::Matrix{Float64}
# 	dim = length(basislist)
# 	projection_mat = zeros(Float64, dim, dim)

# 	for (j, lc_j) in enumerate(basislist)
# 		lc_j_transformed = operate_symop_linearcombo(
# 			basislist,
# 			lc_j,
# 			symop,
# 			map_sym_per_symop,
# 			map_sym_inv,
# 			map_s2p,
# 			atoms_in_prim,
# 			symnum_translation,
# 			time_rev_sym,
# 		)
# 		for (i, lc_i) in enumerate(basislist)
# 			projection_mat[i, j] = dot(lc_i, lc_j_transformed)
# 		end
# 	end

# 	return projection_mat
# end

# """
# 	projection_matrix_linearcombo(
# 		basisdict::Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}},
# 		symmetry::Symmetry,
# 	) -> Vector{Matrix{Float64}}

# Construct projection matrices for each classification label in the `basisdict`.

# This function computes the average projection matrix over all symmetry operations
# (including time-reversal symmetry) for each group of `LinearCombo` objects.

# # Arguments
# - `basisdict::Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}}`: Dictionary mapping classification labels to groups of `LinearCombo` objects
# - `symmetry::Symmetry`: Symmetry information containing all symmetry operations

# # Returns
# - `Vector{Matrix{Float64}}`: List of projection matrices, one for each classification label
# """
# function projection_matrix_linearcombo(
# 	basisdict::Dict{Int, SortedCountingUniqueVector{Basis.LinearCombo}},
# 	symmetry::Symmetry,
# )::Vector{Matrix{Float64}}
# 	result_projections = Vector{Matrix{Float64}}(undef, length(basisdict))

# 	idx_list = sort(collect(keys(basisdict)))
# 	for idx in idx_list
# 		basislist = basisdict[idx]
# 		dim = length(basislist)
# 		local_projection_mat = zeros(Float64, dim, dim)

# 		for (n, symop) in enumerate(symmetry.symdata), time_rev_sym in [false, true]
# 			projection_mat_per_symop = proj_matrix_a_symop_linearcombo(
# 				basislist,
# 				symop,
# 				@view(symmetry.map_sym[:, n]),
# 				symmetry.map_sym_inv,
# 				symmetry.map_s2p,
# 				symmetry.atoms_in_prim,
# 				symmetry.symnum_translation,
# 				time_rev_sym,
# 			)
# 			local_projection_mat += projection_mat_per_symop
# 		end

# 		local_projection_mat = local_projection_mat ./ (2 * symmetry.nsym)
# 		local_projection_mat = hermitianpart(local_projection_mat)

# 		if !is_symmetric(local_projection_mat, tol = 1e-10)
# 			error("Projection matrix is not symmetric. index: $idx")
# 		end

# 		result_projections[idx] = local_projection_mat
# 	end
# 	return result_projections
# end

end # module BasisSets


