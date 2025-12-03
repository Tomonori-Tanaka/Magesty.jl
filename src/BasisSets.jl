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
	struct BasisSet

Represents a set of basis functions for atomic interactions in a crystal structure.
This structure is used to store and manage basis functions that are adapted to the symmetry of the crystal.

# Fields
- `salc_list::Vector{SALC}`: List of symmetry-adapted linear combinations

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
- `BasisSet`: A new basis set instance containing symmetry-adapted linear combinations (SALCs)

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
1. Constructs the basis list by considering all possible combinations of atoms and angular momenta
2. Classifies basis functions by symmetry operations
3. Constructs projection matrices for each symmetry label
4. Generates symmetry-adapted linear combinations (SALCs)
"""
struct BasisSet
	salc_list::Vector{SALC}
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

	# Validate input parameters
	# Construct basis list
	# basislist consists of all possible basis functions which is the product of spherical harmonics.
	if verbosity
		print("Constructing basis list...")
	end
	tesseral_basislist::SortedCountingUniqueVector{Basis.LinearCombo} = construct_tesseral_basislist(
		structure,
		symmetry,
		cluster,
		body1_lmax,
		bodyn_lsum,
		nbody,
	)
	for lc in tesseral_basislist
		@show lc, tesseral_basislist.counts[lc]
	end

	basislist::SortedCountingUniqueVector{SHProduct} =
		construct_basislist(structure, symmetry, cluster, body1_lmax, bodyn_lsum, nbody)

	classified_basisdict::Dict{Int, SortedCountingUniqueVector{SHProduct}} =
		classify_basislist(basislist, symmetry.map_sym)

	classified_basisdict = filter_basisdict(classified_basisdict)
	if verbosity
		print(" Done\n")
	end

	if verbosity
		print("Constructing projection matrices...")
	end
	projection_list = projection_matrix(classified_basisdict, symmetry)
	salc_list = Vector{SALC}()
	for idx in eachindex(projection_list)
		eigenvals, eigenvecs = eigen(projection_list[idx])
		eigenvals = real.(round.(eigenvals, digits = 6))
		eigenvecs = round.(eigenvecs .* (abs.(eigenvecs) .≥ 1e-8), digits = 10)
		if !is_proper_eigenvals(eigenvals)
			println(eigenvals)
			@warn "Critical error: Eigenvalues must be either 0 or 1. index: $idx"
		end
		for idx_eigenval in findall(x -> isapprox(x, 1.0, atol = 1e-8), eigenvals)
			eigenvec = eigenvecs[:, idx_eigenval]
			eigenvec = flip_vector_if_negative_sum(eigenvec)
			eigenvec = round.(eigenvec .* (abs.(eigenvec) .≥ 1e-8), digits = 10)
			push!(salc_list, SALC(classified_basisdict[idx], eigenvec / norm(eigenvec)))
		end
	end
	salc_list = filter(salc -> !is_identically_zero(salc), salc_list)

	if verbosity
		print(" Done\n")
	end

	if verbosity
		print_basisset_stdout(salc_list)
		elapsed_time = (time_ns() - start_time) / 1e9  # Convert to seconds
		println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
		println("-------------------------------------------------------------------")
	end

	return BasisSet(salc_list)
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
)::SortedCountingUniqueVector{SHProduct}

	result_basislist = SortedCountingUniqueVector{SHProduct}()
	cluster_dict::Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}} =
		cluster.cluster_dict

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

"""
	listup_tesseral_basislist(atom_list, lsum; normalize=:none, isotropy::Bool=false)

List up all tesseral (real) basis functions as `LinearCombo` objects for a given atom list and maximum angular momentum sum.

This function generates all possible combinations of angular momenta that sum to `lsum` or less,
and constructs tesseral basis functions using `tesseral_linear_combos_from_tesseral_bases`.

# Arguments
- `atom_list::Vector{<:Integer}`: List of atom indices
- `lsum::Integer`: Maximum sum of angular momentum values
- `normalize::Symbol`: Normalization option (`:none` or `:fro`, default: `:none`)
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`

# Returns
- `Vector{LinearCombo}`: List of `LinearCombo` objects representing tesseral basis functions

# Examples
```julia
atom_list = [1, 2]
lsum = 4
basis_list = listup_tesseral_basislist(atom_list, lsum)
```
"""
function listup_tesseral_basislist(
	atom_list::Vector{<:Integer},
	lsum::Integer;
	normalize::Symbol = :none,
	isotropy::Bool = false,
)::Vector{Basis.LinearCombo}
	result_basislist = Vector{Basis.LinearCombo}()
	for l in 2:lsum
		if l < length(atom_list) || isodd(l)
			continue
		end
		l_list = Combinat.compositions(l, length(atom_list); min = 1)
		for l_vec::Vector{Int} in l_list
			lc_list = tesseral_linear_combos_from_tesseral_bases(
				l_vec,
				atom_list;
				normalize = normalize,
				isotropy = isotropy,
			)
			append!(result_basislist, lc_list)
		end
	end
	# Remove physically equivalent LinearCombos
	return Basis.remove_duplicate_linear_combos(result_basislist)
end

"""
	construct_tesseral_basislist(
		structure::Structure,
		symmetry::Symmetry,
		cluster::Cluster,
		body1_lmax::Vector{Int},
		bodyn_lsum::OffsetArray{Int, 1},
		nbody::Integer;
		normalize::Symbol = :none,
		isotropy::Bool = false,
	)

Construct a list of tesseral (real) basis functions as `LinearCombo` objects, similar to `construct_basislist`
but using `tesseral_linear_combos_from_tesseral_bases` instead of `SHProduct`.

This function follows the same structure as `construct_basislist`:
1. Handles 1-body case (skipping l=1 and odd l due to time-reversal symmetry)
2. Processes multi-body cases (body=2 to nbody)
3. Uses `listup_tesseral_basislist` to generate basis functions for each atom list

# Arguments
- `structure::Structure`: Structure information containing atomic positions and species
- `symmetry::Symmetry`: Symmetry information for the crystal structure
- `cluster::Cluster`: Cluster information for atomic interactions
- `body1_lmax::Vector{Int}`: Maximum angular momentum values for 1-body interactions for each atomic species
- `bodyn_lsum::OffsetArray{Int, 1}`: Maximum sum of angular momentum values for multi-body interactions
- `nbody::Integer`: Maximum number of bodies in interactions
- `normalize::Symbol`: Normalization option (`:none` or `:fro`, default: `:none`)
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0), default: `false`

# Returns
- `Vector{LinearCombo}`: List of `LinearCombo` objects representing tesseral basis functions

# Examples
```julia
body1_lmax = [2, 3]
bodyn_lsum = OffsetArray([0, 0, 4, 6], 0:3)
basis_list = construct_tesseral_basislist(structure, symmetry, cluster, body1_lmax, bodyn_lsum, 3)
```
"""
function construct_tesseral_basislist(
	structure::Structure,
	symmetry::Symmetry,
	cluster::Cluster,
	body1_lmax::Vector{Int},
	bodyn_lsum::OffsetArray{Int, 1},
	nbody::Integer;
	normalize::Symbol = :none,
	isotropy::Bool = false,
)::SortedCountingUniqueVector{Basis.LinearCombo}
	result_basislist = SortedCountingUniqueVector{Basis.LinearCombo}()
	cluster_dict::Dict{Int, Dict{Int, CountingUniqueVector{Vector{Int}}}} =
		cluster.cluster_dict

	# Handle 1-body case
	for iat in symmetry.atoms_in_prim
		lmax = body1_lmax[structure.supercell.kd_int_list[iat]]
		for l in 2:lmax[1] # skip l = 1 because it is prohibited by the time-reversal symmetry
			if l % 2 == 1 # skip odd l cases due to the time-reversal symmetry
				continue
			end
			# For 1-body case, create LinearCombo with single atom
			lc_list = tesseral_linear_combos_from_tesseral_bases(
				[l],
				[iat];
				normalize = normalize,
				isotropy = isotropy,
			)
			for lc::Basis.LinearCombo in lc_list
				push!(result_basislist, lc, 1)  # multiplicity = 1 for 1-body case
			end
		end
	end

	# Process multi-body cases
	for body in 2:nbody
		body_basislist = SortedCountingUniqueVector{Basis.LinearCombo}()
		for prim_atom_sc in symmetry.atoms_in_prim
			cuv::CountingUniqueVector{Vector{Int}} = cluster_dict[body][prim_atom_sc]
			for atom_list::Vector{Int} in cuv
				count = cuv.counts[atom_list]  # Get multiplicity from cluster
				lc_list::Vector{Basis.LinearCombo} =
					listup_tesseral_basislist(atom_list, bodyn_lsum[body];
						normalize = normalize,
						isotropy = isotropy,
					)
				for lc::Basis.LinearCombo in lc_list
					push_unique_tesseral_body!(body_basislist, lc, count, symmetry)
				end
			end
		end
		for lc in body_basislist
			push!(result_basislist, lc, body_basislist.counts[lc])
		end
	end

	return result_basislist
end

"""
	is_translationally_equivalent_linear_combo(lc1::Basis.LinearCombo, lc2::Basis.LinearCombo, symmetry::Symmetry) -> Bool

Check if two `LinearCombo` objects are translationally equivalent.

Two `LinearCombo` objects are translationally equivalent if:
- They are physically equivalent (same `Lf`, `Lseq`, and `(atom, l)` pairs)
- Their atom lists are related by a translation operation in the supercell

This function checks if the atom lists can be mapped to each other via translation operations
defined in `symmetry.symnum_translation`.

# Arguments
- `lc1::Basis.LinearCombo`: First `LinearCombo` to compare
- `lc2::Basis.LinearCombo`: Second `LinearCombo` to compare
- `symmetry::Symmetry`: Symmetry information containing translation mappings

# Returns
- `Bool`: `true` if the `LinearCombo` objects are translationally equivalent, `false` otherwise
"""
function is_translationally_equivalent_linear_combo(
	lc1::Basis.LinearCombo{T1, N1},
	lc2::Basis.LinearCombo{T2, N2},
	symmetry::Symmetry,
) where {T1, T2, N1, N2}
	# Different number of sites
	N1 != N2 && return false

	# Different Lf
	lc1.Lf != lc2.Lf && return false

	# Different Lseq
	lc1.Lseq != lc2.Lseq && return false

	# Different coeff_list means different basis functions
	if lc1.coeff_list != lc2.coeff_list
		return false
	end

	# Different ls values (as multisets) means different basis functions
	ls1_sorted = sort(collect(lc1.ls))
	ls2_sorted = sort(collect(lc2.ls))
	if ls1_sorted != ls2_sorted
		return false
	end

	# Check if atom lists are translationally equivalent
	atom_list1 = lc1.atoms
	atom_list2 = lc2.atoms

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
	push_unique_tesseral_body!(target::SortedCountingUniqueVector{Basis.LinearCombo}, lc::Basis.LinearCombo, count::Integer, symmetry::Symmetry)

Add a `LinearCombo` to the target list only if it is not physically or translationally equivalent
to any existing `LinearCombo` in the list. If equivalent, add the count to the existing entry.

This function checks:
1. If the sum of `l` values matches
2. If the `LinearCombo` is physically equivalent (same `Lf`, `Lseq`, `(atom, l)` pairs, and `coeff_list`)
3. If the `LinearCombo` is translationally equivalent (same cluster shifted by translation)

Similar to `push_unique_body!` but for `LinearCombo` objects.
Note: `coeff_list` must also match for two `LinearCombo` objects to be considered equivalent,
since different `coeff_list` values represent different basis functions.
"""
function push_unique_tesseral_body!(
	target::SortedCountingUniqueVector{Basis.LinearCombo},
	lc::Basis.LinearCombo,
	count::Integer,
	symmetry::Symmetry,
)
	# Quick check: sum of l values
	lsum_lc = sum(collect(lc.ls))
	for existing_lc in target
		lsum_existing = sum(collect(existing_lc.ls))
		if lsum_lc != lsum_existing
			continue
		end
		# Check if physically equivalent (same Lf, Lseq, (atom, l) pairs, and coeff_list)
		# This checks for exact matches (same atoms, same ls order)
		if Basis.is_physically_equivalent(lc, existing_lc)
			# If physically equivalent and atoms are the same, they are duplicates
			# Don't add count, just return (already counted)
			if lc.atoms == existing_lc.atoms
				return
			end
			# If physically equivalent but atoms differ, check if translationally equivalent
			if is_translationally_equivalent_linear_combo(lc, existing_lc, symmetry)
				# Translationally equivalent: don't add count, just return (already counted)
				return
			end
		else
			# If not physically equivalent, still check if translationally equivalent
			# (e.g., [1, 10] vs [9, 2] with same ls and coeff_list)
			if is_translationally_equivalent_linear_combo(lc, existing_lc, symmetry)
				# Translationally equivalent: don't add count, just return (already counted)
				return
			end
		end
	end
	# No equivalent LinearCombo found: add with the given count
	push!(target, lc, count)
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
		atom_l_list_base = get_atom_l_list(basis)
		for isym in 1:size(map_sym, 2)
			mapped_list = map_atom_l_list(atom_l_list_base, map_sym, isym)
			sorted_mapped = sort(mapped_list)
			for (idx2, basis2) in enumerate(basislist)
				if label_list[idx2] == 0 && sort(get_atom_l_list(basis2)) == sorted_mapped
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

end # module BasisSets


