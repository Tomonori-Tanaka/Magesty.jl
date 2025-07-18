"""
	module Symmetries

A module for handling symmetry operations in crystal structures using `Spglib`. It defines `SymmetryOperation` and `Symmetry` types to store symmetry-related data and mappings.

# Types
- **`SymmetryOperation`**: Represents a single symmetry operation with rotation and translation components.
- **`Symmetry`**: Represents the symmetry information of a structure, including symmetry operations, mappings, and metadata.

# Functions
- `is_compatible_cart(rot_cart)`: Check if a rotation matrix is compatible with Cartesian axes.
- `Symmetry(structure, tol)`: Generate a `Symmetry` object for a given structure with a specified tolerance.
"""
module Symmetries

using Base.Threads
using Combinatorics: with_replacement_combinations
using LinearAlgebra
using StaticArrays
using Spglib
using Printf

using ..AtomCells
using ..ConfigParser
using ..Structures
import Base: isless, show

export SymmetryOperation, Maps, Symmetry

"""
	struct SymmetryOperation

Represents a symmetry operation, including rotation and translation components.

# Fields
- `rotation_frac::SMatrix{3, 3, Float64, 9}`: The 3x3 rotation matrix in fractional coordinates.
- `rotation_cart::SMatrix{3, 3, Float64, 9}`: The 3x3 rotation matrix in Cartesian coordinates.
- `translation_frac::SVector{3, Float64}`: The 3x1 translation vector in fractional coordinates.
- `is_translation::Bool`: True if the operation is a pure translation.
- `is_translation_included::Bool`: True if the operation includes a translation.
- `is_proper::Bool`: True if the rotation is proper (determinant > 0).

# Examples
```julia
# Create a symmetry operation
symop = SymmetryOperation(
	rotation_frac = [1 0 0; 0 1 0; 0 0 1],
	rotation_cart = [1 0 0; 0 1 0; 0 0 1],
	translation_frac = [0.0, 0.0, 0.0],
	is_translation = true,
	is_translation_included = false,
	is_proper = true
)
```
"""
struct SymmetryOperation
	rotation_frac::SMatrix{3, 3, Float64, 9}
	rotation_cart::SMatrix{3, 3, Float64, 9}
	translation_frac::SVector{3, Float64}
	is_translation::Bool
	is_translation_included::Bool
	is_proper::Bool
end

function isless(symop1::SymmetryOperation, symop2::SymmetryOperation)
	symop1_rot_flatten = vec(transpose(symop1.rotation_frac))
	symop2_rot_flatten = vec(transpose(symop2.rotation_frac))

	# if matrix element < 0, add 1.0
	symop1_translation = [
		symop1.translation_frac[i] < 0.0 ? 1.0 + symop1.translation_frac[i] :
		symop1.translation_frac[i] for i in 1:3
	]
	symop2_translation = [
		symop2.translation_frac[i] < 0.0 ? 1.0 + symop2.translation_frac[i] :
		symop2.translation_frac[i] for i in 1:3
	]

	return vcat(symop1_rot_flatten, symop1_translation) <
		   vcat(symop2_rot_flatten, symop2_translation)
end

function show(io::IO, symop::SymmetryOperation)
	println("rotation_frac: ")
	display(symop.rotation_frac)
	println("rotation_cart: ")
	display(symop.rotation_cart)
	println("translation_frac: ", symop.translation_frac)
	println("is_translation: ", symop.is_translation)
	println("is_translation_included: ", symop.is_translation_included)
	println("is_proper: ", symop.is_proper)
end

""" Maps
Maps structure is used in Symmetry structure to map supercell atom to primitive index.
"""
struct Maps
	atom::Int
	translation::Int
end

"""
	struct Symmetry

Contains the symmetry information of a structure.

# Fields
- `international_symbol::String`: International symbol of the space group.
- `spacegroup_number::Int`: Space group number.
- `nsym::Int`: Number of symmetry operations.
- `ntran::Int`: Number of pure translation operations.
- `nat_prim::Int`: Number of atoms in the primitive cell.
- `tol::Float64`: Tolerance for symmetry detection.
- `atoms_in_prim::Vector{Int}`: Indices of atoms in the primitive cell.
- `symdata::Vector{SymmetryOperation}`: List of symmetry operations.
- `map_sym::Matrix{Int}`: Maps atoms in the supercell to corresponding atoms under symmetry operations.
- `map_p2s::Matrix{Int}`: Maps atoms in the primitive cell to the supercell.
- `map_s2p::Vector{Maps}`: Maps atoms in the supercell to the primitive cell.
- `symnum_translation::Vector{Int}`: Indices of pure translation operations.

# Constructor
	Symmetry(structure, tol::Float64)

Generate symmetry information for the given structure using the specified tolerance.
"""
struct Symmetry
	international_symbol::String
	spacegroup_number::Int
	nsym::Int   # the number of symmetry operations
	ntran::Int  # the number of translational only operations
	nat_prim::Int   # the total number of atoms in a primitive cell
	tol::Float64
	atoms_in_prim::Vector{Int}

	symdata::Vector{SymmetryOperation}
	map_sym::Matrix{Int}    # [num_atoms, nsym] -> corresponding atom index
	map_p2s::Matrix{Int}    # [nat_prim, ntran] -> corresponding atom index
	map_s2p::Vector{Maps}   # [nat] -> corresponding atom index in primitive cel
	symnum_translation::Vector{Int} # contains the indice of translational only operations

	function Symmetry(structure::Structure, tol::Real; verbosity::Bool = true)

		start_time = time_ns()

		if tol <= 0
			throw(ArgumentError("Tolerance must be positive, got $tol"))
		end

		cell = structure.supercell
		# convert x_frac::Matrix{Float64} to x_frac::Vector{Vector{Float64}}
		x_frac_vec = [col for col in eachcol(cell.x_frac)]
		spglib_data::Spglib.Dataset =
			get_dataset(Spglib.Cell(cell.lattice_vectors, x_frac_vec, cell.kd_int_list))

		if spglib_data.n_operations == 0
			error(
				"Error in symmetry search: No symmetry operations found. Please check the input structure or tolerance setting.",
			)
		end

		# construct symnum_translation and ntran
		symnum_translation = construct_symnum_translation(spglib_data, tol)
		ntran = length(symnum_translation)

		# construct symdata
		symdata = construct_symdata(spglib_data, tol, symnum_translation, cell)

		# construct mapping data
		map_sym = construct_map_sym(spglib_data, tol, structure)

		# generate map_p2s (primitive cell --> supercell)
		map_p2s = construct_map_p2s(spglib_data, cell, map_sym, symnum_translation)

		nat_prim = max(spglib_data.mapping_to_primitive...)
		# generate map_s2p (supercell -> primitive cell)
		map_s2p = construct_map_s2p(cell, map_p2s, nat_prim, ntran)

		# collect atom indices in primitive cell from map_p2s
		atoms_in_prim = Int[map_p2s[i, 1] for i in 1:nat_prim]
		atoms_in_prim = sort(atoms_in_prim)

		if verbosity
			print_symmetry_stdout(
				spglib_data.international_symbol,
				spglib_data.spacegroup_number,
				spglib_data.n_operations,
				ntran,
				nat_prim)
			elapsed_time = (time_ns() - start_time) / 1e9
			println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
			println("-------------------------------------------------------------------")
		end


		return new(
			spglib_data.international_symbol,
			spglib_data.spacegroup_number,
			spglib_data.n_operations,
			ntran,
			nat_prim,
			tol,
			atoms_in_prim,
			symdata,
			map_sym,
			map_p2s,
			map_s2p,
			symnum_translation)
	end
end

function Symmetry(structure::Structure, config::Config4System; verbosity::Bool = true)
	return Symmetry(structure, config.tolerance_sym, verbosity = verbosity)
end

function construct_symnum_translation(spglib_data::Spglib.Dataset, tol::Real)::Vector{Int}
	symnum_translation = Int[]
	for i in 1:spglib_data.n_operations
		if isapprox(spglib_data.rotations[i], I, atol = tol)
			append!(symnum_translation, i)
		end
	end
	return symnum_translation
end

function construct_symdata(
	spglib_data::Spglib.Dataset,
	tol::Real,
	symnum_translation::Vector{Int},
	cell::Structures.Cell,
)::Vector{SymmetryOperation}
	symdata = Vector{SymmetryOperation}(undef, spglib_data.n_operations)
	for i in 1:spglib_data.n_operations
		rotation_cart =
			cell.lattice_vectors * spglib_data.rotations[i] * cell.reciprocal_vectors

		if det(rotation_cart) > 0 ? is_proper = true : is_proper = false
		end
		translation_frac =
			(abs.(spglib_data.translations[i]) .>= tol) .* spglib_data.translations[i]
		# check a translation vector is included in the symmetry operation.
		is_translation_included = false
		for itrans in Base.tail(Tuple(symnum_translation))
			for (elem_translation_frac, elem_translation) in
				zip(translation_frac, spglib_data.translations[itrans])
				if isapprox(elem_translation, 0.0, atol = tol)
					continue
				end
				if elem_translation_frac >= elem_translation
					is_translation_included = true
					break
				end
			end
		end

		symdata_elem = SymmetryOperation(
			spglib_data.rotations[i],
			rotation_cart,
			translation_frac,
			isapprox(spglib_data.rotations[i], I, atol = tol),
			is_translation_included,
			is_proper,
		)
		symdata[i] = symdata_elem
	end
	return symdata
end

function construct_map_sym(
	spglib_data::Spglib.Dataset,
	tol::Real,
	structure::Structure,
)::Matrix{Int}
	natomtypes = length(structure.atomtype_group)
	map_sym = zeros(Int, structure.supercell.num_atoms, Int(spglib_data.n_operations))

	# Process symmetry operations in parallel
	@threads for isym in 1:spglib_data.n_operations
		# Create thread-local storage to avoid memory conflicts
		local_map = zeros(Int, structure.supercell.num_atoms)
		local_x_new = MVector{3, Float64}(undef)
		local_tmp = MVector{3, Float64}(undef)

		for itype in 1:natomtypes
			for iat in structure.atomtype_group[itype]
				# Apply rotation and translation
				local_x_new .=
					spglib_data.rotations[isym] * structure.supercell.x_frac[:, iat] +
					spglib_data.translations[isym]

				for jat in structure.atomtype_group[itype]
					# Calculate relative position
					local_tmp .= (abs.(structure.supercell.x_frac[:, jat] - local_x_new)) .% 1.0
					for (i, val) in enumerate(local_tmp)
						local_tmp[i] = min(val, 1 - val)
					end

					if norm(local_tmp) < tol
						local_map[iat] = jat
						break
					end
				end
			end
		end

		# Write results to shared array
		@inbounds map_sym[:, isym] = local_map
	end

	# Verify results
	zero_pos = CartesianIndices(map_sym)[map_sym .== 0]
	if !isempty(zero_pos)
		error("zero is found in map_sym at $zero_pos")
	end

	return map_sym
end

"""
	construct_map_p2s(spglib_data::Spglib.Dataset, cell::Structures.Cell, map_sym::Matrix{Int}, symnum_translation::Vector{Int}) -> Matrix{Int}

Constructs the mapping from primitive cell to supercell.

# Arguments
- `spglib_data::Spglib.Dataset`: Spglib dataset containing symmetry information
- `cell::Structures.Cell`: Cell information
- `map_sym::Matrix{Int}`: Symmetry mapping matrix
- `symnum_translation::Vector{Int}`: Indices of pure translation operations

# Returns
- `Matrix{Int}`: Mapping matrix from primitive cell to supercell

# Throws
- `ErrorException` if the mapping cannot be constructed properly
"""
function construct_map_p2s(
	spglib_data::Spglib.Dataset,
	cell::Structures.Cell,
	map_sym::AbstractMatrix{<:Integer},
	symnum_translation::AbstractVector{<:Integer},
)::Matrix{Int}
	nat_prim = max(spglib_data.mapping_to_primitive...)
	ntran = length(symnum_translation)
	map_p2s = zeros(Int, nat_prim, ntran)
	is_checked = fill(false, cell.num_atoms)

	jat = 1
	for iat in 1:cell.num_atoms
		if is_checked[iat]
			continue
		end
		for i in 1:ntran
			atomnum_translated = map_sym[iat, symnum_translation[i]]
			map_p2s[jat, i] = atomnum_translated
			is_checked[atomnum_translated] = true
		end
		jat += 1
	end
	if 0 in map_p2s
		error("something wrong in generating map_p2s")
	end
	return map_p2s
end

"""
	construct_map_s2p(cell::Structures.Cell, map_p2s::Matrix{Int}, nat_prim::Int, ntran::Int) -> Vector{Maps}

Constructs the mapping from supercell to primitive cell.

# Arguments
- `cell::Structures.Cell`: Cell information
- `map_p2s::Matrix{Int}`: Mapping matrix from primitive cell to supercell
- `nat_prim::Int`: Number of atoms in primitive cell
- `ntran::Int`: Number of pure translation operations

# Returns
- `Vector{Maps}`: Mapping vector from supercell to primitive cell

# Throws
- `ErrorException` if the mapping cannot be constructed properly
"""
function construct_map_s2p(
	cell::Structures.Cell,
	map_p2s::AbstractMatrix{<:Integer},
	nat_prim::Integer,
	ntran::Integer,
)::Vector{Maps}
	map_s2p = Vector{Maps}(undef, cell.num_atoms)
	initialized = falses(size(map_s2p))
	for iat in 1:nat_prim
		for itran in 1:ntran
			atomnum_translated = map_p2s[iat, itran]
			map_s2p[atomnum_translated] = Maps(iat, itran)
			initialized[atomnum_translated] = true
		end
	end
	undef_indices = findall(x -> x == false, initialized)
	if length(undef_indices) >= 1
		error("undef is detected in `map_s2p` variable: $undef_indices")
	end
	return map_s2p
end

"""
	find_matching_image_cell(symop::SymmetryOperation, x_image::AbstractArray{<:Real, 3}, atom::Integer, cell::Integer) -> Int

Finds the matching image cell index for a given symmetry operation applied to an atom within a supercell.

# Arguments
- `symop::SymmetryOperation`: The symmetry operation to be applied
- `x_image::AbstractArray{<:Real, 3}`: Fractional coordinates of atoms in different images
- `atom::Integer`: Index of the atom to apply symmetry operation to
- `cell::Integer`: Image cell index of the atom
- `tol::Real`: Tolerance for floating point comparisons (default: 1e-5)

# Returns
- `Int`: The matching image cell index if found, -1 if not found

# Throws
- `ArgumentError`: If atom or cell indices are out of bounds
- `ErrorException`: If multiple matches are found
"""
function find_matching_image_cell(
	x_new::MVector{3, Float64},
	x_image::AbstractArray{<:Real, 3},
	;
	tol::Real = 1e-5,
)::Int
	# Use views for better performance
	matches = [
		(n, m) for n in 1:size(x_image, 2), m in 1:size(x_image, 3)
		if isapprox(
			SVector{3, Float64}(@view(x_image[:, n, m])),
			SVector{3, Float64}(x_new);
			atol = tol,
		)
	]

	if isempty(matches)
		return -1
	elseif length(matches) == 1
		return matches[1][2]
	else
		error("Multiple matching image cells found. ")
	end
end

is_in_centeringcell(xf, x_image_frac) = any(xf ≈ vec for vec in x_image_frac[:, :, 1])


function print_symmetry_stdout(
	international_symbol::AbstractString,
	spacegroup_number::Integer,
	nsym::Integer,
	ntran::Integer,
	nat_prim::Integer,
)
	println("""

	SYMMETRY
	========
	""")
	str = """
	 Space group:  $(international_symbol)  ($(spacegroup_number))
	 Number of symmetry operations = $(nsym)

	"""
	if ntran == 1
		str_prim = """
		 Given structure is a primitive cell.
		 Primitive cell contains $(nat_prim) atoms.
		"""
		str *= str_prim
	else
		str_supercell = """
		 Given structure is not a primitive cell.
		 There are $(ntran) translation operations.
		 Primitive cell contains $(nat_prim) atoms.
		"""
		str *= str_supercell
	end

	println(str)
	println("")
end

end

