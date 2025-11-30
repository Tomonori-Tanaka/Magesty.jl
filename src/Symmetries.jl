"""
	module Symmetries

A module for handling symmetry operations in crystal structures using `Spglib`. It defines `SymmetryOperation` and `Symmetry` types to store symmetry-related data and mappings.

# Types
- **`SymmetryOperation`**: Represents a single symmetry operation with rotation and translation components.
- **`Symmetry`**: Represents the symmetry information of a structure, including symmetry operations, mappings, and metadata.

# Functions
- `Symmetry(structure, tol)`: Generate a `Symmetry` object for a given structure with a specified tolerance.
"""
module Symmetries

using Base.Threads
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
- `is_proper::Bool`: True if the rotation is proper (determinant > 0).

# Examples
```julia
# Create a symmetry operation
symop = SymmetryOperation(
	rotation_frac = [1 0 0; 0 1 0; 0 0 1],
	rotation_cart = [1 0 0; 0 1 0; 0 0 1],
	translation_frac = [0.0, 0.0, 0.0],
	is_translation = true,
	is_proper = true
)
```
"""
struct SymmetryOperation
	rotation_frac::SMatrix{3, 3, Float64, 9}
	rotation_cart::SMatrix{3, 3, Float64, 9}
	translation_frac::SVector{3, Float64}
	is_translation::Bool
	is_proper::Bool
end

function isless(symop1::SymmetryOperation, symop2::SymmetryOperation)
	symop1_rot_flatten = vec(transpose(symop1.rotation_frac))
	symop2_rot_flatten = vec(transpose(symop2.rotation_frac))

	# Normalize translation vectors: if element < 0, add 1.0
	symop1_translation = [
		val < 0.0 ? 1.0 + val : val for val in symop1.translation_frac
	]
	symop2_translation = [
		val < 0.0 ? 1.0 + val : val for val in symop2.translation_frac
	]

	return vcat(symop1_rot_flatten, symop1_translation) <
		   vcat(symop2_rot_flatten, symop2_translation)
end

function show(io::IO, symop::SymmetryOperation)
	println(io, "rotation_frac: ")
	print(io, symop.rotation_frac)
	println(io, "\nrotation_cart: ")
	print(io, symop.rotation_cart)
	println(io, "\ntranslation_frac: ", symop.translation_frac)
	println(io, "is_translation: ", symop.is_translation)
	println(io, "is_proper: ", symop.is_proper)
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
	map_sym_inv::Matrix{Int}    # [num_atoms, nsym] -> source atom index
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
		x_frac_vec = collect(eachcol(cell.x_frac))
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
		map_sym_inv = construct_map_sym_inv(map_sym)
		# generate map_p2s (primitive cell --> supercell)
		map_p2s = construct_map_p2s(spglib_data, cell, map_sym, symnum_translation)

		nat_prim = maximum(spglib_data.mapping_to_primitive)
		# generate map_s2p (supercell -> primitive cell)
		map_s2p = construct_map_s2p(cell, map_p2s, nat_prim, ntran)

		# collect atom indices in primitive cell from map_p2s
		atoms_in_prim = [map_p2s[i, 1] for i in 1:nat_prim]
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
			map_sym_inv,
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

		is_proper = det(rotation_cart) > 0
		translation_frac =
			(abs.(spglib_data.translations[i]) .>= tol) .* spglib_data.translations[i]

		symdata_elem = SymmetryOperation(
			SMatrix{3, 3, Float64, 9}(spglib_data.rotations[i]),
			SMatrix{3, 3, Float64, 9}(rotation_cart),
			SVector{3, Float64}(translation_frac),
			isapprox(spglib_data.rotations[i], I, atol = tol),
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
					# Calculate periodic distance in fractional coordinates
					diff = structure.supercell.x_frac[:, jat] - local_x_new
					local_tmp .= abs.(diff) .% 1.0
					# Take minimum distance considering periodic boundary
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
	construct_map_sym_inv(map_sym::Matrix{Int}) -> Matrix{Int}

Constructs the inverse mapping of `map_sym`.

For `map_sym[i, isym] = j` (atom `i` maps to atom `j` under symmetry operation `isym`),
`map_sym_inv[j, isym] = i` (atom `j` comes from atom `i` under symmetry operation `isym`).

# Arguments
- `map_sym::Matrix{Int}`: Symmetry mapping matrix where `map_sym[i, isym]` is the atom index
  that atom `i` maps to under symmetry operation `isym`

# Returns
- `Matrix{Int}`: Inverse mapping matrix where `map_sym_inv[i, isym]` is the source atom index
  that maps to atom `i` under symmetry operation `isym`

# Throws
- `ErrorException` if the inverse mapping cannot be constructed properly (e.g., non-bijective mapping)
"""
function construct_map_sym_inv(map_sym::AbstractMatrix{<:Integer})::Matrix{Int}
	num_atoms, nsym = size(map_sym)
	map_sym_inv = zeros(Int, num_atoms, nsym)

	# Process each symmetry operation
	for isym in 1:nsym
		# For each target atom i, find the source atom j such that map_sym[j, isym] == i
		for i in 1:num_atoms
			found = false
			for j in 1:num_atoms
				if map_sym[j, isym] == i
					if found
						error(
							"Non-bijective mapping detected: atom $i has multiple sources under symmetry operation $isym",
						)
					end
					map_sym_inv[i, isym] = j
					found = true
				end
			end
			if !found
				error(
					"Incomplete mapping: atom $i has no source under symmetry operation $isym",
				)
			end
		end
	end

	return map_sym_inv
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
	if any(==(0), map_p2s)
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
	undef_indices = findall(!, initialized)
	if !isempty(undef_indices)
		error("undef is detected in `map_s2p` variable: $undef_indices")
	end
	return map_s2p
end


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

