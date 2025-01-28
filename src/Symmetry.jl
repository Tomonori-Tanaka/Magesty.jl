"""
	module Symmetries

A module for handling symmetry operations in crystal structures using `Spglib`. It defines `SymmetryOperation` and `Symmetry` types to store symmetry-related data and mappings.

# Types
- **`SymmetryOperation`**: Represents a single symmetry operation with rotation and translation components.
- **`Symmetry`**: Represents the symmetry information of a system, including symmetry operations, mappings, and metadata.

# Functions
- `is_compatible_cart(rot_cart)`: Check if a rotation matrix is compatible with Cartesian axes.
- `Symmetry(system, tol)`: Generate a `Symmetry` object for a given system with a specified tolerance.
"""
module Symmetries

using LinearAlgebra
using StaticArrays
using Spglib

using ..AtomCells
using ..Systems

import Base: isless, show

export SymmetryOperation, Maps, Symmetry

"""
	struct SymmetryOperation

Represents a symmetry operation, including rotation and translation components.

# Fields
- `rotation_frac::SMatrix{3, 3, Float64, 9}`: The 3x3 rotation matrix in fractional coordinates.
- `rotation_cart::SMatrix{3, 3, Float64, 9}`: The 3x3 rotation matrix in Cartesian coordinates.
- `translation_frac::SMatrix{3, 1, Float64, 3}`: The 3x1 translation vector in fractional coordinates.
- `is_translation::Bool`: True if the operation is a pure translation.
- `is_translation_included::Bool`: True if the operation includes a translation.
- `is_compatible_cartesian::Bool`: True if the rotation is compatible with Cartesian axes.
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
	symop1_rot_flatten = [symop1.rotation_frac[i, j] for j in 1:3, i in 1:3]
	symop2_rot_flatten = [symop2.rotation_frac[i, j] for j in 1:3, i in 1:3]

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

Contains the symmetry information of a system.

# Fields
- `international_symbol::String`: International symbol of the space group.
- `nsym::Int`: Number of symmetry operations.
- `ntran::Int`: Number of pure translation operations.
- `nat_prim::Int`: Number of atoms in the primitive cell.
- `tol::Float64`: Tolerance for symmetry detection.
- `symdata::Vector{SymmetryOperation}`: List of symmetry operations.
- `map_sym::Matrix{Int}`: Maps atoms in the supercell to corresponding atoms under symmetry operations.
- `map_p2s::Matrix{Int}`: Maps atoms in the primitive cell to the supercell.
- `map_s2p::Vector{Maps}`: Maps atoms in the supercell to the primitive cell.
- `symnum_translation::Vector{Int}`: Indices of pure translation operations.

# Constructor
	Symmetry(system, tol::Float64)

Generate symmetry information for the given system using the specified tolerance.
"""
struct Symmetry
	international_symbol::String
	spacegroup_number::Int
	nsym::Int   # the number of symmetry operations
	ntran::Int  # the number of translational only operations
	nat_prim::Int   # the number of atoms in a primitive cell
	tol::Float64
	atoms_in_prim::Vector{Int}

	symdata::Vector{SymmetryOperation}
	map_sym::Matrix{Int}    # [num_atoms, nsym] -> corresponding atom index
	map_sym_cell::Array{AtomCell}# [atom, cell, isym] -> corresponding AtomCell instance
	map_p2s::Matrix{Int}    # [nat_prim, ntran] -> corresponding atom index
	map_s2p::Vector{Maps}   # [nat] -> corresponding atom index in primitive cel
	symnum_translation::Vector{Int} # contains the indice of translational only operations
end

function Symmetry(system::System, tol::Real)
	println("""
	========
	SYMMETRY
	========

	""")
	cell = system.supercell
	spglib_cell = Spglib.Cell(cell.lattice_vectors, cell.x_frac, cell.kd_int_list)
	spglib_data::Spglib.Dataset = get_dataset(spglib_cell)

	nsym::Int = spglib_data.n_operations
	if nsym == 0
		error(
			"Error in symmetry search: No symmetry operations found. Please check the input structure or tolerance setting.",
		)
	end

	ntran::Int = 0
	symnum_translation = Int[]
	for i in 1:nsym
		if isapprox(spglib_data.rotations[i], I, atol = tol)
			ntran += 1
			append!(symnum_translation, i)
		end
	end

	# construct symdata
	symdata = Vector{SymmetryOperation}(undef, nsym)
	for i in 1:nsym
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

	# construct mapping data
	natomtypes::Int = length(system.atomtype_group)
	map_sym = zeros(Int, cell.num_atoms, spglib_data.n_operations)
	map_sym_cell = Array{AtomCell}(
		undef,
		cell.num_atoms,
		27,# the number of total image cells
		spglib_data.n_operations,
	)
	initialized = falses(size(map_sym_cell))

	x_new::Vector{Float64} = Vector{Float64}(undef, 3)
	tmp::Vector{Float64} = Vector{Float64}(undef, 3)

	for isym in 1:spglib_data.n_operations
		for itype in 1:natomtypes
			for iat in system.atomtype_group[itype]
				x_new = Vector(spglib_data.rotations[isym] * cell.x_frac[:, iat])
				x_new = x_new + Vector(spglib_data.translations[isym])

				for jat in system.atomtype_group[itype]
					tmp = (abs.(cell.x_frac[:, jat] - x_new)) .% 1.0

					for (i, val) in enumerate(tmp)
						tmp[i] = min(val, 1 - val)
					end

					if norm(tmp) < tol
						map_sym[iat, isym] = jat
						for cell in 1:27# 27 is the total number of neighboring imaginary (virtual) cell including the central cell
							matched_cell = find_matching_image_cell(
								symdata[isym],
								system.x_image_frac,
								iat,
								cell,
								tol = tol,
							)
							map_sym_cell[iat, cell, isym] =
								AtomCell(jat, matched_cell)
							initialized[iat, cell, isym] = true
						end
						break
					end
				end
				if map_sym[iat, isym] == 0
					error(
						"gen_mapping_info: cannot find symmetry for operation number $isym, atom index $iat.",
					)
				elseif false in initialized[iat, :, isym]
					error("false is found in map_sym_cell at $iat, :, $cell")
				end
			end
		end
	end

	# generate map_p2s (primitive cell --> supercell)
	nat_prim = max(spglib_data.mapping_to_primitive...)
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


	# generate map_s2p (supercell -> primitive cell)
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


	if 0 in map_s2p
		error("something wrong in generating map_s2p")
	end

	atoms_in_prim = Int[map_p2s[i, 1] for i in 1:nat_prim]
	atoms_in_prim = sort(atoms_in_prim)

	symmetry = Symmetry(
		spglib_data.international_symbol,
		spglib_data.spacegroup_number,
		spglib_data.n_operations,
		ntran,
		nat_prim,
		tol,
		atoms_in_prim,
		symdata,
		map_sym,
		map_sym_cell,
		map_p2s,
		map_s2p,
		symnum_translation)

	print_symminfo_stdout(symmetry)

	return symmetry
end

"""
	find_matching_image_cell(symop::SymmetryOperation, x_image::AbstractArray{<:Real, 3}, atom::Integer, cell::Integer) -> Union{Nothing, Int}

Finds the matching image cell index for a given symmetry operation applied to an atom within a supercell.

# Arguments

- `symop::SymmetryOperation`: The symmetry operation to be applied, containing rotation and translation matrices.
- `x_image::AbstractArray{<:Real, 3}`: A 3D array containing fractional coordinates of atoms in different images.
- `atom::Integer`: The index of the atom to which the symmetry operation is applied.
- `cell::Integer`: The image cell index of the atom.

# Returns

- `Union{Nothing, Int}`: Returns the matching image cell index `m` if a unique match is found. Returns `nothing` if no match is found. Throws an error if multiple matches are found.

# Examples

```julia
# Assuming appropriate definitions for SymmetryOperation and x_image
symop = SymmetryOperation(rotation=eye(3), translation=[0.0, 0.0, 0.0])
x_image = rand(Float64, 3, 5, 5)  # Example 3D array
atom = 1
cell = 1

result = find_matching_image_cell(symop, x_image, atom, cell)
println(result)  # Example output: nothing or an integer index
"""
function find_matching_image_cell(
	symop::SymmetryOperation,
	x_image::AbstractArray{<:Real, 3},
	atom::Integer,
	cell::Integer,
	;
	tol::Real = 1e-5,
)::Int
	# Apply the symmetry operation to the specified atom and image cell
	x_moved = symop.rotation_frac * x_image[:, atom, cell] + symop.translation_frac
	matches = [
		(n, m) for n in 1:size(x_image, 2), m in 1:size(x_image, 3) if
		isapprox(x_image[:, n, m], x_moved; atol = tol)
	]

	if length(matches) == 0
		return -1
	elseif length(matches) == 1
		return matches[1][2]# cell_moved
	else
		error(
			"Multiple matching image cells found for atom $atom in image cell $cell: $matches",
		)
	end
end

is_in_centeringcell(xf, x_image_frac) = any(xf ≈ vec for vec in x_image_frac[:, :, 1])

function print_symminfo_stdout(symmetry::Symmetry)
	str = """
	Space group:  $(symmetry.international_symbol)  ($(symmetry.spacegroup_number))
	Number of symmetry operations = $(symmetry.nsym)

	"""
	if symmetry.ntran == 1
		str_prim = """
		Given system is a primitive cell.
		Primitive cell contains $(symmetry.nat_prim) atoms.

		"""
		str *= str_prim
	else
		str_supercell = """
		Given system is not a primitive cell.
		There are $(symmetry.ntran) translation operations.
		Primitive cell contains $(symmetry.nat_prim) atoms.

		"""
		str *= str_supercell
	end
	println(str)
end

function __write_symdata(
	symdata::AbstractVector,
	dir::AbstractString,
	filename::AbstractString,
)
	mkpath(dir)
	path = joinpath(dir, filename)
	open(path, "w") do io
		for (i, symdata) in enumerate(symdata)
			write(io, "operation: $i\n")
			write(io, "rotation (fractional)\n")
			for line in 1:3
				vec = Vector(symdata.rotation_frac[line, :])
				vec1 = vec[1]
				vec2 = vec[2]
				vec3 = vec[3]
				write(io, "$vec1\t$vec2\t$vec3\n")
			end
			write(io, "translation (fractional)\n")
			write(
				io,
				"$(symdata.translation_frac[1])\t$(symdata.translation_frac[2])\t$(symdata.translation_frac[3])\n",
			)
			write(io, "\n")
		end
	end
end

end
