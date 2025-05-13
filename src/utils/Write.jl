module Write

using ..Version
using ..AtomicIndices
using ..Structures
using ..Symmetries
using ..BasisSets
using ..SALCs
using ..Optimize
using EzXML
using Printf

"""
	write_sce2xml(structure::Structure, basis_set::BasisSet, filename::String)

Write system information to an XML file using EzXML.

# Arguments
- `structure::Structure`: The structure containing system information
- `filename::String`: Output XML file name

# Throws
- `SystemError` if the file cannot be written
"""
function write_sce2xml(structure::Structure,
	symmetry::Symmetry,
	basis_set::BasisSet,
	optimize::SCEOptimizer,
	filename::String
	;
	write_jphi::Bool = true,
)
	# Create XML document
	doc = XMLDocument()
	root = setroot!(doc, ElementNode("Magesty"))

	# Add version information
	addelement!(root, "version", Version.version_string())

	# Add system information
	system_node = addelement!(root, "System")

	# Add number of atoms and elements
	addelement!(system_node, "NumberOfAtoms", string(structure.supercell.num_atoms))
	addelement!(system_node, "NumberOfElements", string(structure.supercell.num_elements))

	# Add lattice vectors
	lattice_node = addelement!(system_node, "LatticeVector")
	for i in 1:3
		vector = structure.supercell.lattice_vectors[i, :]
		vector_str = join([@sprintf("%.8e", x) for x in vector], " ")
		addelement!(lattice_node, "a$i", vector_str)
	end

	# Convert boolean values to 1 and 0
	periodicity_int = Int.(structure.is_periodic)
	addelement!(system_node, "Periodicity", join(string.(periodicity_int), " "))

	# Add atomic coordinates
	positions_node = addelement!(system_node, "Positions")
	for (i, coord) in enumerate(eachcol(structure.supercell.x_frac))
		coord_str = join([@sprintf("%.12f", x) for x in coord], " ")
		atom_node = addelement!(positions_node, "pos", coord_str)
		atom_node["atom_index"] = string(i)
		atom_node["element"] = structure.kd_name[structure.supercell.kd_int_list[i]]
	end

	# Add symmetry information
	symmetry_node = addelement!(system_node, "Symmetry")
	addelement!(symmetry_node, "NumberOfTranslations", string(symmetry.ntran))
	translations_node = addelement!(symmetry_node, "Translations")
	for (i, isym_trans) in enumerate(symmetry.symnum_translation)
		for iat_prim in symmetry.atoms_in_prim
			map_node = addelement!(translations_node, "map", string(symmetry.map_sym[iat_prim, isym_trans]))
			map_node["trans"] = string(i)
			map_node["atom"] = string(iat_prim)
		end
	end

	basis_set_node = addelement!(system_node, "SCEBasisSet")
	for (i, salc::SALC) in enumerate(basis_set.salc_list)
		salc_node = addelement!(basis_set_node, "SALC")
		salc_node["index"] = string(i)
		for (basis::IndicesUniqueList, coeff, multiplicity) in
			zip(salc.basisset, salc.coeffs, salc.multiplicity)
			basis_node = addelement!(salc_node, "basis", string(coeff))
			basis_node["multiplicity"] = string(multiplicity)
			for (j, indices::Indices) in enumerate(basis)
				basis_node["index-$j"] =
					string(indices.atom) * " " * string(indices.l) * " " * string(indices.m) *
					" " * string(indices.cell)
			end
		end
	end

	if write_jphi
		jphi_node = addelement!(root, "JPhi")
		jphi_node["unit"] = "eV"
		addelement!(jphi_node, "ReferenceEnergy", string(optimize.reference_energy))

		for (i, jphi) in enumerate(optimize.SCE)
			each_jphi_node = addelement!(jphi_node, "jphi", string(jphi))
			each_jphi_node["salc_index"] = string(i)
		end
	end

	# Write to file
	open(filename, "w") do f
		EzXML.prettyprint(f, doc)
	end
end

function write_energy_info(optimize::SCEOptimizer, filename::AbstractString = "energy.txt")
	# Input validation
	if isempty(optimize.spinconfig_dataset.spinconfigs)
		@warn "No spin configurations found in optimizer"
		return
	end

	# Prepare data
	observed_energy_list =
		[spinconfig.energy for spinconfig in optimize.spinconfig_dataset.spinconfigs]
	predicted_energy_list = optimize.predicted_energy_list

	# Format header
	digits_index = length(string(length(observed_energy_list)))
	header = "# Index:" * " "^(digits_index - 1) * "Observed_Energy" * " "^4 * "Predicted_Energy"

	write_list_to_file(observed_energy_list, predicted_energy_list, filename, header)
end

function write_lmf_flattened(
	optimize::SCEOptimizer,
	filename::AbstractString = "lmf_flattened.txt",
)
	# Input validation
	if isempty(optimize.spinconfig_dataset.spinconfigs)
		@warn "No spin configurations found in optimizer"
		return
	end

	# Prepare data
	observed_magfield_vertical_list = optimize.observed_magfield_vertical_flattened_list
	predicted_magfield_vertical_list = optimize.predicted_magfield_vertical_flattened_list

	# Format header
	digits_index = length(string(length(observed_magfield_vertical_list)))
	header =
		"# Index:" * " "^(digits_index - 1) * "Observed_local_magnetic_field" * " "^4 *
		"Predicted_local_magnetic_field"

	write_list_to_file(
		observed_magfield_vertical_list,
		predicted_magfield_vertical_list,
		filename,
		header,
	)
end

"""
	write_list_to_file(data_list::AbstractVector, predicted_list::AbstractVector, filename::AbstractString, header::AbstractString)

Write observed and predicted data to a file with a common format.

# Arguments
- `data_list`: Vector of observed data
- `predicted_list`: Vector of predicted data
- `filename`: Output file name
- `header`: Header string for the output file
"""
function write_list_to_file(
	data_list::AbstractVector,
	predicted_list::AbstractVector,
	filename::AbstractString,
	header::AbstractString,
)
	# Check array lengths
	if length(data_list) != length(predicted_list)
		error("Length mismatch between observed and predicted lists")
	end

	# Format settings
	digits_index = length(string(length(data_list)))

	# Write to file
	try
		open(filename, "w") do f
			# Write header
			println(f, header)

			# Write data
			for (i, (obs, pred)) in enumerate(zip(data_list, predicted_list))
				str = @sprintf(
					"%*d    %15.10f    %15.10f\n",
					digits_index,
					i,
					obs,
					pred
				)
				write(f, str)
			end
		end
	catch e
		@error "Failed to write lists to file" exception = (e, catch_backtrace())
		rethrow(e)
	end
end

end # module Write
