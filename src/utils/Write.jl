module Write

using ..Version
using ..AtomicIndices
using ..Structures
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

	basis_set_node = addelement!(system_node, "SCEBasisSet")
	for (i, salc::SALC) in enumerate(basis_set.salc_list)
		salc_node = addelement!(basis_set_node, "SALC")
		salc_node["index"] = string(i)
		for (basis::IndicesUniqueList, coeff, multiplicity) in zip(salc.basisset, salc.coeffs, salc.multiplicity)
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
		addelement!(jphi_node, "ReferenceEnergy", string(optimize.bias_term))

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

end # module Write
