module Write

using ..Version
using ..Structures
using EzXML

"""
    write_scecoeffs2xml(structure::Structure, filename::String)

Write system information to an XML file using EzXML.

# Arguments
- `structure::Structure`: The structure containing system information
- `filename::String`: Output XML file name

# Throws
- `SystemError` if the file cannot be written
"""
function write_scecoeffs2xml(structure::Structure, filename::String)
	try
		# Create XML document
		doc = XMLDocument()
		root = setroot!(doc, ElementNode("Magesty"))

		# Add version information
		version_node = addelement!(root, "version")
		addtext!(version_node, Version.version_string())

		# Add system information
		system_node = addelement!(root, "System")
		
		# Add number of atoms and elements
		addelement!(system_node, "NumberOfAtoms", string(structure.supercell.num_atoms))
		addelement!(system_node, "NumberOfElements", string(structure.supercell.num_elements))

		# Add lattice vectors
		lattice_node = addelement!(system_node, "LatticeVector")
		for i in 1:3
			vector = structure.supercell.lattice_vectors[i, :]
			vector_str = join(string.(vector), " ")
			addelement!(lattice_node, "a$i", vector_str)
		end

		# Convert boolean values to 1 and 0
		periodicity_int = Int.(structure.supercell.periodicity)
		addelement!(system_node, "Periodicity", join(string.(periodicity_int), " "))

		# Add atomic coordinates
		positions_node = addelement!(system_node, "Positions")
		for i in 1:structure.supercell.num_atoms
			atom_node = addelement!(positions_node, "pos")
			setattribute!(atom_node, "index", string(i))
			setattribute!(atom_node, "element", structure.kd_name[structure.supercell.kd_int_list[i]])
			coord_str = join(string.(structure.supercell.positions[i, :]), " ")
			addtext!(atom_node, coord_str)
		end

		# Write to file
		write(filename, doc)
	catch e
		throw(SystemError("Failed to write XML file: $filename"))
	end
end

end
