module Write

using ..Version
using ..Structures
using EzXML
using Printf

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
			atom_node["index"] = string(i)
			atom_node["element"] = structure.kd_name[structure.supercell.kd_int_list[i]]
		end

		# Write to file
		open(filename, "w") do f
			EzXML.prettyprint(f, doc)
		end
	catch e
		throw(SystemError("Failed to write XML file: $filename.\nError: $e"))
	end
end

end # module Write
