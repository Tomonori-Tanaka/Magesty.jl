module Write

using ..Version
using ..Structures

function write_scecoeffs2xml(structure::Structure, filename::String)
	try
		open(filename, "w") do f
			write(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
			write(f, "<version>$(Version.version_string())</version>\n")
			write(f, "<System>\n")
			write(f, "<NumberOfAtoms>$(structure.supercell.num_atoms)</NumberOfAtoms>\n")
			write(f, "<NumberOfElements>$(structure.supercell.num_elements)</NumberOfElements>\n")
			write(f, "    <LatticeVector>\n")
			write(
				f,
				"        <a1>$(structure.supercell.lattice_vectors[1, 1]) $(structure.supercell.lattice_vectors[1, 2]) $(structure.supercell.lattice_vectors[1, 3])</a1>\n",
			)
			write(
				f,
				"        <a2>$(structure.supercell.lattice_vectors[2, 1]) $(structure.supercell.lattice_vectors[2, 2]) $(structure.supercell.lattice_vectors[2, 3])</a2>\n",
			)
			write(
				f,
				"        <a3>$(structure.supercell.lattice_vectors[3, 1]) $(structure.supercell.lattice_vectors[3, 2]) $(structure.supercell.lattice_vectors[3, 3])</a3>\n",
			)
			write(f, "    </LatticeVector>\n")
			write(f, "</System>\n")
		end
	catch e
		throw(SystemError("Failed to write XML file: $filename"))
	end

end

end
