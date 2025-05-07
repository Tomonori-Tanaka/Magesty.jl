module Write

using ..Magesty

function write_scecoeffs2xml(magesty::SpinCluster, filename::String)
	try
		open(filename, "w") do f
			write(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
			write(f, "<version>$(Magesty.VERSION)</version>\n")
			write(f, "<System>\n")
			write(f, "<NumberOfAtoms>$(magesty.structure.num_atoms)</NumberOfAtoms>\n")
			write(f, "<NumberOfElements>$(magesty.structure.num_elements)</NumberOfElements>\n")
			write(f, "    <LatticeVector>\n")
			write(
				f,
				"            <a1>$(magesty.structure.lattice_vectors[1, 1]) $(magesty.structure.lattice_vectors[1, 2]) $(magesty.structure.lattice_vectors[1, 3])</a1>\n",
			)
			write(
				f,
				"            <a2>$(magesty.structure.lattice_vectors[2, 1]) $(magesty.structure.lattice_vectors[2, 2]) $(magesty.structure.lattice_vectors[2, 3])</a2>\n",
			)
			write(
				f,
				"            <a3>$(magesty.structure.lattice_vectors[3, 1]) $(magesty.structure.lattice_vectors[3, 2]) $(magesty.structure.lattice_vectors[3, 3])</a3>\n",
			)
			write(f, "    </LatticeVector>\n")
			write(f, "</System>\n")
		end
	catch e
		throw(SystemError("Failed to write XML file: $filename"))
	end

end

end
