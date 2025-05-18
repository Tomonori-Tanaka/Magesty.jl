module CalcEnergy
using EzXML
using ..SALCs

export calc_energy

function calc_energy(xml_path::AbstractString, spin_config::AbstractMatrix)
	try
		xml_parsed = EzXML.readxml(xml_path)
	catch e
		throw(ArgumentError("Failed to read XML file: $xml_path, Error: $e"))
	end

	root_node = root(xml_parsed) # root node name is "Magesty"
	system_node = findfirst(root_node, "System")
	num_atoms = parse(Int, nodecontent(findfirst(system_node, "NumberOfAtoms")))

	# check dimensions of spin_config
	if size(spin_config, 2) != num_atoms
		msg =
			"""Number of atoms in XML file and spin_config (number of columns) do not match.\n""" *
			"""  spin configuration size: $(size(spin_config))\n""" *
			"""  number of atoms: $num_atoms"""
		throw(ArgumentError(msg))
	end
    sce_basis_set_node = findfirst(root_node, "SCEBasisSet")
	num_salc = parse(Int, nodecontent(findfirst(sce_basis_set_node, "NumberOfSALCs")))

	salc_list = Vector{SALC}()
	for i in 1:num_salc
		salc_node = findfirst(sce_basis_set_node, "SALC[index='$i']")
		if salc_node === nothing
			throw(ArgumentError("SALC with index $i not found"))
		end
		push!(salc_list, SALC(salc_node))
	end

	return salc_list
end

end # module CalcEnergy