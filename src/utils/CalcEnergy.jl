module CalcEnergy
using EzXML
using LinearAlgebra
using ..SALCs
using ..Symmetries
using ..Optimize

export calc_energy

function calc_energy(
	salc_list::AbstractVector{SALC},
	spin_config::AbstractMatrix{<:Real},
	symmetry::Symmetry,
	optimize::SCEOptimizer,
)::Float64
	if size(spin_config, 1) != 3
		throw(ArgumentError("spin_config must be a 3xN matrix"))
	end

	design_list = Vector{Float64}(undef, length(salc_list))

	for i in eachindex(salc_list)
		design_list[i] = Optimize.calc_X_element_energy(
			salc_list[i],
			spin_config,
			symmetry,
		)
	end

	return dot(design_list, optimize.SCE) + optimize.reference_energy
end

# function calc_energy(xml_path::AbstractString, spin_config::AbstractMatrix)
# 	try
# 		xml_parsed = EzXML.readxml(xml_path)
# 	catch e
# 		throw(ArgumentError("Failed to read XML file: $xml_path, Error: $e"))
# 	end

# 	root_node = root(xml_parsed) # root node name is "Magesty"
# 	system_node = findfirst(root_node, "System")
# 	num_atoms = parse(Int, nodecontent(findfirst(system_node, "NumberOfAtoms")))

# 	# check dimensions of spin_config
# 	if size(spin_config, 2) != num_atoms
# 		msg =
# 			"""Number of atoms in XML file and spin_config (number of columns) do not match.\n""" *
# 			"""  spin configuration size: $(size(spin_config))\n""" *
# 			"""  number of atoms: $num_atoms"""
# 		throw(ArgumentError(msg))
# 	end
#     sce_basis_set_node = findfirst(root_node, "SCEBasisSet")
# 	num_salc = parse(Int, nodecontent(findfirst(sce_basis_set_node, "NumberOfSALCs")))

# 	salc_list = Vector{SALC}()
# 	for i in 1:num_salc
# 		salc_node = findfirst(sce_basis_set_node, "SALC[index='$i']")
# 		if salc_node === nothing
# 			throw(ArgumentError("SALC with index $i not found"))
# 		end
# 		push!(salc_list, SALC(salc_node))
# 	end

# 	return salc_list
# end

end # module CalcEnergy