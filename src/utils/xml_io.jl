module XMLIO

using ..Version
using ..AtomicIndices
using ..Structures
using ..Symmetries
using ..BasisSets
using ..Optimize
using EzXML
using Printf

"""
	write_xml(structure::Structure, basis_set::BasisSet, filename::String)

Write system information to an XML file using EzXML.

# Arguments
- `structure::Structure`: The structure containing system information
- `filename::String`: Output XML file name

# Throws
- `SystemError` if the file cannot be written
"""
function write_xml(structure::Structure,
	symmetry::Symmetry,
	basis_set::BasisSet,
	optimize::Optimizer,
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
			map_node = addelement!(
				translations_node,
				"map",
				string(symmetry.map_sym[iat_prim, isym_trans]),
			)
			map_node["trans"] = string(i)
			map_node["atom"] = string(iat_prim)
		end
	end

	# Add SCE basis set information (CoupledBasis_with_coefficient)
	basis_set_node = addelement!(system_node, "SCEBasisSet")
	num_salc = length(basis_set.salc_list)
	basis_set_node["num_salc"] = string(num_salc)
	for (i, key_group) in enumerate(basis_set.salc_list)
		# Skip empty groups for safety
		if isempty(key_group)
			continue
		end

		# All cbc in the same key_group share Lf/body by construction
		first_cbc = key_group[1]

		salc_node = addelement!(basis_set_node, "SALC")
		salc_node["index"] = string(i)
		salc_node["num_basis"] = string(length(key_group))
		salc_node["body"] = string(length(first_cbc.atoms))
		salc_node["Lf"] = string(first_cbc.Lf)

		for (basis_idx, cbc) in enumerate(key_group)
			# coefficient is a vector over Mf; store as space-separated string
			coeff_str = join(string.(cbc.coefficient), " ")
			basis_node = addelement!(salc_node, "basis", coeff_str)
			basis_node["multiplicity"] = string(cbc.multiplicity)
			basis_node["atoms"] = join(string.(cbc.atoms), " ")
			basis_node["ls"] = join(string.(cbc.ls), " ")
			basis_node["Lseq"] = join(string.(cbc.Lseq), " ")
			basis_node["index"] = string(basis_idx)
		end
	end

	# Add angular momentum coupling results
	if !isempty(basis_set.angular_momentum_couplings)
		amc_node = addelement!(system_node, "AngularMomentumCouplings")
		amc_node["num_couplings"] = string(length(basis_set.angular_momentum_couplings))
		for (i, amc) in enumerate(basis_set.angular_momentum_couplings)
			coupling_node = addelement!(amc_node, "Coupling")
			coupling_node["index"] = string(i)
			coupling_node["ls"] = join(string.(amc.ls), " ")
			coupling_node["Lseq"] = join(string.(amc.Lseq), " ")
			coupling_node["Lf"] = string(amc.Lf)
			# Store tensor dimensions (excluding Mf dimension)
			tensor_shape = size(amc.coeff_tensor)
			site_dims = tensor_shape[1:(end-1)]  # All dimensions except the last (Mf)
			Mf_size = tensor_shape[end]  # Last dimension is Mf
			coupling_node["tensor_shape"] = join(string.(site_dims), " ")
			coupling_node["Mf_size"] = string(Mf_size)
			
			# Output tensor for each Mf value separately
			# Mf values range from -Lf to +Lf in tesseral basis
			# Tesseral basis indexing: m_idx = 1 corresponds to Mf = -Lf, m_idx = 2 corresponds to Mf = -Lf+1, etc.
			# So Mf = m_idx - Lf - 1 for tesseral basis
			for mf_idx in 1:Mf_size
				# Convert tesseral index to Mf value: Mf = m_idx - Lf - 1
				Mf_value = mf_idx - amc.Lf - 1
				mf_node = addelement!(coupling_node, "Mf")
				mf_node["value"] = string(Mf_value)
				mf_node["index"] = string(mf_idx)
				
				# Extract tensor slice for this Mf value
				# Use selectdim to get the last dimension at index mf_idx
				tensor_slice = selectdim(amc.coeff_tensor, ndims(amc.coeff_tensor), mf_idx)
				# Flatten the slice and store as space-separated string
				tensor_flat = vec(tensor_slice)
				tensor_str = join([@sprintf("%.15e", x) for x in tensor_flat], " ")
				addelement!(mf_node, "coeff_tensor", tensor_str)
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



end # module XMLIO

