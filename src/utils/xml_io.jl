module XMLIO

using ..Version
using ..Structures
using ..Symmetries
using ..BasisSets
using ..Optimize
using EzXML
using Printf
using LinearAlgebra
using ..Basis
using ..SortedContainer

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
		addelement!(jphi_node, "ReferenceEnergy", @sprintf("%.15e", optimize.reference_energy))

		for (i, jphi) in enumerate(optimize.SCE)
			each_jphi_node = addelement!(jphi_node, "jphi", @sprintf("%.15e", jphi))
			each_jphi_node["salc_index"] = string(i)
		end
	end

	# Write to file
	open(filename, "w") do f
		EzXML.prettyprint(f, doc)
	end
end

"""
	write_xml(structure::Structure, symmetry::Symmetry, basis_set::BasisSet, filename::String)

Write System information to an XML file without optimization results (JPhi).
This is used for writing System objects.

# Arguments
- `structure::Structure`: The structure containing system information
- `symmetry::Symmetry`: The symmetry information
- `basis_set::BasisSet`: The basis set information
- `filename::String`: Output XML file name

# Throws
- `SystemError` if the file cannot be written
"""
function write_xml(
	structure::Structure,
	symmetry::Symmetry,
	basis_set::BasisSet,
	filename::String,
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

	# Write to file
	open(filename, "w") do f
		EzXML.prettyprint(f, doc)
	end
end

"""
	read_basisset_from_xml(xml_file::AbstractString) -> BasisSet

Read BasisSet from XML file. This reconstructs the basis set from saved SALC information,
avoiding the expensive SALC computation.

# Arguments
- `xml_file::AbstractString`: Path to XML file containing basis set information

# Returns
- `BasisSet`: Reconstructed basis set from XML file

# Throws
- `ErrorException` if the XML file format is invalid or missing required information
"""
function read_basisset_from_xml(
	xml_file::AbstractString,
)::BasisSet
	doc = readxml(xml_file)
	system_node = findfirst("//System", doc)
	if isnothing(system_node)
		throw(ArgumentError("<System> node not found in XML file: $xml_file"))
	end
	
	# Read angular momentum couplings (for coeff_tensor reconstruction)
	amc_dict = Dict{Tuple{Vector{Int}, Vector{Int}, Int}, Basis.AngularMomentumCouplingResult}()
	amc_node = findfirst("AngularMomentumCouplings", system_node)
	if !isnothing(amc_node)
		for coupling_node in findall("Coupling", amc_node)
			ls_str = coupling_node["ls"]
			Lseq_str = coupling_node["Lseq"]
			Lf = parse(Int, coupling_node["Lf"])
			
			ls = parse.(Int, split(ls_str))
			Lseq = isempty(Lseq_str) ? Int[] : parse.(Int, split(Lseq_str))
			
			# Read tensor shape
			tensor_shape_str = coupling_node["tensor_shape"]
			Mf_size = parse(Int, coupling_node["Mf_size"])
			site_dims = parse.(Int, split(tensor_shape_str))
			
			# Reconstruct coeff_tensor
			full_shape = vcat(site_dims, [Mf_size])
			coeff_tensor = zeros(Float64, full_shape...)
			
			for mf_node in findall("Mf", coupling_node)
				mf_idx = parse(Int, mf_node["index"])
				tensor_str_node = findfirst("coeff_tensor", mf_node)
				if isnothing(tensor_str_node)
					error("coeff_tensor node not found for Mf index $mf_idx")
				end
				tensor_str = nodecontent(tensor_str_node)
				tensor_values = parse.(Float64, split(tensor_str))
				
				# Reshape to match site dimensions
				tensor_slice = reshape(tensor_values, site_dims...)
				
				# Set the slice in coeff_tensor using proper indexing
				# Create indices: [:, :, ..., :, mf_idx] for N site dims + 1 Mf dim
				indices = Vector{Any}([Colon() for _ in 1:length(site_dims)])
				push!(indices, mf_idx)
				coeff_tensor[indices...] = tensor_slice
			end
			
			amc = Basis.AngularMomentumCouplingResult(ls, Lseq, Lf, coeff_tensor)
			amc_dict[(ls, Lseq, Lf)] = amc
		end
	end
	
	# Read SCE basis set
	basis_set_node = findfirst("SCEBasisSet", system_node)
	if isnothing(basis_set_node)
		throw(ArgumentError("<SCEBasisSet> node not found in XML file: $xml_file"))
	end
	
	salc_list = Vector{Vector{Basis.CoupledBasis_with_coefficient}}()
	
	for salc_node in findall("SALC", basis_set_node)
		key_group = Vector{Basis.CoupledBasis_with_coefficient}()
		
		for basis_node in findall("basis", salc_node)
			# Parse attributes
			atoms = parse.(Int, split(basis_node["atoms"]))
			ls = parse.(Int, split(basis_node["ls"]))
			Lseq_str = basis_node["Lseq"]
			Lseq = isempty(Lseq_str) ? Int[] : parse.(Int, split(Lseq_str))
			multiplicity = parse(Int, basis_node["multiplicity"])
			
			# Get Lf from SALC node
			Lf = parse(Int, salc_node["Lf"])
			
			# Parse coefficient vector
			coeff_str = nodecontent(basis_node)
			coefficient = parse.(Float64, split(coeff_str))
			
			# Find corresponding coeff_tensor from angular momentum couplings
			key = (ls, Lseq, Lf)
			if haskey(amc_dict, key)
				# Make a copy of coeff_tensor to avoid sharing references
				coeff_tensor = copy(amc_dict[key].coeff_tensor)
			else
				# If coeff_tensor not found, try to reconstruct from tesseral_coupled_bases_from_tesseral_bases
				# This is a fallback - ideally all coeff_tensors should be in XML
				coupled_basis_list = Basis.tesseral_coupled_bases_from_tesseral_bases(ls, atoms)
				# Find matching CoupledBasis
				coeff_tensor = nothing
				for cb in coupled_basis_list
					if cb.Lf == Lf && cb.Lseq == Lseq
						coeff_tensor = copy(cb.coeff_tensor)
						break
					end
				end
				if isnothing(coeff_tensor)
					error("Could not find coeff_tensor for ls=$ls, Lseq=$Lseq, Lf=$Lf")
				end
			end
			
			# Create CoupledBasis_with_coefficient
			cbc = Basis.CoupledBasis_with_coefficient(
				ls,
				Lf,
				Lseq,
				atoms,
				coeff_tensor,
				coefficient,
				multiplicity,
			)
			push!(key_group, cbc)
		end
		
		if !isempty(key_group)
			push!(salc_list, key_group)
		end
	end
	
	# Reconstruct coupled_basislist from salc_list
	coupled_basislist = SortedCountingUniqueVector{Basis.CoupledBasis}()
	for key_group in salc_list
		for cbc in key_group
			cb = Basis.CoupledBasis(
				cbc.ls,
				cbc.Lf,
				cbc.Lseq,
				cbc.atoms,
				cbc.coeff_tensor,
			)
			push!(coupled_basislist, cb)
		end
	end
	
	# Reconstruct angular_momentum_couplings
	angular_momentum_couplings = collect(values(amc_dict))
	
	return BasisSet(coupled_basislist, salc_list, angular_momentum_couplings)
end



end # module XMLIO

