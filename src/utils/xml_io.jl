module XMLIO

using ..Version
using ..Structures
using ..Symmetries
using ..SALCBases
using EzXML
using Printf
using LinearAlgebra
using ..Basis
using ..SortedCounters: SortedCounter

# XML schema constants. Shared between the write_* and read_* paths so
# a schema change here cannot drift between the two sides. Renaming a
# const is the canonical way to evolve the schema.
const TAG_ROOT             = "Magesty"
const TAG_VERSION          = "version"

const TAG_SYSTEM           = "System"
const TAG_NUM_ATOMS        = "NumberOfAtoms"
const TAG_NUM_ELEMENTS     = "NumberOfElements"
const TAG_LATTICE          = "LatticeVector"
const TAG_PERIODICITY      = "Periodicity"
const TAG_POSITIONS        = "Positions"

const TAG_SYMMETRY         = "Symmetry"
const TAG_NUM_TRANSLATIONS = "NumberOfTranslations"
const TAG_TRANSLATIONS     = "Translations"
const ATTR_TOLERANCE_SYM   = "tolerance_sym"

const TAG_SCE_BASIS        = "SCEBasis"
const TAG_SALC             = "SALC"
const TAG_BASIS            = "basis"
const ATTR_ISOTROPY        = "isotropy"

const TAG_AMC              = "AngularMomentumCouplings"
const TAG_COUPLING         = "Coupling"
const TAG_MF               = "Mf"
const TAG_COEFF_TENSOR     = "coeff_tensor"

const TAG_JPHI_BLOCK       = "JPhi"
const TAG_REF_ENERGY       = "ReferenceEnergy"
const TAG_JPHI             = "jphi"

const ATTR_LF              = "Lf"
const ATTR_LSEQ            = "Lseq"
const ATTR_MF_SIZE         = "Mf_size"

# Format helpers. `@sprintf` requires a literal format string at parse
# time, so each well-known format gets a one-liner wrapper.
fmt_lattice(x::Real)    = @sprintf("%.8e", x)
fmt_fractional(x::Real) = @sprintf("%.12f", x)
fmt_tensor(x::Real)     = @sprintf("%.15e", x)

"""
	_write_system_subtree!(root, structure, symmetry, salc_basis, isotropy)

Append the `<System>` subtree (structure, symmetry, SALC basis, angular
momentum couplings) under `root`. Shared by the basis-only and
basis+coefficient writer paths so the schema cannot drift between them.
"""
function _write_system_subtree!(
	root::EzXML.Node,
	structure::Structure,
	symmetry::Symmetry,
	salc_basis::SALCBasis,
	isotropy::Bool,
)
	# Add system information
	system_node = addelement!(root, TAG_SYSTEM)

	# Add number of atoms and elements
	addelement!(system_node, TAG_NUM_ATOMS, string(structure.supercell.num_atoms))
	addelement!(system_node, TAG_NUM_ELEMENTS, string(structure.supercell.num_elements))

	# Add lattice vectors
	lattice_node = addelement!(system_node, TAG_LATTICE)
	for i in 1:3
		vector = structure.supercell.lattice_vectors[i, :]
		vector_str = join([fmt_lattice(x) for x in vector], " ")
		addelement!(lattice_node, "a$i", vector_str)
	end

	# Convert boolean values to 1 and 0
	periodicity_int = Int.(structure.is_periodic)
	addelement!(system_node, TAG_PERIODICITY, join(string.(periodicity_int), " "))

	# Add atomic coordinates
	positions_node = addelement!(system_node, TAG_POSITIONS)
	for (i, coord) in enumerate(eachcol(structure.supercell.x_frac))
		coord_str = join([fmt_fractional(x) for x in coord], " ")
		atom_node = addelement!(positions_node, "pos", coord_str)
		atom_node["atom_index"] = string(i)
		atom_node["element"] = structure.kd_name[structure.supercell.kd_int_list[i]]
	end

	# Add symmetry information. `tolerance_sym` is the only symmetry datum
	# consumed on load — the rest of the subtree (translation map) is
	# emitted for human / external-tool consumption.
	symmetry_node = addelement!(system_node, TAG_SYMMETRY)
	symmetry_node[ATTR_TOLERANCE_SYM] = string(symmetry.tol)
	addelement!(symmetry_node, TAG_NUM_TRANSLATIONS, string(symmetry.ntran))
	translations_node = addelement!(symmetry_node, TAG_TRANSLATIONS)
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
	salc_basis_node = addelement!(system_node, TAG_SCE_BASIS)
	salc_basis_node[ATTR_ISOTROPY] = string(isotropy)
	num_salc = length(salc_basis.salc_list)
	salc_basis_node["num_salc"] = string(num_salc)
	for (i, key_group) in enumerate(salc_basis.salc_list)
		# Skip empty groups for safety
		if isempty(key_group)
			continue
		end

		# All cbc in the same key_group share Lf/body by construction
		first_cbc = key_group[1]

		salc_node = addelement!(salc_basis_node, TAG_SALC)
		salc_node["index"] = string(i)
		salc_node["num_basis"] = string(length(key_group))
		salc_node["body"] = string(length(first_cbc.atoms))
		salc_node[ATTR_LF] = string(first_cbc.Lf)

		for (basis_idx, cbc) in enumerate(key_group)
			# Coefficient is a vector over Mf; store as space-separated text.
			# Normalize negative zero to +0.0 so the textual output is
			# bit-stable across BLAS / LAPACK implementations — eigensolvers
			# on different platforms produce -0.0 vs +0.0 for the same
			# numerical eigenvector (the two are IEEE-equal but `string(...)`
			# preserves the sign bit).
			coeff_str = join(
				(string(iszero(x) ? 0.0 : x) for x in cbc.coefficient),
				" ",
			)
			basis_node = addelement!(salc_node, TAG_BASIS, coeff_str)
			basis_node["multiplicity"] = string(cbc.multiplicity)
			basis_node["atoms"] = join(string.(cbc.atoms), " ")
			basis_node["ls"] = join(string.(cbc.ls), " ")
			basis_node[ATTR_LSEQ] = join(string.(cbc.Lseq), " ")
			basis_node["index"] = string(basis_idx)
		end
	end

	# Add angular momentum coupling results
	if !isempty(salc_basis.angular_momentum_couplings)
		amc_node = addelement!(system_node, TAG_AMC)
		amc_node["num_couplings"] = string(length(salc_basis.angular_momentum_couplings))
		for (i, amc) in enumerate(salc_basis.angular_momentum_couplings)
			coupling_node = addelement!(amc_node, TAG_COUPLING)
			coupling_node["index"] = string(i)
			coupling_node["ls"] = join(string.(amc.ls), " ")
			coupling_node[ATTR_LSEQ] = join(string.(amc.Lseq), " ")
			coupling_node[ATTR_LF] = string(amc.Lf)
			# Store tensor dimensions (excluding Mf dimension)
			tensor_shape = size(amc.coeff_tensor)
			site_dims = tensor_shape[1:(end-1)]  # All dimensions except the last (Mf)
			Mf_size = tensor_shape[end]  # Last dimension is Mf
			coupling_node["tensor_shape"] = join(string.(site_dims), " ")
			coupling_node[ATTR_MF_SIZE] = string(Mf_size)

			# Output tensor for each Mf value separately
			# Mf values range from -Lf to +Lf in tesseral basis
			# Tesseral basis indexing: m_idx = 1 corresponds to Mf = -Lf, m_idx = 2 corresponds to Mf = -Lf+1, etc.
			# So Mf = m_idx - Lf - 1 for tesseral basis
			for mf_idx in 1:Mf_size
				# Convert tesseral index to Mf value: Mf = m_idx - Lf - 1
				Mf_value = mf_idx - amc.Lf - 1
				mf_node = addelement!(coupling_node, TAG_MF)
				mf_node["value"] = string(Mf_value)
				mf_node["index"] = string(mf_idx)

				# Extract tensor slice for this Mf value
				# Use selectdim to get the last dimension at index mf_idx
				tensor_slice = selectdim(amc.coeff_tensor, ndims(amc.coeff_tensor), mf_idx)
				# Flatten the slice and store as space-separated string
				tensor_flat = vec(tensor_slice)
				tensor_str::String = string(join([fmt_tensor(x) for x in tensor_flat], " "))
				addelement!(mf_node, TAG_COEFF_TENSOR, tensor_str)
			end
		end
	end

	return nothing
end

"""
	_write_jphi_block!(root, j0, jphi)

Append the `<JPhi>` block (reference energy + per-SALC coefficients)
under `root`.
"""
function _write_jphi_block!(root::EzXML.Node, j0::Real, jphi::AbstractVector{<:Real})
	jphi_node = addelement!(root, TAG_JPHI_BLOCK)
	jphi_node["unit"] = "eV"
	addelement!(jphi_node, TAG_REF_ENERGY, fmt_tensor(j0))

	for (i, coeff) in enumerate(jphi)
		each_jphi_node = addelement!(jphi_node, TAG_JPHI, fmt_tensor(coeff))
		each_jphi_node["salc_index"] = string(i)
	end
	return nothing
end

"""
	write_basis_xml(structure, symmetry, salc_basis, isotropy, filename)

Write an `SCEBasis` to XML: structure, symmetry parameters, and SALC
basis. No `<JPhi>` block. `XMLIO` is included before the `SCEBasis` type
is defined, so this takes loose component arguments; the `SCEBasis`
unpacking happens in `Magesty.save`.

# Throws
- `SystemError` if the file cannot be written.
"""
function write_basis_xml(
	structure::Structure,
	symmetry::Symmetry,
	salc_basis::SALCBasis,
	isotropy::Bool,
	filename::AbstractString,
)
	doc = XMLDocument()
	root = setroot!(doc, ElementNode(TAG_ROOT))
	addelement!(root, TAG_VERSION, Version.version_string())
	_write_system_subtree!(root, structure, symmetry, salc_basis, isotropy)
	open(filename, "w") do f
		EzXML.prettyprint(f, doc)
	end
end

"""
	write_model_xml(structure, symmetry, salc_basis, isotropy, j0, jphi, filename)

Write an `SCEModel` to XML: the `SCEBasis` subtree plus the fitted
`<JPhi>` block (reference energy `j0` and SCE coefficients `jphi`). Takes
loose component arguments for the same include-order reason as
`write_basis_xml`; the `SCEModel` unpacking happens in `Magesty.save`.

# Throws
- `SystemError` if the file cannot be written.
"""
function write_model_xml(
	structure::Structure,
	symmetry::Symmetry,
	salc_basis::SALCBasis,
	isotropy::Bool,
	j0::Real,
	jphi::AbstractVector{<:Real},
	filename::AbstractString,
)
	doc = XMLDocument()
	root = setroot!(doc, ElementNode(TAG_ROOT))
	addelement!(root, TAG_VERSION, Version.version_string())
	_write_system_subtree!(root, structure, symmetry, salc_basis, isotropy)
	_write_jphi_block!(root, j0, jphi)
	open(filename, "w") do f
		EzXML.prettyprint(f, doc)
	end
end

"""
	read_salcbasis_from_xml(xml_file::AbstractString) -> SALCBasis

Read SALCBasis from XML file. This reconstructs the basis set from saved SALC information,
avoiding the expensive SALC computation.

# Arguments
- `xml_file::AbstractString`: Path to XML file containing basis set information

# Returns
- `SALCBasis`: Reconstructed basis set from XML file

# Throws
- `ErrorException` if the XML file format is invalid or missing required information
"""
function read_salcbasis_from_xml(
	xml_file::AbstractString,
)::SALCBasis
	doc = readxml(xml_file)
	system_node = findfirst("//" * TAG_SYSTEM, doc)
	if isnothing(system_node)
		throw(ArgumentError("<$(TAG_SYSTEM)> node not found in XML file: $xml_file"))
	end

	# Read angular momentum couplings (for coeff_tensor reconstruction).
	# `amc_dict` is the keyed lookup used during SALC reconstruction;
	# `amc_ordered` preserves the file order so a write -> read -> write
	# cycle reproduces the `<AngularMomentumCouplings>` block byte-for-byte.
	amc_dict = Dict{Tuple{Vector{Int}, Vector{Int}, Int}, Basis.AngularMomentumCouplingResult}()
	amc_ordered = Vector{Basis.AngularMomentumCouplingResult}()
	amc_node = findfirst(TAG_AMC, system_node)
	if !isnothing(amc_node)
		for coupling_node in findall(TAG_COUPLING, amc_node)
			ls_str = coupling_node["ls"]
			Lseq_str = coupling_node[ATTR_LSEQ]
			Lf = parse(Int, coupling_node[ATTR_LF])

			ls = parse.(Int, split(ls_str))
			Lseq = isempty(Lseq_str) ? Int[] : parse.(Int, split(Lseq_str))

			# Read tensor shape
			tensor_shape_str = coupling_node["tensor_shape"]
			Mf_size = parse(Int, coupling_node[ATTR_MF_SIZE])
			site_dims = parse.(Int, split(tensor_shape_str))

			# Reconstruct coeff_tensor
			full_shape = vcat(site_dims, [Mf_size])
			coeff_tensor = zeros(Float64, full_shape...)

			for mf_node in findall(TAG_MF, coupling_node)
				mf_idx = parse(Int, mf_node["index"])
				tensor_str_node = findfirst(TAG_COEFF_TENSOR, mf_node)
				if isnothing(tensor_str_node)
					error("<$(TAG_COEFF_TENSOR)> node not found for Mf index $mf_idx")
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
			push!(amc_ordered, amc)
		end
	end

	# Read SCE basis set
	salc_basis_node = findfirst(TAG_SCE_BASIS, system_node)
	if isnothing(salc_basis_node)
		throw(ArgumentError("<$(TAG_SCE_BASIS)> node not found in XML file: $xml_file"))
	end

	salc_list = Vector{Vector{Basis.CoupledBasis_with_coefficient}}()

	for salc_node in findall(TAG_SALC, salc_basis_node)
		key_group = Vector{Basis.CoupledBasis_with_coefficient}()

		for basis_node in findall(TAG_BASIS, salc_node)
			# Parse attributes
			atoms = parse.(Int, split(basis_node["atoms"]))
			ls = parse.(Int, split(basis_node["ls"]))
			Lseq_str = basis_node[ATTR_LSEQ]
			Lseq = isempty(Lseq_str) ? Int[] : parse.(Int, split(Lseq_str))
			multiplicity = parse(Int, basis_node["multiplicity"])

			# Get Lf from SALC node
			Lf = parse(Int, salc_node[ATTR_LF])

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
	coupled_basislist = SortedCounter{Basis.CoupledBasis}()
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

	# Reconstruct angular_momentum_couplings, preserving file order.
	return SALCBasis(coupled_basislist, salc_list, amc_ordered)
end

"""
	_read_tolerance_sym(xml_file::AbstractString) -> Float64

Read the `tolerance_sym` attribute from the `<Symmetry>` node. The XML
reader recomputes `Symmetry` from `structure` + `tolerance_sym` rather
than deserializing it, so this attribute is required.
"""
function _read_tolerance_sym(xml_file::AbstractString)::Float64
	doc = readxml(xml_file)
	symmetry_node = findfirst("//" * TAG_SYMMETRY, doc)
	if isnothing(symmetry_node)
		throw(ArgumentError("<$(TAG_SYMMETRY)> node not found in XML file: $xml_file"))
	end
	if !haskey(symmetry_node, ATTR_TOLERANCE_SYM)
		throw(ArgumentError(
			"<$(TAG_SYMMETRY)> node is missing the '$(ATTR_TOLERANCE_SYM)' " *
			"attribute in XML file: $xml_file (the file predates the SCE " *
			"public-API schema)",
		))
	end
	return parse(Float64, symmetry_node[ATTR_TOLERANCE_SYM])
end

"""
	_read_isotropy(xml_file::AbstractString) -> Bool

Read the `isotropy` attribute from the `<SCEBasis>` node.
"""
function _read_isotropy(xml_file::AbstractString)::Bool
	doc = readxml(xml_file)
	salc_basis_node = findfirst("//" * TAG_SCE_BASIS, doc)
	if isnothing(salc_basis_node)
		throw(ArgumentError("<$(TAG_SCE_BASIS)> node not found in XML file: $xml_file"))
	end
	if !haskey(salc_basis_node, ATTR_ISOTROPY)
		throw(ArgumentError(
			"<$(TAG_SCE_BASIS)> node is missing the '$(ATTR_ISOTROPY)' " *
			"attribute in XML file: $xml_file (the file predates the SCE " *
			"public-API schema)",
		))
	end
	return parse(Bool, salc_basis_node[ATTR_ISOTROPY])
end

"""
	read_basis_components_from_xml(xml_file) -> (structure, symmetry, salcbasis, isotropy)

Read the components of an `SCEBasis` from XML. `structure` is parsed
directly; `symmetry` is *recomputed* via `Symmetry(structure,
tolerance_sym)` (it is a deterministic function of structure + tolerance,
not serialized in full); `salcbasis` is reconstructed by
`read_salcbasis_from_xml`. Accepts both `SCEBasis` and `SCEModel` XML
(the `<JPhi>` block, if present, is ignored). `XMLIO` predates the
`SCEBasis` type, so the tuple is assembled into an `SCEBasis` by
`Magesty.load`.
"""
function read_basis_components_from_xml(xml_file::AbstractString)
	structure = Structure(xml_file; verbosity = false)
	tolerance_sym = _read_tolerance_sym(xml_file)
	symmetry = Symmetry(structure, tolerance_sym; verbosity = false)
	salcbasis = read_salcbasis_from_xml(xml_file)
	isotropy = _read_isotropy(xml_file)
	return structure, symmetry, salcbasis, isotropy
end

"""
	read_model_components_from_xml(xml_file) -> (structure, symmetry, salcbasis, isotropy, j0, jphi)

Read the components of an `SCEModel` from XML: the `SCEBasis` components
(see `read_basis_components_from_xml`) plus the fitted reference energy
`j0` and SCE coefficients `jphi` from the `<JPhi>` block.

# Throws
- `ArgumentError` if the `<JPhi>` block is absent (the file is a
  basis-only `SCEBasis` XML, not an `SCEModel` XML).
"""
function read_model_components_from_xml(xml_file::AbstractString)
	structure, symmetry, salcbasis, isotropy = read_basis_components_from_xml(xml_file)

	doc = readxml(xml_file)
	jphi_block = findfirst("//" * TAG_JPHI_BLOCK, doc)
	if isnothing(jphi_block)
		throw(ArgumentError(
			"<$(TAG_JPHI_BLOCK)> block not found in XML file: $xml_file " *
			"(this is a basis-only SCEBasis XML, not an SCEModel XML)",
		))
	end

	ref_energy_node = findfirst(TAG_REF_ENERGY, jphi_block)
	if isnothing(ref_energy_node)
		throw(ArgumentError(
			"<$(TAG_REF_ENERGY)> node not found in <$(TAG_JPHI_BLOCK)> block: $xml_file",
		))
	end
	j0 = parse(Float64, nodecontent(ref_energy_node))

	jphi_nodes = findall(TAG_JPHI, jphi_block)
	jphi = Vector{Float64}(undef, length(jphi_nodes))
	for node in jphi_nodes
		idx = parse(Int, node["salc_index"])
		jphi[idx] = parse(Float64, nodecontent(node))
	end

	return structure, symmetry, salcbasis, isotropy, j0, jphi
end

end # module XMLIO
