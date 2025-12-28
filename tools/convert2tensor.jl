module ExchangeTensor
using EzXML
using LinearAlgebra
using Printf

if !@isdefined(Magesty)
	include("../src/Magesty.jl")
end
using .Magesty

export ExchangeTensorData, print_full, convert2tensor


"""
	ExchangeTensorData

Exchange interaction tensor for a pair of atoms. Unit is meV.

# Fields
- `jij_tensor::Matrix{Float64}`: 3×3 exchange interaction tensor in Cartesian coordinates (meV)
- `isotropic_jij::Float64`: Isotropic exchange parameter (trace/3 of jij_tensor) (meV)
- `dm_vector::Vector{Float64}`: Dzyaloshinskii-Moriya interaction vector (3 elements) (meV)
- `gamma::Matrix{Float64}`: Anisotropic symmetric exchange tensor (3×3) (meV)
"""
struct ExchangeTensorData<:AbstractMatrix{Float64}
	jij_tensor::Matrix{Float64}
	isotropic_jij::Float64
	dm_vector::Vector{Float64}
	gamma::Matrix{Float64}

	function ExchangeTensorData(jij_tensor::Matrix{Float64})
		iso_jij = tr(jij_tensor) / 3
		# DMI vector from antisymmetric part: W = (1/2)(J - J^T)
		# D_x = W[2,3] = (1/2)(J[2,3] - J[3,2])
		# D_y = W[3,1] = (1/2)(J[3,1] - J[1,3])
		# D_z = W[1,2] = (1/2)(J[1,2] - J[2,1])
		# This matches calc_dm_vector in micromagnetics.jl
		dm_vec =
			[jij_tensor[2, 3] - jij_tensor[3, 2],
				jij_tensor[3, 1] - jij_tensor[1, 3],
				jij_tensor[1, 2] - jij_tensor[2, 1]] / 2
		gamma_mat = (jij_tensor + jij_tensor') / 2 - iso_jij * I(3)
		return new(jij_tensor, iso_jij, dm_vec, gamma_mat)
	end
end

function Base.size(tensor::ExchangeTensorData)
	return size(tensor.jij_tensor)
end

function Base.getindex(tensor::ExchangeTensorData, i::Int, j::Int)
	return tensor.jij_tensor[i, j]
end

function Base.setindex!(tensor::ExchangeTensorData, value::Float64, i::Int, j::Int)
	tensor.jij_tensor[i, j] = value
end

function Base.adjoint(tensor::ExchangeTensorData)
	return ExchangeTensorData(Matrix(tensor.jij_tensor'))
end

function Base.transpose(tensor::ExchangeTensorData)
	return ExchangeTensorData(Matrix(transpose(tensor.jij_tensor)))
end

function Base.:*(tensor::ExchangeTensorData, scalar::Number)
	return ExchangeTensorData(tensor.jij_tensor .* scalar)
end

function Base.:*(scalar::Number, tensor::ExchangeTensorData)
	return tensor * scalar
end

function Base.show(io::IO, tensor::ExchangeTensorData)
	show(io, tensor.jij_tensor)
end

"""
	print_full(io::IO, tensor::ExchangeTensorData)
	print_full(tensor::ExchangeTensorData)

Print all information of ExchangeTensorData in a formatted way.
"""
function print_full(io::IO, tensor::ExchangeTensorData)
	println(io, "Exchange Tensor (meV):")
	println(io, "=" ^ 50)
	println(io, "\nJ_ij tensor (3×3):")
	for i in 1:3
		println(
			io,
			@sprintf("  % 12.6f  % 12.6f  % 12.6f",
				tensor.jij_tensor[i, 1],
				tensor.jij_tensor[i, 2],
				tensor.jij_tensor[i, 3])
		)
	end
	println(io, "\nIsotropic exchange J_ij: ", @sprintf("% 12.6f meV", tensor.isotropic_jij))
	println(io, "\nDzyaloshinskii-Moriya vector (meV):")
	println(io, @sprintf("  D_x: % 12.6f", tensor.dm_vector[1]))
	println(io, @sprintf("  D_y: % 12.6f", tensor.dm_vector[2]))
	println(io, @sprintf("  D_z: % 12.6f", tensor.dm_vector[3]))
	println(io, "\nAnisotropic symmetric exchange Γ (3×3):")
	for i in 1:3
		println(
			io,
			@sprintf("  % 12.6f  % 12.6f  % 12.6f",
				tensor.gamma[i, 1],
				tensor.gamma[i, 2],
				tensor.gamma[i, 3])
		)
	end
	println(io, "=" ^ 50)
end

function print_full(tensor::ExchangeTensorData)
	print_full(stdout, tensor)
end

function convert2tensor(input::AbstractString, atoms::AbstractVector{Int})::ExchangeTensorData
	atom1 = atoms[1]
	atom2 = atoms[2]

	structure = Magesty.Structure(input, verbosity = false)
	symmetry = Magesty.Symmetry(structure, 1e-5, verbosity = false)

	# Read XML document
	doc = readxml(input)

	# Strategy: First try mapping atom1 to primitive cell, then atom2 if not found
	# This ensures 1-2 and 16-1 give the same result

	# Try 1: Map atom1 to primitive cell, apply translation to atom2
	atom1_in_prim = symmetry.atoms_in_prim[symmetry.map_s2p[atom1].atom]
	translation_global = symmetry.symnum_translation[symmetry.map_s2p[atom1].translation]
	atom2_translated = symmetry.map_sym_inv[atom2, translation_global]

	# Normalize pair order: always use smaller index as first atom
	atom1_prim = atom1_in_prim
	atom2_prim = atom2_translated
	is_permuted = false

	if atom1_prim > atom2_prim
		is_permuted = true
		atom1_prim, atom2_prim = atom2_prim, atom1_prim
	end

	# Check if this pair exists
	pair_found = false
	try
		pair_found = check_pair_exists(doc, atom1_prim, atom2_prim)
	catch
		pair_found = false
	end

	result_debug = calculate_tensor_for_pair_debug(doc, atom1_prim, atom2_prim)
	println("result_debug:")
	display(result_debug)

	if pair_found
		result = calculate_tensor_for_pair(doc, atom1_prim, atom2_prim)
		if is_permuted
			# When order is swapped, need to transpose the result
			# J(j,i) = J(i,j)^T
			result = result'
		end
		# Convert from Hartree to meV (multiply by 1000)
		return result * 1000
	end

	# Try 2: Map atom2 to primitive cell, apply translation to atom1
	atom2_in_prim = symmetry.atoms_in_prim[symmetry.map_s2p[atom2].atom]
	translation_global = symmetry.symnum_translation[symmetry.map_s2p[atom2].translation]
	atom1_translated = symmetry.map_sym_inv[atom1, translation_global]

	# Normalize pair order again
	atom1_prim = atom1_translated
	atom2_prim = atom2_in_prim
	is_permuted = false

	if atom1_prim > atom2_prim
		is_permuted = true
		atom1_prim, atom2_prim = atom2_prim, atom1_prim
	end

	# Check if this pair exists
	pair_found = false
	try
		pair_found = check_pair_exists(doc, atom1_prim, atom2_prim)
	catch
		pair_found = false
	end


	if pair_found
		result = calculate_tensor_for_pair(doc, atom1_prim, atom2_prim)
		if is_permuted
			# When order is swapped, need to transpose the result
			# J(j,i) = J(i,j)^T
			result = result'
		end
		# Convert from Hartree to meV (multiply by 1000)
		return result * 1000
	end
	# If no pair found, return zero tensor
	return ExchangeTensorData(zeros(3, 3)) * 1000
end

function calculate_tensor_for_pair(doc, atom1::Int, atom2::Int)::ExchangeTensorData
	# Exchange interaction tensor J_ij (3x3 matrix)
	# Separate contributions by Lf:
	# - Lf=0: Isotropic exchange Jij (scalar)
	# - Lf=1: Dzyaloshinskii-Moriya interaction (DMI, antisymmetric)
	# - Lf=2: Anisotropic symmetric exchange (symmetric traceless)

	# Pre-fetch XML nodes
	sce_basis_set = findfirst("//SCEBasisSet", doc)
	if isnothing(sce_basis_set)
		throw(ArgumentError("<SCEBasisSet> node not found in the XML file."))
	end

	JPhi_node = findfirst("//JPhi", doc)
	if isnothing(JPhi_node)
		throw(ArgumentError("<JPhi> node not found in the XML file."))
	end

	# Pre-parse JPhi values
	jphi_dict = Dict{String, Float64}()
	for jphi in EzXML.findall("jphi", JPhi_node)
		jphi_dict[jphi["salc_index"]] = parse(Float64, nodecontent(jphi))
	end

	# Reconstruct coeff_tensor for each Lf (atom-pair independent quantities)
	# These are computed once and reused for all atom pairs
	coeff_tensor_Lf0, coeff_tensor_Lf1, coeff_tensor_Lf2 = reconstruct_coeff_tensors(doc)


	# Scaling factor: (4π)^(n_C/2) where n_C=2 for pair
	# This matches the scaling in design matrix (see Optimize.jl:266)
	# Paper Eq. (3): Φ = (√(4π))^{n_C} ∏ Y_{lm}

	# Store SALC information (salc_index, j_phi, coefficient) for each Lf
	# Structure: salc_info_Lf[Lf] = Vector of (salc_index, j_phi, coefficient)
	salc_info_Lf0 = Vector{Tuple{Int, Float64, Vector{Float64}, Int}}()  # Lf=0: Isotropic exchange
	salc_info_Lf1 = Vector{Tuple{Int, Float64, Vector{Float64}, Int}}()  # Lf=1: DMI
	salc_info_Lf2 = Vector{Tuple{Int, Float64, Vector{Float64}, Int}}()  # Lf=2: Anisotropic symmetric

	# First pass: Collect SALC information for matching atom pairs
	for salc_node in EzXML.findall("SALC", sce_basis_set)
		salc_index = parse(Int, salc_node["index"])
		index_str = string(salc_index)
		j_phi = get(jphi_dict, index_str, 0.0)

		# Skip if j_phi is zero or not found
		if j_phi == 0.0 || !haskey(jphi_dict, index_str)
			continue
		end

		Lf = parse(Int, salc_node["Lf"])
		body = parse(Int, salc_node["body"])

		# Only process 2-body interactions with ls=[1,1]
		if body != 2
			continue
		end

		# Only process Lf=0, 1, 2 for pair interactions
		if !(Lf in [0, 1, 2])
			continue
		end

		# Find the basis in this SALC group that matches the target atom pair
		# Note: One SALC group never contains the same atom pair twice
		for basis_node in EzXML.findall("basis", salc_node)
			atoms = parse.(Int, split(basis_node["atoms"]))
			ls = parse.(Int, split(basis_node["ls"]))
			multiplicity = parse(Int, basis_node["multiplicity"])

			# Check if this basis matches the target atom pair
			if ls != [1, 1]
				continue
			end
			if atoms != [atom1, atom2]
				continue
			end

			# Found matching basis - store SALC information
			coeff_str = nodecontent(basis_node)
			coefficient = parse.(Float64, split(coeff_str))

			# Store in the appropriate Lf list
			if Lf == 0
				push!(salc_info_Lf0, (salc_index, j_phi, coefficient, multiplicity))
			elseif Lf == 1
				push!(salc_info_Lf1, (salc_index, j_phi, coefficient, multiplicity))
			elseif Lf == 2
				push!(salc_info_Lf2, (salc_index, j_phi, coefficient, multiplicity))
			end

		end
	end

	coeff_total_l0 = zeros(3, 3)
	coeff_total_l1 = zeros(3, 3)
	coeff_total_l2 = zeros(3, 3)

	if !isnothing(coeff_tensor_Lf0)
		for (salc_index, j_phi, coefficient, multiplicity) in salc_info_Lf0
			for mf_idx in 1:1
				coeff_total_l0 +=
					j_phi * coeff_tensor_Lf0[:, :, mf_idx] .* coefficient[mf_idx] * 3
			end
			# coeff_total_l0 += j_phi * coefficient[1] * I(3) * sqrt(3)
		end
	end
	# Lf=1: Dzyaloshinskii-Moriya interaction (DMI, antisymmetric)
	# DMI is antisymmetric part of exchange tensor: D = (1/2) * (J - J^T)
	# In tesseral coordinates: m = -1 (y), 0 (z), 1 (x)
	# Mf ranges from -1 to 1, so mf_idx = 1, 2, 3 corresponds to Mf = -1, 0, 1
	if !isnothing(coeff_tensor_Lf1)
		for (salc_index, j_phi, coefficient, multiplicity) in salc_info_Lf1
			for mf_idx in 1:3
				# coefficient[mf_idx] is the SALC coefficient for Mf = mf_idx - 2 (i.e., -1, 0, 1)
				# Scaling factor: 3 * multiplicity (from design matrix scaling)
				coeff_total_l1 +=
					j_phi * coeff_tensor_Lf1[:, :, mf_idx] .* coefficient[mf_idx] * 3
			end
		end
	end
	if !isnothing(coeff_tensor_Lf2)
		for (salc_index, j_phi, coefficient, multiplicity) in salc_info_Lf2
			for mf_idx in 1:5
				coeff_total_l2 +=
					j_phi * coeff_tensor_Lf2[:, :, mf_idx] .* coefficient[mf_idx] * 3
			end
		end
	end

	tensor = coeff_total_l0 + coeff_total_l1 + coeff_total_l2
	# Convert from tesseral to Cartesian coordinates
	# Tesseral order: m = -1 (y), 0 (z), 1 (x) → Cartesian order: x, y, z
	# Conversion matrix: [y, z, x] → [x, y, z]
	convert_matrix = [0 0 1; 1 0 0; 0 1 0]  # Maps: [tesseral_y, tesseral_z, tesseral_x] → [cart_x, cart_y, cart_z]
	tensor_cartesian = convert_matrix * tensor * convert_matrix'

	# Second pass: Compute linear combinations for each Lf using stored SALC information
	return ExchangeTensorData(tensor_cartesian)
end

"""
	check_pair_exists(doc, atom1::Int, atom2::Int) -> Bool

Check if the atom pair exists in the XML SALC basis set.
"""
function check_pair_exists(doc, atom1::Int, atom2::Int)::Bool
	sce_basis_set = findfirst("//SCEBasisSet", doc)
	if isnothing(sce_basis_set)
		return false
	end

	# Check if pair exists in any SALC
	for salc_node in EzXML.findall("SALC", sce_basis_set)
		body = parse(Int, salc_node["body"])
		if body != 2
			continue
		end

		for basis_node in EzXML.findall("basis", salc_node)
			atoms = parse.(Int, split(basis_node["atoms"]))
			ls = parse.(Int, split(basis_node["ls"]))

			if ls == [1, 1] && atoms == [atom1, atom2]
				return true
			end
		end
	end

	return false
end

"""
	reconstruct_coeff_tensors(doc) -> Tuple{Union{Array{Float64,3},Nothing}, ...}

Reconstruct coeff_tensor for Lf=0, 1, 2 from XML document.
These are atom-pair independent quantities.
Returns (coeff_tensor_Lf0, coeff_tensor_Lf1, coeff_tensor_Lf2).
"""
function reconstruct_coeff_tensors(
	doc,
)::Tuple{
	Union{Array{Float64, 3}, Nothing},
	Union{Array{Float64, 3}, Nothing},
	Union{Array{Float64, 3}, Nothing},
}
	coeff_tensor_Lf0 = nothing
	coeff_tensor_Lf1 = nothing
	coeff_tensor_Lf2 = nothing

	system_node = findfirst("//System", doc)
	if !isnothing(system_node)
		amc_node = findfirst("AngularMomentumCouplings", system_node)
		if !isnothing(amc_node)
			for coupling_node in findall("Coupling", amc_node)
				ls_str = coupling_node["ls"]
				Lseq_str = coupling_node["Lseq"]
				Lf = parse(Int, coupling_node["Lf"])

				ls = parse.(Int, split(ls_str))
				Lseq = isempty(Lseq_str) ? Int[] : parse.(Int, split(Lseq_str))

				# Only process ls=[1,1] for pair interactions
				if length(ls) != 2 || ls != [1, 1]
					continue
				end

				# Only process Lf=0, 1, 2
				if !(Lf in [0, 1, 2])
					continue
				end

				# Reconstruct coeff_tensor from XML
				coeff_tensor = reconstruct_coeff_tensor_from_node(coupling_node)

				# Store in the appropriate Lf variable
				if Lf == 0
					coeff_tensor_Lf0 = coeff_tensor
				elseif Lf == 1
					coeff_tensor_Lf1 = coeff_tensor
				elseif Lf == 2
					coeff_tensor_Lf2 = coeff_tensor
				end
			end
		end
	end

	return (coeff_tensor_Lf0, coeff_tensor_Lf1, coeff_tensor_Lf2)
end

"""
	reconstruct_coeff_tensor_from_node(coupling_node) -> Array{Float64, 3}

Reconstruct coeff_tensor from XML coupling node.
Returns a tensor of shape (site_dims..., Mf_size).
"""
function reconstruct_coeff_tensor_from_node(coupling_node)::Array{Float64, 3}
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

	return coeff_tensor
end

"""
	calculate_tensor_for_pair_debug(doc, atom1::Int, atom2::Int)::ExchangeTensorData

Calculate the exchange tensor for a pair of atoms using the debug mode.
Directly calculate the tensor brute force.
"""
function calculate_tensor_for_pair_debug(doc, atom1::Int, atom2::Int)::ExchangeTensorData
	# Pre-fetch XML nodes
	sce_basis_set = findfirst("//SCEBasisSet", doc)
	if isnothing(sce_basis_set)
		throw(ArgumentError("<SCEBasisSet> node not found in the XML file."))
	end

	JPhi_node = findfirst("//JPhi", doc)
	if isnothing(JPhi_node)
		throw(ArgumentError("<JPhi> node not found in the XML file."))
	end

	# Pre-parse JPhi values
	jphi_dict = Dict{String, Float64}()
	for jphi in EzXML.findall("jphi", JPhi_node)
		jphi_dict[jphi["salc_index"]] = parse(Float64, nodecontent(jphi))
	end

	coeff_tensor_Lf0, coeff_tensor_Lf1, coeff_tensor_Lf2 = reconstruct_coeff_tensors(doc)

	# Calculate the tensor brute force.
	# alpha and beta means m and m' in the tesseral basis.
	# Convert to Cartesian basis (x, y, z) later.
	tensor = zeros(Float64, 3, 3)
	for salc_node in EzXML.findall("SALC", sce_basis_set)
		salc_index = parse(Int, salc_node["index"])
		index_str = string(salc_index)
		j_phi = get(jphi_dict, index_str, 0.0)

		if j_phi == 0.0 || !haskey(jphi_dict, index_str)
			continue
		end

		# Lf = 0: Isotropic exchange
		Lf = parse(Int, salc_node["Lf"])
		if Lf == 0
			if isnothing(coeff_tensor_Lf0)
				continue
			end
			
			for basis_node in EzXML.findall("basis", salc_node)
				atoms = parse.(Int, split(basis_node["atoms"]))
				ls = parse.(Int, split(basis_node["ls"]))
				multiplicity = parse(Int, basis_node["multiplicity"])

				if ls != [1, 1]
					continue
				end
				# Check if atoms match (considering order: [atom1, atom2] or [atom2, atom1])
				if sort(atoms) != sort([atom1, atom2])
					continue
				end

				coeff_str = nodecontent(basis_node)
				coefficient::Vector{Float64} = parse.(Float64, split(coeff_str))
				# because Lf = 0, there is only one basis function.
				coeff = coefficient[1]
				coupling_tensor = coeff_tensor_Lf0[:, :, 1]

				for alpha in 1:3, beta in 1:3
					tensor[alpha, beta] += j_phi * coupling_tensor[alpha, beta] * coeff * multiplicity 
				end
			end
		end

	end
	tensor = tensor * 1000 * 3

	return ExchangeTensorData(tensor)
end

end # module ExchangeTensor

# Script execution part (outside module)
using ArgParse

function main()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"input"
		help = "Input XML file"
		required = true
		arg_type = String

		"--atoms", "-a"
		help = "Two atom indices"
		required = true
		nargs = 2
		arg_type = Int

	end

	args = parse_args(s)
	tensor_matrix = ExchangeTensor.convert2tensor(args["input"], args["atoms"])
	ExchangeTensor.print_full(tensor_matrix)
	println("")
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end