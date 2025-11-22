using ArgParse
using EzXML
using LinearAlgebra
using Printf
using StaticArrays

# Cache for XML documents to avoid repeated file reads
const XML_CACHE = Dict{String, EzXML.Document}()

function get_cached_xml(input::AbstractString)::EzXML.Document
	if !haskey(XML_CACHE, input)
		XML_CACHE[input] = readxml(input)
	end
	return XML_CACHE[input]
end

function convert2tensor(input::AbstractString, atoms::Vector{Int})::Matrix{Float64}
	atom1 = atoms[1]
	atom2 = atoms[2]

	# tensor matrix composing jxx, jxy, jxz, jyx, jyy, jyz, jzx, jzy, jzz
	result = zeros(3, 3)

	# Use cached XML document
	doc = get_cached_xml(input)

	# Calculate tensor for atom1 -> atom2
	result_forward = calculate_tensor_for_pair(doc, atom1, atom2)

	# Calculate tensor for atom2 -> atom1 (swapped)
	result_backward = calculate_tensor_for_pair(doc, atom2, atom1)
	antisymmetric_part = 0.5*(result_backward - result_backward')
	symmetric_part = 0.5*(result_backward + result_backward')
	antisymmetric_part = antisymmetric_part'
	result_backward = antisymmetric_part + symmetric_part

	# Add both tensors together
	result = result_forward + result_backward

	return result
end

function calculate_tensor_for_pair(doc, atom1::Int, atom2::Int)::Matrix{Float64}
	# tensor matrix composing j-1-1, j-10, j-11, j0-1, j00, j01, j1-1, j10, j11
	result_tmp = zeros(3, 3)
	is_found = false# Flag to check the target interaction is found

	# Pre-fetch XML nodes to avoid repeated searches
	sce_basis_set = findfirst("//SCEBasisSet", doc)
	if isnothing(sce_basis_set)
		throw(ArgumentError("<SCEBasisSet> node not found in the XML file."))
	end

	JPhi_node = findfirst("//JPhi", doc)
	if isnothing(JPhi_node)
		throw(ArgumentError("<JPhi> node not found in the XML file."))
	end

	# Pre-parse JPhi values for faster lookup
	jphi_dict = Dict{String, Float64}()
	for jphi in EzXML.findall("jphi", JPhi_node)
		jphi_dict[jphi["salc_index"]] = parse(Float64, nodecontent(jphi))
	end

	for alpha in -1:1, beta in -1:1
		for salc in EzXML.findall("SALC", sce_basis_set)
			index = parse(Int, salc["index"])
			index_str = string(index)

			# Use pre-parsed JPhi values
			j_phi = get(jphi_dict, index_str, -1000.0)

			num_basis = parse(Int, salc["num_basis"])
			# println("SALC index: $index, num_basis: $num_basis")
			for basis in EzXML.findall("basis", salc)
				# Check if the basis consists of the two atoms
				num_index_attrs = count(attr ->
						startswith(EzXML.nodename(attr), "index-"),
					EzXML.eachattribute(basis),
				)
				if num_index_attrs != 2
					break
				end

				m_vec = MVector{2, Int}(-10, -10)	# Initialize with invalid m values
				for i in 1:2
					name = "index-$i"
					# "1 1 -1 1" → ["1","1","-1","1"] → [1,1,-1,1]
					idx = parse.(Int, split(basis[name]))
					if idx[2] != 1  # Check if the angular momentum is 1
						break
					elseif i == 1 && idx[1] == atom1
						m_vec[1] = idx[3]  # Collect the m value
					elseif i == 2 && idx[1] == atom2
						is_found = true
						m_vec[2] = idx[3]  # Collect the m value
					end
				end
				if any(x -> abs(x) > 1, m_vec)	# Skip if the m value is invalid
					continue
				end
				if m_vec[1] == alpha && m_vec[2] == beta
					result_tmp[alpha+2, beta+2] +=
						j_phi * (parse(Float64, nodecontent(basis))) * (3/(4π))
				end
			end
		end
	end

	# if !is_found
	# 	throw(ArgumentError("The target interaction between atoms $atom1 and $atom2 in cell $cell is not found."))
	# end

	# Convert the result_tmp to the tensor matrix
	result = zeros(3, 3)
	result[1, 1] = result_tmp[3, 3]  # jxx
	result[1, 2] = result_tmp[3, 1]  # jxy
	result[1, 3] = result_tmp[3, 2]  # jxz
	result[2, 1] = result_tmp[1, 3]  # jyx
	result[2, 2] = result_tmp[1, 1]  # jyy
	result[2, 3] = result_tmp[1, 2]  # jyz
	result[3, 1] = result_tmp[2, 3]  # jzx
	result[3, 2] = result_tmp[2, 1]  # jzy
	result[3, 3] = result_tmp[2, 2]  # jzz

	result = result .* 1000

	return result
end

function calculate_biquadratic_term(
	input::AbstractString,
	atoms::Vector{Int},
)::Float64
	atom1 = atoms[1]
	atom2 = atoms[2]

	doc = get_cached_xml(input)

	JPhi_node = findfirst("//JPhi", doc)
	if isnothing(JPhi_node)
		throw(ArgumentError("<JPhi> node not found in the XML file."))
	end

	result_forward = calculate_biquadratic_term_for_pair(doc, atom1, atom2)
	result_backward = calculate_biquadratic_term_for_pair(doc, atom2, atom1)

	return result_forward + result_backward
end

function calculate_biquadratic_term_for_pair(doc, atom1::Int, atom2::Int)::Float64
	sce_basis_set = findfirst("//SCEBasisSet", doc)
	if isnothing(sce_basis_set)
		throw(ArgumentError("<SCEBasisSet> node not found in the XML file."))
	end

	JPhi_node = findfirst("//JPhi", doc)
	if isnothing(JPhi_node)
		throw(ArgumentError("<JPhi> node not found in the XML file."))
	end

	jphi_dict = Dict{String, Float64}()
	for jphi in EzXML.findall("jphi", JPhi_node)
		jphi_dict[jphi["salc_index"]] = parse(Float64, nodecontent(jphi))
	end

	result = 0.0

	for salc in EzXML.findall("SALC", sce_basis_set)
		index = parse(Int, salc["index"])
		index_str = string(index)

		j_phi = get(jphi_dict, index_str, Inf)

		num_basis = parse(Int, salc["num_basis"])
		for basis in EzXML.findall("basis", salc)
			num_index_attrs = count(attr ->
					startswith(EzXML.nodename(attr), "index-"),
				EzXML.eachattribute(basis),
			)
			if num_index_attrs != 2
				break
			end

			m_vec = MVector{2, Int}(-10, -10)
			for i in 1:2
				name = "index-$i"
				idx = parse.(Int, split(basis[name]))
				if idx[2] != 2 # Check if the angular momentum is 2
					break
				elseif i == 1 && idx[1] == atom1
					m_vec[1] = idx[3]  # Collect the angular momentum m value
				elseif i == 2 && idx[1] == atom2
					m_vec[2] = idx[3]  # Collect the angular momentum m value
				end
			end
			if any(x -> abs(x) > 2, m_vec)
				continue
			end
			if m_vec[1] == m_vec[2]
				result += j_phi * (parse(Float64, nodecontent(basis))) * (15/(8π))
			end
		end
	end

	return result
end

function main()
	s = ArgParseSettings()
	@add_arg_table s begin
		"input"
		help = "Input XML file"
		required = true
		arg_type = String

		"--atoms", "-a"
		help = "Two atom indices"
		required = true
		nargs = 2
		arg_type = Int


		"--biquadratic", "-b"
		help = "Calculate biquadratic term"
		action = :store_true
	end

	args = parse_args(s)
	if args["biquadratic"]
		biquadratic_term = calculate_biquadratic_term(args["input"], args["atoms"])
		println("Biquadratic term (meV):")
		println(biquadratic_term)
	else
		tensor_matrix::Matrix{Float64} = convert2tensor(args["input"], args["atoms"])
		println("Tensor matrix (meV):")
		display(tensor_matrix)
		println("")
		println("--------------------------------")
		println(@sprintf("Isotropic Jij (meV): %.6f", 1/3 * tr(tensor_matrix)))
		println("--------------------------------")
		println("Anisotropic symmetric part (meV): ")
		display(1/2 * (tensor_matrix + tensor_matrix') - 1/3 * tr(tensor_matrix) * I)
		println("--------------------------------")
		println("Anisotropic antisymmetric part (meV): ")
		display(1/2 * (tensor_matrix - tensor_matrix'))
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
