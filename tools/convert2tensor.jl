using ArgParse
using EzXML

function convert2tensor(input::AbstractString, atoms::Vector{Int}, cell::Integer)::Matrix{Float64}
	atom1 = atoms[1]
	atom2 = atoms[2]


	# tensor matrix composing jxx, jxy, jxz, jyx, jyy, jyz, jzx, jzy, jzz
	result = zeros(3, 3)

	# Read the XML file
	doc = readxml(input)

	# tensor matrix composing j-1-1, j-10, j-11, j0-1, j00, j01, j1-1, j10, j11
	result_tmp = zeros(3, 3)
	is_found = false	# Flag to check the target interaction is found
	for alpha in -1:1, beta in -1:1
		# Extract the <SCEBasisSet> node
		sce_basis_set = findfirst("//SCEBasisSet", doc)
		if isnothing(sce_basis_set)
			throw(ArgumentError("<SCEBasisSet> node not found in the XML file."))
		end

		for salc in EzXML.findall("SALC", sce_basis_set)
			index = parse(Int, salc["index"])

			# Get the SCE coefficient for the current SALC
			JPhi_node = findfirst("//JPhi", doc)
			j_phi = -1000
			for jphi in EzXML.findall("jphi", JPhi_node)
				if jphi["salc_index"] == string(index)
					j_phi = parse(Float64, nodecontent(jphi))
					break
				end
			end


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

				m_vec = fill(-10, 2)
				for i in 1:2
					name = "index-$i"
					# "1 1 -1 1" → ["1","1","-1","1"] → [1,1,-1,1]
					idx = parse.(Int, split(basis[name]))
					if idx[2] != 1  # Check if the angular momentum is 1
						break
					elseif i == 1 && idx[1] == atom1
						m_vec[1] = idx[3]  # Collect the angular momentum l value
					elseif i == 2 && idx[1] == atom2 && idx[4] == cell
						is_found = true
						m_vec[2] = idx[3]  # Collect the angular momentum l value
					end
				end
				if any(x -> abs(x) > 1, m_vec)
					continue
				end
				if m_vec[1] == alpha && m_vec[2] == beta
					result_tmp[alpha+2, beta+2] += j_phi * (parse(Float64, nodecontent(basis))) * (3/4π)
				end
			end
		end
	end

	if !is_found
		throw(ArgumentError("The target interaction between atoms $atom1 and $atom2 in cell $cell is not found."))
	end


	# Convert the result_tmp to the tensor matrix
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

		"--cell", "-c"
		help = "Cell Index of the second atom (the cell of first atom is always 1). "
		required = true
		arg_type = Int
	end

	args = parse_args(s)
	tensor_matrix::Matrix{Float64} = convert2tensor(args["input"], args["atoms"], args["cell"])
	display(tensor_matrix)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
