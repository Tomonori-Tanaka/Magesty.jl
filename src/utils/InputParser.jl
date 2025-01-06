module InputParser

export Parser

mutable struct Parser
	# general
	name::String
	num_atoms::Int
	kd_name::Vector{String}
	is_periodic::Vector{Bool}
	j_zero_thr::Float64

	# symmetry
	tolerance_sym::Float64

	# interaction
	model::Int
	nbody::Int
	lmax::Matrix{Int}# [≤kd, ≤nbody]
	cutoff_radii::Array{Float64, 3}

	# regression
	weight::Float64
	datafile::String

	# structure
	scale::Float64
	lattice_vectors::Matrix{Float64}
	kd_int_list::Vector{Int}
	x_fractional::Matrix{Float64}
end

function Parser(input_dict::Dict{String, Any})

	if !haskey(input_dict, "general")
		error("\"general\" field is not found in the input.")
	elseif !haskey(input_dict, "symmetry")
		error("\"symmetry\" field is not found in the input.")
	elseif !haskey(input_dict, "interaction")
		error("\"interaction\" field is not found in the input.")
	elseif !haskey(input_dict, "regression")
		error("\"regression\" field is not found in the input.")
	elseif !haskey(input_dict, "structure")
		error("\"structure\" field is not found in the input.")
	end

	# general variables
	general_dict = input_dict["general"]
	name::String = general_dict["name"]
	num_atoms::Int = general_dict["nat"]
	if length(general_dict["kd"]) != length(Set(general_dict["kd"]))
		println("Duplication is detected in \"kd\" tag.")
		exit()
	end
	kd_name::Vector{String} = general_dict["kd"]
	is_periodic::Vector{Bool} = general_dict["periodicity"]
	if haskey(general_dict, "j_zero_thr")
		j_zero_thr::Float64 = general_dict["j_zero_thr"]
	else
		j_zero_thr = 1e-8
	end

	# symmetry variables
	symmetry_dict = input_dict["symmetry"]
	if haskey(symmetry_dict, "tolerance")
		tolerance_sym::Float64 = symmetry_dict["tolerance"]
	else
		tolerance_sym = 1e-3
	end

	# interaction variables
	interaction_dict = input_dict["interaction"]
	check_interaction_field(interaction_dict, kd_name)
	model::Int = interaction_dict["model"]
	if model == 1
		nbody::Int = 2
	else
		nbody = interaction_dict["nbody"]
	end
	lmax::Matrix{Int} = parse_lmax(interaction_dict["lmax"], kd_name, nbody, model)
	cutoff_radii::Array{Float64} = parse_cutoff(interaction_dict["cutoff"], kd_name, nbody)

	# regression
	regression_dict = input_dict["regression"]
	weight::Float64 = regression_dict["weight"]
	datafile::String = regression_dict["datafile"]

	# structure
	structure_dict = input_dict["structure"]
	scale = structure_dict["scale"]
	lattice_vectors = scale * hcat(structure_dict["lattice"]...)
	kd_int_list, x_fractional = parse_position(structure_dict["position"], num_atoms)

	return Parser(
		name,
		num_atoms,
		kd_name,
		is_periodic,
		j_zero_thr,
		tolerance_sym,
		model,
		nbody,
		lmax,
		cutoff_radii,
		weight,
		datafile,
		scale,
		lattice_vectors,
		kd_int_list,
		x_fractional,
	)
end

function check_interaction_field(
	interaction_dict::Dict{String, Any},
	kd_name::Vector{String},
)

	if interaction_dict["model"] == 1
		nbody = 2
	else
		nbody = interaction_dict["nbody"]
	end
	lmax_dict::Dict{String, Any} = interaction_dict["lmax"]
	cutoff_dict::Dict{String, Any} = interaction_dict["cutoff"]

	# length check
	for (key, value) ∈ lmax_dict
		if length(value) != nbody
			error("The size of $key in \"lmax\" tag is different from \"nbody\" tag.")
		end
		if !(key in kd_name)
			error("Specified element $key in \"lmax\" tag is not in \"kd\"")
		end
	end
	for (key, value) ∈ cutoff_dict
		if length(value) != nbody
			error("The size of $key in \"cutoff\" tag is different from \"nbody\" tag.")
		end
		key1, key2 = split(key, "-")
		if !(key1 in kd_name) && key1 != "*"
			error("Specified element $key1 in \"cutoff\" tag is not in \"kd\"")
		elseif !(key2 in kd_name) && key2 != "*"
			error("Specified element $key2 in \"cutoff\" tag is not in \"kd\"")
		end
	end

	return nothing
end

function parse_lmax(
	lmax_dict::Dict{String, Any},
	kd_name::Vector{String},
	nbody::Int,
	model::Int,
)
	kd_num = length(kd_name)
	lmax_tmp = fill(-1, kd_num, nbody)
	lmax_check = fill(false, kd_num, nbody)

	for (kd_index, kd) ∈ enumerate(kd_name)
		for j ∈ 1:nbody
			lmax_tmp[kd_index, j] = lmax_dict[kd][j]
			lmax_check[kd_index, j] = true
		end
	end

	# check 
	indices = findall(x -> x == false, lmax_check)
	if length(indices) != 0
		println("ERROR found in \"lmax\" tag")
		for index ∈ indices
			println(
				"element: ",
				kd_name[index[1]],
				", l = ",
				index[2],
				" is not specified.",
			)
		end
		error()
	end

	if model == 1
		indices = findall(x -> !(x in [0, 1]), lmax_tmp[:, 1])
		if !isempty(indices)
			error("lmax on the 1-body interaction should be 0 or 1 when model = 1.")
		end
		indices = findall(x -> !(x == 1), lmax_tmp[:, 2])
		if !isempty(indices)
			error("lmax on the 2-body interaction should be 1 when model = 1.")
		end
	end

	return Matrix{Int}(lmax_tmp)
end

function parse_cutoff(cutoff_dict::Dict{String, Any}, kd_name::Vector{String}, nbody::Int)
	kd_num = length(kd_name)
	cutoff_tmp = fill(0.0, kd_num, kd_num, nbody)
	cutoff_check = fill(false, kd_num, kd_num, nbody)

	for n ∈ 1:nbody
		if n == 1
			cutoff_tmp[:, :, n] .= 0.0
			cutoff_check[:, :, n] .= true
			continue
		end
		for (key, value) ∈ cutoff_dict
			elem1 = strip(split(key, "-")[1])
			elem2 = strip(split(key, "-")[2])

			index_elem1 = findall(x -> x == elem1, kd_name)[1]
			index_elem2 = findall(x -> x == elem2, kd_name)[1]
			cutoff_tmp[index_elem1, index_elem2, n] = value[n]
			cutoff_tmp[index_elem2, index_elem1, n] = value[n]
			cutoff_check[index_elem1, index_elem2, n] = true
			cutoff_check[index_elem2, index_elem1, n] = true

		end
	end

	indices = findall(x -> x == false, cutoff_check)
	if length(indices) != 0
		println("ERROR found in \"cutoff\" tag")
		for index ∈ indices
			println(
				"element1: ",
				kd_name[index[1]],
				", element2: ",
				kd_name[index[2]],
				", nbody: ",
				index[3],
				" is not specified.",
			)
		end
		exit()
	end
	# println(cutoff_tmp)
	return Array{Float64}(cutoff_tmp)
end

function parse_position(
	position_list::Vector{Dict{String, Any}},
	num_atoms::Int,
)
	kd_int_list_tmp = fill(-1, num_atoms)
	kd_int_list_check = fill(false, num_atoms)
	position_tmp = fill(-1.0, 3, num_atoms)
	position_check = fill(false, 3, num_atoms)

	for (i, dict) ∈ enumerate(position_list)
		kd_int_list_tmp[i] = dict["kd"]
		kd_int_list_check[i] = true
		position_tmp[:, i] = dict["coords"]
		position_check[:, i] .= true
	end

	# check missing
	if false in kd_int_list_check
		error("Error in \"kd\" list in position field.")
	end
	if false in position_check
		error("Error in \"coords\" list in position field.")
	end

	kd_int_list = Vector{Int}(kd_int_list_tmp)
	position = Matrix{Float64}(position_tmp)

	return kd_int_list, position
end

end # module
