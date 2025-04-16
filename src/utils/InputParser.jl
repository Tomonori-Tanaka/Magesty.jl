"""
    Module InputParser

This module provides functionality for parsing input data and constructing a Parser object.
The Parser object contains all necessary information for the simulation, including
structure structure, symmetry, interaction parameters, and regression settings.
"""
module InputParser

export Parser

"""
    struct Parser

A structure that holds all input parameters for the simulation.

# Fields
- `name::String`: Name of the structure
- `mode::String`: Operation mode (default: "optimize")
- `num_atoms::Int`: Number of atoms in the structure
- `kd_name::Vector{String}`: List of chemical species names
- `is_periodic::Vector{Bool}`: Periodicity flags for each direction
- `j_zero_thr::Float64`: Threshold for zero interaction (default: 1e-8)
- `tolerance_sym::Float64`: Symmetry tolerance (default: 1e-3)
- `nbody::Int`: Maximum number of bodies in interactions
- `lmax::Matrix{Int}`: Maximum angular momentum values [nkd × nbody]
- `cutoff_radii::Array{Float64, 3}`: Cutoff radii for interactions [nkd × nkd × nbody]
- `weight::Float64`: Weight parameter for regression
- `datafile::String`: Path to data file
- `scale::Float64`: Scale factor for lattice vectors
- `lattice_vectors::Matrix{Float64}`: Lattice vectors (3×3 matrix)
- `kd_int_list::Vector{Int}`: List of chemical species indices
- `x_fractional::Matrix{Float64}`: Fractional coordinates of atoms
"""
mutable struct Parser
	# general parameters
	name::String
	mode::String
	num_atoms::Int
	kd_name::Vector{String}
	is_periodic::Vector{Bool}
	j_zero_thr::Float64

	# symmetry parameters
	tolerance_sym::Float64

	# interaction parameters
	nbody::Int
	lmax::Matrix{Int}# [≤kd, ≤nbody]
	cutoff_radii::Array{Float64, 3}

	# regression parameters
	weight::Float64
	datafile::String

	# structure parameters
	scale::Float64
	lattice_vectors::Matrix{Float64}
	kd_int_list::Vector{Int}
	x_fractional::Matrix{Float64}
end

"""
    Parser(input_dict::AbstractDict{<:AbstractString, Any})

Construct a Parser object from an input dictionary.

# Arguments
- `input_dict::AbstractDict{<:AbstractString, Any}`: Dictionary containing input parameters

# Returns
- `Parser`: A new Parser instance

# Throws
- `ErrorException` if required fields are missing or invalid
"""
function Parser(input_dict::AbstractDict{<:AbstractString, Any})

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
	if haskey(general_dict, "mode")
		mode::String = general_dict["mode"]
	else
		mode = "optimize"
	end
	num_atoms::Int = general_dict["nat"]
	if length(general_dict["kd"]) != length(Set(general_dict["kd"]))
		error("Duplication is detected in \"kd\" tag.")
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
	nbody = interaction_dict["nbody"]
	lmax::Matrix{Int} = parse_lmax(interaction_dict["lmax"], kd_name, nbody)
	cutoff_radii::Array{Float64} = parse_cutoff(interaction_dict["cutoff"], kd_name, nbody)

	# regression
	regression_dict = input_dict["regression"]
	weight::Float64 = regression_dict["weight"]
	datafile::String = regression_dict["datafile"]

	# structure
	structure_dict = input_dict["structure"]
	scale = structure_dict["scale"]
	lattice_vectors::Matrix{Float64} = scale * hcat(structure_dict["lattice"]...)
	if size(lattice_vectors) != (3, 3)
		error(
			"The size of \"lattice\":$(size(lattice_vectors)) is different with the correct size:(3, 3).",
		)
	end
	kd_int_list::Vector{Int} = structure_dict["kd_list"]
	if length(kd_int_list) != num_atoms
		error(
			"The length of \"kd_list\":$(length(kd_int_list)) is diffrent with that of \"nat\":$(num_atoms).",
		)
	end
	kd_int_set = Set(kd_int_list)# temporaly variable
	nkd = length(kd_name)
	if nkd != length(kd_int_set)
		error(
			"The number of elements: $nkd is different with that obtained from \"kd_list\".",
		)
	elseif any(x -> (x < 1 || nkd < x), kd_int_set)
		error("The integers in \"kd_list\" must be consecutive numbers starting from 1.")
	end
	if length(structure_dict["position"]) != num_atoms
		error(
			"The length of \"position\":$(length(structure_dict["position"])) is different with that of \"nat\":$(num_atoms).",
		)
	end

	x_fractional::Matrix{Float64} = parse_position(structure_dict["position"], num_atoms)

	return Parser(
		name,
		mode,
		num_atoms,
		kd_name,
		is_periodic,
		j_zero_thr,
		tolerance_sym,
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

"""
    check_interaction_field(interaction_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString})

Validate the interaction field parameters.

# Arguments
- `interaction_dict::AbstractDict`: Dictionary containing interaction parameters
- `kd_name::AbstractVector`: List of chemical species names

# Throws
- `ErrorException` if parameters are invalid
"""
function check_interaction_field(
	interaction_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
)

	nbody = interaction_dict["nbody"]
	lmax_dict::Dict{String, Any} = interaction_dict["lmax"]
	cutoff_dict::Dict{String, Any} = interaction_dict["cutoff"]

	# length check
	for (key, value) in lmax_dict
		if length(value) != nbody
			error("The size of $key in \"lmax\" tag is different from \"nbody\" tag.")
		end
		if !(key in kd_name)
			error("Specified element $key in \"lmax\" tag is not in \"kd\"")
		end
	end
	for (key, value) in cutoff_dict
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

"""
    parse_lmax(lmax_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString}, nbody::Integer)

Parse and validate maximum angular momentum values.

# Arguments
- `lmax_dict::AbstractDict`: Dictionary of lmax values
- `kd_name::AbstractVector`: List of chemical species names
- `nbody::Integer`: Maximum number of bodies

# Returns
- `Matrix{Int}`: Matrix of lmax values [nkd × nbody]

# Throws
- `ErrorException` if parameters are invalid
"""
function parse_lmax(
	lmax_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
	nbody::Integer,
)
	kd_num = length(kd_name)
	lmax_tmp = fill(-1, kd_num, nbody)
	lmax_check = fill(false, kd_num, nbody)

	for (kd_index, kd) in enumerate(kd_name)
		# Check if the length of lmax values matches nbody
		if length(lmax_dict[kd]) != nbody
			error("""
			Invalid lmax values for element $kd:
			- Expected length: $nbody (matching nbody parameter)
			- Actual length: $(length(lmax_dict[kd]))
			Please ensure the number of lmax values matches the nbody parameter.
			""")
		end
		for j in 1:nbody
			lmax_tmp[kd_index, j] = lmax_dict[kd][j]
			lmax_check[kd_index, j] = true
		end
	end

	# Check for missing values
	missing_indices = findall(x -> x == false, lmax_check)
	if !isempty(missing_indices)
		error("""
		Missing lmax values in the input:
		$(join([
			"Element $(kd_name[idx[1]]), l = $(idx[2])" 
			for idx in missing_indices
		], "\n"))
		Please ensure all lmax values are properly specified.
		""")
	end

	return Matrix{Int}(lmax_tmp)
end

"""
    parse_cutoff(cutoff_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString}, nbody::Integer)

Parse and validate cutoff radii.

# Arguments
- `cutoff_dict::AbstractDict`: Dictionary of cutoff values
- `kd_name::AbstractVector`: List of chemical species names
- `nbody::Integer`: Maximum number of bodies

# Returns
- `Array{Float64, 3}`: Array of cutoff radii [nkd × nkd × nbody]

# Throws
- `ErrorException` if parameters are invalid
"""
function parse_cutoff(
	cutoff_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
	nbody::Integer,
)
	kd_num = length(kd_name)
	cutoff_tmp = fill(0.0, kd_num, kd_num, nbody)
	cutoff_check = fill(false, kd_num, kd_num, nbody)

	for n in 1:nbody
		if n == 1
			cutoff_tmp[:, :, n] .= 0.0
			cutoff_check[:, :, n] .= true
			continue
		end
		for (key, value) in cutoff_dict
			# Check if the length of cutoff values matches nbody
			if length(value) != nbody
				error("""
				Invalid cutoff values for pair $key:
				- Expected length: $nbody (matching nbody parameter)
				- Actual length: $(length(value))
				Please ensure the number of cutoff values matches the nbody parameter.
				""")
			end
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

	# Check for missing values
	missing_indices = findall(x -> x == false, cutoff_check)
	if !isempty(missing_indices)
		error("""
		Missing cutoff values in the input:
		$(join([
			"Element1: $(kd_name[idx[1]]), Element2: $(kd_name[idx[2]]), nbody: $(idx[3])" 
			for idx in missing_indices
		], "\n"))
		Please ensure all cutoff values are properly specified.
		""")
	end

	return Array{Float64}(cutoff_tmp)
end

"""
    parse_position(position_list::AbstractVector{<:AbstractVector{<:Real}}, num_atoms::Integer)

Parse and validate atomic positions.

# Arguments
- `position_list::AbstractVector`: List of position vectors
- `num_atoms::Integer`: Number of atoms

# Returns
- `Matrix{Float64}`: Matrix of atomic positions [3 × num_atoms]

# Throws
- `ErrorException` if positions are invalid
"""
function parse_position(
	position_list::AbstractVector{<:AbstractVector{<:Real}},
	num_atoms::Integer,
)::Matrix{Float64}
	position_tmp = fill(-1.0, 3, num_atoms)
	position_check = fill(false, 3, num_atoms)

	for (i, vec) in enumerate(position_list)
		position_tmp[:, i] = vec
		position_check[:, i] .= true
	end

	# check missing
	if false in position_check
		error("Error in \"coords\" list in position field.")
	end

	position = Matrix{Float64}(position_tmp)

	return position
end

end # module
