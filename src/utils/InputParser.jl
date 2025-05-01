"""
    Module InputParser

This module provides functionality for parsing input data and constructing a Parser object.
The Parser object contains all necessary information for the simulation, including
structure, symmetry, interaction parameters, and regression settings.

# Overview
The module handles the parsing of input parameters from a dictionary format into a structured
Parser object. It performs validation of input parameters and ensures all required fields
are present and properly formatted.

# Key Components
- `Parser` struct: Holds all simulation parameters
- Input validation functions
- Parameter parsing utilities
"""
module InputParser

using Base.Threads
using LinearAlgebra
using Printf
using Statistics

export Parser

"""
    struct Parser

A structure that holds all input parameters for the simulation.

# Fields
- `name::String`: Name of the structure
- `mode::String`: Operation mode ("optimize" by default)
- `num_atoms::Int`: Number of atoms in the structure
- `kd_name::Vector{String}`: List of chemical species names
- `is_periodic::Vector{Bool}`: Periodicity flags for each direction
- `j_zero_thr::Float64`: Threshold for zero interaction (1e-8 by default)
- `tolerance_sym::Float64`: Symmetry tolerance (1e-3 by default)
- `nbody::Int`: Maximum number of bodies in interactions
- `lmax::Matrix{Int}`: Maximum angular momentum values [nkd × nbody]
- `cutoff_radii::Array{Float64, 3}`: Cutoff radii for interactions [nkd × nkd × nbody]
- `weight::Union{Float64, String}`: Weight parameter for regression or "auto" for automatic tuning
- `training_ratio::Float64`: Ratio of data to use for training (1.0 by default)
- `datafile::String`: Path to data file
- `scale::Float64`: Scale factor for lattice vectors
- `lattice_vectors::Matrix{Float64}`: Lattice vectors (3×3 matrix)
- `kd_int_list::Vector{Int}`: List of chemical species indices
- `x_fractional::Matrix{Float64}`: Fractional coordinates of atoms

# Examples
```julia
parser = Parser(input_dict)
```
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
	weight::Union{Float64, String}
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
  with the following required sections:
  - "general": General simulation parameters
  - "symmetry": Symmetry-related parameters
  - "interaction": Interaction parameters
  - "regression": Regression parameters
  - "structure": Structure parameters

# Returns
- `Parser`: A new Parser instance with validated parameters

# Throws
- `ErrorException` if:
  - Required fields are missing
  - Parameters are invalid
  - Data format is incorrect

# Examples
```julia
input_dict = Dict(
    "general" => Dict(
        "name" => "Fe",
        "nat" => 2,
        "kd" => ["Fe"],
        "periodicity" => [true, true, true]
    ),
    # ... other sections ...
)
parser = Parser(input_dict)
```
"""
function Parser(input_dict::AbstractDict{<:AbstractString, Any})
	# Check required sections
	required_sections = ["general", "symmetry", "interaction", "regression", "structure"]
	for section in required_sections
		if !haskey(input_dict, section)
			throw(ArgumentError("Required section \"$section\" is missing in the input dictionary."))
		end
	end

	# Parse general parameters
	general_dict = input_dict["general"]
	name = get(general_dict, "name", "")::String
	if isempty(name)
		throw(ArgumentError("Structure name is required in the \"general\" section."))
	end

	mode = get(general_dict, "mode", "optimize")::String
	if !(mode in ["optimize", "predict"])
		throw(ArgumentError("Invalid mode: $mode. Must be either \"optimize\" or \"predict\"."))
	end

	num_atoms = get(general_dict, "nat", 0)::Int
	if num_atoms <= 0
		throw(ArgumentError("Number of atoms (nat) must be positive, got $num_atoms."))
	end

	kd_name = get(general_dict, "kd", String[])::Vector{String}
	if isempty(kd_name)
		throw(ArgumentError("Chemical species list (kd) is required in the \"general\" section."))
	end
	if length(kd_name) != length(Set(kd_name))
		throw(ArgumentError("Duplicate chemical species found in \"kd\" list."))
	end

	is_periodic = get(general_dict, "periodicity", [true, true, true])::Vector{Bool}
	if length(is_periodic) != 3
		throw(ArgumentError("Periodicity must be specified for all three directions."))
	end

	j_zero_thr = get(general_dict, "j_zero_thr", 1e-8)::Float64
	if j_zero_thr < 0
		throw(ArgumentError("j_zero_thr must be non-negative, got $j_zero_thr."))
	end

	# Parse symmetry parameters
	symmetry_dict = input_dict["symmetry"]
	tolerance_sym = get(symmetry_dict, "tolerance", 1e-3)::Float64
	if tolerance_sym <= 0
		throw(ArgumentError("Symmetry tolerance must be positive, got $tolerance_sym."))
	end

	# Parse interaction parameters
	interaction_dict = input_dict["interaction"]
	check_interaction_field(interaction_dict, kd_name)
	nbody = get(interaction_dict, "nbody", 0)::Int
	if nbody <= 0
		throw(ArgumentError("nbody must be positive, got $nbody."))
	end

	lmax = parse_lmax(interaction_dict["lmax"], kd_name, nbody)
	cutoff_radii = parse_cutoff(interaction_dict["cutoff"], kd_name, nbody)

	# Parse regression parameters
	regression_dict = input_dict["regression"]
	weight = get(regression_dict, "weight", 0.0)::Real
	if weight < 0.0
		throw(ArgumentError("Weight must be non-negative, got $weight."))
	elseif weight > 1.0
		weight = ceil(weight)
	end

	cv = get(regression_dict, "cv", 0)::Int
	if cv < 0
		throw(ArgumentError("cv must be non-negative, got $cv."))
	end

	datafile = get(regression_dict, "datafile", "")::String
	if isempty(datafile)
		throw(ArgumentError("Data file path is required in the \"regression\" section."))
	end

	# Parse structure parameters
	structure_dict = input_dict["structure"]
	scale = get(structure_dict, "scale", 1.0)::Real
	if scale <= 0
		throw(ArgumentError("Scale factor must be positive, got $scale."))
	end

	lattice_vectors = scale * hcat(structure_dict["lattice"]...)
	if size(lattice_vectors) != (3, 3)
		throw(ArgumentError("Lattice vectors must be a 3×3 matrix, got $(size(lattice_vectors))."))
	end

	kd_int_list = get(structure_dict, "kd_list", Int[])::Vector{Int}
	if length(kd_int_list) != num_atoms
		throw(ArgumentError("Length of kd_list ($(length(kd_int_list))) must match number of atoms ($num_atoms)."))
	end

	kd_int_set = Set(kd_int_list)
	nkd = length(kd_name)
	if nkd != length(kd_int_set)
		throw(ArgumentError("Number of chemical species ($nkd) doesn't match kd_list entries."))
	elseif any(x -> (x < 1 || nkd < x), kd_int_set)
		throw(ArgumentError("Chemical species indices must be consecutive numbers starting from 1."))
	end

	positions = get(structure_dict, "position", Vector{Float64}[])::Vector{Vector{Float64}}
	if length(positions) != num_atoms
		throw(ArgumentError("Number of positions ($(length(positions))) must match number of atoms ($num_atoms)."))
	end

	x_fractional = parse_position(positions, num_atoms)

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
  - "nbody": Maximum number of bodies
  - "lmax": Dictionary of maximum angular momentum values
  - "cutoff": Dictionary of cutoff radii
- `kd_name::AbstractVector`: List of chemical species names

# Throws
- `ErrorException` if:
  - Length of lmax values doesn't match nbody
  - Length of cutoff values doesn't match nbody
  - Invalid chemical species names are used
  - Duplicate entries are found

# Examples
```julia
interaction_dict = Dict(
    "nbody" => 2,
    "lmax" => Dict("Fe" => [2, 2]),
    "cutoff" => Dict("Fe-Fe" => [0.0, 5.0])
)
check_interaction_field(interaction_dict, ["Fe"])
```
"""
function check_interaction_field(
	interaction_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
)
	# Check required fields
	required_fields = ["nbody", "lmax", "cutoff"]
	for field in required_fields
		if !haskey(interaction_dict, field)
			throw(ArgumentError("Required field \"$field\" is missing in interaction parameters."))
		end
	end

	nbody = interaction_dict["nbody"]
	lmax_dict = interaction_dict["lmax"]
	cutoff_dict = interaction_dict["cutoff"]

	# Validate lmax values
	for (key, value) in lmax_dict
		if !(key in kd_name)
			throw(ArgumentError("Invalid chemical species \"$key\" in lmax dictionary."))
		end
		if length(value) != nbody
			throw(ArgumentError("Length of lmax values for $key ($(length(value))) doesn't match nbody ($nbody)."))
		end
		if any(x -> x < 0, value)
			throw(ArgumentError("Negative lmax values are not allowed for $key."))
		end
	end

	# Validate cutoff values
	for (key, value) in cutoff_dict
		key1, key2 = split(key, "-")
		if !(key1 in kd_name) && key1 != "*"
			throw(ArgumentError("Invalid chemical species \"$key1\" in cutoff dictionary."))
		end
		if !(key2 in kd_name) && key2 != "*"
			throw(ArgumentError("Invalid chemical species \"$key2\" in cutoff dictionary."))
		end
		if length(value) != nbody
			throw(ArgumentError("Length of cutoff values for $key ($(length(value))) doesn't match nbody ($nbody)."))
		end
	end

	return nothing
end

"""
    parse_lmax(lmax_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString}, nbody::Integer)

Parse and validate maximum angular momentum values.

# Arguments
- `lmax_dict::AbstractDict`: Dictionary of lmax values for each chemical species
- `kd_name::AbstractVector`: List of chemical species names
- `nbody::Integer`: Maximum number of bodies

# Returns
- `Matrix{Int}`: Matrix of lmax values [nkd × nbody]

# Throws
- `ErrorException` if:
  - Length of lmax values doesn't match nbody
  - Missing values for any chemical species
  - Invalid lmax values

# Examples
```julia
lmax_dict = Dict("Fe" => [2, 2], "Ni" => [2, 2])
lmax = parse_lmax(lmax_dict, ["Fe", "Ni"], 2)
```
"""
function parse_lmax(
	lmax_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
	nbody::Integer,
)::Matrix{Int}
	kd_num = length(kd_name)
	lmax_tmp = fill(-1, kd_num, nbody)
	lmax_check = fill(false, kd_num, nbody)

	# Parse lmax values
	for (kd_index, kd) in enumerate(kd_name)
		if !haskey(lmax_dict, kd)
			rror("Missing lmax values for chemical species \"$kd\".")
		end

		values = lmax_dict[kd]
		if length(values) != nbody
			error("Length of lmax values for $kd ($(length(values))) doesn't match nbody ($nbody).")
		end

		for j in 1:nbody
			value = values[j]
			if value < 0
				error("Negative lmax value ($value) is not allowed for $kd.")
			end
			lmax_tmp[kd_index, j] = value
			lmax_check[kd_index, j] = true
		end
	end

	# Check for missing values
	missing_indices = findall(x -> x == false, lmax_check)
	if !isempty(missing_indices)
		error("Missing lmax values for: " * join([
			"$(kd_name[idx[1]]) (l = $(idx[2]))"
			for idx in missing_indices
		], ", "))
	end

	return Matrix{Int}(lmax_tmp)
end

"""
    parse_cutoff(cutoff_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString}, nbody::Integer)

Parse and validate cutoff radii for interactions.

# Arguments
- `cutoff_dict::AbstractDict`: Dictionary of cutoff values for each pair of chemical species
- `kd_name::AbstractVector`: List of chemical species names
- `nbody::Integer`: Maximum number of bodies

# Returns
- `Array{Float64, 3}`: Array of cutoff radii [nkd × nkd × nbody]

# Throws
- `ErrorException` if:
  - Length of cutoff values doesn't match nbody
  - Missing values for any pair of chemical species
  - Invalid cutoff values

# Examples
```julia
cutoff_dict = Dict("Fe-Fe" => [0.0, 5.0], "Fe-Ni" => [0.0, 4.5])
cutoff = parse_cutoff(cutoff_dict, ["Fe", "Ni"], 2)
```
"""
function parse_cutoff(
	cutoff_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
	nbody::Integer,
)::Array{Float64, 3}
	kd_num = length(kd_name)
	cutoff_tmp = fill(0.0, kd_num, kd_num, nbody)
	cutoff_check = fill(false, kd_num, kd_num, nbody)

	# Parse cutoff values
	for n in 1:nbody
		if n == 1
			cutoff_tmp[:, :, n] .= 0.0
			cutoff_check[:, :, n] .= true
			continue
		end

		for (key, value) in cutoff_dict
			if length(value) != nbody
				throw(ArgumentError("Length of cutoff values for $key ($(length(value))) doesn't match nbody ($nbody)."))
			end

			elem1, elem2 = split(key, "-")
			index_elem1 = elem1 == "*" ? 1 : findfirst(x -> x == elem1, kd_name)
			index_elem2 = elem2 == "*" ? 1 : findfirst(x -> x == elem2, kd_name)

			if isnothing(index_elem1) || isnothing(index_elem2)
				throw(ArgumentError("Invalid chemical species pair \"$key\" in cutoff dictionary."))
			end

			cutoff_value = value[n]
			cutoff_tmp[index_elem1, index_elem2, n] = cutoff_value
			cutoff_tmp[index_elem2, index_elem1, n] = cutoff_value
			cutoff_check[index_elem1, index_elem2, n] = true
			cutoff_check[index_elem2, index_elem1, n] = true
		end
	end

	# Check for missing values
	missing_indices = findall(x -> x == false, cutoff_check)
	if !isempty(missing_indices)
		throw(ArgumentError("Missing cutoff values for: " * join([
			"$(kd_name[idx[1]])-$(kd_name[idx[2]]) (nbody = $(idx[3]))"
			for idx in missing_indices
		], ", ")))
	end

	return Array{Float64}(cutoff_tmp)
end

"""
    parse_position(position_list::AbstractVector{<:AbstractVector{<:Real}}, num_atoms::Integer)

Parse and validate atomic positions.

# Arguments
- `position_list::AbstractVector`: List of position vectors for each atom
- `num_atoms::Integer`: Number of atoms

# Returns
- `Matrix{Float64}`: Matrix of atomic positions [3 × num_atoms]

# Throws
- `ErrorException` if:
  - Number of positions doesn't match num_atoms
  - Invalid position values
  - Missing position data

# Examples
```julia
positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
x_fractional = parse_position(positions, 2)
```
"""
function parse_position(
	position_list::AbstractVector{<:AbstractVector{<:Real}},
	num_atoms::Integer,
)::Matrix{Float64}
	position_tmp = fill(-1.0, 3, num_atoms)
	position_check = fill(false, 3, num_atoms)

	# Parse positions
	for (i, vec) in enumerate(position_list)
		if length(vec) != 3
			throw(ArgumentError("Position vector for atom $i must have 3 components, got $(length(vec))."))
		end
		position_tmp[:, i] = vec
		position_check[:, i] .= true
	end

	# Check for missing values
	if any(x -> x == false, position_check)
		throw(ArgumentError("Missing or invalid position data for some atoms."))
	end

	return Matrix{Float64}(position_tmp)
end

end # module
