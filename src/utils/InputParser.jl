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

export Parser

"""
	struct Parser

A structure that holds all input parameters for the simulation.

# Required Fields
- `name::String`: Name of the structure
- `num_atoms::Int`: Number of atoms in the structure
- `kd_name::Vector{String}`: List of chemical species names
- `nbody::Int`: Maximum number of bodies in interactions
- `lmax::Matrix{Int}`: Maximum angular momentum values [nkd x nbody]
- `cutoff_radii::Array{Float64, 3}`: Cutoff radii for interactions [nkd x nkd x nbody]
- `datafile::String`: Path to data file
- `lattice_vectors::Matrix{Float64}`: Lattice vectors (3x3 matrix)
- `kd_int_list::Vector{Int}`: List of chemical species indices
- `x_fractional::Matrix{Float64}`: Fractional coordinates of atoms

# Optional Fields
- `is_periodic::Vector{Bool}`: Periodicity flags for each direction (default: [true, true, true])
- `tolerance_sym::Float64`: Symmetry tolerance (default: 1e-3)
- `ndata::Int`: Number of data points (default: -1)
- `weight::Float64`: Weight parameter for regression (default: 0.0)

# Examples
```julia
parser = Parser(input_dict)
```
"""
struct Parser
	# Required parameters
	name::String
	num_atoms::Int
	kd_name::Vector{String}
	nbody::Int
	lmax::Matrix{Int}
	cutoff_radii::Array{Float64, 3}
	datafile::String
	lattice_vectors::Matrix{Float64}
	kd_int_list::Vector{Int}
	x_fractional::Matrix{Float64}

	# Optional parameters
	is_periodic::Vector{Bool}
	tolerance_sym::Float64
	ndata::Int
	weight::Float64
end

# Default values for optional parameters
const DEFAULT_VALUES = Dict{Symbol, Any}(
	:is_periodic => [true, true, true],
	:tolerance_sym => 1e-3,
	:ndata => -1,
	:weight => 0.0,
)

"""
	struct ValidationRule

A structure that defines a validation rule for a field.

# Fields
- `field::Symbol`: The field name to validate
- `validator::Function`: The validation function
- `error_message::Union{String, Function}`: The error message to display if validation fails
"""
struct ValidationRule
	field::Symbol
	validator::Function
	error_message::Union{String, Function}
end

"""
	const VALIDATION_RULES

Collection of validation rules for Parser parameters.
"""
const VALIDATION_RULES = [
	# Required parameters
	ValidationRule(:name, x -> !isempty(x), "Structure name cannot be empty"),
	ValidationRule(:num_atoms, x -> x > 0, x -> "Number of atoms must be positive, got $(x)"),
	ValidationRule(:kd_name, x -> !isempty(x), "Chemical species list cannot be empty"),
	ValidationRule(:nbody, x -> x > 0, x -> "nbody must be positive, got $(x)"),
	ValidationRule(:datafile, x -> !isempty(x), "Data file path cannot be empty"),
	ValidationRule(
		:lattice_vectors,
		x -> size(x) == (3, 3),
		x -> "Lattice vectors must be a 3x3 matrix, got $(size(x))",
	),
	ValidationRule(
		:kd_int_list,
		x -> x isa Vector{Int},
		"kd_int_list must be a vector of integers",
	),

	# Optional parameters
	ValidationRule(
		:is_periodic,
		x -> length(x) == 3,
		"Periodicity must be specified for all three directions",
	),
	ValidationRule(
		:tolerance_sym,
		x -> x > 0,
		x -> "Symmetry tolerance must be positive, got $(x)",
	),
	ValidationRule(:ndata, x -> x isa Int, x -> "Number of data must be an integer, got $(x)"),
	ValidationRule(:weight, x -> x >= 0, x -> "Weight must be non-negative, got $(x)"),
]

"""
	validate_kd_list(params::NamedTuple)::Nothing

Validate the chemical species list and indices.

# Arguments
- `params::NamedTuple`: The parameters to validate

# Returns
- `Nothing`

# Throws
- `ArgumentError` if validation fails
"""
function validate_kd_list(params::NamedTuple)::Nothing
	if length(params.kd_int_list) != params.num_atoms
		throw(
			ArgumentError(
				"Length of kd_list ($(length(params.kd_int_list))) must match number of atoms ($(params.num_atoms))",
			),
		)
	end

	kd_int_set = Set(params.kd_int_list)
	nkd = length(params.kd_name)
	if nkd != length(kd_int_set)
		throw(ArgumentError("Number of chemical species ($nkd) doesn't match kd_list entries"))
	elseif any(x -> (x < 1 || nkd < x), kd_int_set)
		throw(ArgumentError("Chemical species indices must be consecutive numbers starting from 1"))
	end
	return nothing
end

"""
	validate_parser_parameters(params::NamedTuple)::Nothing

Validate all parameters for Parser construction.

# Arguments
- `params::NamedTuple`: The parameters to validate

# Returns
- `Nothing`

# Throws
- `ArgumentError` if any parameter is invalid
"""
function validate_parser_parameters(params::NamedTuple)::Nothing
	# Apply basic validation rules
	for rule in VALIDATION_RULES
		value = getfield(params, rule.field)
		if !rule.validator(value)
			error_message =
				rule.error_message isa Function ? rule.error_message(value) : rule.error_message
			throw(ArgumentError(error_message))
		end
	end

	# Apply custom validation
	validate_kd_list(params)
	return nothing
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
			throw(
				ArgumentError("Required section \"$section\" is missing in the input dictionary."),
			)
		end
	end

	# Parse general parameters
	general_dict = input_dict["general"]
	name = get(general_dict, "name", "")::String
	num_atoms = get(general_dict, "nat", 0)::Int
	kd_name = get(general_dict, "kd", String[])::Vector{String}
	is_periodic = get(general_dict, "periodicity", DEFAULT_VALUES[:is_periodic])::Vector{Bool}

	# Parse symmetry parameters
	symmetry_dict = input_dict["symmetry"]
	tolerance_sym = get(symmetry_dict, "tolerance", DEFAULT_VALUES[:tolerance_sym])::Float64

	# Parse interaction parameters
	interaction_dict = input_dict["interaction"]
	nbody = get(interaction_dict, "nbody", 0)::Int
	lmax = parse_lmax(interaction_dict["lmax"], kd_name, nbody)
	cutoff_radii = parse_cutoff(interaction_dict["cutoff"], kd_name, nbody)

	# Parse regression parameters
	regression_dict = input_dict["regression"]
	weight = get(regression_dict, "weight", DEFAULT_VALUES[:weight])::Real
	datafile = get(regression_dict, "datafile", "")::String
	ndata = get(regression_dict, "ndata", DEFAULT_VALUES[:ndata])::Int

	# Parse structure parameters
	structure_dict = input_dict["structure"]
	lattice_vectors = hcat(structure_dict["lattice"]...)
	kd_int_list = get(structure_dict, "kd_list", Int[])::Vector{Int}
	positions = get(structure_dict, "position", Vector{Float64}[])::Vector{Vector{Float64}}
	x_fractional = parse_position(positions, num_atoms)

	# Validate all parameters
	params = (
		# Required parameters
		name = name,
		num_atoms = num_atoms,
		kd_name = kd_name,
		nbody = nbody,
		lmax = lmax,
		cutoff_radii = cutoff_radii,
		datafile = datafile,
		lattice_vectors = lattice_vectors,
		kd_int_list = kd_int_list,
		x_fractional = x_fractional,

		# Optional parameters
		is_periodic = is_periodic,
		tolerance_sym = tolerance_sym,
		ndata = ndata,
		weight = weight,
	)
	validate_parser_parameters(params)

	return Parser(
		params.name,
		params.num_atoms,
		params.kd_name,
		params.nbody,
		params.lmax,
		params.cutoff_radii,
		params.datafile,
		params.lattice_vectors,
		params.kd_int_list,
		params.x_fractional,
		params.is_periodic,
		params.tolerance_sym,
		params.ndata,
		params.weight,	
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
			throw(
				ArgumentError(
					"Length of lmax values for $key ($(length(value))) doesn't match nbody ($nbody).",
				),
			)
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
			throw(
				ArgumentError(
					"Length of cutoff values for $key ($(length(value))) doesn't match nbody ($nbody).",
				),
			)
		end
	end

	return nothing
end

"""
	parse_lmax(lmax_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString}, nbody::Integer)::Matrix{Int}

Parse and validate maximum angular momentum values.

# Arguments
- `lmax_dict::AbstractDict`: Dictionary of lmax values for each chemical species
- `kd_name::AbstractVector`: List of chemical species names
- `nbody::Integer`: Maximum number of bodies

# Returns
- `Matrix{Int}`: Matrix of lmax values [nkd × nbody]

# Throws
- `ArgumentError` if:
  - Length of lmax values doesn't match nbody
  - Missing values for any chemical species
  - Invalid lmax values
"""
function parse_lmax(
	lmax_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
	nbody::Integer,
)::Matrix{Int}
	kd_num = length(kd_name)
	lmax_tmp = fill(-1, kd_num, nbody)

	# Parse lmax values
	for (kd_index, kd) in enumerate(kd_name)
		if !haskey(lmax_dict, kd)
			throw(ArgumentError("Missing lmax values for chemical species \"$kd\"."))
		end

		values = lmax_dict[kd]
		if length(values) != nbody
			throw(
				ArgumentError(
					"Length of lmax values for $kd ($(length(values))) doesn't match nbody ($nbody).",
				),
			)
		end

		for j in 1:nbody
			value = values[j]
			lmax_tmp[kd_index, j] = value
		end
	end

	return Matrix{Int}(lmax_tmp)
end

"""
	parse_cutoff(cutoff_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString}, nbody::Integer)::Array{Float64, 3}

Parse and validate cutoff radii for interactions.

# Arguments
- `cutoff_dict::AbstractDict`: Dictionary of cutoff values for each pair of chemical species
- `kd_name::AbstractVector`: List of chemical species names
- `nbody::Integer`: Maximum number of bodies

# Returns
- `Array{Float64, 3}`: Array of cutoff radii [nkd × nkd × nbody]

# Throws
- `ArgumentError` if:
  - Length of cutoff values doesn't match nbody
  - Missing values for any pair of chemical species
  - Invalid cutoff values
"""
function parse_cutoff(
	cutoff_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
	nbody::Integer,
)::Array{Float64, 3}
	kd_num = length(kd_name)
	cutoff_tmp = fill(0.0, kd_num, kd_num, nbody)

	# Parse cutoff values
	for n in 1:nbody
		if n == 1
			cutoff_tmp[:, :, n] .= 0.0
			continue
		end

		for (key, value) in cutoff_dict
			if length(value) != nbody
				throw(
					ArgumentError(
						"Length of cutoff values for $key ($(length(value))) doesn't match nbody ($nbody).",
					),
				)
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
		end
	end

	return Array{Float64}(cutoff_tmp)
end

"""
	parse_position(position_list::AbstractVector{<:AbstractVector{<:Real}}, num_atoms::Integer)::Matrix{Float64}

Parse and validate atomic positions.

# Arguments
- `position_list::AbstractVector`: List of position vectors for each atom
- `num_atoms::Integer`: Number of atoms

# Returns
- `Matrix{Float64}`: Matrix of atomic positions [3 × num_atoms]

# Throws
- `ArgumentError` if:
  - Number of positions doesn't match num_atoms
  - Invalid position values
  - Missing position data
"""
function parse_position(
	position_list::AbstractVector{<:AbstractVector{<:Real}},
	num_atoms::Integer,
)::Matrix{Float64}
	position_tmp = fill(-1.0, 3, num_atoms)

	# Parse positions
	for (i, vec) in enumerate(position_list)
		if length(vec) != 3
			throw(
				ArgumentError(
					"Position vector for atom $i must have 3 components, got $(length(vec)).",
				),
			)
		end
		position_tmp[:, i] = vec
	end

	return Matrix{Float64}(position_tmp)
end

end # module
