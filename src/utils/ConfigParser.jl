module ConfigParser

export Config4System, Config4Optimize

struct ValidationRule
	field::Symbol
	validator::Function
	error_message::Union{String, Function}
end

"""
    Config4System

A structure that holds system configuration parameters for molecular dynamics simulations.

# Required Parameters
- `name::String`: Name of the system
- `num_atoms::Int`: Number of atoms in the system
- `kd_name::Vector{String}`: List of chemical species names
- `nbody::Int`: Maximum order of many-body interactions
- `lmax::Matrix{Int}`: Maximum angular momentum for each chemical species and interaction order
- `cutoff_radii::Array{Float64, 3}`: Cutoff radii for interactions between chemical species
- `lattice_vectors::Matrix{Float64}`: Lattice vectors of the system
- `kd_int_list::Vector{Int}`: List of chemical species indices for each atom
- `x_fractional::Matrix{Float64}`: Fractional coordinates of atoms

# Optional Parameters
- `is_periodic::Vector{Bool}`: Periodicity flags for each direction [default: [true, true, true]]
- `tolerance_sym::Float64`: Tolerance for symmetry operations [default: 1e-3]
"""
struct Config4System
	# required parameters
	name::String
	num_atoms::Int
	kd_name::Vector{String}
	nbody::Int
	lmax::Matrix{Int}
	cutoff_radii::Array{Float64, 3}
	lattice_vectors::Matrix{Float64}
	kd_int_list::Vector{Int}
	x_fractional::Matrix{Float64}

	# optional parameters
	is_periodic::Vector{Bool}
	tolerance_sym::Float64
end

# Default values for optional parameters
const DEFAULT_VALUES_SYSTEM = Dict{Symbol, Any}(
	:is_periodic => [true, true, true],
	:tolerance_sym => 1e-3,
)

# Validation rules for system configuration
const VALIDATION_RULES_SYSTEM = [
	# Required parameters
	ValidationRule(:name, x -> !isempty(x), "Structure name cannot be empty"),
	ValidationRule(:num_atoms, x -> x > 0, x -> "Number of atoms must be positive, got $(x)"),
	ValidationRule(:kd_name, x -> !isempty(x), "Chemical species list cannot be empty"),
	ValidationRule(:nbody, x -> x > 0, x -> "nbody must be positive, got $(x)"),
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
]

"""
    Config4System(input_dict::AbstractDict{<:AbstractString, Any})

Construct a Config4System from a dictionary of parameters.

# Arguments
- `input_dict::AbstractDict{<:AbstractString, Any}`: Dictionary containing system configuration parameters

# Returns
- `Config4System`: A new system configuration object

# Throws
- `ArgumentError`: If required parameters are missing or invalid
"""
function Config4System(input_dict::AbstractDict{<:AbstractString, Any})
	# Check required sections
	required_sections = ["general", "symmetry", "interaction", "structure"]
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
	is_periodic = get(general_dict, "periodicity", DEFAULT_VALUES_SYSTEM[:is_periodic])::Vector{Bool}

	# Parse symmetry parameters
	symmetry_dict = input_dict["symmetry"]
	tolerance_sym = get(symmetry_dict, "tolerance", DEFAULT_VALUES_SYSTEM[:tolerance_sym])::Float64

	# Parse interaction parameters
	interaction_dict = input_dict["interaction"]
	check_interaction_field(interaction_dict, kd_name)
	nbody = get(interaction_dict, "nbody", 0)::Int
	lmax = parse_lmax(interaction_dict["lmax"], kd_name, nbody)
	cutoff_radii = parse_cutoff(interaction_dict["cutoff"], kd_name, nbody)

	# Parse structure parameters
	structure_dict = input_dict["structure"]
	lattice_vectors = hcat(structure_dict["lattice"]...)
	kd_int_list = get(structure_dict, "kd_list", Int[])::Vector{Int}
	positions = get(structure_dict, "position", Vector{Float64}[])::Vector{Vector{Float64}}
	x_fractional = parse_position(positions, num_atoms)

	# Validate parameters
	params = (
		name = name,
		num_atoms = num_atoms,
		is_periodic = is_periodic,
		kd_name = kd_name,
		nbody = nbody,
		lmax = lmax,
		cutoff_radii = cutoff_radii,
		lattice_vectors = lattice_vectors,
		kd_int_list = kd_int_list,
		x_fractional = x_fractional,
		tolerance_sym = tolerance_sym,
	)
	validate_system_parameters(params)

	return Config4System(params...)
end

# Validate system parameters using defined rules
function validate_system_parameters(params::NamedTuple)::Nothing
	for rule in VALIDATION_RULES_SYSTEM
		value = getfield(params, rule.field)
		if !rule.validator(value)
			error_message =
				rule.error_message isa Function ? rule.error_message(value) : rule.error_message
			throw(ArgumentError(error_message))
		end
	end
	validate_kd_list(params)
	return nothing
end

# Validate kd_list parameters
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
    check_interaction_field(interaction_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString})

Validate interaction parameters in the configuration dictionary.

# Arguments
- `interaction_dict::AbstractDict{<:AbstractString, Any}`: Dictionary containing interaction parameters
- `kd_name::AbstractVector{<:AbstractString}`: List of chemical species names

# Throws
- `ArgumentError`: If required fields are missing or parameters are invalid
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
    parse_lmax(lmax_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString}, nbody::Integer)

Parse and validate maximum angular momentum values.

# Arguments
- `lmax_dict::AbstractDict{<:AbstractString, Any}`: Dictionary containing lmax values for each chemical species
- `kd_name::AbstractVector{<:AbstractString}`: List of chemical species names
- `nbody::Integer`: Maximum order of many-body interactions

# Returns
- `Matrix{Int}`: Matrix of lmax values for each chemical species and interaction order

# Throws
- `ArgumentError`: If lmax values are missing or invalid
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
    parse_cutoff(cutoff_dict::AbstractDict{<:AbstractString, Any}, kd_name::AbstractVector{<:AbstractString}, nbody::Integer)

Parse and validate cutoff radius values.

# Arguments
- `cutoff_dict::AbstractDict{<:AbstractString, Any}`: Dictionary containing cutoff values for each chemical species pair
- `kd_name::AbstractVector{<:AbstractString}`: List of chemical species names
- `nbody::Integer`: Maximum order of many-body interactions

# Returns
- `Array{Float64, 3}`: Array of cutoff radii for each chemical species pair and interaction order

# Throws
- `ArgumentError`: If cutoff values are missing or invalid
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
    parse_position(position_list::AbstractVector{<:AbstractVector{<:Real}}, num_atoms::Integer)

Parse and validate atomic positions.

# Arguments
- `position_list::AbstractVector{<:AbstractVector{<:Real}}`: List of position vectors for each atom
- `num_atoms::Integer`: Number of atoms in the system

# Returns
- `Matrix{Float64}`: Matrix of fractional coordinates for each atom

# Throws
- `ArgumentError`: If position vectors are invalid
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

"""
    Config4Optimize

A structure that holds optimization parameters.

# Required Parameters
- `datafile::String`: Path to the data file

# Optional Parameters
- `j_zero_thr::Float64`: Threshold for zero values [default: 1e-8]
- `ndata::Int`: Number of data points to use [default: -1]
- `weight::Float64`: Weight for optimization [default: 0.0]
"""
struct Config4Optimize
	# required parameters
	datafile::String

	# optional parameters
	j_zero_thr::Float64
	ndata::Int
	weight::Float64
end

# Default values for optional parameters
const DEFAULT_VALUES_OPTIMIZE = Dict{Symbol, Any}(
	:j_zero_thr => 1e-8,
	:ndata => -1,
	:weight => 0.0,
)

# Validation rules for optimization parameters
const VALIDATION_RULES_OPTIMIZE = [
	ValidationRule(:datafile, x -> !isempty(x), "Data file path cannot be empty"),
	ValidationRule(:j_zero_thr, x -> x > 0, x -> "j_zero_thr must be positive, got $(x)"),
	ValidationRule(:ndata, x -> x isa Int, x -> "ndata must be an integer, got $(x)"),
	ValidationRule(:weight, x -> 0 <= x <= 1, x -> "weight must be between 0 and 1, got $(x)"),
]

"""
    Config4Optimize(input_dict::AbstractDict{<:AbstractString, Any})

Construct a Config4Optimize from a dictionary of parameters.

# Arguments
- `input_dict::AbstractDict{<:AbstractString, Any}`: Dictionary containing optimization parameters

# Returns
- `Config4Optimize`: A new optimization configuration object

# Throws
- `ArgumentError`: If required parameters are missing or invalid
"""
function Config4Optimize(input_dict::AbstractDict{<:AbstractString, Any})
	# Check required sections
	required_sections = ["general", "regression"]
	for section in required_sections
		if !haskey(input_dict, section)
			throw(ArgumentError("Required section \"$section\" is missing in the input dictionary."))
		end
	end

	datafile = input_dict["general"]["datafile"]
	j_zero_thr = get(input_dict["regression"], "j_zero_thr", DEFAULT_VALUES_OPTIMIZE[:j_zero_thr])::Float64
	ndata = get(input_dict["regression"], "ndata", DEFAULT_VALUES_OPTIMIZE[:ndata])::Int
	weight = get(input_dict["regression"], "weight", DEFAULT_VALUES_OPTIMIZE[:weight])::Float64

	params = (
		datafile = datafile,
		j_zero_thr = j_zero_thr,
		ndata = ndata,
		weight = weight,
	)
	validate_optimize_parameters(params)

	return Config4Optimize(params...)
end

function validate_optimize_parameters(params::NamedTuple)::Nothing
	for rule in VALIDATION_RULES_OPTIMIZE
		value = getfield(params, rule.field)
		if !rule.validator(value)
			error_message =
				rule.error_message isa Function ? rule.error_message(value) : rule.error_message
			throw(ArgumentError(error_message))
		end
	end
	return nothing
end

end # module
