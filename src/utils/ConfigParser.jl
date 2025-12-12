module ConfigParser
using OffsetArrays

export Config4System, Config4Optimize

struct ValidationRule
	field::Symbol
	validator::Function
	error_message::Union{String, Function}
end

# Default values for system configuration
const DEFAULT_VALUES_SYSTEM = Dict{Symbol, Any}(
	:is_periodic => [true, true, true],
	:tolerance_sym => 1e-3,
	:isotropy => false,
)

# Validation rules for system configuration
const VALIDATION_RULES_SYSTEM = [
	# Required parameters
	ValidationRule(:name, x -> !isempty(x), "Structure name cannot be empty"),
	ValidationRule(
		:num_atoms,
		x -> x > 0,
		x -> "Number of atoms must be positive, got $(x)",
	),
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
	Config4System

A structure that holds system configuration parameters for molecular dynamics simulations.

# Required Parameters
- `name::String`: Name of the system
- `num_atoms::Int`: Number of atoms in the system
- `kd_name::Vector{String}`: List of chemical species names
- `nbody::Int`: Maximum order of many-body interactions
- `body1_lmax::Vector{Int}`: Maximum angular momentum for each body 1 and interaction order
- `bodyn_lmax::OffsetArray{Int, 1}`: Maximum angular momentum for each body n and interaction order
- `bodyn_cutoff::OffsetArray{Float64, 3}`: Cutoff radii for interactions between elements
- `lattice_vectors::Matrix{Float64}`: Lattice vectors of the system
- `kd_int_list::Vector{Int}`: List of chemical species indices for each atom
- `x_fractional::Matrix{Float64}`: Fractional coordinates of atoms

# Optional Parameters
- `is_periodic::Vector{Bool}`: Periodicity flags for each direction [default: [true, true, true]]
- `tolerance_sym::Float64`: Tolerance for symmetry operations [default: 1e-3]
- `isotropy::Bool`: If `true`, only include isotropic terms (Lf=0) [default: false]
"""
struct Config4System
	# required parameters
	name::String
	num_atoms::Int
	kd_name::Vector{String}
	nbody::Int
	body1_lmax::Vector{Int}
	bodyn_lsum::OffsetArray{Int, 1}
	bodyn_cutoff::OffsetArray{Float64, 3}
	lattice_vectors::Matrix{Float64}
	kd_int_list::Vector{Int}
	x_fractional::Matrix{Float64}

	# optional parameters
	is_periodic::Vector{Bool}
	tolerance_sym::Float64
	isotropy::Bool

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
					ArgumentError(
						"Required section \"$section\" is missing in the input dictionary.",
					),
				)
			end
		end

		# Parse general parameters
		general_dict = input_dict["general"]
		name = general_dict["name"]::String
		num_atoms = general_dict["nat"]::Int
		kd_name = general_dict["kd"]
		is_periodic = get(
			general_dict,
			"periodicity",
			DEFAULT_VALUES_SYSTEM[:is_periodic],
		)::Vector{Bool}

		# Parse symmetry parameters
		symmetry_dict = input_dict["symmetry"]
		tolerance_sym =
			get(symmetry_dict, "tolerance", DEFAULT_VALUES_SYSTEM[:tolerance_sym])::Float64
		isotropy = get(symmetry_dict, "isotropy", DEFAULT_VALUES_SYSTEM[:isotropy])::Bool

		# Parse interaction parameters
		interaction_dict = input_dict["interaction"]
		nbody = interaction_dict["nbody"]::Int
		if nbody < 1
			throw(ArgumentError("nbody must be positive, got $(nbody)"))
		end
		body1_lmax::Vector{Int} =
			parse_interaction_body1(interaction_dict, kd_name)

		bodyn_lsum::OffsetArray{Int, 1}, bodyn_cutoff::OffsetArray{Float64, 3} =
			parse_interaction_bodyn(interaction_dict, kd_name, nbody)
		# lmax = parse_lmax(interaction_dict["lmax"], kd_name, nbody)
		# cutoff_radii = parse_cutoff(interaction_dict["cutoff"], kd_name, nbody)

		# Parse structure parameters
		structure_dict = input_dict["structure"]
		lattice_vectors = hcat(structure_dict["lattice"]...)
		kd_int_list = structure_dict["kd_list"]::Vector{Int}
		positions_any = structure_dict["position"]
		positions = [Float64.(vec) for vec in positions_any]
		x_fractional = parse_position(positions, num_atoms)

		# Validate parameters
		params = (
			# required parameters
			name = name,
			num_atoms = num_atoms,
			kd_name = kd_name,
			nbody = nbody,
			body1_lmax = body1_lmax,
			bodyn_lsum = bodyn_lsum,
			bodyn_cutoff = bodyn_cutoff,
			lattice_vectors = lattice_vectors,
			kd_int_list = kd_int_list,
			x_fractional = x_fractional,
			# optional parameters
			is_periodic = is_periodic,
			tolerance_sym = tolerance_sym,
			isotropy = isotropy,
		)


		return new(params...)
	end
end

function parse_interaction_body1(
	interaction_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
)::Vector{Int}
	# Check if body1 key exists, if not return default values
	if !haskey(interaction_dict, "body1")
		# Return default values (all zeros) if body1 section is missing
		return fill(0, length(kd_name))
	end

	body1_dict = interaction_dict["body1"]
	result = fill(-1, length(kd_name))

	for (key, value) in body1_dict["lmax"]
		index = findfirst(x -> x == key, kd_name)
		if isnothing(index)
			throw(ArgumentError("Invalid chemical species \"$key\" in interaction.body1 field."))
		end
		result[index] = value
	end

	# Check for missing species only if body1 section exists
	missed_kd = kd_name[result .== -1]
	if !isempty(missed_kd)
		throw(
			ArgumentError(
				"lmax value (interaction.body1) is missing for some chemical species: $(join(missed_kd, ", "))",
			),
		)
	end

	return result
end

function parse_interaction_bodyn(
	interaction_dict::AbstractDict{<:AbstractString, Any},
	kd_name::AbstractVector{<:AbstractString},
	nbody::Integer,
)::Tuple{OffsetArray{Int, 1}, OffsetArray{Float64, 3}}
	# Handle nbody = 1 case: return empty arrays with appropriate ranges
	if nbody == 1
		# Create empty arrays with range 2:1 (empty range)
		result_lsum = OffsetArray(Int[], 2:1)
		nkd = length(kd_name)
		# Create empty 3D array with proper dimensions
		result_cutoff = OffsetArray(zeros(Float64, 0, nkd, nkd), 2:1, 1:nkd, 1:nkd)
		return result_lsum, result_cutoff
	end

	# Parse lsum values for nbody >= 2
	result_lsum = OffsetArray(fill(Int(-1), 2:nbody), 2:nbody)
	for n in 2:nbody
		bodyn_dict = interaction_dict["body$n"]
		result_lsum[n] = bodyn_dict["lsum"]::Int
	end

	missed_lsum::Vector{Int} = findall(x -> x == -1, result_lsum)
	if !isempty(missed_lsum)
		throw(ArgumentError("lsum value is missing for some body: $(join(missed_lsum, ", "))"))
	end

	# Parse cutoff values for nbody >= 2
	nkd = length(kd_name)
	result_cutoff = OffsetArray(zeros(Float64, nbody-1, nkd, nkd), 2:nbody, 1:nkd, 1:nkd)
	initialized_check = OffsetArray(fill(false, nbody-1, nkd, nkd), 2:nbody, 1:nkd, 1:nkd)
	for n in 2:nbody
		bodyn_dict = interaction_dict["body$n"]
		for (key, value) in bodyn_dict["cutoff"]
			# key is "<elem1>-<elem2>"
			elem1, elem2 = split(key, "-")
			index_elem1 = findfirst(x -> x == elem1, kd_name)
			index_elem2 = findfirst(x -> x == elem2, kd_name)
			if isnothing(index_elem1) || isnothing(index_elem2)
				throw(
					ArgumentError(
						"Invalid chemical species pair \"$key\" in cutoff dictionary.",
					),
				)
			end
			result_cutoff[n, index_elem1, index_elem2] = value
			result_cutoff[n, index_elem2, index_elem1] = value
			initialized_check[n, index_elem1, index_elem2] = true
			initialized_check[n, index_elem2, index_elem1] = true
		end
	end

	for n in 2:nbody
		if !all(initialized_check[n, :, :])
			throw(
				ArgumentError(
					"cutoff value is missing for body $(n)",
				),
			)
		end
	end

	return result_lsum, result_cutoff
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
	if isempty(position_list)
		throw(ArgumentError("Position list cannot be empty"))
	end

	if length(position_list) != num_atoms
		throw(
			ArgumentError(
				"Number of positions ($(length(position_list))) must match number of atoms ($num_atoms)",
			),
		)
	end

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


# Default values for optimization parameters
const DEFAULT_VALUES_OPTIMIZE = Dict{Symbol, Any}(
	:ndata => -1,
	:weight => 0.0, # 0 means magnetic field only, 1 means energy only
	:alpha => 0.0,# 0 means ridge regression, 1 means LASSO regression
	:lambda => 0.0, # 0 means no regularization
)

# Validation rules for optimization parameters
const VALIDATION_RULES_OPTIMIZE = [
	ValidationRule(:datafile, x -> !isempty(x), "Data file path cannot be empty"),
	ValidationRule(:ndata, x -> x isa Int, x -> "ndata must be an integer, got $(x)"),
	ValidationRule(
		:weight,
		x -> 0 <= x <= 1,
		x -> "weight must be between 0 and 1, got $(x)",
	),
	ValidationRule(
		:alpha,
		x -> 0.0 <= x <= 1.0,
		x -> "alpha must be between 0 and 1, got $(x)",
	),
	ValidationRule(
		:lambda,
		x -> x isa Real && x >= 0.0,
		x -> "lambda must be a non-negative real number, got $(x)",
	),
]
"""
	Config4Optimize

A structure that holds optimization parameters.

# Required Parameters
- `datafile::String`: Path to the data file

# Optional Parameters
- `ndata::Int`: Number of data points to use [default: -1]
- `weight::Float64`: Weight for optimization [default: 0.0]
"""
struct Config4Optimize
	# required parameters
	datafile::String

	# optional parameters
	ndata::Int
	weight::Float64
	alpha::Float64
	lambda::Float64

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
		required_sections = ["regression"]
		for section in required_sections
			if !haskey(input_dict, section)
				throw(
					ArgumentError(
						"Required section \"$section\" is missing in the input dictionary.",
					),
				)
			end
		end

		datafile = input_dict["regression"]["datafile"]::String
		ndata = get(input_dict["regression"], "ndata", DEFAULT_VALUES_OPTIMIZE[:ndata])::Int
		weight = Float64(get(
			input_dict["regression"],
			"weight",
			DEFAULT_VALUES_OPTIMIZE[:weight],
		))

		alpha = Float64(get(
			input_dict["regression"],
			"alpha",
			DEFAULT_VALUES_OPTIMIZE[:alpha],
		))

		lambda = Float64(get(
			input_dict["regression"],
			"lambda",
			DEFAULT_VALUES_OPTIMIZE[:lambda],
		))

		params = (
			datafile = datafile,
			ndata = ndata,
			weight = weight,
			alpha = alpha,
			lambda = lambda,
		)
		validate_optimize_parameters(params)

		return new(params...)
	end
end

function validate_optimize_parameters(params::NamedTuple)::Nothing
	for rule in VALIDATION_RULES_OPTIMIZE
		value = getfield(params, rule.field)
		if !rule.validator(value)
			error_message =
				rule.error_message isa Function ? rule.error_message(value) :
				rule.error_message
			throw(ArgumentError(error_message))
		end
	end
	return nothing
end

end # module
