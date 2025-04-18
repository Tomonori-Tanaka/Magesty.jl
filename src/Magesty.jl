"""
	Magesty

The main module of the Magesty package, providing the entry point for spin cluster expansion analysis and optimization.

This module provides the following main features:
- Magnetic structure setup and management
- Symmetry analysis
- Cluster expansion
- Basis function generation
- Spin configuration optimization

# Usage
```julia
using Magesty

# Load configuration from a TOML file
sc = SpinCluster("config.toml")

# Display system information
print_info(sc)

# Output energy lists
write_energy_lists(sc)
```

# Main Types
- `System`: Collection of structure, symmetry, cluster, and basis set
- `SpinCluster`: Extension of System with optimization capabilities

# Main Functions
- `print_info`: Display detailed system information
- `write_energy_lists`: Output energy lists to a file
- `write_magfield_vertical_list`: Output magnetic field vertical components to a file

# Submodules
- `Structures`: Crystal structure processing
- `Symmetries`: Symmetry operations processing
- `Clusters`: Cluster expansion processing
- `BasisSets`: Basis function generation
- `Optimize`: Spin configuration optimization
"""
module Magesty

using Printf
using TOML

include("common/SortedContainer.jl")

include("types/AtomCell.jl")
include("types/AtomicIndices.jl")
include("types/UnitaryMatrixCl.jl")
include("types/SALC.jl")
include("types/SpinConfig.jl")

include("utils/InputParser.jl")
include("utils/InputSetter.jl")
include("utils/RotationMatrix.jl")
include("utils/MySphericalHarmonics.jl")
using .InputParser

include("Structure.jl")
include("Symmetry.jl")
include("Cluster.jl")
include("BasisSet.jl")
include("Optimize.jl")

using .Structures
using .Symmetries
using .Clusters
using .BasisSets
using .Optimize

export System, SpinCluster, print_info, write_energy_lists, write_magfield_vertical_list

"""
	System

A collection of structure, symmetry, cluster, and basis set.

# Fields
- `config::Parser`: Configuration parser
- `structure::Structure`: Crystal structure information
- `symmetry::Symmetry`: Symmetry operations
- `cluster::Cluster`: Cluster information
- `basisset::BasisSet`: Basis set information
"""
struct System
	config::Parser
	structure::Structure
	symmetry::Symmetry
	cluster::Cluster
	basisset::BasisSet
end

"""
	System(input_dict::Dict{<:AbstractString, <:Any})

Create a System from a dictionary of input parameters.

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing input parameters

# Returns
- `System`: A new System instance

# Throws
- `ErrorException` if required parameters are missing or invalid
"""
function System(input_dict::Dict{<:AbstractString, <:Any})
	parser::Parser = Parser(input_dict)
	structure::Structure = set_system(parser)
	symmetry::Symmetry = set_symmetry(parser, structure)
	cluster::Cluster = set_cluster(parser, structure, symmetry)
	basisset::BasisSet = set_basisset(parser, structure, symmetry, cluster)

	return System(parser, structure, symmetry, cluster, basisset)
end

"""
	System(toml_file::AbstractString)

Create a System from a TOML configuration file.

# Arguments
- `toml_file::AbstractString`: Path to the TOML configuration file

# Returns
- `System`: A new System instance

# Throws
- `SystemError` if the file cannot be read
- `ErrorException` if the TOML parsing fails
"""
function System(toml_file::AbstractString)
	try
		open(toml_file) do io
			toml = read(io, String)
			config = TOML.parse(toml)
			return System(config)
		end
	catch e
		if isa(e, SystemError)
			throw(SystemError("Failed to read file: $toml_file"))
		else
			throw(ErrorException("Failed to parse TOML file: $toml_file"))
		end
	end
end

"""
	SpinCluster

An extension of System with optimization capabilities.

# Fields
- `config::Parser`: Configuration parser
- `structure::Structure`: Crystal structure information
- `symmetry::Symmetry`: Symmetry operations
- `cluster::Cluster`: Cluster information
- `basisset::BasisSet`: Basis set information
- `optimize::Union{SCEOptimizer, Nothing}`: Optimizer instance or nothing
"""
struct SpinCluster
	config::Parser
	structure::Structure
	symmetry::Symmetry
	cluster::Cluster
	basisset::BasisSet
	optimize::Union{SCEOptimizer, Nothing}
end

"""
	SpinCluster(input_dict::Dict{<:AbstractString, <:Any})

Create a SpinCluster from a dictionary of input parameters.
Differ from System in that it also sets the optimizer.

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing input parameters

# Returns
- `SpinCluster`: A new SpinCluster instance

# Throws
- `ErrorException` if required parameters are missing or invalid
"""
function SpinCluster(input_dict::Dict{<:AbstractString, <:Any})
	parser = Parser(input_dict)
	structure::Structure = set_system(parser)
	symmetry::Symmetry = set_symmetry(parser, structure)
	cluster::Cluster = set_cluster(parser, structure, symmetry)
	basisset::BasisSet = set_basisset(parser, structure, symmetry, cluster)
	optimize = if parser.mode == "optimize"
		set_optimize(parser, structure, symmetry, basisset)
	else
		nothing
	end

	return SpinCluster(parser, structure, symmetry, cluster, basisset, optimize)
end

"""
	SpinCluster(toml_file::AbstractString)

Create a SpinCluster from a TOML configuration file.

# Arguments
- `toml_file::AbstractString`: Path to the TOML configuration file

# Returns
- `SpinCluster`: A new SpinCluster instance

# Throws
- `SystemError` if the file cannot be read
- `ErrorException` if the TOML parsing fails
"""
function SpinCluster(toml_file::AbstractString)
	try
		open(toml_file) do io
			toml = read(io, String)
			config = TOML.parse(toml)
			return SpinCluster(config)
		end
	catch e
		if isa(e, SystemError)
			throw(SystemError("Failed to read file: $toml_file"))
		else
			throw(ErrorException("Failed to parse TOML file: $toml_file"))
		end
	end
end

"""
	SpinCluster(system::System)

Create a SpinCluster from an existing System.

# Arguments
- `system::System`: An existing System instance

# Returns
- `SpinCluster`: A new SpinCluster instance
"""
function SpinCluster(system::System)
	optimize = if system.config.mode == "optimize"
		set_optimize(system.config, system.structure, system.symmetry, system.basisset)
	else
		nothing
	end

	return SpinCluster(system.config, system.structure, system.symmetry, system.cluster, system.basisset, optimize)
end

"""
	print_info(sc::SpinCluster)

Print detailed information about the SpinCluster.

# Arguments
- `sc::SpinCluster`: The SpinCluster to display information about
"""
function print_info(sc::SpinCluster)
	println(
		"""
		+-----------------------------------+
		|          Magesty v0.1.0           |
		+-----------------------------------+

		""",
	)

	Structures.print_info(sc.structure)
	Symmetries.print_info(sc.symmetry)
	BasisSets.print_info(sc.basisset)
	Optimize.print_info(sc.optimize)
end

"""
	write_energy_lists(sc::SpinCluster, filename::AbstractString = "energy_lists.txt")

Write energy lists to a file.

# Arguments
- `sc::SpinCluster`: The SpinCluster containing energy data
- `filename::AbstractString`: Output file name (default: "energy_lists.txt")

# Throws
- `ErrorException` if the optimizer is not set
"""
function write_energy_lists(sc::SpinCluster, filename::AbstractString = "energy_lists.txt")
	if sc.optimize === nothing
		throw(ErrorException("Optimizer is not set"))
	end
	Optimize.write_energy_lists(sc.optimize, filename)
end

"""
	write_magfield_vertical_list(sc::SpinCluster, filename::AbstractString = "magfield_vertical_list.txt")

Write magnetic field vertical components to a file.

# Arguments
- `sc::SpinCluster`: The SpinCluster containing magnetic field data
- `filename::AbstractString`: Output file name (default: "magfield_vertical_list.txt")

# Throws
- `ErrorException` if the optimizer is not set
"""
function write_magfield_vertical_list(sc::SpinCluster, filename::AbstractString = "magfield_vertical_list.txt")
	if sc.optimize === nothing
		throw(ErrorException("Optimizer is not set"))
	end
	Optimize.write_magfield_vertical_list(sc.optimize, filename)
end

end
