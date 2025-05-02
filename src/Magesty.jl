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

include("utils/ConfigParser.jl")
include("utils/InputParser.jl")
include("utils/InputSetter.jl")
include("utils/RotationMatrix.jl")
include("utils/MySphericalHarmonics.jl")
using .ConfigParser
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
function System(input_dict::Dict{<:AbstractString, <:Any}, verbosity::Bool = true)
	config::Config4System = Config4System(input_dict)

	structure::Structure = Structure(config)
	verbosity && Structures.print_info(structure)

	symmetry::Symmetry = Symmetry(structure, config)
	verbosity && Symmetries.print_info(symmetry)

	cluster::Cluster = Cluster(structure, symmetry, config)
	verbosity && Clusters.print_info(cluster)

	basisset::BasisSet = BasisSet(structure, symmetry, cluster, config)
	verbosity && BasisSets.print_info(basisset)

	return System(structure, symmetry, cluster, basisset)
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
function System(toml_file::AbstractString, verbosity::Bool = true)
	try
		open(toml_file) do io
			toml = read(io, String)
			input_dict = TOML.parse(toml)
			return System(input_dict, verbosity)
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
	structure::Structure
	symmetry::Symmetry
	cluster::Cluster
	basisset::BasisSet
	optimize::SCEOptimizer
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
function SpinCluster(input_dict::Dict{<:AbstractString, <:Any}, verbosity::Bool = true)
	parser = Parser(input_dict)
	structure::Structure = set_system(parser)
	if verbosity
		Structures.print_info(structure)
	end
	symmetry::Symmetry = set_symmetry(parser, structure)
	if verbosity
		Symmetries.print_info(symmetry)
	end
	cluster::Cluster = set_cluster(parser, structure, symmetry)
	if verbosity
		Clusters.print_info(cluster)
	end
	basisset::BasisSet = set_basisset(parser, structure, symmetry, cluster)
	if verbosity
		BasisSets.print_info(basisset)
	end
	optimize::SCEOptimizer = set_optimize(parser, structure, symmetry, basisset)
	if verbosity
		Optimize.print_info(optimize)
	end

	return SpinCluster(structure, symmetry, cluster, basisset, optimize)
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
function SpinCluster(toml_file::AbstractString, verbosity::Bool = true)
	try
		open(toml_file) do io
			toml = read(io, String)
			input_dict = TOML.parse(toml)
			return SpinCluster(input_dict, verbosity)
		end
	catch e
		if isa(e, SystemError)
			throw(SystemError("Failed to read file: $toml_file"))
		else
			throw(ErrorException("Failed to parse TOML file: $toml_file"))
		end
	end
end

function SpinCluster(
	system::System,
	input_dict::Dict{<:AbstractString, <:Any},
	verbosity::Bool = true,
)
	parser = Parser(input_dict)
	optimize = set_optimize(parser, system.structure, system.symmetry, system.basisset)
	if verbosity
		Optimize.print_info(optimize)
	end
	return SpinCluster(
		system.structure,
		system.symmetry,
		system.cluster,
		system.basisset,
		optimize,
	)
end

function SpinCluster(system::System, toml_file::AbstractString, verbosity::Bool = true)
	try
		open(toml_file) do io
			toml = read(io, String)
			input_dict = TOML.parse(toml)
			return SpinCluster(system, input_dict, verbosity)
		end
	catch e
		if isa(e, SystemError)
			throw(SystemError("Failed to read file: $toml_file"))
		else
			throw(ErrorException("Failed to parse TOML file: $toml_file"))
		end
	end
end

function SpinCluster(
	spincluster::SpinCluster,
	input_dict::AbstractDict{<:AbstractString, <:Any},
	verbosity::Bool = true,
)
	parser = Parser(input_dict)
	sce_with_bias = vcat(spincluster.optimize.bias_term, spincluster.optimize.SCE)
	optimize = SCEOptimizer(
		spincluster.structure,
		spincluster.symmetry,
		spincluster.basisset,
		parser.weight,
		spincluster.optimize.spinconfig_dataset,
		sce_with_bias,
	)
	if verbosity
		Optimize.print_info(optimize)
	end
	return SpinCluster(
		spincluster.structure,
		spincluster.symmetry,
		spincluster.cluster,
		spincluster.basisset,
		optimize,
	)
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
	Clusters.print_info(sc.cluster)
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
function write_magfield_vertical_list(
	sc::SpinCluster,
	filename::AbstractString = "magfield_vertical_list.txt",
)
	if sc.optimize === nothing
		throw(ErrorException("Optimizer is not set"))
	end
	Optimize.write_magfield_vertical_list(sc.optimize, filename)
end

end
