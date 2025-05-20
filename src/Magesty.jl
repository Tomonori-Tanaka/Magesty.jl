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

include("common/version.jl")
using .Version

include("common/SortedContainer.jl")

include("types/AtomCell.jl")
include("types/AtomicIndices.jl")
include("types/UnitaryMatrixCl.jl")
include("types/SALC.jl")
include("SpinConfig.jl")

include("utils/ConfigParser.jl")
include("utils/RotationMatrix.jl")
include("utils/MySphericalHarmonics.jl")
using .ConfigParser

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

include("utils/Write.jl")
using .Write
include("utils/CalcEnergy.jl")
using .CalcEnergy

export System,
	SpinCluster, print_info, write_energy_lists, write_magfield_vertical_list, VERSION, Write

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
	config::Config4System = Config4System(input_dict)
	structure::Structure = Structure(config)
	verbosity && Structures.print_info(structure)

	symmetry::Symmetry = Symmetry(structure, config)
	verbosity && Symmetries.print_info(symmetry)

	cluster::Cluster = Cluster(structure, symmetry, config)
	verbosity && Clusters.print_info(cluster)

	basisset::BasisSet = BasisSet(structure, symmetry, cluster, config)
	verbosity && BasisSets.print_info(basisset)

	optimize::SCEOptimizer = SCEOptimizer(structure, symmetry, basisset, config)
	verbosity && Optimize.print_info(optimize)

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
	config::Config4Optimize = Config4Optimize(input_dict)
	optimize = SCEOptimizer(system.structure, system.symmetry, system.basisset, config)
	verbosity && Optimize.print_info(optimize)
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
	config::Config4Optimize = Config4Optimize(input_dict)
	sce_with_bias = vcat(spincluster.optimize.reference_energy, spincluster.optimize.SCE)
	optimize = SCEOptimizer(
		spincluster.structure,
		spincluster.symmetry,
		spincluster.basisset,
		config.weight,
		spincluster.optimize.spinconfig_list,
		sce_with_bias,
	)
	verbosity && Optimize.print_info(optimize)

	return SpinCluster(
		spincluster.structure,
		spincluster.symmetry,
		spincluster.cluster,
		spincluster.basisset,
		optimize,
	)
end

function SpinCluster(
	spincluster::SpinCluster,
	weight::Real,
	spinconfig_list::AbstractVector{SpinConfig},
	verbosity::Bool = true,
)
	optimize = SCEOptimizer(
		spincluster.structure,
		spincluster.symmetry,
		spincluster.basisset,
		weight,
		spinconfig_list,
		vcat(spincluster.optimize.reference_energy, spincluster.optimize.SCE),
	)
	verbosity && Optimize.print_info(optimize)
	return SpinCluster(
		spincluster.structure,
		spincluster.symmetry,
		spincluster.cluster,
		spincluster.basisset,
		optimize,
	)
end

function calc_energy(sc::SpinCluster, spin_config::AbstractMatrix{<:Real})
	if sc.structure.supercell.num_atoms != size(spin_config, 2)
		num_atoms = sc.structure.supercell.num_atoms
		throw(
			ArgumentError(
				"spin_config must be 3xN matrix where N is the number of atoms in the supercell. $num_atoms",
			),
		)
	end
	return CalcEnergy.calc_energy(sc.basisset.salc_list, spin_config, sc.symmetry, sc.optimize)
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
		|          Magesty v$(VERSION)      |
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
	write_sce2xml(structure::Structure, basis_set::BasisSet, optimize::SCEOptimizer, filename::AbstractString="jphi.xml"; write_jphi::Bool=true)
	write_sce2xml(sc::SpinCluster, filename::AbstractString="jphi.xml"; write_jphi::Bool=true)

Write the structure, basis set, and optimization results to an XML file in SCE format.

# Arguments
- `structure::Structure`: Crystal structure information
- `basis_set::BasisSet`: Basis set information
- `optimize::SCEOptimizer`: Optimization results
- `sc::SpinCluster`: Spin cluster object containing structure, basis set, and optimization results
- `filename::AbstractString="jphi.xml"`: Output XML file name
- `write_jphi::Bool=true`: Whether to write the J_ij parameters

# Examples
```julia
# Using individual components
write_sce2xml(structure, basis_set, optimizer)
write_sce2xml(structure, basis_set, optimizer, "output.xml", write_jphi=false)

# Using SpinCluster object
write_sce2xml(spin_cluster)
write_sce2xml(spin_cluster, "output.xml", write_jphi=false)
```
"""
function write_sce2xml(
	sc::SpinCluster,
	filename::AbstractString = "jphi.xml";
	write_jphi::Bool = true,
)
	write_sce2xml(
		sc.structure,
		sc.symmetry,
		sc.basisset,
		sc.optimize,
		filename;
		write_jphi = write_jphi,
	)
end
function write_sce2xml(
	structure::Structure,
	symmetry::Symmetry,
	basis_set::BasisSet,
	optimize::SCEOptimizer,
	filename::AbstractString = "jphi.xml";
	write_jphi::Bool = true,
)
	Write.write_sce2xml(structure, symmetry, basis_set, optimize, filename; write_jphi = write_jphi)
end

function write_energy_info(sc::SpinCluster, filename::AbstractString = "energy.txt")
	Write.write_energy_info(sc.optimize, filename)
end

function write_lmf_flattened(sc::SpinCluster, filename::AbstractString = "lmf_flattened.txt")
	Write.write_lmf_flattened(sc.optimize, filename)
end

end # module Magesty
