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

```

# Main Types
- `System`: Collection of structure, symmetry, cluster, and basis set
- `SpinCluster`: Extension of System with optimization capabilities

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
using Dates

include("common/version.jl")
using .Version

include("common/SortedContainer.jl")

include("types/AtomCells.jl")
include("types/AtomicIndices.jl")
include("types/UnitaryMatrixCl.jl")
include("types/SALCs.jl")
include("SpinConfigs.jl")
using .SpinConfigs

include("utils/ConfigParser.jl")
include("utils/RotationMatrix.jl")
include("utils/MySphericalHarmonics.jl")
using .ConfigParser

include("Structures.jl")
include("Symmetries.jl")
include("Clusters.jl")
include("BasisSets.jl")
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

export System, SpinCluster, VERSION

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
	System

Create a `System` instance from either a dictionary of input parameters or a TOML configuration file.

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing input parameters.
- `toml_file::AbstractString`: Path to the TOML configuration file.
- `verbosity::Bool=true`: Whether to print detailed information during initialization.

# Returns
- `System`: A new `System` instance containing structure, symmetry, cluster, and basis set.

# Throws
- `SystemError`: If the TOML file cannot be read.
- `ErrorException`: If required parameters are missing, invalid, or the TOML parsing fails.

# Examples
```julia
# Create a System from a dictionary
input_dict = Dict("key" => "value")
system = System(input_dict)

# Create a System from a TOML file
system = System("config.toml")
```
"""
function System(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)

	if verbosity
		print_header()
	end
	config::Config4System = Config4System(input_dict)

	structure::Structure = Structure(config, verbosity = verbosity)

	symmetry::Symmetry = Symmetry(structure, config, verbosity = verbosity)

	cluster::Cluster = Cluster(structure, symmetry, config, verbosity = verbosity)

	basisset::BasisSet = BasisSet(structure, symmetry, cluster, config, verbosity = verbosity)

	return System(structure, symmetry, cluster, basisset)
end

function System(toml_file::AbstractString; verbosity::Bool = true)
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
- `optimize::Optimizer`: Optimizer instance
"""
struct SpinCluster
	structure::Structure
	symmetry::Symmetry
	cluster::Cluster
	basisset::BasisSet
	optimize::Optimizer
end

"""
	SpinCluster

Create a `SpinCluster` instance from either a dictionary of input parameters, a TOML configuration file, or an existing `System` instance. This is an extension of `System` that includes optimization capabilities.

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing input parameters.
- `toml_file::AbstractString`: Path to the TOML configuration file.
- `system::System`: An existing `System` instance.
- `verbosity::Bool=true`: Whether to print detailed information during initialization.

# Returns
- `SpinCluster`: A new `SpinCluster` instance containing structure, symmetry, cluster, basis set, and optimizer.

# Throws
- `SystemError`: If the TOML file cannot be read.
- `ErrorException`: If required parameters are missing, invalid, or the TOML parsing fails.

# Examples
```julia
# Create a SpinCluster from a dictionary
input_dict = Dict("key" => "value")
spin_cluster = SpinCluster(input_dict)

# Create a SpinCluster from a TOML file
spin_cluster = SpinCluster("config.toml")

# Create a SpinCluster from an existing System
system = System("config.toml")
spin_cluster = SpinCluster(system)
```
"""
function SpinCluster(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)

	if verbosity
		print_header()
	end

	config_system::Config4System = Config4System(input_dict)
	structure::Structure = Structure(config_system, verbosity = verbosity)

	symmetry::Symmetry = Symmetry(structure, config_system, verbosity = verbosity)

	cluster::Cluster = Cluster(structure, symmetry, config_system, verbosity = verbosity)

	basisset::BasisSet =
		BasisSet(structure, symmetry, cluster, config_system, verbosity = verbosity)

	config_optimize::Config4Optimize = Config4Optimize(input_dict)
	optimize::Optimizer =
		Optimizer(structure, symmetry, basisset, config_optimize, verbosity = verbosity)

	return SpinCluster(structure, symmetry, cluster, basisset, optimize)
end

function SpinCluster(toml_file::AbstractString; verbosity::Bool = true)
	try
		open(toml_file) do io
			toml = read(io, String)
			input_dict = TOML.parse(toml)
			return SpinCluster(input_dict, verbosity = verbosity)
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
	SpinCluster(system::System, input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)

Creates a `SpinCluster` instance by extending an existing `System` object with optimization capabilities. 
This constructor uses a dictionary of input parameters to configure the optimization process.

# Arguments
- `system::System`: An existing `System` instance containing structure, symmetry, cluster, and basis set information.
- `input_dict::Dict{<:AbstractString, <:Any}`: A dictionary containing input parameters for optimization.
- `verbosity::Bool=true`: Whether to print detailed information during initialization.

# Returns
- `SpinCluster`: A new `SpinCluster` instance containing structure, symmetry, cluster, basis set, and optimizer.

# Throws
- `ErrorException`: If required parameters are missing or invalid.

# Example
```julia
input_dict = Dict("key" => "value")
system = System(input_dict)
spin_cluster = SpinCluster(system, input_dict)
```
"""
function SpinCluster(
	system::System,
	input_dict::Dict{<:AbstractString, <:Any},
	;
	verbosity::Bool = true,
)
	if verbosity
		print_header()
	end
	config::Config4Optimize = Config4Optimize(input_dict)
	optimize =
		Optimizer(system.structure, system.symmetry, system.basisset, config, verbosity = verbosity)
	return SpinCluster(
		system.structure,
		system.symmetry,
		system.cluster,
		system.basisset,
		optimize,
	)
end

"""
	SpinCluster
Create a `SpinCluster` instance from a `System` and a list of `SpinConfig` objects.
This constructor is used when the optimization process is based on predefined spin configurations.
"""
function SpinCluster(
	system::System,
	input_dict::AbstractDict{<:AbstractString, <:Any},
	spinconfig_list::AbstractVector{SpinConfig},
	;
	verbosity::Bool = true,
)
	if verbosity
		print_header()
	end
	config::Config4Optimize = Config4Optimize(input_dict)
	optimize = Optimizer(
		system.structure,
		system.symmetry,
		system.basisset,
		config.alpha,
		config.lambda,
		config.weight,
		spinconfig_list,
		verbosity = verbosity,
	)
	return SpinCluster(
		system.structure,
		system.symmetry,
		system.cluster,
		system.basisset,
		optimize,
	)

end


"""
	calc_energy(sc::SpinCluster, spin_config::AbstractMatrix{<:Real})

Calculate the energy of a spin configuration using the spin cluster expansion.

# Arguments
- `sc::SpinCluster`: A `SpinCluster` instance containing structure, symmetry, basis set, and optimization information.
- `spin_config::AbstractMatrix{<:Real}`: A 3xN matrix representing the spin configuration, where N is the number of atoms in the supercell.

# Returns
- `Float64`: The calculated energy of the given spin configuration.

# Throws
- `ArgumentError`: If the number of columns in `spin_config` does not match the number of atoms in the supercell.

# Example
```julia
spin_config = rand(3, sc.structure.supercell.num_atoms) # Random spin configuration
energy = calc_energy(sc, spin_config)
```
"""
function calc_energy(spincluster::SpinCluster, spin_config::AbstractMatrix{<:Real})::Float64
	if spincluster.structure.supercell.num_atoms != size(spin_config, 2)
		num_atoms = spincluster.structure.supercell.num_atoms
		throw(
			ArgumentError(
				"spin_config must be 3xN matrix where N is the number of atoms in the supercell. $num_atoms",
			),
		)
	end
	return CalcEnergy.calc_energy(
		spincluster.basisset.salc_list,
		spin_config,
		spincluster.symmetry,
		spincluster.optimize,
	)
end

function calc_torque(spincluster::SpinCluster, spin_config::AbstractMatrix{<:Real})::Matrix{Float64}
	if spincluster.structure.supercell.num_atoms != size(spin_config, 2)
		num_atoms = spincluster.structure.supercell.num_atoms
		throw(
			ArgumentError(
				"spin_config must be 3xN matrix where N is the number of atoms in the supercell. $num_atoms",
			),
		)
	end
	return CalcEnergy.calc_torque(
		spincluster.basisset.salc_list,
		spin_config,
		spincluster.symmetry,
		spincluster.optimize,
	)
end

"""
	write_xml(structure::Structure, basis_set::BasisSet, optimize::Optimizer, filename::AbstractString="jphi.xml"; write_jphi::Bool=true)
	write_xml(sc::SpinCluster, filename::AbstractString="jphi.xml"; write_jphi::Bool=true)

Write the structure, basis set, and optimization results to an XML file in SCE format.

# Arguments
- `structure::Structure`: Crystal structure information
- `basis_set::BasisSet`: Basis set information
- `optimize::Optimizer`: Optimization results
- `sc::SpinCluster`: Spin cluster object containing structure, basis set, and optimization results
- `filename::AbstractString="jphi.xml"`: Output XML file name
- `write_jphi::Bool=true`: Whether to write the J_ij parameters

# Examples
```julia
# Using individual components
write_xml(structure, basis_set, optimizer)
write_xml(structure, basis_set, optimizer, "output.xml", write_jphi=false)

# Using SpinCluster object
write_xml(spin_cluster)
write_xml(spin_cluster, "output.xml", write_jphi=false)
```
"""
function write_xml(
	sc::SpinCluster,
	filename::AbstractString = "jphi.xml";
	write_jphi::Bool = true,
)
	write_xml(
		sc.structure,
		sc.symmetry,
		sc.basisset,
		sc.optimize,
		filename;
		write_jphi = write_jphi,
	)
end
function write_xml(
	structure::Structure,
	symmetry::Symmetry,
	basis_set::BasisSet,
	optimize::Optimizer,
	filename::AbstractString = "jphi.xml";
	write_jphi::Bool = true,
)
	Write.write_xml(
		structure,
		symmetry,
		basis_set,
		optimize,
		filename;
		write_jphi = write_jphi,
	)
end

function write_energies(
	sc::SpinCluster,
	filename::AbstractString = "energy_list.txt",
)

	observed_energy_list = [spinconfig.energy for spinconfig in sc.optimize.spinconfig_list]
	predicted_energy_list = sc.optimize.predicted_energy_list
	# Write to file
	try
		open(filename, "w") do f
			# Write header
			println(f, "# data index,    DFT_Energy,    SCE_Energy\n# unit of energy is eV")

			# Write data
			idx_width = ndigits(length(observed_energy_list))
			for i in eachindex(sc.optimize.spinconfig_list)
				str = @sprintf(
					" %*d    % 15.10e    % 15.10e\n",
					idx_width,
					i,
					observed_energy_list[i],
					predicted_energy_list[i],
				)
				write(f, str)
			end
		end
	catch e
		@error "Failed to write lists to file" exception = (e, catch_backtrace())
		rethrow(e)
	end
end

function write_torque_list(
	sc::SpinCluster,
	filename::AbstractString = "torque_list.txt",
)
    predicted_torque_list::Vector{Matrix{Float64}} = sc.optimize.predicted_torque_list
    observed_torque_list::Vector{Matrix{Float64}} = [spinconfig.torques for spinconfig in sc.optimize.spinconfig_list]

	# Write to file
	open(filename, "w") do f
		# Write header
		println(
			f,
			"# atom index,    element,   DFT_torque_x,    DFT_torque_y,    DFT_torque_z,    SCE_torque_x,    SCE_torque_y,    SCE_torque_z\n# unit of torque is eV",
		)

		# Write data
		idx_width = ndigits(length(sc.optimize.spinconfig_list))
		element_string_list = [sc.structure.kd_name[elm_idx] for elm_idx in sc.structure.supercell.kd_int_list]
		element_width = maximum(length.(element_string_list))

    for (ndata, (obs_torque_matrix, pred_torque_matrix)) in
        enumerate(zip(observed_torque_list, predicted_torque_list))
        println(f, "# data index: $ndata")
        for (iatom, (obs_torque, pred_torque)) in
            enumerate(zip(eachcol(obs_torque_matrix), eachcol(pred_torque_matrix)))
            # obs_torque and pred_torque are length-3 vectors (x, y, z)
            str = @sprintf(
                " %*d %*s  % 15.10e   % 15.10e   % 15.10e    % 15.10e   % 15.10e   % 15.10e\n",
                idx_width,
                iatom,
                element_width,
                element_string_list[iatom],
                obs_torque[1],
                obs_torque[2],
                obs_torque[3],
                pred_torque[1],
                pred_torque[2],
                pred_torque[3],
            )
            write(f, str)
        end
    end
	end
end

"""
	get the reference energy
"""
function get_j0(sc::SpinCluster)::Float64
	return sc.optimize.reference_energy
end

"""
	get the spin-cluster coefficients
"""
function get_jphi(sc::SpinCluster)::Vector{Float64}
	return sc.optimize.SCE
end

"""
	get the reference energy and spin-cluster coefficients
"""
function get_j0_jphi(sc::SpinCluster)::Tuple{Float64, Vector{Float64}}
	return sc.optimize.reference_energy, sc.optimize.SCE
end

function print_header()
	println(
		"""
		+-----------------------------------+
				  Magesty v$(Version.version_string())      
		+-----------------------------------+

		Julia version: $(VERSION)

		Number of threads: $(Threads.nthreads())

		Job started at $(now())

		""")
end

end # module Magesty
