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
include("common/CountingContainer.jl")
include("types/AtomCells.jl")
include("SpinConfigs.jl")
using .SpinConfigs

include("utils/SphericalHarmonicsTransforms.jl")
include("utils/AngularMomentumCoupling.jl")
include("types/Basis.jl")
using .AngularMomentumCoupling
using .Basis

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

include("utils/xml_io.jl")
using .XMLIO
include("utils/EnergyTorque.jl")
using .EnergyTorque

export System, SpinCluster, VERSION, install_tools
export SCEModel, fit_sce_model, predict_energy, AbstractEstimator, OLS, Ridge
export build_sce_basis, build_sce_basis_from_xml
export write_xml

# Re-export read_embset from SpinConfigs for user convenience
const read_embset = SpinConfigs.read_embset
export read_embset

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

# Shared skeleton for all input-driven constructors below.
# Returns the (structure, symmetry, cluster) triplet; callers append
# the BasisSet — either computed via `BasisSet(...)` or loaded from XML.
function _build_structure_skeleton(
	config::Config4System;
	verbosity::Bool = true,
)
	structure::Structure = Structure(config, verbosity = verbosity)
	symmetry::Symmetry = Symmetry(structure, config, verbosity = verbosity)
	cluster::Cluster = Cluster(structure, symmetry, config, verbosity = verbosity)
	return structure, symmetry, cluster
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
	structure, symmetry, cluster = _build_structure_skeleton(config; verbosity = verbosity)
	basisset::BasisSet = BasisSet(structure, symmetry, cluster, config, verbosity = verbosity)
	return System(structure, symmetry, cluster, basisset)
end

function System(toml_file::AbstractString; verbosity::Bool = true)
	try
		open(toml_file) do io
			toml = read(io, String)
			input_dict = TOML.parse(toml)
			return System(input_dict; verbosity = verbosity)
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
	build_sce_basis(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool=false) -> System

Build a `System` (structure + symmetry + cluster + basis set) from a parsed input
dictionary. This is the headless alternative to `System(input_dict)` and skips the
header banner by default; the SALC basis is computed from scratch (use
`build_sce_basis_from_xml` to load a precomputed basis instead).

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing the parsed
  TOML input parameters (typically the result of `TOML.parsefile`).
- `verbosity::Bool=false`: Whether to print detailed progress information.

# Returns
- `System`: A `System` instance containing structure, symmetry, cluster, and basis set.

# Throws
- `ErrorException`: If required parameters are missing or invalid.

# Examples
```julia
using TOML
input_dict = TOML.parsefile("input.toml")
system = build_sce_basis(input_dict)
```
"""
function build_sce_basis(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = false)::System
	config::Config4System = Config4System(input_dict)
	structure, symmetry, cluster = _build_structure_skeleton(config; verbosity = verbosity)
	basisset::BasisSet = BasisSet(structure, symmetry, cluster, config, verbosity = verbosity)
	return System(structure, symmetry, cluster, basisset)
end

"""
	build_sce_basis_from_xml(input_dict::Dict{<:AbstractString, <:Any}, xml_file::AbstractString; verbosity::Bool = true)::System

Build System from input.toml dictionary and XML file. This function constructs structure, symmetry, and cluster
from the input dictionary, but loads the basis set from the XML file to avoid expensive SALC computation.

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing input parameters (parsed from input.toml)
- `xml_file::AbstractString`: Path to XML file containing basis set information
- `verbosity::Bool=true`: Whether to print detailed information during initialization

# Returns
- `System`: A new `System` instance with basis set loaded from XML file

# Throws
- `ErrorException` if required parameters are missing or XML file format is invalid

# Examples
```julia
using TOML
input_dict = TOML.parsefile("input.toml")
system = build_sce_basis_from_xml(input_dict, "scecoeffs.xml")
```
"""
function build_sce_basis_from_xml(
	input_dict::Dict{<:AbstractString, <:Any},
	xml_file::AbstractString;
	verbosity::Bool = true,
)::System
	if verbosity
		print_header()
	end

	config::Config4System = Config4System(input_dict)
	structure, symmetry, cluster = _build_structure_skeleton(config; verbosity = verbosity)

	# Load basis set from XML file instead of computing it
	if verbosity
		println("Loading basis set from XML file: $xml_file")
	end
	basisset::BasisSet = XMLIO.read_basisset_from_xml(xml_file)
	if verbosity
		println("Successfully loaded basis set from XML file")
	end

	return System(structure, symmetry, cluster, basisset)
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
	structure, symmetry, cluster = _build_structure_skeleton(config_system; verbosity = verbosity)
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
	SpinCluster(system::System,
	            input_dict::AbstractDict{<:AbstractString, <:Any},
	            spinconfig_list::AbstractVector{SpinConfig};
	            verbosity::Bool=true)

Create a `SpinCluster` from an existing `System`, an input dictionary, and a
caller-supplied list of `SpinConfig` objects. Use this overload when the spin
configurations come from somewhere other than the EMBSET path declared in
`input_dict` (for example, configurations generated programmatically or pulled
from a non-standard format).

# Arguments
- `system::System`: An existing `System` instance.
- `input_dict::AbstractDict{<:AbstractString, <:Any}`: Dictionary with
  optimization parameters (`[regression]` section).
- `spinconfig_list::AbstractVector{SpinConfig}`: Training spin configurations
  used in place of those that would be loaded from disk.
- `verbosity::Bool=true`: Whether to print detailed information during fitting.

# Returns
- `SpinCluster`: A fitted `SpinCluster` instance.

# Throws
- `ErrorException`: If required `[regression]` parameters are missing or invalid.

# Examples
```julia
system = System("input.toml")
configs = [SpinConfig(...), SpinConfig(...)]
sc = SpinCluster(system, input_dict, configs)
```
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
	return EnergyTorque.calc_energy(
		spincluster.basisset.salc_list,
		spin_config,
		spincluster.symmetry,
		spincluster.optimize,
	)
end

"""
	calc_torque(sc::SpinCluster, spin_config::AbstractMatrix{<:Real})

Calculate the torque (local magnetic field) for each atom in a spin configuration using the spin cluster expansion.

# Arguments
- `sc::SpinCluster`: A `SpinCluster` instance containing structure, symmetry, basis set, and optimization information.
- `spin_config::AbstractMatrix{<:Real}`: A 3xN matrix representing the spin configuration, where N is the number of atoms in the supercell.

# Returns
- `Matrix{Float64}`: A 3xN matrix representing the torque for each atom, where each column corresponds to the torque vector (in eV) for one atom.

# Throws
- `ArgumentError`: If the number of columns in `spin_config` does not match the number of atoms in the supercell.

# Example
```julia
spin_config = rand(3, sc.structure.supercell.num_atoms) # Random spin configuration
torque = calc_torque(sc, spin_config)
```
"""
function calc_torque(spincluster::SpinCluster, spin_config::AbstractMatrix{<:Real})::Matrix{Float64}
	if spincluster.structure.supercell.num_atoms != size(spin_config, 2)
		num_atoms = spincluster.structure.supercell.num_atoms
		throw(
			ArgumentError(
				"spin_config must be 3xN matrix where N is the number of atoms in the supercell. $num_atoms",
			),
		)
	end
	return EnergyTorque.calc_torque(
		spincluster.basisset.salc_list,
		spin_config,
		spincluster.symmetry,
		spincluster.optimize,
	)
end

"""
	write_xml(sc::SpinCluster, filename::AbstractString="jphi.xml"; write_jphi::Bool=true)

Write the structure, symmetry, basis set, and fitted SCE coefficients
held by a `SpinCluster` to an XML file in SCE format.

# Arguments
- `sc::SpinCluster`: Spin cluster (structure + symmetry + basis set + optimizer)
- `filename::AbstractString="jphi.xml"`: Output XML file name
- `write_jphi::Bool=true`: Whether to embed the J_ij parameters
  (set to `false` to write only the structure/symmetry/basis section)

# Examples
```julia
write_xml(spin_cluster)
write_xml(spin_cluster, "output.xml", write_jphi=false)
```

To export a `System` (no fitted coefficients), use the
`write_xml(::System, ...)` overload instead. To export from
standalone components (`structure`, `symmetry`, `basis_set`,
`optimizer`), wrap them with `SpinCluster(structure, symmetry,
cluster, basis_set, optimizer)` first.
"""
function write_xml(
	sc::SpinCluster,
	filename::AbstractString = "jphi.xml";
	write_jphi::Bool = true,
)
	XMLIO.write_xml(
		sc.structure,
		sc.symmetry,
		sc.basisset,
		sc.optimize,
		filename;
		write_jphi = write_jphi,
	)
end

"""
	write_xml(system::System, filename::AbstractString="jphi.xml")

Write System information to an XML file. This saves the structure, symmetry, and basis set
information without optimization results (JPhi).

# Arguments
- `system::System`: The System to write
- `filename::AbstractString`: Output XML file name (default: "jphi.xml")

# Examples
```julia
system = build_sce_basis(input_dict)
write_xml(system, "system.xml")
```
"""
function write_xml(
	system::System,
	filename::AbstractString = "jphi.xml",
)
	XMLIO.write_xml(system.structure, system.symmetry, system.basisset, filename)
end

"""
	write_energies(sc::SpinCluster, filename::AbstractString="energy_list.txt")

Write the observed (DFT) and predicted (SCE) energies to a text file.

# Arguments
- `sc::SpinCluster`: A `SpinCluster` instance containing optimization results.
- `filename::AbstractString="energy_list.txt"`: Output file name.

# Output Format
The file contains:
- Header line: `# data index,    DFT_Energy,    SCE_Energy`
- Data lines: index, observed energy (eV), predicted energy (eV)

# Example
```julia
write_energies(sc, "my_energies.txt")
```
"""
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

"""
	write_torques(sc::SpinCluster, filename::AbstractString="torque_list.txt")

Write the observed (DFT) and predicted (SCE) torques for each atom to a text file.

# Arguments
- `sc::SpinCluster`: A `SpinCluster` instance containing optimization results.
- `filename::AbstractString="torque_list.txt"`: Output file name.

# Output Format
The file contains:
- Header line: `# atom index,    element,   DFT_torque_x,    DFT_torque_y,    DFT_torque_z,    SCE_torque_x,    SCE_torque_y,    SCE_torque_z`
- Data lines: atom index, element symbol, observed torque components (eV), predicted torque components (eV)
- Data is grouped by configuration index

# Example
```julia
write_torques(sc, "my_torques.txt")
```
"""
function write_torques(
	sc::SpinCluster,
	filename::AbstractString = "torque_list.txt",
)
	predicted_torque_list::Vector{Matrix{Float64}} = sc.optimize.predicted_torque_list
	observed_torque_list::Vector{Matrix{Float64}} =
		[spinconfig.torques for spinconfig in sc.optimize.spinconfig_list]

	# Write to file
	open(filename, "w") do f
		# Write header
		println(
			f,
			"# atom index,    element,   DFT_torque_x,    DFT_torque_y,    DFT_torque_z,    SCE_torque_x,    SCE_torque_y,    SCE_torque_z\n# unit of torque is eV",
		)

		# Write data
		idx_width = ndigits(length(sc.optimize.spinconfig_list))
		element_string_list =
			[sc.structure.kd_name[elm_idx] for elm_idx in sc.structure.supercell.kd_int_list]
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
	get_j0(sc::SpinCluster) -> Float64

Return the fitted reference energy (bias term `j0`) from a `SpinCluster`.

# Arguments
- `sc::SpinCluster`: Spin cluster with a fitted optimizer.

# Returns
- `Float64`: Reference energy in eV (the constant offset of the SCE expansion).

# Examples
```julia
sc = SpinCluster("input.toml")
j0 = get_j0(sc)
```
"""
function get_j0(sc::SpinCluster)::Float64
	return sc.optimize.reference_energy
end

"""
	get_jphi(sc::SpinCluster) -> Vector{Float64}

Return the fitted SCE coefficients (`jphi`) from a `SpinCluster`.

# Arguments
- `sc::SpinCluster`: Spin cluster with a fitted optimizer.

# Returns
- `Vector{Float64}`: SCE coefficients, one entry per SALC basis function.

# Examples
```julia
sc = SpinCluster("input.toml")
jphi = get_jphi(sc)
```
"""
function get_jphi(sc::SpinCluster)::Vector{Float64}
	return sc.optimize.SCE
end

"""
	get_j0_jphi(sc::SpinCluster) -> Tuple{Float64, Vector{Float64}}

Return the reference energy and SCE coefficients as a tuple.

# Arguments
- `sc::SpinCluster`: Spin cluster with a fitted optimizer.

# Returns
- `Tuple{Float64, Vector{Float64}}`: `(j0, jphi)` where `j0` is the bias term in eV
  and `jphi` is the vector of SCE coefficients.

# Examples
```julia
sc = SpinCluster("input.toml")
j0, jphi = get_j0_jphi(sc)
```
"""
function get_j0_jphi(sc::SpinCluster)::Tuple{Float64, Vector{Float64}}
	return sc.optimize.reference_energy, sc.optimize.SCE
end

"""
    install_tools(; bindir::AbstractString=joinpath(homedir(), ".julia", "bin"))

Install CLI wrapper scripts for Magesty tools into `bindir` (default: `~/.julia/bin`).
Each wrapper calls `julia /path/to/script.jl "\$@"`, so the correct package version
is always used regardless of where the package is installed.

After running this once, add `~/.julia/bin` to your PATH:
    export PATH="\$HOME/.julia/bin:\$PATH"

# Arguments
- `bindir::AbstractString=joinpath(homedir(), ".julia", "bin")`: Target directory
  for the wrapper scripts. The directory is created if it does not exist.

# Returns
- `Nothing`: Writes wrapper scripts as a side effect and prints their paths.

# Throws
- `ErrorException`: If the Magesty package directory cannot be determined.

# Available commands after installation
- `vasp2extxyz`           — convert a single VASP output to extxyz
- `vasp2extxyz_recursive` — recursively convert VASP outputs under a directory

# Examples
```julia
using Magesty
Magesty.install_tools()
# or to install into a custom location
Magesty.install_tools(bindir="/usr/local/bin")
```
"""
function install_tools(; bindir::AbstractString = joinpath(homedir(), ".julia", "bin"))
    pkg_dir = pkgdir(Magesty)
    if isnothing(pkg_dir)
        error("Cannot determine Magesty package directory")
    end

    mkpath(bindir)

    tools = [
        "vasp2extxyz"           => joinpath("tools", "vasp", "vasp2extxyz.jl"),
        "vasp2extxyz_recursive" => joinpath("tools", "vasp", "vasp2extxyz_recursive.jl"),
    ]

    global_env = "@v$(Sys.VERSION.major).$(Sys.VERSION.minor)"
    for (name, rel_path) in tools
        script = joinpath(pkg_dir, rel_path)
        wrapper = joinpath(bindir, name)
        write(wrapper, "#!/bin/sh\nexec julia --project=$global_env \"$script\" \"\$@\"\n")
        chmod(wrapper, 0o755)
        println("Installed: $wrapper")
    end

    println("\nMake sure ~/.julia/bin is in your PATH:")
    println("  export PATH=\"\$HOME/.julia/bin:\$PATH\"")
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
