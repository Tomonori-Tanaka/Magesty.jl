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

include("System.jl")
include("Symmetry.jl")
include("Cluster.jl")
include("BasisSet.jl")
include("Optimize.jl")

using .Systems
using .Symmetries
using .Clusters
using .BasisSets
using .Optimize

export SpinCluster

struct SpinCluster
	config::Parser
	system::System
	symmetry::Symmetry
	cluster::Cluster
	basisset::BasisSet
	optimize::Union{SCEOptimizer, Nothing}
end

function SpinCluster(input_dict::Dict{<:AbstractString, <:Any})
	parser = Parser(input_dict)
	system::System = set_system(parser)
	symmetry::Symmetry = set_symmetry(parser, system)
	cluster::Cluster = set_cluster(parser, system, symmetry)
	basisset::BasisSet = set_basisset(parser, system, symmetry, cluster)
	optimize = if parser.mode == "optimize"
		set_optimize(parser, system, symmetry, basisset)
	else
		nothing
	end

	return SpinCluster(parser, system, symmetry, cluster, basisset, optimize)
end

function SpinCluster(toml_file::AbstractString)
	open(toml_file) do io
		toml = read(io, String)
		config = TOML.parse(toml)
		return SpinCluster(config)
	end
end

function write_xml(sc::SpinCluster)
	system = sc.system

end

function print_info(sc::SpinCluster)
	println(
		"""
		+-----------------------------------+
		|          Magesty v0.1.0           |
		+-----------------------------------+

		""",
	)

	Systems.print_info(sc.system)
	Symmetries.print_info(sc.symmetry)
	BasisSets.print_info(sc.basisset)
	Optimize.print_info(sc.optimize)
end

function write_energy_lists(sc::SpinCluster, filename::AbstractString = "energy_lists.txt")
	Optimize.write_energy_lists(sc.optimize, filename)
end

function write_magfield_vertical_list(sc::SpinCluster, filename::AbstractString = "magfield_vertical_list.txt")
	Optimize.write_magfield_vertical_list(sc.optimize, filename)
end

end
