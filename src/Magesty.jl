module Magesty

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
	optimize::SCEOptimizer
end

function SpinCluster(input_dict::Dict{<:AbstractString, <:Any})
	parser = Parser(input_dict)
	system::System = set_system(parser)
	symmetry::Symmetry = set_symmetry(parser, system)
	cluster::Cluster = set_cluster(parser, system, symmetry)
	basisset::BasisSet = set_basisset(parser, system, symmetry, cluster)
	optimize::SCEOptimizer = set_optimize(parser, system, symmetry, cluster, basisset)

	return SpinCluster(parser, system, symmetry, cluster, basisset, optimize)
end

function SpinCluster(toml_file::AbstractString)
	open(toml_file) do io
		toml = read(io, String)
		config = TOML.parse(toml)
		return SpinCluster(config)
	end
end

end
