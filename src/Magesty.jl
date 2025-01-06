module Magesty

using TOML

include("common/SortedContainer.jl")

include("types/AtomCell.jl")
include("types/AtomicIndices.jl")
include("types/UnitaryMatrixCl.jl")

include("utils/InputParser.jl")
include("utils/InputSetter.jl")
using .InputParser

include("System.jl")
include("Symmetry.jl")
include("Cluster.jl")
using .Systems
using .Symmetries
using .Clusters

export SpinCluster

struct SpinCluster
	config::Parser
	system::System
	symmetry::Symmetry
	cluster::Cluster
end

function SpinCluster(input_dict::Dict{<:AbstractString, <:Any})
	parser = Parser(input_dict)
	system::System = set_system(parser)
	symmetry::Symmetry = set_symmetry(parser, system)
	cluster::Cluster = set_cluster(parser, system, symmetry)

	return SpinCluster(parser, system, symmetry, cluster)
end

function SpinCluster(toml_file::AbstractString)
	open(toml_file) do io
		toml = read(io, String)
		config = TOML.parse(toml)
		return SpinCluster(config)
	end
end

end
