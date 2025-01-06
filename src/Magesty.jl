module Magesty

using TOML

include("common/SortedContainer.jl")

include("types/AtomCell.jl")
include("types/AtomicIndices.jl")
include("types/UnitaryMatrixCl.jl")

include("System.jl")
include("Symmetry.jl")
include("Cluster.jl")

using .Systems
using .Symmetries
using .Clusters

struct SpinCluster
	system::System
	symmetry::Symmetry
	cluster::Cluster
end

function SpinCluster(input_dict::Dict{AbstractString, Any})
    
end

function SpinCluster(toml_file::AbstractString)
	open(toml_file) do io
		toml = read(io, String)
		config = TOML.parse(toml)
		return SpinCluster(config)
	end
end

end
