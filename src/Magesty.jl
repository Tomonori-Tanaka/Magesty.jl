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
	optimizer = sc.optimize
	# Input validation
	if isempty(optimizer.spinconfig_list)
		@warn "No spin configurations found in optimizer"
		return
	end

	# Prepare data
	spinconfig_list = optimizer.spinconfig_list
	observed_energy_list = [spinconfig.energy for spinconfig in spinconfig_list]
	predicted_energy_list = optimizer.predicted_energy_list

	# Check array lengths
	if length(observed_energy_list) != length(predicted_energy_list)
		error("Length mismatch between observed and predicted energy lists")
	end

	# Format settings
	digits_index = length(string(length(observed_energy_list)))
	header = "# Index:" * " "^(digits_index-1) * "Observed_Energy" * " "^4 * "Predicted_Energy"

	# Write to file
	try
		open(filename, "w") do f
			# Write header
			println(f, header)

			# Write data
			for (i, (obs_energy, pred_energy)) in enumerate(zip(observed_energy_list, predicted_energy_list))
				str = @sprintf("%*d    %15.10f    %15.10f\n", digits_index, i, obs_energy, pred_energy)
				write(f, str)
			end
		end
	catch e
		@error "Failed to write energy lists to file" exception=(e, catch_backtrace())
		rethrow(e)
	end
end

function write_torque_list(sc::SpinCluster, filename::AbstractString = "torque_list.txt")
	optimizer = sc.optimize
end

function write_sce(sc::SpinCluster)
end

end
