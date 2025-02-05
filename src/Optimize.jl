"""
	Optimize.jl

This module contains functions for optimizing the SCE coefficients.
"""
module Optimize

using ..Systems 
using ..Symmetries
using ..Clusters
using ..BasisSets
using ..SpinConfigs

export SCEOptimizer

struct SCEOptimizer
	SCE::Vector{Float64}
end

function SCEOptimizer(
	system::System,
	symmetry::Symmetry,
	cluster::Cluster,
	basisset::BasisSet,
	j_zero_thr::Real,
	weight::Real,
    datafile::AbstractString,
)
    # read datafile
    spinconfig_list::Vector{SpinConfig} = read_embset(datafile, system.supercell.num_atoms)

    tmp = [1.0]
    return SCEOptimizer(tmp)
end

function ols_energy(salc_list::AbstractVector{SALC}, spinconfig_list::AbstractVector{SpinConfig})
end

end
