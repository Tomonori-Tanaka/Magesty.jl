"""
    compare_energy.jl

This script compares the energy of the spin configurations from two different EMBSET files.
"""

using LinearAlgebra
include("../src/SpinConfig.jl")
using .SpinConfigs

using ArgParse


function compare_energy(embset1::AbstractString, embset2::AbstractString, denominator::Real, num_atoms::Integer)
    spin_configs1::Vector{SpinConfig} = read_embset(embset1, num_atoms)
    spin_configs2::Vector{SpinConfig} = read_embset(embset2, num_atoms)
    
    # shorter length is the number of datasets
    num_datasets = length(spin_configs1) < length(spin_configs2) ? length(spin_configs1) : length(spin_configs2)

    println("energy1 - energy2")
    for i in 1:num_datasets
        energy1 = spin_configs1[i].energy
        energy2 = spin_configs2[i].energy
        println("$(i): $((energy1 - energy2) / denominator)")
    end
end

s = ArgParseSettings()
@add_arg_table s begin
    "--embset1"
    help = "The path to the first EMBSET file"
    required = true

    "--embset2"
    help = "The path to the second EMBSET file"
    required = true

    "--denominator", "-d"
    help = "The denominator of the energy"
    required = true
    arg_type = Float64
    default = 1.0

    "--num_atoms", "-n"
    help = "The number of atoms"
    required = true
    arg_type = Int
    
end

parsed_args = parse_args(s)

compare_energy(parsed_args["embset1"], parsed_args["embset2"], parsed_args["denominator"], parsed_args["num_atoms"])
