include("/Users/tomorin/Packages/Magesty.jl/src/Magesty.jl")
using .Magesty
using TOML
using Printf
using JLD2
using Plots
using LinearAlgebra
using Statistics
using Logging


input = TOML.parse(open("input.toml", "r"))
system = System(input, verbosity = true)

@save "system.jld2" system
Magesty.write_xml(system)
@load "system.jld2" system

sclus = SpinCluster(system, input, verbosity = true)
#Magesty.write_energies(sclus)
#Magesty.write_torques(sclus)
Magesty.write_xml(sclus, "jphi.xml")
