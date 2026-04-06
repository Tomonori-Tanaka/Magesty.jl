using TOML
using Printf
using Test
using StaticArrays
using Random
using Magesty


input = TOML.parse(open(joinpath(@__DIR__, "input.toml"), "r"))
# input["regression"]["datafile"] = joinpath(@__DIR__, "EMBSET.dat")

@testset "Square Lattice Tests" begin
	system = System(input, verbosity = false)
	Magesty.write_xml(system, joinpath(@__DIR__, "system.xml"))

	@testset "Hypothetical SCE model" begin
		# Assign Jij = -1.0 eV (ferromagnetic interaction) to all pairs.
		# Jij is the Heisenberg exchange interaction.
	end

end
