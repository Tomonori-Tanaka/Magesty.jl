using TOML
using Printf
using Test
using StaticArrays
using Random
using Magesty


input = TOML.parse(open(joinpath(@__DIR__, "input.toml"), "r"))
input["regression"]["datafile"] = joinpath(@__DIR__, "EMBSET.dat")

@testset "Dimer Tests" begin
	system = System(input, verbosity = false)

	spinconfig_list = Magesty.SpinConfigs.SpinConfig[]
	fm_energy = -1.0
	afm_energy = 1.0
	fm_spin_directions = [                                    0.0 0.0;
		0.0 0.0;
		1.0 1.0]
	afm_spin_directions = [                                      0.0 0.0;
		0.0 0.0;
		-1.0 1.0]
	# magfield is dummy
	fm_local_magfield =
		[                                                                        0.0 0.0;
			0.0 0.0;
			0.0 0.0]
	afm_local_magfield =
		[                                                                            0.0 0.0;
			0.0 0.0;
			0.0 0.0]
	push!(
		spinconfig_list,
		Magesty.SpinConfigs.SpinConfig(
			fm_energy,
			[1.0, 1.0], # magmom_size_list
			fm_spin_directions,
			fm_local_magfield,
		),
	)
	push!(
		spinconfig_list,
		Magesty.SpinConfigs.SpinConfig(
			afm_energy,
			[1.0, 1.0],
			afm_spin_directions,
			afm_local_magfield,
		),
	)


	# regression with w = 0
	optimize = Magesty.Optimize.Optimizer(
		system.structure,
		system.symmetry,
		system.basisset,
		0.0,
		0.0,
		0.0,
		spinconfig_list,
		verbosity = false,
	)
	spincluster = Magesty.SpinCluster(
		system.structure,
		system.symmetry,
		system.cluster,
		system.basisset,
		optimize,
	)
	path = joinpath(@__DIR__, "dimer.xml")
	Magesty.write_xml(spincluster, path)


	@test Magesty.MySphericalHarmonics.Zₗₘ(1, -1, [0.0, 0.0, 1.0]) ≈ 0.0 atol = 1e-6
	@test Magesty.MySphericalHarmonics.Zₗₘ(1, 0, [0.0, 0.0, 1.0]) ≈ √(3 / 4π) atol = 1e-6
	@test Magesty.MySphericalHarmonics.Zₗₘ(1, 1, [0.0, 0.0, 1.0]) ≈ 0.0 atol = 1e-6

	@test spincluster.optimize.reference_energy ≈ 0.0 atol = 1e-6
	@test length(spincluster.optimize.SCE) == 1
	@test spincluster.optimize.SCE[1]*√(3) ≈ -1.0 atol = 1e-6


end
