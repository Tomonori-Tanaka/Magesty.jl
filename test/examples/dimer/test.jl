using TOML
using Printf
using Test
using StaticArrays
using Random
using Magesty
include("../../../tools/convert2tensor.jl")
using .ExchangeTensor


input = TOML.parse(open(joinpath(@__DIR__, "input.toml"), "r"))
input["regression"]["datafile"] = joinpath(@__DIR__, "EMBSET.dat")

@testset "Dimer Tests" begin
	system = System(input, verbosity = false)

	spinconfig_list = Magesty.SpinConfigs.SpinConfig[]
	fm_energy = -1.0
	afm_energy = 1.0
	fm_spin_directions =
		[                                                               0.0 0.0;
			0.0 0.0;
			1.0 1.0]
	afm_spin_directions =
		[                                                                   0.0 0.0;
			0.0 0.0;
			-1.0 1.0]
	# magfield is dummy
	fm_local_magfield =
		[                                                         0.0 0.0;
			0.0 0.0;
			0.0 0.0]
	afm_local_magfield =
		[                                                             0.0 0.0;
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

	@testset "DMI energy" begin
		# H = D⋅(e_1 × e_2), D = (0, 0, -1)
		input["regression"]["weight"] = 0.0
		input["symmetry"]["isotropy"] = false
		system = System(input, verbosity = false)
		spinconfig_list = Magesty.SpinConfigs.SpinConfig[]
		spin_direction1 =
			[1.0 0.0;
					  0.0 1.0;
					  0.0 0.0]
		local_magfield1 =
			[0.0 0.0;
					  0.0 0.0;
					  0.0 0.0]
		energy1 = -1.0
		spin_direction2 =
			[0.0 1.0;
					  1.0 0.0;
					  0.0 0.0]
		local_magfield2 =
			[0.0 0.0;
					  0.0 0.0;
					  0.0 0.0]
		energy2 = 1.0
		push!(
			spinconfig_list,
			Magesty.SpinConfigs.SpinConfig(energy1, [1.0, 1.0], spin_direction1, local_magfield1),
			Magesty.SpinConfigs.SpinConfig(energy2, [1.0, 1.0], spin_direction2, local_magfield2),
		)

		# adjast basis set
		basisset = system.basisset
		coeff_tensor = basisset.angular_momentum_couplings[2].coeff_tensor
		coefficient = [0.0, 1.0, 0.0]
		coupled_basis_with_coefficient = Magesty.Basis.CoupledBasis_with_coefficient(
			[1, 1],
			1,
			Int[],
			[1, 2],
			coeff_tensor,
			coefficient,
			1,
		)
		salc_list = Vector{Vector{Magesty.Basis.CoupledBasis_with_coefficient}}()
		push!(salc_list, [coupled_basis_with_coefficient])
		coupled_basislist = basisset.coupled_basislist
		angular_momentum_couplings = basisset.angular_momentum_couplings
		basisset = Magesty.BasisSet(coupled_basislist, salc_list, angular_momentum_couplings)

		optimize = Magesty.Optimize.Optimizer(
			system.structure,
			system.symmetry,
			basisset,
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
			basisset,
			optimize,
		)
		@test spincluster.optimize.reference_energy ≈ 0.0 atol = 1e-6
		@test length(spincluster.optimize.SCE) == 1
		@test spincluster.optimize.SCE[1]*3/√(2) ≈ -1.0 atol = 1e-6
		Magesty.write_xml(spincluster, joinpath(@__DIR__, "dimer_dmi.xml"))
		exchange_tensor = convert2tensor(
			joinpath(@__DIR__, "dimer_dmi.xml"),
			[1, 2],
		)
		display(exchange_tensor)
		display(exchange_tensor.isotropic_jij)
		display(exchange_tensor.dm_vector)
		display(exchange_tensor.gamma)

	end

end
