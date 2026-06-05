using TOML
using Printf
using Test
using StaticArrays
using Random
using Magesty
include("../../../tools/convert2tensor.jl")
using .ExchangeTensor


input = TOML.parse(open(joinpath(@__DIR__, "input.toml"), "r"))

@testset "Dimer Tests" begin
	basis = SCEBasis(input; verbosity = false)

	fm = Magesty.SpinConfigs.SpinConfig(
		-1.0,
		[1.0, 1.0],
		[0.0 0.0;
			0.0 0.0;
			1.0 1.0],
		[0.0 0.0;
			0.0 0.0;
			0.0 0.0],
	)
	afm = Magesty.SpinConfigs.SpinConfig(
		1.0,
		[1.0, 1.0],
		[0.0 0.0;
			0.0 0.0;
			-1.0 1.0],
		[0.0 0.0;
			0.0 0.0;
			0.0 0.0],
	)
	spinconfig_list = [fm, afm]

	dataset = SCEDataset(basis, spinconfig_list)
	fitted = fit(SCEFit, dataset, OLS(); torque_weight = 0.0, verbosity = false)
	model = SCEModel(fitted)
	Magesty.save(model, joinpath(@__DIR__, "dimer.xml"))

	@test Magesty.TesseralHarmonics.Zₗₘ(1, -1, [0.0, 0.0, 1.0]) ≈ 0.0 atol = 1e-6
	@test Magesty.TesseralHarmonics.Zₗₘ(1, 0, [0.0, 0.0, 1.0]) ≈ √(3 / 4π) atol = 1e-6
	@test Magesty.TesseralHarmonics.Zₗₘ(1, 1, [0.0, 0.0, 1.0]) ≈ 0.0 atol = 1e-6

	@test intercept(fitted) ≈ 0.0 atol = 1e-6
	@test length(coef(fitted)) == 1
	@test coef(fitted)[1] * √3 ≈ -1.0 atol = 1e-6

	@testset "DMI energy" begin
		# H = D⋅(e_1 × e_2), D = (0, 0, -1)
		input["symmetry"]["isotropy"] = false
		dmi_basis = SCEBasis(input; verbosity = false)
		spin_direction1 =
			[1.0 0.0;
				0.0 1.0;
				0.0 0.0]
		spin_direction2 =
			[0.0 1.0;
				1.0 0.0;
				0.0 0.0]
		local_magfield =
			[0.0 0.0;
				0.0 0.0;
				0.0 0.0]
		sc1 = Magesty.SpinConfigs.SpinConfig(
			-1.0,
			[1.0, 1.0],
			spin_direction1,
			local_magfield,
		)
		sc2 = Magesty.SpinConfigs.SpinConfig(
			1.0,
			[1.0, 1.0],
			spin_direction2,
			local_magfield,
		)
		dmi_configs = [sc1, sc2]

		# Replace the salc_list with a single hand-crafted Lf=1 SALC whose
		# coefficient picks out only the z-component of the antisymmetric
		# l=1 × l=1 coupling — i.e. D = (0, 0, -1) in our convention.
		original_salcbasis = dmi_basis.salcbasis
		coeff_tensor = original_salcbasis.angular_momentum_couplings[2].coeff_tensor
		clusters = Magesty.CoupledBases.enumerate_orbit_clusters(
			[1, 2],
			dmi_basis.symmetry.map_sym,
			dmi_basis.symmetry.symnum_translation,
		)
		cbc = Magesty.CoupledBases.CoupledBasis_with_coefficient(
			[1, 1],
			1,
			Int[],
			[1, 2],
			coeff_tensor,
			[0.0, 1.0, 0.0],
			1,
			clusters,
		)
		salc_list = [[cbc]]
		modified_salcbasis = Magesty.SALCBasis(
			original_salcbasis.coupled_basislist,
			salc_list,
			original_salcbasis.angular_momentum_couplings,
		)
		dmi_basis = SCEBasis(
			dmi_basis.structure,
			dmi_basis.symmetry,
			modified_salcbasis,
			false,
		)

		dmi_dataset = SCEDataset(dmi_basis, dmi_configs)
		dmi_fitted = fit(SCEFit, dmi_dataset, OLS(); torque_weight = 0.0, verbosity = false)
		dmi_model = SCEModel(dmi_fitted)

		@test intercept(dmi_fitted) ≈ 0.0 atol = 1e-6
		@test length(coef(dmi_fitted)) == 1
		@test coef(dmi_fitted)[1] * 3 / √2 ≈ -1.0 atol = 1e-6

		dmi_xml = joinpath(@__DIR__, "dimer_dmi.xml")
		Magesty.save(dmi_model, dmi_xml)
		exchange_tensor = convert2tensor(dmi_xml, [1, 2])
		#display(exchange_tensor)
		#display(exchange_tensor.isotropic_jij)
		#display(exchange_tensor.dm_vector)
		#display(exchange_tensor.gamma)
	end
end
