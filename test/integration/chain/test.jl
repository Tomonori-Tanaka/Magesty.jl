using TOML
using Printf
using Test
using StaticArrays
using Random
using Magesty


input = TOML.parse(open(joinpath(@__DIR__, "input.toml"), "r"))

@testset "Chain Tests" begin
	# E = ∑_{i < j} Jij (e_i ⋅ e_j) , Jij_1NN = -1.0, Jij_2NN = -2.0, Jij_3NN = -3.0
	basis = SCEBasis(input; verbosity = false)

	fm = Magesty.SpinConfigs.SpinConfig(
		-2.0,
		[1.0, 1.0],
		[0.0 0.0;
			0.0 0.0;
			1.0 1.0],
		[0.0 0.0;
			0.0 0.0;
			0.0 0.0],
	)
	afm = Magesty.SpinConfigs.SpinConfig(
		2.0,
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
	# Energy-only fit (torque is the dummy zero field above; weight=0 reproduces the legacy w=0 path).
	fitted = fit(SCEFit, dataset, OLS(); torque_weight = 0.0)
	model = SCEModel(fitted)
	Magesty.save(model, joinpath(@__DIR__, "chain.xml"))

	@test Magesty.TesseralHarmonics.Zₗₘ(1, 0, [0.0, 0.0, 1.0]) ≈ √(3 / 4π) atol = 1e-6

	# The FM / AFM design rows are symmetric (flipping atom 2 flips the l=1 contribution),
	# so the fit recovers j0 = 0 and a single SCE coefficient whose Heisenberg projection
	# `coef * √3` equals the underlying Jij = -1.
	@test intercept(fitted) ≈ 0.0 atol = 1e-6
	@test length(coef(fitted)) == 1
	@test coef(fitted)[1] * √3 ≈ -1.0 atol = 1e-6
end
