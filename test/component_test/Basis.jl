@testset "Basis" begin
	include("../../src/utils/AngularMomentumCoupling.jl")
	using .AngularMomentumCoupling
	include("../../src/types/Basis.jl")
	using .Basis
	using Test

	@testset "LinearCombo" begin
		@testset "Constructor" begin
			atoms = [10, 20]
			ls = [1, 1]
			lc_list = tesseral_linear_combos_from_tesseral_bases(ls, atoms)
			@test length(lc_list) == 9
			# Lf = 0
			lc_lf0 = lc_list[1]
			@test lc_lf0.ls == (1, 1)
			@test lc_lf0.Lf == 0
			@test lc_lf0.Lseq == Int[]
			@test lc_lf0.atoms == [10, 20]
			@test lc_lf0.coeff_list == [1.0]
			@test isapprox(
				lc_lf0.coeff_tensor,
				[1/sqrt(3) 0.0 0.0; 0.0 1/sqrt(3) 0.0; 0.0 0.0 1/sqrt(3)],
				atol = 1e-10,
			)


		end
	end
end
