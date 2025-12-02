@testitem "Basis" begin
	include("../../src/utils/AngularMomentumCoupling.jl")
	using .AngularMomentumCoupling
	include("../../src/types/Basis.jl")
	using .Basis
	using Test

	@testset "LinearCombo" begin
		@testset "Constructor" begin
			atoms = [10, 20]
			ls = [1, 1]
			bases, paths = build_all_complex_bases(ls)

			keys = sort(collect(keys(bases)))   # [0, 1, 2]

			# Lf = 0, Lseq = []
		end
	end
end
