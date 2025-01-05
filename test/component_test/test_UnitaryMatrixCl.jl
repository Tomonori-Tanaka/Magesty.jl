using ..UnitaryMatrixCl
using LinearAlgebra

@testset "UnitaryMatrixCl" begin

	l1 = 1
	C1 = UniMatCl(l1)
	@test getindex_m(C1, -l1, -l1) ≈ im / √2
	@test getindex_m(C1, l1, l1) ≈ -1 / √2
	@test getindex_m(C1, 0, 0) ≈ 1
	@test getindex_m(C1, 0, 0) == C1[2, 2]
	@test getindex_m(C1, l1, 0) ≈ 0
	C1adj = C1'
	D = Diagonal([1 + 0im for _ in 1:3])
	A = Matrix(D)
	@test isapprox((C1 * C1adj), A, atol = 1e-10)

	l2 = 2
	C2 = UniMatCl(l2)
	C2inv = inv(C2)
	C2adj = adjoint(C2)

	@test getindex_m(C2, -l2, l2) ≈ -im / √2
	@test getindex_m(C2, 0, 0) ≈ 1
	@test getindex_m(C2, 0, 0) == C2[3, 3]
	@test getindex_m(C2, l2, 0) ≈ 0

end
