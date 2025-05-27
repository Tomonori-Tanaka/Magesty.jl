using ..UnitaryMatrixCl
using LinearAlgebra
using Test

@testset "UnitaryMatrixCl" begin
	@testset "Constructor" begin
		# Test l=0 case
		C0 = UniMatCl(0)
		@test size(C0) == (1, 1)
		@test C0[1, 1] ≈ 1

		# Test l=1 case
		l1 = 1
		C1 = UniMatCl(l1)
		@test size(C1) == (3, 3)
		@test C1[2, 2] ≈ 1  # m = m' = 0

		# Test l=2 case
		l2 = 2
		C2 = UniMatCl(l2)
		@test size(C2) == (5, 5)
		@test C2[3, 3] ≈ 1  # m = m' = 0

		# Test matrix constructor
		mat = Complex[1 0 0; 0 1 0; 0 0 1]  # Explicitly specify Complex type
		C_from_mat = UniMatCl(mat)
		@test size(C_from_mat) == (3, 3)
		@test C_from_mat.l == 1

		# Test error cases
		@test_throws ArgumentError UniMatCl(-1)  # Negative l
		@test_throws ArgumentError UniMatCl(Complex[1 0; 0 1])  # Non-square matrix
		@test_throws ArgumentError UniMatCl(Complex[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])  # Invalid dimension
	end

	@testset "Matrix Elements (l = 1)" begin
		l1 = 1
		C1 = UniMatCl(l1)

		# Test specific elements
		@test getindex_m(C1, -l1, -l1) ≈ im / √2
		@test getindex_m(C1, l1, l1) ≈ -1 / √2
		@test getindex_m(C1, 0, 0) ≈ 1
		@test getindex_m(C1, 0, 0) == C1[2, 2]
		@test getindex_m(C1, l1, 0) ≈ 0

		# Test bounds checking
		@test_throws BoundsError getindex_m(C1, -2, 0)  # m1 out of bounds
		@test_throws BoundsError getindex_m(C1, 0, 2)   # m2 out of bounds
	end

	@testset "Matrix Elements (l = 2)" begin
		l2 = 2
		C2 = UniMatCl(l2)

		correct_matrix =
			1/√2 *
			[im  0  0  0 -im;
				0 im  0  im 0;
				  0 0 √2 0 0;
				  0 1 0 -1 0;
				  1 0 0 0 1;
			]
		@test isapprox(C2, correct_matrix, atol = 1e-10)
		@test isapprox(transpose(C2), transpose(correct_matrix), atol = 1e-10)
		@test isapprox(conj(C2), conj(correct_matrix), atol = 1e-10)	

	end

	@testset "Matrix Operations" begin
		l1 = 1
		C1 = UniMatCl(l1)
		C1adj = C1'
		D = Diagonal([1 + 0im for _ in 1:3])
		A = Matrix(D)

		# Test unitarity
		@test isapprox((C1 * C1adj), A, atol = 1e-10)
		@test isapprox((C1adj * C1), A, atol = 1e-10)

		# Test inverse
		C1inv = inv(C1)
		@test isapprox((C1 * C1inv), A, atol = 1e-10)
		@test isapprox(C1inv, C1adj, atol = 1e-10)

		# Test multiplication with different l values
		l2 = 2
		C2 = UniMatCl(l2)
		@test_throws ArgumentError C1 * C2
	end

	@testset "Additional Operations" begin
		l1 = 1
		C1 = UniMatCl(l1)

		# Test transpose
		C1t = transpose(C1)
		@test size(C1t) == size(C1)
		@test C1t.l == C1.l

		# Test conjugate
		C1c = conj(C1)
		@test size(C1c) == size(C1)
		@test C1c.l == C1.l

		# Test equality
		C1_copy = UniMatCl(l1)
		@test C1 == C1_copy
		@test isapprox(C1, C1_copy, atol = 1e-10)

		# Test immutability
		@test_throws ArgumentError C1[1, 1] = 0
	end
end
