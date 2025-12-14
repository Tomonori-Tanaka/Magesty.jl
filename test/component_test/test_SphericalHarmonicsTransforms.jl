using ..SphericalHarmonicsTransforms
using LinearAlgebra
using Test

@testset "SphericalHarmonicsTransforms" begin
	@testset "c2r_sph_harm_matrix" begin
		@testset "Basic properties" begin
			# Test l=0 case
			U0 = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(0)
			@test size(U0) == (1, 1)
			@test U0[1, 1] ≈ 1

			# Test l=1 case
			U1 = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(1)
			@test size(U1) == (3, 3)
			@test U1[2, 2] ≈ 1  # m = m' = 0

			# Test l=2 case
			U2 = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(2)
			@test size(U2) == (5, 5)
			@test U2[3, 3] ≈ 1  # m = m' = 0

			# Test error cases
			@test_throws ArgumentError SphericalHarmonicsTransforms.c2r_sph_harm_matrix(-1)  # Negative l
		end

		@testset "Matrix Elements (l = 0)" begin
			U = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(0)
			# Expected matrix: [1]
			U_expected = [1.0 + 0.0im]
			@test isapprox(U, U_expected, atol = 1e-10)
		end

		@testset "Matrix Elements (l = 1)" begin
			U = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(1)
			# Expected matrix (rows/cols: m = -1, 0, +1)
			# Row m = -1: [im/√2, 0, im/√2]
			# Row m =  0: [0, 1, 0]
			# Row m = +1: [1/√2, 0, -1/√2]
			sqrt2_inv = 1 / sqrt(2)
			U_expected = [
				0.0 + im * sqrt2_inv  0.0 + 0.0im  0.0 + im * sqrt2_inv;
				0.0 + 0.0im           1.0 + 0.0im  0.0 + 0.0im;
				sqrt2_inv              0.0 + 0.0im -sqrt2_inv
			]
			@test isapprox(U, U_expected, atol = 1e-10)
		end

		@testset "Matrix Elements (l = 2)" begin
			U = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(2)
			# Expected matrix (rows/cols: m = -2, -1, 0, +1, +2)
			# Row m = -2: [im/√2, 0, 0, 0, -im/√2]
			# Row m = -1: [0, im/√2, 0, im/√2, 0]
			# Row m =  0: [0, 0, 1, 0, 0]
			# Row m = +1: [0, 1/√2, 0, -1/√2, 0]
			# Row m = +2: [1/√2, 0, 0, 0, 1/√2]
			sqrt2_inv = 1 / sqrt(2)
			U_expected = [
				0.0 + im * sqrt2_inv  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 - im * sqrt2_inv;
				0.0 + 0.0im           0.0 + im * sqrt2_inv  0.0 + 0.0im  0.0 + im * sqrt2_inv  0.0 + 0.0im;
				0.0 + 0.0im           0.0 + 0.0im  1.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im;
				0.0 + 0.0im           sqrt2_inv        0.0 + 0.0im -sqrt2_inv       0.0 + 0.0im;
				sqrt2_inv             0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  sqrt2_inv
			]
			@test isapprox(U, U_expected, atol = 1e-10)
		end
	end

	@testset "r2c_sph_harm_matrix" begin
		@testset "Basic properties" begin
			# Test l=0 case
			U0 = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(0)
			@test size(U0) == (1, 1)
			@test U0[1, 1] ≈ 1

			# Test l=1 case
			U1 = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(1)
			@test size(U1) == (3, 3)
			@test U1[2, 2] ≈ 1  # m = m' = 0

			# Test l=2 case
			U2 = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(2)
			@test size(U2) == (5, 5)
			@test U2[3, 3] ≈ 1  # m = m' = 0

			# Test error cases
			@test_throws ArgumentError SphericalHarmonicsTransforms.r2c_sph_harm_matrix(-1)  # Negative l
		end

		@testset "Matrix Elements (l = 0)" begin
			U = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(0)
			# Expected matrix: [1]
			U_expected = [1.0 + 0.0im]
			@test isapprox(U, U_expected, atol = 1e-10)
		end

		@testset "Matrix Elements (l = 1)" begin
			U = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(1)
			# Expected matrix (rows/cols: m = -1, 0, +1)
			# Row m = -1: [-im/√2, 0, 1/√2]
			# Row m =  0: [0, 1, 0]
			# Row m = +1: [-im/√2, 0, -1/√2]
			sqrt2_inv = 1 / sqrt(2)
			U_expected = [
				-im * sqrt2_inv   0.0 + 0.0im  sqrt2_inv;
				0.0 + 0.0im       1.0 + 0.0im  0.0 + 0.0im;
				-im * sqrt2_inv   0.0 + 0.0im -sqrt2_inv
			]
			@test isapprox(U, U_expected, atol = 1e-10)
		end

		@testset "Matrix Elements (l = 2)" begin
			U = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(2)
			# Expected matrix (rows/cols: m = -2, -1, 0, +1, +2)
			# Row m = -2: [-im/√2, 0, 0, 0, 1/√2]
			# Row m = -1: [0, -im/√2, 0, 1/√2, 0]
			# Row m =  0: [0, 0, 1, 0, 0]
			# Row m = +1: [0, -im/√2, 0, -1/√2, 0]
			# Row m = +2: [im/√2, 0, 0, 0, 1/√2]
			sqrt2_inv = 1 / sqrt(2)
			U_expected = [
				-im * sqrt2_inv   0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  sqrt2_inv;
				0.0 + 0.0im      -im * sqrt2_inv  0.0 + 0.0im  sqrt2_inv        0.0 + 0.0im;
				0.0 + 0.0im      0.0 + 0.0im  1.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im;
				0.0 + 0.0im      -im * sqrt2_inv  0.0 + 0.0im -sqrt2_inv       0.0 + 0.0im;
				im * sqrt2_inv   0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  sqrt2_inv
			]
			@test isapprox(U, U_expected, atol = 1e-10)
		end
	end

	@testset "Consistency: Inverse relationship" begin
		@testset "l = 0" begin
			U_c2r = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(0)
			U_r2c = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(0)
			I_expected = Matrix{ComplexF64}(I, 1, 1)

			@test isapprox(U_c2r * U_r2c, I_expected, atol = 1e-10)
			@test isapprox(U_r2c * U_c2r, I_expected, atol = 1e-10)
		end

		@testset "l = 1" begin
			U_c2r = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(1)
			U_r2c = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(1)
			I_expected = Matrix{ComplexF64}(I, 3, 3)

			@test isapprox(U_c2r * U_r2c, I_expected, atol = 1e-10)
			@test isapprox(U_r2c * U_c2r, I_expected, atol = 1e-10)
		end

		@testset "l = 2" begin
			U_c2r = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(2)
			U_r2c = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(2)
			I_expected = Matrix{ComplexF64}(I, 5, 5)

			@test isapprox(U_c2r * U_r2c, I_expected, atol = 1e-10)
			@test isapprox(U_r2c * U_c2r, I_expected, atol = 1e-10)
		end

		@testset "l = 3" begin
			U_c2r = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(3)
			U_r2c = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(3)
			I_expected = Matrix{ComplexF64}(I, 7, 7)

			@test isapprox(U_c2r * U_r2c, I_expected, atol = 1e-10)
			@test isapprox(U_r2c * U_c2r, I_expected, atol = 1e-10)
		end

		@testset "Multiple l values" begin
			for l in [0, 1, 2, 3, 4, 5]
				U_c2r = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(l)
				U_r2c = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(l)
				I_expected = Matrix{ComplexF64}(I, 2l+1, 2l+1)

				@test isapprox(U_c2r * U_r2c, I_expected, atol = 1e-10)
				@test isapprox(U_r2c * U_c2r, I_expected, atol = 1e-10)
			end
		end
	end
end

