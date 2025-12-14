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

		@testset "Matrix Elements (l = 1)" begin
			U1 = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(1)
			idx(m) = m + 1 + 1  # m + l + 1 for l=1

			# Test m = 0 element
			@test U1[idx(0), idx(0)] ≈ 1

			# Test m = ±1 elements
			@test U1[idx(1), idx(1)] ≈ (-1)^1 / sqrt(2)
			@test U1[idx(1), idx(-1)] ≈ im * (-1)^1 / sqrt(2)
			@test U1[idx(-1), idx(1)] ≈ 1 / sqrt(2)
			@test U1[idx(-1), idx(-1)] ≈ -im / sqrt(2)
		end

		@testset "Matrix Elements (l = 2)" begin
			U2 = SphericalHarmonicsTransforms.c2r_sph_harm_matrix(2)
			idx(m) = m + 2 + 1  # m + l + 1 for l=2

			# Test m = 0 element
			@test U2[idx(0), idx(0)] ≈ 1

			# Test m = ±1 elements
			@test U2[idx(1), idx(1)] ≈ (-1)^1 / sqrt(2)
			@test U2[idx(1), idx(-1)] ≈ im * (-1)^1 / sqrt(2)
			@test U2[idx(-1), idx(1)] ≈ 1 / sqrt(2)
			@test U2[idx(-1), idx(-1)] ≈ -im / sqrt(2)

			# Test m = ±2 elements
			@test U2[idx(2), idx(2)] ≈ (-1)^2 / sqrt(2)
			@test U2[idx(2), idx(-2)] ≈ im * (-1)^2 / sqrt(2)
			@test U2[idx(-2), idx(2)] ≈ 1 / sqrt(2)
			@test U2[idx(-2), idx(-2)] ≈ -im / sqrt(2)
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

		@testset "Matrix Elements (l = 1)" begin
			U1 = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(1)
			idx(m) = m + 1 + 1  # m + l + 1 for l=1

			# Test m = 0 element
			@test U1[idx(0), idx(0)] ≈ 1

			# Test m = ±1 elements
			@test U1[idx(1), idx(1)] ≈ (-1)^1 / sqrt(2)
			@test U1[idx(1), idx(-1)] ≈ 1 / sqrt(2)
			@test U1[idx(-1), idx(1)] ≈ -im * (-1)^1 / sqrt(2)
			@test U1[idx(-1), idx(-1)] ≈ im / sqrt(2)
		end

		@testset "Matrix Elements (l = 2)" begin
			U2 = SphericalHarmonicsTransforms.r2c_sph_harm_matrix(2)
			idx(m) = m + 2 + 1  # m + l + 1 for l=2

			# Test m = 0 element
			@test U2[idx(0), idx(0)] ≈ 1

			# Test m = ±1 elements
			@test U2[idx(1), idx(1)] ≈ (-1)^1 / sqrt(2)
			@test U2[idx(1), idx(-1)] ≈ 1 / sqrt(2)
			@test U2[idx(-1), idx(1)] ≈ -im * (-1)^1 / sqrt(2)
			@test U2[idx(-1), idx(-1)] ≈ im / sqrt(2)

			# Test m = ±2 elements
			@test U2[idx(2), idx(2)] ≈ (-1)^2 / sqrt(2)
			@test U2[idx(2), idx(-2)] ≈ 1 / sqrt(2)
			@test U2[idx(-2), idx(2)] ≈ -im * (-1)^2 / sqrt(2)
			@test U2[idx(-2), idx(-2)] ≈ im / sqrt(2)
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

