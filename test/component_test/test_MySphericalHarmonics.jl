using Test
using LinearAlgebra
using LegendrePolynomials

using ..MySphericalHarmonics
using ..MySphericalHarmonics: P̄ₗₘ, dP̄ₗₘ, Yₗₘ, ∂Yₗₘ_∂r̂x, ∂Yₗₘ_∂r̂y, ∂Yₗₘ_∂r̂z, yₗₘ,
	∂Sₗₘ_∂r̂x, ∂Sₗₘ_∂r̂y, ∂Sₗₘ_∂r̂z, ssₗₘ, ∂Sₗₘ_∂x, ∂Sₗₘ_∂y, ∂Sₗₘ_∂z

@testset "Legendre polynomials" begin
	# (l, m, r̂z)
	@testset "P̄ₗₘ" begin
		# Error handling
		@test_throws ArgumentError P̄ₗₘ(0, 0, 1.1)
		@test_throws ArgumentError P̄ₗₘ(0, 0, [1, 1, 1])
		@test_throws ArgumentError P̄ₗₘ(0, 0, [1, 0])
		@test_throws ArgumentError P̄ₗₘ(1, -1, [1, 0, 0])
		@test_throws ArgumentError P̄ₗₘ(1, 2, [1, 0, 0])
		@test_throws ArgumentError P̄ₗₘ(-1, 0, [1, 0, 0])

		# l = 0, m = 0, √(1/4π)
		@test isapprox(P̄ₗₘ(0, 0, 0), √(1 / 4π))
		@test isapprox(P̄ₗₘ(0, 0, 1), √(1 / 4π))
		@test isapprox(P̄ₗₘ(0, 0, -1), √(1 / 4π))

		# l = 1, m = 0, √(3/4π)*r̂z
		@test isapprox(P̄ₗₘ(1, 0, 0), 0)
		@test isapprox(P̄ₗₘ(1, 0, 1), √(3 / 4π) * 1)
		@test isapprox(P̄ₗₘ(1, 0, -1), √(3 / 4π) * -1)

		# l = 1, m = 1, -√(3/8π)
		@test isapprox(P̄ₗₘ(1, 1, 0), -√(3 / 8π))
		@test isapprox(P̄ₗₘ(1, 1, 1), -√(3 / 8π))
		@test isapprox(P̄ₗₘ(1, 1, -1), -√(3 / 8π))

		# l = 2, m = 0, √(5/16π)*(3r̂z^2 - 1)
		@test isapprox(P̄ₗₘ(2, 0, 0), -√(5 / 16π))
		@test isapprox(P̄ₗₘ(2, 0, 1), √(5 / 16π) * (3 * 1^2 - 1))
		@test isapprox(P̄ₗₘ(2, 0, -1), √(5 / 16π) * (3 * (-1)^2 - 1))

		# l = 2, m = 1, -√(15/8π)*r̂z
		@test isapprox(P̄ₗₘ(2, 1, 0), 0)
		@test isapprox(P̄ₗₘ(2, 1, 1), -√(15 / 8π) * 1)
		@test isapprox(P̄ₗₘ(2, 1, -1), -√(15 / 8π) * -1)

		# l = 2, m = 2, √(45/96π)
		@test isapprox(P̄ₗₘ(2, 2, 0), √(45 / 96π))
		@test isapprox(P̄ₗₘ(2, 2, 1), √(45 / 96π))
		@test isapprox(P̄ₗₘ(2, 2, -1), √(45 / 96π))
	end

	@testset "dP̄ₗₘ" begin
		# l = 0, m = 0, 0
		@test isapprox(dP̄ₗₘ(0, 0, 0), 0)
		@test isapprox(dP̄ₗₘ(0, 0, 1), 0)
		@test isapprox(dP̄ₗₘ(0, 0, -1), 0)

		# l = 1, m = 0, √(3/4π)
		@test isapprox(dP̄ₗₘ(1, 0, 0), √(3 / 4π))
		@test isapprox(dP̄ₗₘ(1, 0, 1), √(3 / 4π))
		@test isapprox(dP̄ₗₘ(1, 0, -1), √(3 / 4π))

		# l = 1, m = 1, 0
		@test isapprox(dP̄ₗₘ(1, 1, 0), 0)
		@test isapprox(dP̄ₗₘ(1, 1, 1), 0)
		@test isapprox(dP̄ₗₘ(1, 1, -1), 0)

		# l = 2, m = 0, √(5/16π)*(6r̂z)
		@test isapprox(dP̄ₗₘ(2, 0, 0), 0)
		@test isapprox(dP̄ₗₘ(2, 0, 1), √(5 / 16π) * 6 * 1)
		@test isapprox(dP̄ₗₘ(2, 0, -1), √(5 / 16π) * 6 * -1)

		# l = 2, m = 1, -√(15/8π)
		@test isapprox(dP̄ₗₘ(2, 1, 0), -√(15 / 8π))
		@test isapprox(dP̄ₗₘ(2, 1, 1), -√(15 / 8π))
		@test isapprox(dP̄ₗₘ(2, 1, -1), -√(15 / 8π))

		# l = 2, m = 2, 0
		@test isapprox(dP̄ₗₘ(2, 2, 0), 0)
		@test isapprox(dP̄ₗₘ(2, 2, 1), 0)
		@test isapprox(dP̄ₗₘ(2, 2, -1), 0)
	end
end

@testset "complex spherical harmonics" begin
	@testset "Yₗₘ" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# Error handling
		@test_throws ArgumentError Yₗₘ(0, 0, [1, 1, 0])
		@test_throws ArgumentError Yₗₘ(0, 0, [1, 0])
		@test_throws ArgumentError Yₗₘ(0, 1, [1, 0, 0])
		@test_throws ArgumentError Yₗₘ(0, -1, [1, 0, 0])

		# l = 0, m = 0, √(1/4π)
		@test isapprox(Yₗₘ(0, 0, [1, 0, 0]), 1 / √(4π))
		@test isapprox(Yₗₘ(0, 0, [0, 1, 0]), 1 / √(4π))
		@test isapprox(Yₗₘ(0, 0, [0, 0, 1]), 1 / √(4π))
		@test isapprox(Yₗₘ(0, 0, [0, 0, -1]), 1 / √(4π))

		# l = 1, m = 0, √(3/4π)*r̂z
		@test isapprox(Yₗₘ(1, 0, [1, 0, 0]), 0)
		@test isapprox(Yₗₘ(1, 0, [0, 1, 0]), 0)
		@test isapprox(Yₗₘ(1, 0, [0, 0, 1]), √(3 / 4π))
		@test isapprox(Yₗₘ(1, 0, [0, 0, -1]), -√(3 / 4π))

		# l = 1, m = 1, -√(3/8π)*(r̂x + ir̂y)
		@test isapprox(Yₗₘ(1, 1, [1, 0, 0]), -√(3 / 8π))
		@test isapprox(Yₗₘ(1, 1, [-1, 0, 0]), √(3 / 8π))
		@test isapprox(Yₗₘ(1, 1, [0, 1, 0]), -√(3 / 8π)im)
		@test isapprox(Yₗₘ(1, 1, [0, -1, 0]), √(3 / 8π)im)
		@test isapprox(Yₗₘ(1, 1, [0, 0, 1]), 0)
		@test isapprox(Yₗₘ(1, 1, [0, 0, -1]), 0)

		# l = 1, m = -1, √(3/8π)*(r̂x - ir̂y)
		@test isapprox(Yₗₘ(1, -1, [1, 0, 0]), √(3 / 8π))
		@test isapprox(Yₗₘ(1, -1, [-1, 0, 0]), -√(3 / 8π))
		@test isapprox(Yₗₘ(1, -1, [0, 1, 0]), -√(3 / 8π)im)
		@test isapprox(Yₗₘ(1, -1, [0, -1, 0]), √(3 / 8π)im)
		@test isapprox(Yₗₘ(1, -1, [0, 0, 1]), 0)
		@test isapprox(Yₗₘ(1, -1, [0, 0, -1]), 0)

		# l = 2, m = 0, √(5/16π)*(3r̂z^2 - 1)
		@test isapprox(Yₗₘ(2, 0, [1, 0, 0]), -√(5 / 16π))
		@test isapprox(Yₗₘ(2, 0, [0, 1, 0]), -√(5 / 16π))
		@test isapprox(Yₗₘ(2, 0, [0, 0, 1]), √(5 / 16π) * (3 * 1^2 - 1))
		@test isapprox(Yₗₘ(2, 0, [0, 0, -1]), √(5 / 16π) * (3 * (-1)^2 - 1))

		# l = 2, m = 1, -√(15/8π)*(r̂x + ir̂y)*r̂z
		@test isapprox(Yₗₘ(2, 1, [1, 0, 0]), 0)
		@test isapprox(Yₗₘ(2, 1, [0, 1, 0]), 0)
		@test isapprox(Yₗₘ(2, 1, [0, 0, 1]), 0)
		@test isapprox(Yₗₘ(2, 1, [0, 0, -1]), 0)
		@test isapprox(Yₗₘ(2, 1, vec111), -√(15 / 8π) * (1 / 3) * (1 + im))
		@test isapprox(Yₗₘ(2, 1, -vec111), -√(15 / 8π) * (1 / 3) * (1 + im))

		# l = 2, m = -1, √(15/8π)*(r̂x - ir̂y)*r̂z
		@test isapprox(Yₗₘ(2, -1, [1, 0, 0]), 0)
		@test isapprox(Yₗₘ(2, -1, [0, 1, 0]), 0)
		@test isapprox(Yₗₘ(2, -1, [0, 0, 1]), 0)
		@test isapprox(Yₗₘ(2, -1, [0, 0, -1]), 0)
		@test isapprox(Yₗₘ(2, -1, vec111), √(15 / 8π) * (1 / 3) * (1 - im))
		@test isapprox(Yₗₘ(2, -1, -vec111), √(15 / 8π) * (1 / 3) * (1 - im))


		# l = 2, m = 2, √(45/96π)*(r̂x + ir̂y)^2
		@test isapprox(Yₗₘ(2, 2, [1, 0, 0]), √(45 / 96π))
		@test isapprox(Yₗₘ(2, 2, [-1, 0, 0]), √(45 / 96π))
		@test isapprox(Yₗₘ(2, 2, [0, 1, 0]), -√(45 / 96π))
		@test isapprox(Yₗₘ(2, 2, [0, -1, 0]), -√(45 / 96π))
		@test isapprox(Yₗₘ(2, 2, [0, 0, 1]), 0)
		@test isapprox(Yₗₘ(2, 2, [0, 0, -1]), 0)
		@test isapprox(Yₗₘ(2, 2, vec111), √(45 / 96π) * (1 / 3) * 2im)
		@test isapprox(Yₗₘ(2, 2, -vec111), √(45 / 96π) * (1 / 3) * (2im))

		# l = 2, m = -2, √(45/96π)*(r̂x - ir̂y)^2
		@test isapprox(Yₗₘ(2, -2, [1, 0, 0]), √(45 / 96π))
		@test isapprox(Yₗₘ(2, -2, [-1, 0, 0]), √(45 / 96π))
		@test isapprox(Yₗₘ(2, -2, [0, 1, 0]), -√(45 / 96π))
		@test isapprox(Yₗₘ(2, -2, [0, -1, 0]), -√(45 / 96π))
		@test isapprox(Yₗₘ(2, -2, [0, 0, 1]), 0)
		@test isapprox(Yₗₘ(2, -2, [0, 0, -1]), 0)
		@test isapprox(
			Yₗₘ(2, -2, vec111),
			√(45 / 96π) * (1 / 3) * (-2im),
		)
		@test isapprox(
			Yₗₘ(2, -2, -vec111),
			√(45 / 96π) * (1 / 3) * (-2im),
		)
	end

	@testset "∂Yₗₘ_∂r̂x" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Yₗₘ_∂r̂x(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(0, 0, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(0, 0, [0, 0, -1]), 0)

		# l = 1, m = 0, 0
		@test isapprox(∂Yₗₘ_∂r̂x(1, 0, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(1, 0, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(1, 0, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(1, 0, [0, 0, -1]), 0)
		# l = 1, m = 1, -√(3/8π)
		@test isapprox(∂Yₗₘ_∂r̂x(1, 1, [1, 0, 0]), -√(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, 1, [-1, 0, 0]), -√(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, 1, [0, 1, 0]), -√(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, 1, [0, -1, 0]), -√(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, 1, [0, 0, 1]), -√(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, 1, [0, 0, -1]), -√(3 / 8π))
		# l = 1, m = -1, √(3/8π)
		@test isapprox(∂Yₗₘ_∂r̂x(1, -1, [1, 0, 0]), √(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, -1, [-1, 0, 0]), √(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, -1, [0, 1, 0]), √(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, -1, [0, -1, 0]), √(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, -1, [0, 0, 1]), √(3 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(1, -1, [0, 0, -1]), √(3 / 8π))

		# l = 2, m = 0, 0
		@test isapprox(∂Yₗₘ_∂r̂x(2, 0, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 0, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 0, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 0, [0, 0, -1]), 0)
		# l = 2, m = 1, -√(15/8π)*r̂z
		@test isapprox(∂Yₗₘ_∂r̂x(2, 1, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 1, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 1, [0, 0, 1]), -√(15 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(2, 1, [0, 0, -1]), √(15 / 8π))
		# l = 2, m = -1, √(15/8π)*r̂z
		@test isapprox(∂Yₗₘ_∂r̂x(2, -1, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, -1, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, -1, [0, 0, 1]), √(15 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂x(2, -1, [0, 0, -1]), -√(15 / 8π))
		# l = 2, m = 2, √(45/96π)*2*(r̂x + ir̂y)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 2, [1, 0, 0]), 2 * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂x(2, 2, [-1, 0, 0]), -2 * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂x(2, 2, [0, 1, 0]), 2 * √(45 / 96π)im)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 2, [0, -1, 0]), -2 * √(45 / 96π)im)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 2, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, 2, [0, 0, -1]), 0)
		@test isapprox(
			∂Yₗₘ_∂r̂x(2, 2, vec111),
			2 * √(45 / 96π) * (1 / √(3)) * (1 + im),
		)
		@test isapprox(
			∂Yₗₘ_∂r̂x(2, 2, -vec111),
			-2 * √(45 / 96π) * (1 / √(3)) * (1 + im),
		)
		# l = 2, m = -2, √(45/96π)*2*(r̂x - ir̂y)
		@test isapprox(∂Yₗₘ_∂r̂x(2, -2, [1, 0, 0]), 2 * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂x(2, -2, [-1, 0, 0]), -2 * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂x(2, -2, [0, 1, 0]), -2 * √(45 / 96π)im)
		@test isapprox(∂Yₗₘ_∂r̂x(2, -2, [0, -1, 0]), 2 * √(45 / 96π)im)
		@test isapprox(∂Yₗₘ_∂r̂x(2, -2, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂x(2, -2, [0, 0, -1]), 0)
		@test isapprox(
			∂Yₗₘ_∂r̂x(2, -2, vec111),
			2 * √(45 / 96π) * (1 / √(3)) * (1 - im),
		)
		@test isapprox(
			∂Yₗₘ_∂r̂x(2, -2, -vec111),
			-2 * √(45 / 96π) * (1 / √(3)) * (1 - im),
		)
	end

	@testset "∂Yₗₘ_∂r̂y" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Yₗₘ_∂r̂y(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(0, 0, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(0, 0, [0, 0, -1]), 0)

		# l = 1, m = 0, 0
		@test isapprox(∂Yₗₘ_∂r̂y(1, 0, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(1, 0, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(1, 0, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(1, 0, [0, 0, -1]), 0)
		# l = 1, m = 1, -√(3/8π)im
		@test isapprox(∂Yₗₘ_∂r̂y(1, 1, [1, 0, 0]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, 1, [-1, 0, 0]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, 1, [0, 1, 0]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, 1, [0, -1, 0]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, 1, [0, 0, 1]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, 1, [0, 0, -1]), -√(3 / 8π)im)
		# l = 1, m = -1, -√(3/8π)im
		@test isapprox(∂Yₗₘ_∂r̂y(1, -1, [1, 0, 0]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, -1, [-1, 0, 0]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, -1, [0, 1, 0]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, -1, [0, -1, 0]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, -1, [0, 0, 1]), -√(3 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(1, -1, [0, 0, -1]), -√(3 / 8π)im)

		# l = 2, m = 0, 0
		@test isapprox(∂Yₗₘ_∂r̂y(2, 0, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 0, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 0, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 0, [0, 0, -1]), 0)
		# l = 2, m = 1, -im*1*√(15/8π)*r̂z
		@test isapprox(∂Yₗₘ_∂r̂y(2, 1, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 1, [-1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 1, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 1, [0, -1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 1, [0, 0, 1]), -√(15 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 1, [0, 0, -1]), √(15 / 8π)im)
		# l = 2, m = -1, -im*1*√(15/8π)*r̂z
		@test isapprox(∂Yₗₘ_∂r̂y(2, -1, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, -1, [-1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, -1, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, -1, [0, -1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, -1, [0, 0, 1]), -√(15 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂y(2, -1, [0, 0, -1]), √(15 / 8π)im)
		# l = 2, m = 2, im*2(r̂x + im*r̂y)√(45/96π)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 2, [1, 0, 0]), 2im * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂y(2, 2, [-1, 0, 0]), -2im * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂y(2, 2, [0, 1, 0]), -2 * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂y(2, 2, [0, -1, 0]), 2 * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂y(2, 2, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, 2, [0, 0, -1]), 0)
		@test isapprox(
			∂Yₗₘ_∂r̂y(2, 2, vec111),
			2im * (1 / √(3)) * √(45 / 96π) * (1 + im),
		)
		@test isapprox(
			∂Yₗₘ_∂r̂y(2, 2, -vec111),
			-2im * (1 / √(3)) * √(45 / 96π) * (1 + im),
		)
		# l = 2, m = -2, -im*2(r̂x - im*r̂y)√(45/96π)
		@test isapprox(∂Yₗₘ_∂r̂y(2, -2, [1, 0, 0]), -2im * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂y(2, -2, [-1, 0, 0]), 2im * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂y(2, -2, [0, 1, 0]), -2 * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂y(2, -2, [0, -1, 0]), 2 * √(45 / 96π))
		@test isapprox(∂Yₗₘ_∂r̂y(2, -2, [0, 0, 1]), 0)
		@test isapprox(∂Yₗₘ_∂r̂y(2, -2, [0, 0, -1]), 0)
		@test isapprox(
			∂Yₗₘ_∂r̂y(2, -2, vec111),
			-2im * (1 / √(3)) * √(45 / 96π) * (1 - im),
		)
		@test isapprox(
			∂Yₗₘ_∂r̂y(2, -2, -vec111),
			2im * (1 / √(3)) * √(45 / 96π) * (1 - im),
		)
	end

	@testset "∂Yₗₘ_∂r̂z" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Yₗₘ_∂r̂z(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(0, 0, [0, 0, 1]), 0)

		# l = 1, m = 0, √(3/4π)
		@test isapprox(∂Yₗₘ_∂r̂z(1, 0, [1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Yₗₘ_∂r̂z(1, 0, [0, 1, 0]), √(3 / 4π))
		@test isapprox(∂Yₗₘ_∂r̂z(1, 0, [0, 0, 1]), √(3 / 4π))
		# l = 1, m = 1, 0
		@test isapprox(∂Yₗₘ_∂r̂z(1, 1, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(1, 1, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(1, 1, [0, 0, 1]), 0)
		# l = 1, m = -1, 0
		@test isapprox(∂Yₗₘ_∂r̂z(1, -1, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(1, -1, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(1, -1, [0, 0, 1]), 0)

		# l = 2, m = 0, 6√(5/16π)r̂z
		@test isapprox(∂Yₗₘ_∂r̂z(2, 0, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(2, 0, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(2, 0, [0, 0, 1]), 6 * √(5 / 16π))
		@test isapprox(∂Yₗₘ_∂r̂z(2, 0, [0, 0, -1]), -6 * √(5 / 16π))
		# l = 2, m = 1, -√(15/8π)(r̂x + im*r̂y)
		@test isapprox(∂Yₗₘ_∂r̂z(2, 1, [1, 0, 0]), -√(15 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂z(2, 1, [0, 1, 0]), -√(15 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂z(2, 1, [0, 0, 1]), 0)
		@test isapprox(
			∂Yₗₘ_∂r̂z(2, 1, vec111),
			-√(15 / 8π) * (1 / √(3)) * (1 + im),
		)
		@test isapprox(
			∂Yₗₘ_∂r̂z(2, 1, -vec111),
			√(15 / 8π) * (1 / √(3)) * (1 + im),
		)
		# l = 2, m = -1, √(15/8π)(r̂x - im*r̂y)
		@test isapprox(∂Yₗₘ_∂r̂z(2, -1, [1, 0, 0]), √(15 / 8π))
		@test isapprox(∂Yₗₘ_∂r̂z(2, -1, [0, 1, 0]), -√(15 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂z(2, -1, [0, -1, 0]), √(15 / 8π)im)
		@test isapprox(∂Yₗₘ_∂r̂z(2, -1, [0, 0, 1]), 0)
		@test isapprox(
			∂Yₗₘ_∂r̂z(2, -1, vec111),
			√(15 / 8π) * (1 / √(3)) * (1 - im),
		)
		@test isapprox(
			∂Yₗₘ_∂r̂z(2, -1, -vec111),
			-√(15 / 8π) * (1 / √(3)) * (1 - im),
		)
		# l = 2, m = 2, 0
		@test isapprox(∂Yₗₘ_∂r̂z(2, 2, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(2, 2, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(2, 2, [0, 0, 1]), 0)
		# l = 2, m = -2, 0
		@test isapprox(∂Yₗₘ_∂r̂z(2, -2, [1, 0, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(2, -2, [0, 1, 0]), 0)
		@test isapprox(∂Yₗₘ_∂r̂z(2, -2, [0, 0, 1]), 0)
	end

	@testset "yₗₘ" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(yₗₘ(0, 0, [1, 0, 0]), 0)
		@test isapprox(yₗₘ(0, 0, [0, 1, 0]), 0)
		@test isapprox(yₗₘ(0, 0, [0, 0, 1]), 0)

		# l = 1, m = 0, √(3/4π)*r̂z
		@test isapprox(yₗₘ(1, 0, [1, 0, 0]), 0)
		@test isapprox(yₗₘ(1, 0, [0, 1, 0]), 0)
		@test isapprox(yₗₘ(1, 0, [0, 0, 1]), √(3 / 4π))
		@test isapprox(yₗₘ(1, 0, vec111), √(1 / 4π))
		@test isapprox(yₗₘ(1, 0, -vec111), -√(1 / 4π))
		# l = 1, m = 1, -√(3/8π)*(r̂x + im*r̂y)
		@test isapprox(yₗₘ(1, 1, [1, 0, 0]), -√(3 / 8π))
		@test isapprox(yₗₘ(1, 1, [-1, 0, 0]), √(3 / 8π))
		@test isapprox(yₗₘ(1, 1, [0, 1, 0]), -√(3 / 8π)im)
		@test isapprox(yₗₘ(1, 1, [0, -1, 0]), √(3 / 8π)im)
		@test isapprox(yₗₘ(1, 1, [0, 0, 1]), 0)
		@test isapprox(yₗₘ(1, 1, vec111), -√(1 / 8π) * (1 + im))
		@test isapprox(yₗₘ(1, 1, -vec111), √(1 / 8π) * (1 + im))
		# l = 1, m = -1, √(3/8π)*(r̂x - im*r̂y)
		@test isapprox(yₗₘ(1, -1, [1, 0, 0]), √(3 / 8π))
		@test isapprox(yₗₘ(1, -1, [-1, 0, 0]), -√(3 / 8π))
		@test isapprox(yₗₘ(1, -1, [0, 1, 0]), -√(3 / 8π)im)
		@test isapprox(yₗₘ(1, -1, [0, -1, 0]), √(3 / 8π)im)
		@test isapprox(yₗₘ(1, -1, [0, 0, 1]), 0)
		@test isapprox(yₗₘ(1, -1, vec111), √(1 / 8π) * (1 - im))
		@test isapprox(yₗₘ(1, -1, -vec111), -√(1 / 8π) * (1 - im))

		# l = 2, m = 0, 6√(5/16π)*r̂z^2
		@test isapprox(yₗₘ(2, 0, [1, 0, 0]), 0)
		@test isapprox(yₗₘ(2, 0, [0, 1, 0]), 0)
		@test isapprox(yₗₘ(2, 0, [0, 0, 1]), 6 * √(5 / 16π))
		@test isapprox(yₗₘ(2, 0, [0, 0, -1]), 6 * √(5 / 16π))
		@test isapprox(yₗₘ(2, 0, vec111), 2 * √(5 / 16π))
		@test isapprox(yₗₘ(2, 0, -vec111), 2 * √(5 / 16π))
		# l = 2, m = 1, -√(15/2π)*(r̂x + im*r̂y)*r̂z
		@test isapprox(yₗₘ(2, 1, [1, 0, 0]), 0)
		@test isapprox(yₗₘ(2, 1, [0, 1, 0]), 0)
		@test isapprox(yₗₘ(2, 1, [0, 0, 1]), 0)
		@test isapprox(yₗₘ(2, 1, vec111), -(1 / 3) * √(15 / 2π) * (1 + im))
		@test isapprox(yₗₘ(2, 1, -vec111), -(1 / 3) * √(15 / 2π) * (1 + im))
		# l = 2, m = -1, √(15/2π)*(r̂x - im*r̂y)*r̂z
		@test isapprox(yₗₘ(2, -1, [1, 0, 0]), 0)
		@test isapprox(yₗₘ(2, -1, [0, 1, 0]), 0)
		@test isapprox(yₗₘ(2, -1, [0, 0, 1]), 0)
		@test isapprox(yₗₘ(2, -1, vec111), (1 / 3) * √(15 / 2π) * (1 - im))
		@test isapprox(yₗₘ(2, -1, -vec111), (1 / 3) * √(15 / 2π) * (1 - im))
		# l = 2, m = 2, 2*√(45/96π)*(r̂x + im*r̂y)^2
		@test isapprox(yₗₘ(2, 2, [1, 0, 0]), 2 * √(45 / 96π))
		@test isapprox(yₗₘ(2, 2, [-1, 0, 0]), 2 * √(45 / 96π))
		@test isapprox(yₗₘ(2, 2, [0, 1, 0]), -2 * √(45 / 96π))
		@test isapprox(yₗₘ(2, 2, [0, -1, 0]), -2 * √(45 / 96π))
		@test isapprox(yₗₘ(2, 2, [0, 0, 1]), 0)
		@test isapprox(yₗₘ(2, 2, [0, 0, -1]), 0)
		@test isapprox(yₗₘ(2, 2, vec111), 2 / 3 * √(45 / 96π) * (1 + im)^2)
		@test isapprox(yₗₘ(2, 2, -vec111), 2 / 3 * √(45 / 96π) * (1 + im)^2)
		# l = 2, m = -2, 2*√(45/96π)*(r̂x - im*r̂y)^2
		@test isapprox(yₗₘ(2, -2, [1, 0, 0]), 2 * √(45 / 96π))
		@test isapprox(yₗₘ(2, -2, [-1, 0, 0]), 2 * √(45 / 96π))
		@test isapprox(yₗₘ(2, -2, [0, 1, 0]), -2 * √(45 / 96π))
		@test isapprox(yₗₘ(2, -2, [0, -1, 0]), -2 * √(45 / 96π))
		@test isapprox(yₗₘ(2, -2, [0, 0, 1]), 0)
		@test isapprox(yₗₘ(2, -2, [0, 0, -1]), 0)
		@test isapprox(yₗₘ(2, -2, vec111), 2 / 3 * √(45 / 96π) * (1 - im)^2)
		@test isapprox(yₗₘ(2, -2, -vec111), 2 / 3 * √(45 / 96π) * (1 - im)^2)
	end
end

@testset "real spherical harmonics" begin

	@testset "Sₗₘ" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 1/√(4π)
		@test isapprox(Sₗₘ(0, 0, [1, 0, 0]), 1 / √(4π))
		@test isapprox(Sₗₘ(0, 0, [0, 1, 0]), 1 / √(4π))
		@test isapprox(Sₗₘ(0, 0, [0, 0, 1]), 1 / √(4π))
		@test isapprox(Sₗₘ(0, 0, vec111), 1 / √(4π))
		@test isapprox(Sₗₘ(0, 0, -vec111), 1 / √(4π))

		# l = 1, m = 0, √(3/4π)*r̂z
		@test isapprox(Sₗₘ(1, 0, [1, 0, 0]), 0)
		@test isapprox(Sₗₘ(1, 0, [0, 1, 0]), 0)
		@test isapprox(Sₗₘ(1, 0, [0, 0, 1]), √(3 / 4π))
		@test isapprox(Sₗₘ(1, 0, [0, 0, -1]), -√(3 / 4π))
		@test isapprox(Sₗₘ(1, 0, vec111), √(1 / 4π))
		@test isapprox(Sₗₘ(1, 0, -vec111), -√(1 / 4π))
		# l = 1, m = 1, √(3/4π)*r̂x
		@test isapprox(Sₗₘ(1, 1, [1, 0, 0]), √(3 / 4π))
		@test isapprox(Sₗₘ(1, 1, [-1, 0, 0]), -√(3 / 4π))
		@test isapprox(Sₗₘ(1, 1, [0, 1, 0]), 0)
		@test isapprox(Sₗₘ(1, 1, [0, 0, 1]), 0)
		@test isapprox(Sₗₘ(1, 1, vec111), √(1 / 4π))
		@test isapprox(Sₗₘ(1, 1, -vec111), -√(1 / 4π))
		# l = 1, m = -1, √(3/4π)*r̂y
		@test isapprox(Sₗₘ(1, -1, [1, 0, 0]), 0)
		@test isapprox(Sₗₘ(1, -1, [-1, 0, 0]), 0)
		@test isapprox(Sₗₘ(1, -1, [0, 1, 0]), √(3 / 4π))
		@test isapprox(Sₗₘ(1, -1, [0, -1, 0]), -√(3 / 4π))
		@test isapprox(Sₗₘ(1, -1, [0, 0, 1]), 0)
		@test isapprox(Sₗₘ(1, -1, [0, 0, -1]), 0)
		@test isapprox(Sₗₘ(1, -1, vec111), √(1 / 4π))
		@test isapprox(Sₗₘ(1, -1, -vec111), -√(1 / 4π))

		# l = 2, m = 0, √(5/16π)*(3r̂z^2 - 1)
		@test isapprox(Sₗₘ(2, 0, [1, 0, 0]), -√(5 / 16π))
		@test isapprox(Sₗₘ(2, 0, [0, 1, 0]), -√(5 / 16π))
		@test isapprox(Sₗₘ(2, 0, [0, 0, 1]), √(5 / 16π) * 2)
		@test isapprox(Sₗₘ(2, 0, [0, 0, -1]), √(5 / 16π) * 2)
		@test isapprox(Sₗₘ(2, 0, vec111), 0, atol = 1e-15)
		@test isapprox(Sₗₘ(2, 0, -vec111), 0, atol = 1e-15)
		# l = 2, m = 1, √(15/4π)*r̂z*r̂x
		@test isapprox(Sₗₘ(2, 1, [1, 0, 0]), 0)
		@test isapprox(Sₗₘ(2, 1, [-1, 0, 0]), 0)
		@test isapprox(Sₗₘ(2, 1, [0, 1, 0]), 0)
		@test isapprox(Sₗₘ(2, 1, [0, -1, 0]), 0)
		@test isapprox(Sₗₘ(2, 1, [0, 0, 1]), 0)
		@test isapprox(Sₗₘ(2, 1, [0, 0, -1]), 0)
		@test isapprox(Sₗₘ(2, 1, vec111), √(15 / 4π) * (1 / 3))
		@test isapprox(Sₗₘ(2, 1, -vec111), √(15 / 4π) * (1 / 3))
		@test isapprox(Sₗₘ(2, 1, [1 / √(2), 0, 1 / √(2)]), √(15 / 4π) * (1 / 2))
		@test isapprox(Sₗₘ(2, 1, [1 / √(2), 0, -1 / √(2)]), -√(15 / 4π) * (1 / 2))
		@test isapprox(Sₗₘ(2, 1, [-1 / √(2), 0, 1 / √(2)]), -√(15 / 4π) * (1 / 2))
		@test isapprox(Sₗₘ(2, 1, [-1 / √(2), 0, -1 / √(2)]), √(15 / 4π) * (1 / 2))
		# l = 2, m = -1, √(15/4π)*r̂z*r̂y
		@test isapprox(Sₗₘ(2, -1, [1, 0, 0]), 0)
		@test isapprox(Sₗₘ(2, -1, [-1, 0, 0]), 0)
		@test isapprox(Sₗₘ(2, -1, [0, 1, 0]), 0)
		@test isapprox(Sₗₘ(2, -1, [0, -1, 0]), 0)
		@test isapprox(Sₗₘ(2, -1, [0, 0, 1]), 0)
		@test isapprox(Sₗₘ(2, -1, [0, 0, -1]), 0)
		@test isapprox(Sₗₘ(2, -1, vec111), √(15 / 4π) * (1 / 3))
		@test isapprox(Sₗₘ(2, -1, -vec111), √(15 / 4π) * (1 / 3))
		@test isapprox(Sₗₘ(2, -1, [0, 1 / √(2), 1 / √(2)]), √(15 / 4π) * (1 / 2))
		@test isapprox(Sₗₘ(2, -1, [0, 1 / √(2), -1 / √(2)]), -√(15 / 4π) * (1 / 2))
		@test isapprox(Sₗₘ(2, -1, [0, -1 / √(2), 1 / √(2)]), -√(15 / 4π) * (1 / 2))
		@test isapprox(Sₗₘ(2, -1, [0, -1 / √(2), -1 / √(2)]), √(15 / 4π) * (1 / 2))
		# l = 2, m = 2, √(2)*√(45/96π)*(r̂x^2-r̂y^2)
		@test isapprox(Sₗₘ(2, 2, [1, 0, 0]), √(2) * √(45 / 96π))
		@test isapprox(Sₗₘ(2, 2, [-1, 0, 0]), √(2) * √(45 / 96π))
		@test isapprox(Sₗₘ(2, 2, [0, 1, 0]), -√(2) * √(45 / 96π))
		@test isapprox(Sₗₘ(2, 2, [0, -1, 0]), -√(2) * √(45 / 96π))
		@test isapprox(Sₗₘ(2, 2, [0, 0, 1]), 0)
		@test isapprox(Sₗₘ(2, 2, [0, 0, -1]), 0)
		@test isapprox(Sₗₘ(2, 2, vec111), 0)
		@test isapprox(Sₗₘ(2, 2, -vec111), 0)
		@test isapprox(Sₗₘ(2, 2, [1 / √(2), 0, 1 / √(2)]), √(2) * √(45 / 96π) * (1 / 2))
		@test isapprox(Sₗₘ(2, 2, [-1 / √(2), 0, 1 / √(2)]), √(2) * √(45 / 96π) * (1 / 2))
		@test isapprox(Sₗₘ(2, 2, [1 / √(2), 0, -1 / √(2)]), √(2) * √(45 / 96π) * (1 / 2))
		@test isapprox(Sₗₘ(2, 2, [0, 1 / √(2), 1 / √(2)]), -√(2) * √(45 / 96π) * (1 / 2))
		@test isapprox(Sₗₘ(2, 2, [0, 1 / √(2), -1 / √(2)]), -√(2) * √(45 / 96π) * (1 / 2))
		@test isapprox(Sₗₘ(2, 2, [1 / √(2), -1 / √(2), 0]), 0)
		@test isapprox(Sₗₘ(2, 2, [-1 / √(2), -1 / √(2), 0]), 0)
		# l = 2, m = -2, √(15/4π)*r̂x*r̂y
		@test isapprox(Sₗₘ(2, -2, [1, 0, 0]), 0)
		@test isapprox(Sₗₘ(2, -2, [-1, 0, 0]), 0)
		@test isapprox(Sₗₘ(2, -2, [0, 1, 0]), 0)
		@test isapprox(Sₗₘ(2, -2, [0, -1, 0]), 0)
		@test isapprox(Sₗₘ(2, -2, [0, 0, 1]), 0)
		@test isapprox(Sₗₘ(2, -2, [0, 0, -1]), 0)
		@test isapprox(Sₗₘ(2, -2, vec111), √(15 / 4π) / 3)
		@test isapprox(Sₗₘ(2, -2, -vec111), √(15 / 4π) / 3)
		@test isapprox(Sₗₘ(2, -2, [0, 1 / √(2), 1 / √(2)]), 0)
		@test isapprox(Sₗₘ(2, -2, [0, -1 / √(2), -1 / √(2)]), 0)
		@test isapprox(Sₗₘ(2, -2, [1 / √(2), 1 / √(2), 0]), √(15 / 4π) / 2)
		@test isapprox(Sₗₘ(2, -2, [-1 / √(2), 1 / √(2), 0]), -√(15 / 4π) / 2)
		@test isapprox(Sₗₘ(2, -2, [1 / √(2), -1 / √(2), 0]), -√(15 / 4π) / 2)
		@test isapprox(Sₗₘ(2, -2, [-1 / √(2), -1 / √(2), 0]), √(15 / 4π) / 2)
	end

	@testset "∂Sₗₘ_∂r̂x" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Sₗₘ_∂r̂x(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(0, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(0, 0, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(0, 0, -vec111), 0)

		# l = 1, m = 0, 0
		@test isapprox(∂Sₗₘ_∂r̂x(1, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, 0, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, 0, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, 0, -vec111), 0)
		# l = 1, m = 1, √(3/4π)
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, [1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, [-1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, [0, 1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, [0, -1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, [0, 0, 1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, [0, 0, -1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, vec111), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, -vec111), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(1, 1, [1 / √(2), 0, 1 / √(2)]), √(3 / 4π))
		# l = 1, m = -1, 0
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, -vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [1 / √(2), 0, 1 / √(2)]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [1 / √(2), 0, -1 / √(2)]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [-1 / √(2), 0, 1 / √(2)]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(1, -1, [-1 / √(2), 0, -1 / √(2)]), 0)

		# l = 2, m = 0, 0
		@test isapprox(∂Sₗₘ_∂r̂x(2, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 0, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 0, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 0, -vec111), 0)
		# l = 2, m = 1, √(15/4π)*r̂z
		@test isapprox(∂Sₗₘ_∂r̂x(2, 1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 1, [0, 0, 1]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(2, 1, [0, 0, -1]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(2, 1, vec111), √(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂x(2, 1, -vec111), -√(15 / 4π) / √(3))
		# l = 2, m = -1, 0
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, -vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [1 / √(2), 0, 1 / √(2)]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [1 / √(2), 0, -1 / √(2)]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [-1 / √(2), 0, 1 / √(2)]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -1, [-1 / √(2), 0, -1 / √(2)]), 0)
		# l = 2, m = 2, √(15/4π)*r̂x
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, [1, 0, 0]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, [-1, 0, 0]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, vec111), √(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, -vec111), -√(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, [1 / √(2), 0, 1 / √(2)]), √(15 / 4π) / √(2))
		@test isapprox(∂Sₗₘ_∂r̂x(2, 2, [-1 / √(2), 0, 1 / √(2)]), -√(15 / 4π) / √(2))
		# l = 2, m = -2, √(15/4π)*r̂y
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [0, 1, 0]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [0, -1, 0]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, vec111), √(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, -vec111), -√(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [1 / √(2), 0, 1 / √(2)]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [-1 / √(2), 0, 1 / √(2)]), 0)
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [1 / √(2), 1 / √(2), 0]), √(15 / 4π) / √(2))
		@test isapprox(∂Sₗₘ_∂r̂x(2, -2, [-1 / √(2), -1 / √(2), 0]), -√(15 / 4π) / √(2))
	end

	@testset "∂Sₗₘ_∂r̂y" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Sₗₘ_∂r̂y(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(0, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(0, 0, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(0, 0, -vec111), 0)

		# l = 1, m = 0, 0
		@test isapprox(∂Sₗₘ_∂r̂y(1, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 0, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 0, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 0, -vec111), 0)
		# l = 1, m = 1, 0
		@test isapprox(∂Sₗₘ_∂r̂y(1, 1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 1, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 1, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(1, 1, -vec111), 0)
		# l = 1, m = -1, √(3/4π)
		@test isapprox(∂Sₗₘ_∂r̂y(1, -1, [1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(1, -1, [-1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(1, -1, [0, 1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(1, -1, [0, -1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(1, -1, [0, 0, 1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(1, -1, [0, 0, -1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(1, -1, vec111), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(1, -1, -vec111), √(3 / 4π))

		# l = 2, m = 0, 0
		@test isapprox(∂Sₗₘ_∂r̂y(2, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 0, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 0, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 0, -vec111), 0)
		# l = 2, m = 1, 0
		@test isapprox(∂Sₗₘ_∂r̂y(2, 1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 1, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 1, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 1, -vec111), 0)
		# l = 2, m = -1, √((15/4π)*r̂z
		@test isapprox(∂Sₗₘ_∂r̂y(2, -1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, -1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, -1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, -1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, -1, [0, 0, 1]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(2, -1, [0, 0, -1]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(2, -1, vec111), √(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂y(2, -1, -vec111), -√(15 / 4π) / √(3))
		# l = 2, m = 2, -√((15/4π)*r̂y
		@test isapprox(∂Sₗₘ_∂r̂y(2, 2, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 2, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 2, [0, 1, 0]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(2, 2, [0, -1, 0]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(2, 2, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 2, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, 2, vec111), -√(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂y(2, 2, -vec111), √(15 / 4π) / √(3))
		# l = 2, m = -2, √((15/4π)*r̂x
		@test isapprox(∂Sₗₘ_∂r̂y(2, -2, [1, 0, 0]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(2, -2, [-1, 0, 0]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂y(2, -2, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, -2, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, -2, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, -2, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂y(2, -2, vec111), √(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂y(2, -2, -vec111), -√(15 / 4π) / √(3))
	end

	@testset "∂Sₗₘ_∂r̂z" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Sₗₘ_∂r̂z(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(0, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(0, 0, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(0, 0, -vec111), 0)

		# l = 1, m = 0, √(3/4π)
		@test isapprox(∂Sₗₘ_∂r̂z(1, 0, [1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(1, 0, [-1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(1, 0, [0, 1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(1, 0, [0, -1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(1, 0, [0, 0, 1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(1, 0, [0, 0, -1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(1, 0, vec111), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(1, 0, -vec111), √(3 / 4π))
		# l = 1, m = 1, 0
		@test isapprox(∂Sₗₘ_∂r̂z(1, 1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, 1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, 1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, 1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, 1, [0, 0, -1]), 0)
		# l = 1, m = -1, 0
		@test isapprox(∂Sₗₘ_∂r̂z(1, -1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, -1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, -1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, -1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, -1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(1, -1, [0, 0, -1]), 0)
		# l = 2, m = 0,  √(5/16π)*6r̂z
		@test isapprox(∂Sₗₘ_∂r̂z(2, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 0, [0, 0, 1]), √(5 / 16π) * 6)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 0, [0, 0, -1]), -√(5 / 16π) * 6)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 0, vec111), √(5 / 16π) * 6 / √(3))
		@test isapprox(∂Sₗₘ_∂r̂z(2, 0, -vec111), -√(5 / 16π) * 6 / √(3))
		# l = 2, m = 1, √(15/4π)*r̂x
		@test isapprox(∂Sₗₘ_∂r̂z(2, 1, [1, 0, 0]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(2, 1, [-1, 0, 0]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(2, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 1, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 1, vec111), √(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂z(2, 1, -vec111), -√(15 / 4π) / √(3))
		# l = 2, m = -1, √(15/4π)*r̂y
		@test isapprox(∂Sₗₘ_∂r̂z(2, -1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -1, [0, 1, 0]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(2, -1, [0, -1, 0]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂r̂z(2, -1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -1, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -1, vec111), √(15 / 4π) / √(3))
		@test isapprox(∂Sₗₘ_∂r̂z(2, -1, -vec111), -√(15 / 4π) / √(3))
		# l = 2, m = 2, 0
		@test isapprox(∂Sₗₘ_∂r̂z(2, 2, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 2, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 2, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 2, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 2, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 2, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 2, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, 2, -vec111), 0)
		# l = 2, m = -2, 0
		@test isapprox(∂Sₗₘ_∂r̂z(2, -2, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -2, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -2, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -2, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -2, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -2, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -2, vec111), 0)
		@test isapprox(∂Sₗₘ_∂r̂z(2, -2, -vec111), 0)
	end

	@testset "ssₗₘ" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(ssₗₘ(0, 0, [1, 0, 0]), 0)
		@test isapprox(ssₗₘ(0, 0, [0, 1, 0]), 0)
		@test isapprox(ssₗₘ(0, 0, [0, 0, 1]), 0)
		@test isapprox(ssₗₘ(0, 0, vec111), 0)
		@test isapprox(ssₗₘ(0, 0, -vec111), 0)

		# l = 1, m = 0, √(3/4π)*r̂z
		@test isapprox(ssₗₘ(1, 0, [1, 0, 0]), 0)
		@test isapprox(ssₗₘ(1, 0, [-1, 0, 0]), 0)
		@test isapprox(ssₗₘ(1, 0, [0, 1, 0]), 0)
		@test isapprox(ssₗₘ(1, 0, [0, -1, 0]), 0)
		@test isapprox(ssₗₘ(1, 0, [0, 0, 1]), √(3 / 4π))
		@test isapprox(ssₗₘ(1, 0, [0, 0, -1]), -√(3 / 4π))
		@test isapprox(ssₗₘ(1, 0, vec111), √(1 / 4π))
		@test isapprox(ssₗₘ(1, 0, -vec111), -√(1 / 4π))
		# l = 1, m = 1, √(3/4π)*r̂x
		@test isapprox(ssₗₘ(1, 1, [1, 0, 0]), √(3 / 4π))
		@test isapprox(ssₗₘ(1, 1, [-1, 0, 0]), -√(3 / 4π))
		@test isapprox(ssₗₘ(1, 1, [0, 1, 0]), 0)
		@test isapprox(ssₗₘ(1, 1, [0, -1, 0]), 0)
		@test isapprox(ssₗₘ(1, 1, [0, 0, 1]), 0)
		@test isapprox(ssₗₘ(1, 1, [0, 0, -1]), 0)
		@test isapprox(ssₗₘ(1, 1, vec111), √(1 / 4π))
		@test isapprox(ssₗₘ(1, 1, -vec111), -√(1 / 4π))
		# l = 1, m = -1, √(3/4π)*r̂y
		@test isapprox(ssₗₘ(1, -1, [1, 0, 0]), 0)
		@test isapprox(ssₗₘ(1, -1, [-1, 0, 0]), 0)
		@test isapprox(ssₗₘ(1, -1, [0, 1, 0]), √(3 / 4π))
		@test isapprox(ssₗₘ(1, -1, [0, -1, 0]), -√(3 / 4π))
		@test isapprox(ssₗₘ(1, -1, [0, 0, 1]), 0)
		@test isapprox(ssₗₘ(1, -1, [0, 0, -1]), 0)
		@test isapprox(ssₗₘ(1, -1, vec111), √(1 / 4π))
		@test isapprox(ssₗₘ(1, -1, -vec111), -√(1 / 4π))
		# l = 2, m = 0, √(5/16π)*6r̂z^2
		@test isapprox(ssₗₘ(2, 0, [1, 0, 0]), 0)
		@test isapprox(ssₗₘ(2, 0, [0, 1, 0]), 0)
		@test isapprox(ssₗₘ(2, 0, [0, 0, 1]), √(5 / 16π) * 6)
		@test isapprox(ssₗₘ(2, 0, [0, 0, -1]), √(5 / 16π) * 6)
		@test isapprox(ssₗₘ(2, 0, vec111), √(5 / 16π) * 6 / 3)
		@test isapprox(ssₗₘ(2, 0, -vec111), √(5 / 16π) * 6 / 3)
		# l = 2, m = 1, √(15/π)*r̂x*r̂z
		@test isapprox(ssₗₘ(2, 1, [1, 0, 0]), 0)
		@test isapprox(ssₗₘ(2, 1, [-1, 0, 0]), 0)
		@test isapprox(ssₗₘ(2, 1, [0, 1, 0]), 0)
		@test isapprox(ssₗₘ(2, 1, [0, -1, 0]), 0)
		@test isapprox(ssₗₘ(2, 1, [0, 0, 1]), 0)
		@test isapprox(ssₗₘ(2, 1, [0, 0, -1]), 0)
		@test isapprox(ssₗₘ(2, 1, vec111), √(15 / π) / 3)
		@test isapprox(ssₗₘ(2, 1, -vec111), √(15 / π) / 3)
		@test isapprox(ssₗₘ(2, 1, [1 / √2, 0, 1 / √2]), √(15 / π) / 2)
		@test isapprox(ssₗₘ(2, 1, [-1 / √2, 0, 1 / √2]), -√(15 / π) / 2)
		@test isapprox(ssₗₘ(2, 1, [1 / √2, 0, -1 / √2]), -√(15 / π) / 2)
		@test isapprox(ssₗₘ(2, 1, [-1 / √2, 0, -1 / √2]), √(15 / π) / 2)
		# l = 2, m = -1, √(15/π)*r̂y*r̂z
		@test isapprox(ssₗₘ(2, -1, [1, 0, 0]), 0)
		@test isapprox(ssₗₘ(2, -1, [-1, 0, 0]), 0)
		@test isapprox(ssₗₘ(2, -1, [0, 1, 0]), 0)
		@test isapprox(ssₗₘ(2, -1, [0, -1, 0]), 0)
		@test isapprox(ssₗₘ(2, -1, [0, 0, 1]), 0)
		@test isapprox(ssₗₘ(2, -1, [0, 0, -1]), 0)
		@test isapprox(ssₗₘ(2, -1, vec111), √(15 / π) / 3)
		@test isapprox(ssₗₘ(2, -1, -vec111), √(15 / π) / 3)
		@test isapprox(ssₗₘ(2, -1, [0, 1 / √2, -1 / √2]), -√(15 / π) / 2)
		# l = 2, m = 2, √(15/4π)*(r̂x^2-r̂y^2)
		@test isapprox(ssₗₘ(2, 2, [1, 0, 0]), √(15 / 4π))
		@test isapprox(ssₗₘ(2, 2, [-1, 0, 0]), √(15 / 4π))
		@test isapprox(ssₗₘ(2, 2, [0, 1, 0]), -√(15 / 4π))
		@test isapprox(ssₗₘ(2, 2, [0, -1, 0]), -√(15 / 4π))
		@test isapprox(ssₗₘ(2, 2, [0, 0, 1]), 0)
		@test isapprox(ssₗₘ(2, 2, [0, 0, -1]), 0)
		@test isapprox(ssₗₘ(2, 2, vec111), 0)
		@test isapprox(ssₗₘ(2, 2, -vec111), 0)
		@test isapprox(ssₗₘ(2, 2, [1 / √2, 0, 1 / √2]), √(15 / 4π) / 2)
		# l = 2, m = -2, √(15/π)*(r̂x*r̂y)
		@test isapprox(ssₗₘ(2, -2, [1, 0, 0]), 0)
		@test isapprox(ssₗₘ(2, -2, [0, 1, 0]), 0)
		@test isapprox(ssₗₘ(2, -2, [0, 0, 1]), 0)
		@test isapprox(ssₗₘ(2, -2, vec111), √(15 / π) / 3)
		@test isapprox(ssₗₘ(2, -2, -vec111), √(15 / π) / 3)
		@test isapprox(ssₗₘ(2, -2, [1 / √2, 1 / √2, 0]), √(15 / π) / 2)
		@test isapprox(ssₗₘ(2, -2, [-1 / √2, 1 / √2, 0]), -√(15 / π) / 2)
	end

	@testset "∂Sₗₘ_∂x" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Sₗₘ_∂x(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(0, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂x(0, 0, vec111), 0)

		# l = 1, m = 0, -√(3/4π)*r̂z*r̂x
		@test isapprox(∂Sₗₘ_∂x(1, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, 0, vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂x(1, 0, -vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [1 / √2, 0, 1 / √2]), -√(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [-1 / √2, 0, 1 / √2]), √(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [1 / √2, 0, -1 / √2]), √(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂x(1, 0, [-1 / √2, 0, -1 / √2]), -√(3 / 4π) / 2)
		# l = 1, m = 1, √(3/4π)*(1-r̂x^2)
		@test isapprox(∂Sₗₘ_∂x(1, 1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, 1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, 1, [0, 1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂x(1, 1, [0, -1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂x(1, 1, [0, 0, 1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂x(1, 1, [0, 0, -1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂x(1, 1, vec111), √(3 / 4π) * (1 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂x(1, 1, -vec111), √(3 / 4π) * (1 - 1 / 3))
		# l = 1, m = -1, -√(3/4π)*r̂x*r̂y
		@test isapprox(∂Sₗₘ_∂x(1, -1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [0, 0, -1]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, -1, vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂x(1, -1, -vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [1 / √2, 0, 1 / √2]), 0)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [1 / √2, 1 / √2, 0]), -√(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [-1 / √2, -1 / √2, 0]), -√(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂x(1, -1, [-1 / √2, 1 / √2, 0]), √(3 / 4π) / 2)
		# l = 2, m = 0, -6√(5/16π)*r̂x*r̂z^2
		@test isapprox(∂Sₗₘ_∂x(2, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(2, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(2, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂x(2, 0, vec111), -6 * √(5 / 16π) * (1 / 3√3))
		@test isapprox(∂Sₗₘ_∂x(2, 0, -vec111), 6 * √(5 / 16π) * (1 / 3√3))
		@test isapprox(∂Sₗₘ_∂x(2, 0, [1 / √2, 0, 1 / √2]), -6 * √(5 / 16π) * (1 / 2√2))
		@test isapprox(∂Sₗₘ_∂x(2, 0, [-1 / √2, 0, 1 / √2]), 6 * √(5 / 16π) * (1 / 2√2))
		# l = 2, m = 1,  √(15 / π)*r̂z*(1/√4-r̂x^2)
		@test isapprox(∂Sₗₘ_∂x(2, 1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(2, 1, [-1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(2, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(2, 1, [0, -1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂x(2, 1, [0, 0, 1]), √(15 / 4π))
		@test isapprox(∂Sₗₘ_∂x(2, 1, [0, 0, -1]), -√(15 / 4π))
		@test isapprox(∂Sₗₘ_∂x(2, 1, vec111), √(15 / π) / √3 * (1 / √4 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂x(2, 1, -vec111), -√(15 / π) / √3 * (1 / √4 - 1 / 3))
	end

	@testset "∂Sₗₘ_∂y" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Sₗₘ_∂y(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂y(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂y(0, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂y(0, 0, vec111), 0)
		# l = 1, m = 0, -√(3/4π)*r̂z*r̂y
		@test isapprox(∂Sₗₘ_∂y(1, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂y(1, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂y(1, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂y(1, 0, vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂y(1, 0, -vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂y(1, 0, [0, 1 / √2, 1 / √2]), -√(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂y(1, 0, [0, -1 / √2, 1 / √2]), √(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂y(1, 0, [0, 1 / √2, -1 / √2]), √(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂y(1, 0, [0, -1 / √2, -1 / √2]), -√(3 / 4π) / 2)
		# l = 1, m = 1, -√(3/4π)*r̂x*r̂y
		@test isapprox(∂Sₗₘ_∂y(1, 1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂y(1, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂y(1, 1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂y(1, 1, vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂y(1, 1, -vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂y(1, 1, [1 / √2, 1 / √2, 0]), -√(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂y(1, 1, [-1 / √2, -1 / √2, 0]), -√(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂y(1, 1, [1 / √2, -1 / √2, 0]), √(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂y(1, 1, [-1 / √2, 1 / √2, 0]), √(3 / 4π) / 2)
		# l = 1, m = -1, √(3/4π)*(1-r̂y^2)
		@test isapprox(∂Sₗₘ_∂y(1, -1, [1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂y(1, -1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂y(1, -1, [0, 0, 1]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂y(1, -1, vec111), √(3 / 4π) * (1 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂y(1, -1, -vec111), √(3 / 4π) * (1 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂y(1, -1, [1 / √2, 1 / √2, 0]), √(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂y(1, -1, [-1 / √2, -1 / √2, 0]), √(3 / 4π) / 2)
	end

	@testset "∂Sₗₘ_∂z" begin
		vec111 = 1 / √(3) * [1, 1, 1]
		# l = 0, m = 0, 0
		@test isapprox(∂Sₗₘ_∂z(0, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(0, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(0, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂z(0, 0, vec111), 0)

		# l = 1, m = 0, √(3/4π)*(1-r̂z^2)
		@test isapprox(∂Sₗₘ_∂z(1, 0, [1, 0, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂z(1, 0, [0, 1, 0]), √(3 / 4π))
		@test isapprox(∂Sₗₘ_∂z(1, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, 0, vec111), √(3 / 4π) * (1 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂z(1, 0, -vec111), √(3 / 4π) * (1 - 1 / 3))
		# l = 1, m = 1, -√(3/4π)*r̂x*r̂z
		@test isapprox(∂Sₗₘ_∂z(1, 1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, 1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, 1, vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂z(1, 1, -vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂z(1, 1, [1 / √2, 1 / √2, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, 1, [-1 / √2, 0, 1 / √2]), √(3 / 4π) / 2)
		# l = 1, m = -1, -√(3/4π)*r̂y*r̂z
		@test isapprox(∂Sₗₘ_∂z(1, -1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, -1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, -1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, -1, vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂z(1, -1, -vec111), -√(3 / 4π) / 3)
		@test isapprox(∂Sₗₘ_∂z(1, -1, [1 / √2, 1 / √2, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(1, -1, [0, 1 / √2, 1 / √2]), -√(3 / 4π) / 2)
		@test isapprox(∂Sₗₘ_∂z(1, -1, [0, -1 / √2, 1 / √2]), √(3 / 4π) / 2)
		# l = 2, m = 0, 3/2*√(5/π)*r̂z*(1-r̂z^2)
		@test isapprox(∂Sₗₘ_∂z(2, 0, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 0, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 0, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 0, vec111), 3 / 2 * √(5 / π) * (1 / √3) * (1 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂z(2, 0, -vec111), -3 / 2 * √(5 / π) * (1 / √3) * (1 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂z(2, 0, [0, 1 / √2, 1 / √2]), 3 / 2 * √(5 / π) * (1 / √2) * (1 - 1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, 0, [0, 1 / √2, -1 / √2]), -3 / 2 * √(5 / π) * (1 / √2) * (1 - 1 / 2))
		# l = 2, m = 1, √(15/π)*r̂x*(1/2-r̂z^2)
		@test isapprox(∂Sₗₘ_∂z(2, 1, [1, 0, 0]), √(15 / π) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, 1, [-1, 0, 0]), -√(15 / π) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, 1, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 1, vec111), √(15 / π) * (1 / √3) * (1/2 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂z(2, 1, -vec111), -√(15 / π) * (1 / √3) * (1/2 - 1 / 3))
		# l = 2, m = -1, √(15/π)*r̂y*(1/2-r̂z^2)
		@test isapprox(∂Sₗₘ_∂z(2, -1, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, -1, [0, 1, 0]), √(15 / π) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, -1, [0, -1, 0]), -√(15 / π) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, -1, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, -1, vec111), √(15 / π) * (1 / √3) * (1/2 - 1 / 3))
		@test isapprox(∂Sₗₘ_∂z(2, -1, -vec111), -√(15 / π) * (1 / √3) * (1/2 - 1 / 3))
		# l = 2, m = 2, -√(15/4π)*r̂z*(r̂x^2-r̂y^2)
		@test isapprox(∂Sₗₘ_∂z(2, 2, [1, 0, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 2, [0, 1, 0]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 2, [0, 0, 1]), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 2, vec111), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 2, -vec111), 0)
		@test isapprox(∂Sₗₘ_∂z(2, 2, [1 / √2, 0, 1 / √2]), -√(15 / 4π) * (1 / √2) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, 2, [1 / √2, 0, -1 / √2]), √(15 / 4π) * (1 / √2) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, 2, [-1 / √2, 0, 1 / √2]), -√(15 / 4π) * (1 / √2) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, 2, [0, 1 / √2, 1 / √2]), √(15 / 4π) * (1 / √2) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, 2, [0, -1 / √2, 1 / √2]), √(15 / 4π) * (1 / √2) * (1 / 2))
		@test isapprox(∂Sₗₘ_∂z(2, 2, [0, 1 / √2, -1 / √2]), -√(15 / 4π) * (1 / √2) * (1 / 2))
	end
end
