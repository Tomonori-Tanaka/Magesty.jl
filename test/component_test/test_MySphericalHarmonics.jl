using Test
using LinearAlgebra

using ..MySphericalHarmonics

@testset "Y_lm analytical" begin
	@test isapprox(MySphericalHarmonics.Yₗₘ(0, 0, [1, 0, 0]), 1 / sqrt(4π))
	@test isapprox(MySphericalHarmonics.Yₗₘ(1, 1, [-1, 0, 0]), sqrt(3 / 8π))
	@test isapprox(MySphericalHarmonics.Yₗₘ(1, -1, [-1, 0, 0]), -sqrt(3 / 8π))
	@test isapprox(MySphericalHarmonics.Yₗₘ(1, 1, [0, 0, 1]), 0)
	@test isapprox(MySphericalHarmonics.Yₗₘ(2, 2, [0, 1, 0]), -3 * sqrt(5 / 96π))
	@test isapprox(MySphericalHarmonics.Yₗₘ(2, 0, [0, 0, 1]), sqrt(5 / 4π))
end

@testset "S_lm analytical" begin
	@test isapprox(MySphericalHarmonics.Sₗₘ(0, 0, [1, 0, 0]), 1 / sqrt(4π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, -1, [0, 0, 1]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, -1, [0, 0, -1]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, -1, [0, 1, 0]), 1 / 2 * sqrt(3 / π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, -1, [0, -1, 0]), -1 / 2 * sqrt(3 / π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, -1, [1, 0, 0]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, -1, [-1, 0, 0]), 0)
	@test isapprox(
		MySphericalHarmonics.Sₗₘ(1, -1, 1 / sqrt(2) * [0, 1, 1]),
		1 / 2 * sqrt(3 / (2π)),
	)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, -1, 1 / sqrt(2) * [1, 0, 1]), 0)
	@test isapprox(
		MySphericalHarmonics.Sₗₘ(1, -1, 1 / sqrt(2) * [1, 1, 0]),
		1 / 2 * sqrt(3 / (2π)),
	)

	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, [0, 0, 1]), sqrt(3 / 4π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, [0, 0, -1]), -sqrt(3 / 4π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, [0, 1, 0]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, [0, -1, 0]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, [1, 0, 0]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, [-1, 0, 0]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, 1 / sqrt(2) * [0, 1, 1]), sqrt(3 / 8π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, 1 / sqrt(2) * [0, 1, -1]), -sqrt(3 / 8π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, 1 / sqrt(2) * [1, 0, 1]), sqrt(3 / 8π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 0, 1 / sqrt(2) * [1, 1, 0]), 0)

	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, [0, 0, 1]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, [0, 0, -1]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, [0, 1, 0]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, [0, -1, 0]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, [1, 0, 0]), 1 / 2 * sqrt(3 / π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, [-1, 0, 0]), -1 / 2 * sqrt(3 / π))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, 1 / sqrt(2) * [0, 1, 1]), 0)
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, 1 / sqrt(2) * [1, 0, 1]), 1 / 2 * sqrt(3 / (2π)))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, 1 / sqrt(2) * [1, 0, -1]), 1 / 2 * sqrt(3 / (2π)))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, 1 / sqrt(2) * [1, 1, 0]), 1 / 2 * sqrt(3 / (2π)))
	@test isapprox(MySphericalHarmonics.Sₗₘ(1, 1, 1 / sqrt(2) * [-1, 1, 0]), -1 / 2 * sqrt(3 / (2π)))
	
end

@testset "Real Spherical Harmonics (Slm)" begin

	# uvec_list = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1 / sqrt(2), 1 / sqrt(2), 0], [0, 0, -1]]
	uvec_list = [randn(3) for _ in 1:10]
	for (i, uvec) in enumerate(uvec_list)
		uvec_list[i] = uvec / norm(uvec)
	end
	result = 0 + 0im
	for l in 1:10
		for m in -l:l
			for uvec in uvec_list
				if m < 0
					n = abs(m)
					result =
						(-1)^n / (sqrt(2) * im) * (
							MySphericalHarmonics.Yₗₘ(l, n, uvec) -
							conj(MySphericalHarmonics.Yₗₘ(l, n, uvec))
						)
				elseif m == 0
					result = MySphericalHarmonics.Yₗₘ(l, m, uvec)
				else
					result =
						(-1)^m / sqrt(2) * (
							MySphericalHarmonics.Yₗₘ(l, m, uvec) +
							conj(MySphericalHarmonics.Yₗₘ(l, m, uvec))
						)
				end
				# println(l, m, uvec)
				@test isapprox(Sₗₘ(l, m, uvec), result, atol = 1e-10)
			end
		end
	end

	# test for derivatives
	for uvec in uvec_list
		# _∂r̂x
		# l = 1
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂x(1, -1, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂x(1, 0, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂x(1, 1, uvec) ≈ √(3 / 4π)
		# l = 2
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂x(2, -2, uvec) ≈ 6 * √(5 / 48π) * uvec[2]
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂x(2, -1, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂x(2, 0, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂x(2, 1, uvec) ≈ 3 * √(5 / 12π) * uvec[3]
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂x(2, 2, uvec) ≈ 6 * √(5 / 48π) * uvec[1]

		# _∂r̂y
		# l = 1
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂y(1, -1, uvec) ≈ √(3 / 4π)
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂y(1, 0, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂y(1, 1, uvec) ≈ 0
		# l = 2
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂y(2, -1, uvec) ≈ 3 * √(5 / 12π) * uvec[3]
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂y(2, 0, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂y(2, 1, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂y(2, 2, uvec) ≈ -6 * √(5 / 48π) * uvec[2]

		# _∂r̂z
		# l = 1
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂z(1, -1, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂z(1, 0, uvec) ≈ √(3 / 4π)
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂z(1, 1, uvec) ≈ 0
		# l = 2
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂z(2, -2, uvec) ≈ 0
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂z(2, -1, uvec) ≈ 3 * √(5 / 12π) * uvec[2]
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂z(2, 0, uvec) ≈ 3 * √(5 / 4π) * uvec[3]
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂z(2, 1, uvec) ≈ 3 * √(5 / 12π) * uvec[1]
		@test MySphericalHarmonics.∂Sₗₘ_∂r̂z(2, 2, uvec) ≈ 0

		@test d_Slm[1](1, 1, uvec) ≈ √(3 / 4π) - uvec[1] * (√(3 / 4π) * uvec[1] + 0 + 0)
		@test d_Slm[1](1, 0, uvec) ≈ 0 - uvec[1] * (0 + 0 + √(3 / 4π) * uvec[3])

	end
end
