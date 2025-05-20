using Rotations
using WignerD
using ..UnitaryMatrixCl
using ..RotationMatrix


@testset "RotationMatrix" begin

	Δ1(a, b, c) =
		[
			cos(a)*cos(c)-sin(a)*sin(c)*cos(b)   sin(a)*sin(b) cos(a)*sin(c)+sin(a)*cos(c)*cos(b)
											  sin(c)*sin(b) cos(b) -cos(c)*sin(b)
											  -cos(a)*sin(c)*cos(b)-sin(a)*cos(c) cos(a)*sin(b) cos(a)*cos(c)*cos(b)-sin(a)*sin(c)
		]

	for _ in 1:20
		α = 2π * rand() - π# -π ≤ α ≤ π 
		β = π * rand()
		γ = 2π * rand() - π

		#= 		R1 = RotZ(α)
				R2 = RotY(β)
				R3 = RotZ(γ)
				rotmat = R1 * R2 * R3 =#

		rotmat = RotZYZ(α, β, γ)
		@test all(rotmat2euler(rotmat) .≈ (α, β, γ))
		@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)
	end

	# 90 deg. rotation around z-axis
	rotmat = RotZ(deg2rad(90))
	@test all(rad2deg.(rotmat2euler(rotmat, false)) .≈ (90, 0, 0))
	@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)

	# -90 deg. rotation around z-axis
	rotmat = RotZ(deg2rad(-90))
	@test all(rad2deg.(rotmat2euler(rotmat, false)) .≈ (-90, 0, 0))
	@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)

	# 90 deg. rotation around y-axis
	rotmat = RotY(deg2rad(90))
	@test all(rad2deg.(rotmat2euler(rotmat, false)) .≈ (0, 90, 0))
	@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)

	# -90 deg. rotation around y-axis
	# Note that returned β is positive because of the definition of the Eular angles (0 ≤ β ≤ π).
	rotmat = RotY(deg2rad(-90))
	@test all(rad2deg.(rotmat2euler(rotmat, false)) .≈ (180, 90, 180))
	@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)

	# 180 deg. rotation around y-axis
	rotmat = RotY(deg2rad(180))
	@test all(rad2deg.(rotmat2euler(rotmat, false)) .≈ (0, 180, 0))
	@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)

	# 90 deg. rotation around x-axis
	rotmat = RotX(deg2rad(90))
	@test all(rad2deg.(rotmat2euler(rotmat, false)) .≈ (-90, 90, 90))
	@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)

	# 180 deg. rotation around x-axis
	rotmat = RotX(deg2rad(180))
	@test all(rad2deg.(rotmat2euler(rotmat)) .≈ (-180, 180, 0))
	@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)

	# -180 deg. rotation around x-axis
	rotmat = RotX(deg2rad(-180))
	@test all(rad2deg.(rotmat2euler(rotmat)) .≈ (-180, 180, 0))
	@test Δ1(rotmat2euler(rotmat)...) ≈ Δl(1, rotmat2euler(rotmat)...)
end
