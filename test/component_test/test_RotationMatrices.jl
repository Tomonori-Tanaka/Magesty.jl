using Rotations
using WignerD
using ..UnitaryMatrixCl
using ..RotationMatrices


@testset "RotationMatrices" begin
	Rz(α) = [cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
	Ry(β) = [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)]

	Δ1(a, b, c) =
		[
			cos(a)*cos(c)-sin(a)*sin(c)*cos(b)   sin(a)*sin(b) cos(a)*sin(c)+sin(a)*cos(c)*cos(b)
											  sin(c)*sin(b) cos(b) -cos(c)sin(b)
											  -cos(a)*sin(c)*cos(b)-sin(a)*cos(c) cos(a)*sin(b) cos(a)*cos(c)*cos(b)-sin(a)*sin(c)
		]

	for _ in 1:20
		α = 2π * rand()
		β = π * rand()
		γ = 2π * rand()

		R1 = Rz(α)
		R2 = Ry(β)
		R3 = Rz(γ)
		rotmat = R1 * R2 * R3
		rotmat_ = RotZYZ(α, β, γ)
		@test α ≈ rotmat2euler(rotmat)[1]
		@test β ≈ rotmat2euler(rotmat)[2]
		@test γ ≈ rotmat2euler(rotmat)[3]
	end

	rotmat = [1 0 0; 0 1 0; 0 0 -1]
	rotmat = Rz(1)

	rotmat2 = [cos(π / 2) -sin(π / 2) 0; sin(π / 2) cos(π / 2) 0; 0 0 1]
	# rotmat2 = [cos(π) -sin(π) 0; sin(π) cos(π) 0; 0 0 1]
	# rotmat2 = [cos(π / 4) -sin(π / 4) 0; sin(π / 4) cos(π / 4) 0; 0 0 1]
	# rotmat2 = Rz(0)
	# display(Rz(0))

	#delta = Δl(1, π/2, 0, 0)
	#@show delta * [1, 0, 0]
	#@show delta * [0, 1, 0]
	#@show delta * [0, 0, 1]
end
