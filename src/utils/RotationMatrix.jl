module RotationMatrices

using LinearAlgebra
using WignerD

using ..UnitaryMatrixCl

export rotmat2euler, Δl

"""
	rotmat2euler(m::AbstractMatrix{<:Real}, mod_positive::Bool = true) -> Tuple{Float64, Float64, Float64}

Converts a 3x3 rotation matrix `m` to Euler angles `(α, β, γ)`.

# Arguments
- `m::AbstractMatrix{<:Real}`: A 3x3 rotation matrix.
- `mod_positive::Bool`: If `true`, the angles are adjusted to be within `[0, 2π)` for `α` and `γ`. Defaults to `true`.

# Returns
- `Tuple{Float64, Float64, Float64}`: The Euler angles `(α, β, γ)`.

# Raises
- `ArgumentError`: If the input matrix `m` is not 3x3.

# Examples
```julia
m = [1 0 0; 0 cos(pi/4) -sin(pi/4); 0 sin(pi/4) cos(pi/4)]
angles = rotmat2euler(m)
println(angles)  # (0.0, 0.7853981633974483, 0.0)
"""
function rotmat2euler(
	m::AbstractMatrix{<:Real},
	mod_positive::Bool = true;
	tol::Real = 1e-12,
)::Tuple{Float64, Float64, Float64}
	if size(m) != (3, 3)
		throw(ArgumentError("Only 3x3 matrices are allowed"))
	end

	β = acos(m[3, 3])# 0 ≤ β ≤ π

	if isapprox(β, 0.0, atol = tol)
		γ = 0
		α = acos(m[1, 1])
	elseif isapprox(β, π, atol = tol)
		γ = 0
		α = acos(-m[1, 1])
	else
		α = atan(m[2, 3], m[1, 3])
		γ = atan(m[3, 2], -m[3, 1])
	end

	if mod_positive
		α = mod(α, 2π)
		γ = mod(γ, 2π)
	end

	return (α, β, γ)
end


"""
	Δl(l::Int, α::Float64, β::Float64, γ::Float64)::Matrix{Float64}

Compute the Δ matrix for a given `l`, `α`, `β`, and `γ`.

# Arguments
- `l::Int`: Positive integer representing the angular momentum quantum number.
- `α::Float64`, `β::Float64`, `γ::Float64`: Euler angles in radians.

# Returns
- A real-valued matrix `Δ::Matrix{Float64}`.

# Throws
- `ArgumentError` if `l` is not positive.
- `ArgumentError` if the resulting matrix contains imaginary parts exceeding the threshold `1e-12`.

# Notes
- The resulting matrix is converted to real if the imaginary parts are negligible.

# References
- M.A. Blanco et al., Journal of Molecular Structure (Theochem) 419 19-27 (1997).
"""
function Δl(l::Integer, α::Real, β::Real, γ::Real; tol::Real = 1e-12)::Matrix{Float64}
	# Ensure l is positive
	if l < 0
		throw(ArgumentError("Only positive l is allowed. Given: $l"))
	end

	cl_mat = UniMatCl(l)# Assuming UniMatCl returns a Matrix{Complex}
	wigD::Matrix{Complex} = wignerD(l, α, β, γ)# Assuming wignerD returns a Matrix{Complex}

	# Compute Δ matrix (eq. (35) in the reference)
	Δ = conj(cl_mat) * wigD * transpose(cl_mat)

	# Check for excessive imaginary parts
	if any(abs.(imag.(Δ)) .> tol)
		throw(
			ArgumentError("Matrix contains imaginary parts exceeding the threshold 1e-12"),
		)
	end

	realΔ = real(Δ)

	if !(is_orthogonal(realΔ))
		error("the rotation matrix is not orthogonal")
	end

	return real(Δ)
end

function is_orthogonal(Q::AbstractMatrix{<:Real}; tol::Float64 = 1e-10)::Bool
	n = size(Q, 1)
	return (norm(Q' * Q - I(n)) < tol)
end

end
