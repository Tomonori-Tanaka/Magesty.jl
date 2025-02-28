"""
	module MySphericalHarmonics

This module provides functions to compute spherical harmonics ( Y_{l,m} ) and related derivatives, following the formalism described in Drautz (Phys. Rev. B 102, 024104, 2020). It includes:

- Normalized associated Legendre polynomials ( bar{P}_{l,m} ).
- Spherical harmonics ( Y_{l,m} ) and real-valued ( S_{l,m} ).
- Partial derivatives with respect to Cartesian coordinates.

# Functions
- `P̄ₗₘ(l, m, r̂z)`: Compute the normalized associated Legendre polynomial.
- `Yₗₘ(l, m, uvec)`: Compute the spherical harmonic.
- `Sₗₘ(l, m, uvec)`: Compute the real-valued spherical harmonic.
- `∂Sₗₘ_∂r̂x(l, m, uvec)`, `∂Sₗₘ_∂r̂y(l, m, uvec)`, `∂Sₗₘ_∂r̂z(l, m, uvec)`: Partial derivatives of ( S_{l,m} ).
- `∂ᵢSlm(l, m, uvec)`: Compute the gradient of ( S_{l,m} ) as a vector.
"""
module MySphericalHarmonics

using LegendrePolynomials
using LinearAlgebra

# abstract type SphericalHarmonicsProduct end
export Sₗₘ, d_Slm, ∂ᵢSlm

"""
Note that dnPl in LegendrePolynomials package includes the Condon-Shortley Phase ( (-1)^m ).

ref)
R. Drautz, Phys. Rev. B 102, 024104 (2020).
"""
function P̄ₗₘ(l::Int, m::Int, r̂z::Real)::Float64
	return (-1)^m * √((2l + 1) / (4π) * factorial(l - m) / factorial(l + m)) *
		   dnPl(r̂z, l, m)
end

function P̄ₗₘ(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	return P̄ₗₘ(l, m, uvec[3])
end

"""
dP̄ₗₘ/dr̂z
"""
function dP̄ₗₘ(l::Int, m::Int, r̂z::Real)::Float64
	return (-1)^m * √((2l + 1) / (4π) * factorial(l - m) / factorial(l + m)) *
		   dnPl(r̂z, l, m + 1)
end

function dP̄ₗₘ(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return (-1)^m * √((2l + 1) / (4π) * factorial(l - m) / factorial(l + m)) *
		   dnPl(uvec[3], l, m + 1)
end

"""
equation (D18) in R. Drautz, Phys. Rev. B 102, 024104 (2020).
"""
function ∂Yₗₘ_∂r̂x(l::Int, m::Int, uvec::Vector{<:Real})::Complex
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return m * (uvec[1] + uvec[2] * im)^(m - 1) * P̄ₗₘ(l, m, uvec[3])
end

"""
equation (D19) in R. Drautz, Phys. Rev. B 102, 024104 (2020).
"""
function ∂Yₗₘ_∂r̂y(l::Int, m::Int, uvec::Vector{<:Real})::Complex
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return 1.0im * ∂Yₗₘ_∂r̂x(l, m, uvec)
end

"""
equation (D20) in R. Drautz, Phys. Rev. B 102, 024104 (2020).
"""
function ∂Yₗₘ_∂r̂z(l::Int, m::Int, uvec::Vector{<:Real})::Complex
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return (uvec[1] + uvec[2] * im)^m * dnPl(l, m + 1, uvec[3])
end

"""
ref)
equation (D22) in R. Drautz, Phys. Rev. B 102, 024104 (2020).
"""
function yₗₘ(l::Int, m::Int, uvec::Vector{<:Real})::Complex
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return uvec[1] * ∂Yₗₘ_∂r̂x(l, m, uvec) + uvec[2] * ∂Yₗₘ_∂r̂y(l, m, uvec) +
		   uvec[3] * ∂Yₗₘ_∂r̂z(l, m, uvec)
end

function Yₗₘ(l::Int, m::Int, uvec::Vector{<:Real})::Complex
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	if m < 0
		return (-1)^(abs(m)) * conj(uvec[1] + uvec[2] * im)^abs(m) *
			   P̄ₗₘ(l, abs(m), uvec[3])
	else
		return (uvec[1] + uvec[2] * im)^m * P̄ₗₘ(l, m, uvec[3])
	end
end

function Sₗₘ(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if l < 0
		error("Invalid value for l: l must be positive.")
	end

	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	if m == 0
		return Real(Yₗₘ(l, m, uvec))
	elseif m > 0
		return (-1)^m * √2 * P̄ₗₘ(l, m, uvec[3]) *
			   sum(
				   (-1)^k * binomial(m, 2k) * uvec[1]^(m - 2k) * uvec[2]^(2k) for
				   k in 0:floor(Int, m / 2)
			   )
	elseif m < 0
		n = abs(m)
		return (-1)^n * √2 * P̄ₗₘ(l, n, uvec[3]) *
			   sum(
				   (-1)^k * binomial(n, 2k + 1) * uvec[1]^(n - (2k + 1)) *
				   uvec[2]^(2k + 1) for k in 0:floor(Int, (n - 1) / 2)
			   )
	end
end

function ∂Sₗₘ_∂r̂x(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	if m == 0
		return 0.0
	elseif m > 0
		return (-1)^m * √2 * m * P̄ₗₘ(l, m, uvec[3]) *
			   sum(
				   (-1)^k * binomial(m - 1, 2k) * uvec[1]^(m - 1 - 2k) * uvec[2]^(2k) for
				   k in 0:floor(Int, (m - 1) / 2)
			   )
	else # m < 0
		n = abs(m)
		if n == 1
			return 0.0
		else
			return (-1)^n * √2 * n * P̄ₗₘ(l, n, uvec[3]) *
				   sum(
					   (-1)^k * binomial(n - 1, 2k + 1) * uvec[1]^(n - 2 - 2k) *
					   uvec[2]^(2k + 1) for
					   k in 0:floor(Int, (n - 2) / 2)
				   )
		end
	end
end

function ∂Sₗₘ_∂r̂y(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	if m == 0
		return 0.0
	elseif m > 0
		return -∂Sₗₘ_∂r̂x(l, -m, uvec)
	else # m < 0
		return ∂Sₗₘ_∂r̂x(l, abs(m), uvec)
	end
end

function ∂Sₗₘ_∂r̂z(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	if m == 0
		return dP̄ₗₘ(l, m, uvec[3])
	elseif m > 0
		return (-1)^m * √2 * dP̄ₗₘ(l, m, uvec[3]) *
			   sum(
				   (-1)^k * binomial(m, 2k) * uvec[1]^(m - 2k) * uvec[2]^(2k) for
				   k in 0:floor(Int, m / 2)
			   )
	else# m < 0
		n = abs(m)
		return (-1)^n * √2 * dP̄ₗₘ(l, n, uvec[3]) *
			   sum(
				   (-1)^k * binomial(n, 2k + 1) * uvec[1]^(n - (2k + 1)) *
				   uvec[2]^(2k + 1) for k in 0:floor(Int, (n - 1) / 2)
			   )
	end
end

# ss imply "small s"
function ssₗₘ(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return uvec[1] * ∂Sₗₘ_∂r̂x(l, m, uvec) +
		   uvec[2] * ∂Sₗₘ_∂r̂y(l, m, uvec) +
		   uvec[3] * ∂Sₗₘ_∂r̂z(l, m, uvec)
end

function ∂Sₗₘ_∂x(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return ∂Sₗₘ_∂r̂x(l, m, uvec) - uvec[1] * ssₗₘ(l, m, uvec)
end

function ∂Sₗₘ_∂y(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return ∂Sₗₘ_∂r̂y(l, m, uvec) - uvec[2] * ssₗₘ(l, m, uvec)
end

function ∂Sₗₘ_∂z(l::Int, m::Int, uvec::Vector{<:Real})::Float64
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return ∂Sₗₘ_∂r̂z(l, m, uvec) - uvec[3] * ssₗₘ(l, m, uvec)
end

d_Slm = [∂Sₗₘ_∂x, ∂Sₗₘ_∂y, ∂Sₗₘ_∂z]

function ∂ᵢSlm(l::Int, m::Int, uvec::Vector{<:Real})::Vector{Float64}
	if !is_normalized(uvec)
		error("uvec is not normalized")
	end

	return [∂Sₗₘ_∂x(l, m, uvec), ∂Sₗₘ_∂y(l, m, uvec), ∂Sₗₘ_∂z(l, m, uvec)]
end

function is_normalized(vec::AbstractVector{<:Real}, tol::Real = 1e-8)::Bool
	if length(vec) != 3
		error("vec must be a 3-element vector")
	end

	norm_vec = norm(vec)
	if isapprox(norm_vec, 1.0, atol = tol)
		return true
	else
		return false
	end
end

end