"""
	module MySphericalHarmonics

This module provides functions to compute spherical harmonics ( Y_{l,m} ) and related derivatives, following the formalism described in Drautz (Phys. Rev. B 102, 024104, 2020). It includes:

- Normalized associated Legendre polynomials ( bar{P}_{l,m} ).
- Spherical harmonics ( Y_{l,m} ) and tesseral harmonics ( Z_{l,m} ).
- Partial derivatives with respect to Cartesian coordinates.

# Functions
- `P̄ₗₘ(l, m, r̂z)`: Compute the normalized associated Legendre polynomial.
- `Yₗₘ(l, m, uvec)`: Compute the spherical harmonic.
- `Zₗₘ(l, m, uvec)`: Compute the tesseral harmonic.
- `∂Zₗₘ_∂r̂x(l, m, uvec)`, `∂Zₗₘ_∂r̂y(l, m, uvec)`, `∂Zₗₘ_∂r̂z(l, m, uvec)`: Partial derivatives of ( Z_{l,m} ).
- `∂ᵢZlm(l, m, uvec)`: Compute the gradient of ( Z_{l,m} ) as a vector.
"""
module MySphericalHarmonics

using LegendrePolynomials
using LinearAlgebra

# abstract type SphericalHarmonicsProduct end
export Zₗₘ, d_Zlm, ∂ᵢZlm

"""
	P̄ₗₘ(l::Integer, m::Integer, r̂z::Real) -> Float64
	P̄ₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Compute P̄ₗₘ (Drautz's convention).

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `r̂z`: Cosine of polar angle, cos(θ) ∈ [-1,1]
- `uvec`: Alternative input as normalized 3D direction vector

# Mathematical Details
Drautz's convention:
P̄ₗₘ = (-1)ᵐ √((2l+1)/(4π) * (l-m)!/(l+m)!) * d^m/d(r̂z)^m Pₗ

where Pₗₘ includes the Condon-Shortley phase from LegendrePolynomials.jl.

# Notes
- The function automatically handles the Condon-Shortley phase
- For vector input, uses the z-component (uvec[3]) as r̂z
- Caller must ensure valid `(l, m)`, `r̂z ∈ [-1, 1]`, and for the vector method a 3D normalized `uvec` (no validation; invalid input may yield incorrect results or errors from dependencies)

# References
- R. Drautz, Phys. Rev. B 102, 024104 (2020)

# Examples
```julia
P̄ₗₘ(2, 1, 0.5)              # Scalar input
P̄ₗₘ(2, 1, [0.0, 0.0, 1.0])  # Vector input
```
"""
function P̄ₗₘ(l::Integer, m::Integer, r̂z::Real)::Float64
	const_factor = √((2l + 1) / (4π))
	factorial_ratio = √(factorial(l - abs(m)) / factorial(l + abs(m)))
	phase = (-1)^m

	return phase * const_factor * factorial_ratio * dnPl(r̂z, l, m)
end

function P̄ₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	return P̄ₗₘ(l, m, uvec[3])
end

"""
	dP̄ₗₘ(l::Integer, m::Integer, r̂z::Real) -> Float64
	dP̄ₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Compute the derivative of P̄ₗₘ with respect to r̂z (cos θ).

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `r̂z`: Cosine of polar angle, cos(θ) ∈ [-1,1]
- `uvec`: Alternative input as normalized 3D direction vector

# Returns
- Value of dP̄ₗₘ/dr̂z

# Mathematical Details
The derivative is computed using the relationship:
dP̄ₗₘ/dr̂z = (-1)ᵐ √((2l+1)/(4π) * (l-m)!/(l+m)!) * d^(m+1)/d(r̂z)^(m+1) Pₗ

where dPₗₘ/dr̂z is computed using the LegendrePolynomials.jl package.

# Notes
- Uses the recurrence relation for associated Legendre polynomials
- For vector input, uses the z-component (uvec[3]) as r̂z

# Examples
```julia
dP̄ₗₘ(2, 1, 0.5)              # Scalar input
dP̄ₗₘ(2, 1, [0.0, 0.0, 1.0])  # Vector input
```
"""
function dP̄ₗₘ(l::Integer, m::Integer, r̂z::Real)::Float64
	validate_lm(l, m)
	validate_r̂z(r̂z)
	normalization = √((2l + 1) / (4π) * factorial(l - abs(m)) / factorial(l + abs(m)))
	phase = (-1)^m

	return phase * normalization * dnPl(r̂z, l, m + 1)
end

function dP̄ₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	return dP̄ₗₘ(l, m, uvec[3])
end


"""
	Yₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Compute the complex spherical harmonic Yₗₘ(θ,φ).

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `uvec`: Normalized 3D direction vector [x, y, z]

# Returns
- Complex spherical harmonic value Yₗₘ(θ,φ)

# Mathematical Details
Yₗₘ(θ,φ) = P̄ₗₘ(cos θ) exp(imφ)(sin θ)^m, subject to m ≥ 0
where:
- For m ≥ 0: exp(imφ)*(sin θ)^m = (x + iy)ᵐ
- For m < 0: Yₗ,-ₘ = (-1)ᵐ Yₗₘ*

# References
- R. Drautz, Phys. Rev. B 102, 024104 (2020)


# Examples
```julia
# z-axis
Yₗₘ(1, 0, [0.0, 0.0, 1.0])

# x-axis
Yₗₘ(1, 1, [1.0, 0.0, 0.0])
```
"""
function Yₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex
	validate_lm(l, m)
	validate_uvec(uvec)

	if m < 0
		return (-1)^(abs(m)) * conj(uvec[1] + uvec[2] * im)^abs(m) *
			   P̄ₗₘ(l, abs(m), uvec[3])
	else
		return (uvec[1] + uvec[2] * im)^m * P̄ₗₘ(l, m, uvec[3])
	end
end

"""
	∂Yₗₘ_∂r̂x(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Compute the partial derivative of Yₗₘ with respect to r̂x.

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `uvec`: Normalized 3D direction vector [r̂x, r̂y, r̂z]

# Returns
- Complex value of ∂Yₗₘ/∂r̂x

# Mathematical Details
For m > 0:  ∂Yₗₘ/∂r̂x = m(r̂x + iy)^(m-1) P̄ₗₘ(z)
For m < 0:  ∂Yₗₘ/∂r̂x = (-1)^|m| |m|(r̂x - iy)^(|m|-1) P̄ₗₘ(z)
For m = 0:  ∂Yₗₘ/∂r̂x = 0

where:
- n = |m|
- P̄ₗₘ is the normalized associated Legendre polynomial

Reference: Equation (D18) in R. Drautz, Phys. Rev. B 102, 024104 (2020)
"""
function ∂Yₗₘ_∂r̂x(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	validate_lm(l, m)
	validate_uvec(uvec)

	# Early return for m = 0
	m == 0 && return Complex{Float64}(0)

	# Common calculations
	n = abs(m)
	plm = P̄ₗₘ(l, n, uvec[3])

	# Complex coordinate z = x ± iy
	z_xy = if m < 0
		ComplexF64(uvec[1], -uvec[2])  # x - iy
	else
		ComplexF64(uvec[1], uvec[2])   # x + iy
	end

	# Compute derivative with phase factor
	return (m < 0 ? (-1)^n * n : m) * z_xy^(n - 1) * plm
end

"""
	∂Yₗₘ_∂r̂y(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Compute the partial derivative of Yₗₘ with respect to r̂y.

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `uvec`: Normalized 3D direction vector [r̂x, r̂y, r̂z]

# Returns
- Complex value of ∂Yₗₘ/∂r̂y

# Mathematical Details
For m > 0:  ∂Yₗₘ/∂r̂y = +i ∂Yₗₘ/∂r̂x
For m < 0:  ∂Yₗₘ/∂r̂y = -i ∂Yₗₘ/∂r̂x
For m = 0:  ∂Yₗₘ/∂r̂y = 0

Reference: Equation (D19) in R. Drautz, Phys. Rev. B 102, 024104 (2020)
"""
function ∂Yₗₘ_∂r̂y(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	validate_lm(l, m)
	validate_uvec(uvec)

	# Early return for m = 0
	m == 0 && return Complex{Float64}(0)

	# The sign of the imaginary unit depends on the sign of m
	im_factor = m < 0 ? -im : im
	return im_factor * ∂Yₗₘ_∂r̂x(l, m, uvec)
end

"""
	∂Yₗₘ_∂r̂z(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Compute the partial derivative of complex spherical harmonic Yₗₘ with respect to r̂z.

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `uvec`: Normalized 3D direction vector [r̂x, r̂y, r̂z]

# Mathematical Details
For m ≥ 0:  ∂Yₗₘ/∂r̂z = (r̂x + ir̂y)ᵐ dP̄ₗₘ/dr̂z
For m < 0:  ∂Yₗₘ/∂r̂z = (-1)ⁿ(r̂x - ir̂y)ⁿ dP̄ₗₘ/dr̂z

where:
- n = |m|
- dP̄ₗₘ/dr̂z is the derivative of the normalized associated Legendre polynomial

# Notes
The angular dependence on φ is carried by the complex exponential term (r̂x ± ir̂y)ᵐ,
while the θ dependence is in the derivative of P̄ₗₘ.

Reference: Equation (D20) in R. Drautz, Phys. Rev. B 102, 024104 (2020)
"""
function ∂Yₗₘ_∂r̂z(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	validate_lm(l, m)
	validate_uvec(uvec)

	# Common calculations
	n = abs(m)
	dplm = dP̄ₗₘ(l, n, uvec[3])

	# Complex coordinate z = x ± iy
	z_xy = if m < 0
		ComplexF64(uvec[1], -uvec[2])  # x - iy
	else
		ComplexF64(uvec[1], uvec[2])   # x + iy
	end

	# Phase factor for negative m
	phase = m < 0 ? (-1)^n : 1

	return phase * z_xy^n * dplm
end

"""
	yₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Compute the directional derivative (uvec ⋅ ∇)Yₗₘ.

# Mathematical Details
yₗₘ = x ∂Yₗₘ/∂x̂ + y ∂Yₗₘ/∂ŷ + z ∂Yₗₘ/∂ẑ

Reference: Equation (D22) in R. Drautz, Phys. Rev. B 102, 024104 (2020)
"""
function yₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex
	validate_lm(l, m)
	validate_uvec(uvec)

	return uvec[1] * ∂Yₗₘ_∂r̂x(l, m, uvec) + uvec[2] * ∂Yₗₘ_∂r̂y(l, m, uvec) +
		   uvec[3] * ∂Yₗₘ_∂r̂z(l, m, uvec)
end

"""
	Zₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Compute the tesseral harmonic Zₗₘ.

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `uvec`: Normalized 3D direction vector [r̂x, r̂y, r̂z]

# Mathematical Details
For m = 0:  Zₗₘ = P̄ₗₘ(r̂z)
For m > 0:  Zₗₘ = (-1)ᵐ√2 P̄ₗₘ(r̂z) ∑ₖ (-1)ᵏ (m,2k) r̂x^(m-2k) r̂y^(2k)
For m < 0:  Zₗₘ = (-1)ⁿ√2 P̄ₗₘ(r̂z) ∑ₖ (-1)ᵏ (n,2k+1) r̂x^(n-2k-1) r̂y^(2k+1)

where n = |m| and (n,k) denotes binomial coefficient.

Reference: Equation (***) in T. Tanaka and Y. Gohda, ***
"""
function Zₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	# Special case for m = 0
	m == 0 && return real(Yₗₘ(l, 0, uvec))

	# Common factors for m ≠ 0
	n = abs(m)
	plm = P̄ₗₘ(l, n, uvec[3])
	phase = (-1)^n
	common_factor = phase * √2 * plm

	if m > 0
		return common_factor * sum(
			(-1)^k * binomial(m, 2k) * uvec[1]^(m - 2k) * uvec[2]^(2k) for
			k in 0:floor(Int, m / 2)
		)
	else
		return common_factor * sum(
			(-1)^k * binomial(n, 2k + 1) *
			uvec[1]^(n - (2k + 1)) * uvec[2]^(2k + 1)
			for k in 0:floor(Int, (n - 1) / 2)
		)
	end
end

"""
	∂Zₗₘ_∂r̂x(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Compute the partial derivative of tesseral harmonic Zₗₘ with respect to r̂x.

# Mathematical Details
For m = 0:  ∂Zₗₘ/∂r̂x = 0
For m > 0:  ∂Zₗₘ/∂r̂x = (-1)ᵐ√2 m P̄ₗₘ ∑ₖ (-1)ᵏ (m-1,2k) r̂x^(m-1-2k) r̂y^(2k)
For m < 0:  ∂Zₗₘ/∂r̂x = (-1)ⁿ√2 n P̄ₗₘ ∑ₖ (-1)ᵏ (n-1,2k+1) r̂x^(n-2-2k) r̂y^(2k+1)

where n = |m| and (n,k) denotes binomial coefficient.
Special case: For m = -1, ∂Zₗₘ/∂r̂x = 0

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `uvec`: Normalized 3D direction vector [r̂x, r̂y, r̂z]

# Returns
- Value of ∂Zₗₘ/∂r̂x

Reference: Equation (***) in T. Tanaka and Y. Gohda, ***
"""
function ∂Zₗₘ_∂r̂x(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	# Early returns for special cases
	m == 0 && return 0.0
	m == -1 && return 0.0

	# Common factors
	m_abs = abs(m)
	plm = P̄ₗₘ(l, m_abs, uvec[3])
	phase = (-1)^m_abs

	if m > 0
		return phase * √2 * m * plm *
			   sum(
				   (-1)^k * binomial(m - 1, 2k) *
				   uvec[1]^(m - 1 - 2k) * uvec[2]^(2k)
				   for k in 0:floor(Int, (m - 1) / 2)
			   )
	else
		return phase * √2 * m_abs * plm *
			   sum(
				   (-1)^k * binomial(m_abs - 1, 2k + 1) *
				   uvec[1]^(m_abs - 2 - 2k) * uvec[2]^(2k + 1)
				   for k in 0:floor(Int, (m_abs - 2) / 2)
			   )
	end
end

"""
	∂Zₗₘ_∂r̂y(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Compute the partial derivative of tesseral harmonic Zₗₘ with respect to r̂y.

# Mathematical Details
For m = 0:  ∂Zₗₘ/∂r̂y = 0
For m > 0:  ∂Zₗₘ/∂r̂y = -∂Zₗₘ/∂r̂x(l, -m)
For m < 0:  ∂Zₗₘ/∂r̂y = +∂Zₗₘ/∂r̂x(l, |m|)

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `uvec`: Normalized 3D direction vector [r̂x, r̂y, r̂z]

# Returns
- Value of ∂Zₗₘ/∂r̂y

# Notes
The y-derivative is related to the x-derivative through sign changes and 
magnetic quantum number inversion, reducing computational complexity.
"""
function ∂Zₗₘ_∂r̂y(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	# Early return for m = 0 case
	m == 0 && return 0.0

	# Use the relationship between x and y derivatives
	if m > 0
		return -∂Zₗₘ_∂r̂x(l, -m, uvec)
	else
		return ∂Zₗₘ_∂r̂x(l, abs(m), uvec)
	end
end

"""
	∂Zₗₘ_∂r̂z(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Compute the partial derivative of tesseral harmonic Zₗₘ with respect to r̂z.

# Mathematical Details
For m = 0:  ∂Zₗₘ/∂r̂z = dP̄ₗₘ/dr̂z
For m > 0:  ∂Zₗₘ/∂r̂z = (-1)ᵐ√2 (dP̄ₗₘ/dr̂z) ∑ₖ (-1)ᵏ (m,2k) r̂x^(m-2k) r̂y^(2k)
For m < 0:  ∂Zₗₘ/∂r̂z = (-1)ⁿ√2 (dP̄ₗₘ/dr̂z) ∑ₖ (-1)ᵏ (n,2k+1) r̂x^(n-2k-1) r̂y^(2k+1)

where n = |m| and (n,k) denotes binomial coefficient.

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `m`: Magnetic quantum number (-l ≤ m ≤ l)
- `uvec`: Normalized 3D direction vector [r̂x, r̂y, r̂z]

# Returns
- Value of ∂Zₗₘ/∂r̂z
"""
function ∂Zₗₘ_∂r̂z(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	# special case for m = 0
	m == 0 && return dP̄ₗₘ(l, m, uvec[3])

	# common factors for m ≠ 0
	n = abs(m)
	dplm = dP̄ₗₘ(l, n, uvec[3])
	phase = (-1)^m
	common_factor = phase * √2 * dplm

	if m > 0
		return common_factor * sum(
			(-1)^k * binomial(m, 2k) * uvec[1]^(m - 2k) * uvec[2]^(2k) for
			k in 0:floor(Int, m / 2)
		)
	else
		return common_factor * sum(
			(-1)^k * binomial(n, 2k + 1) * uvec[1]^(n - (2k + 1)) * uvec[2]^(2k + 1)
			for k in 0:floor(Int, (n - 1) / 2)
		)
	end
end

# zz imply "small z"
function zzₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	return uvec[1] * ∂Zₗₘ_∂r̂x(l, m, uvec) +
		   uvec[2] * ∂Zₗₘ_∂r̂y(l, m, uvec) +
		   uvec[3] * ∂Zₗₘ_∂r̂z(l, m, uvec)
end

function ∂Zₗₘ_∂x(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	return ∂Zₗₘ_∂r̂x(l, m, uvec) - uvec[1] * zzₗₘ(l, m, uvec)
end

function ∂Zₗₘ_∂y(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	return ∂Zₗₘ_∂r̂y(l, m, uvec) - uvec[2] * zzₗₘ(l, m, uvec)
end

function ∂Zₗₘ_∂z(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)

	return ∂Zₗₘ_∂r̂z(l, m, uvec) - uvec[3] * zzₗₘ(l, m, uvec)
end

d_Zlm = [∂Zₗₘ_∂x, ∂Zₗₘ_∂y, ∂Zₗₘ_∂z]

function ∂ᵢZlm(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Vector{Float64}
	validate_lm(l, m)
	validate_uvec(uvec)

	return [∂Zₗₘ_∂x(l, m, uvec), ∂Zₗₘ_∂y(l, m, uvec), ∂Zₗₘ_∂z(l, m, uvec)]
end

"""
	validate_lm(l::Integer, m::Integer)

Validate the quantum numbers l and m for spherical harmonics.

# Arguments
- `l::Integer`: Angular momentum quantum number (l ≥ 0)
- `m::Integer`: Magnetic quantum number (-l ≤ m ≤ l)

# Throws
- `ArgumentError`: If l is negative or if |m| > l

# Examples
```julia
validate_lm(2, 1)    # OK
validate_lm(2, -2)   # OK
validate_lm(-1, 0)   # Throws ArgumentError
validate_lm(1, 2)    # Throws ArgumentError
```
"""
function validate_lm(l::Integer, m::Integer)
	l < 0 && throw(ArgumentError("l must be non-negative (got l = $l)"))
	abs(m) > l && throw(ArgumentError("abs(m) must be ≤ l (got m = $m, l = $l)"))
	return nothing
end

"""
	validate_uvec(uvec::AbstractVector{<:Real}, tol::Real = 1e-8)

Validate that the input vector is a 3D unit vector within specified tolerance.

# Arguments
- `uvec::AbstractVector{<:Real}`: Vector to validate
- `tol::Real`: Tolerance for normalization check (default: 1e-8)

# Throws
- `ArgumentError`: If the vector is not 3D or not normalized
- `DomainError`: If tolerance is negative

# Examples
```julia
validate_uvec([1.0, 0.0, 0.0])     # OK
validate_uvec([1.0, 1.0, 1.0])     # Throws ArgumentError (not normalized)
validate_uvec([1.0, 0.0])          # Throws ArgumentError (not 3D)
```
"""
function validate_uvec(uvec::AbstractVector{<:Real}, tol::Real = 1e-8)
	tol < 0 && throw(DomainError(tol, "tolerance must be non-negative"))
	length(uvec) == 3 ||
		throw(ArgumentError("vector must have exactly 3 elements (got $(length(uvec)))"))

	norm_uvec = norm(uvec)
	isapprox(norm_uvec, 1.0, atol = tol) ||
		throw(ArgumentError("vector must be normalized: norm = $norm_uvec, vec = $uvec"))
	return nothing
end

"""
	validate_r̂z(r̂z::Real)

Validate that the input is a real number in the range [-1, 1].

# Arguments
- `r̂z::Real`: Real number to validate

# Throws
- `ArgumentError`: If r̂z is not in the range [-1, 1]

# Examples
```julia
validate_r̂z(0.0)    # OK
validate_r̂z(1.0)    # OK
validate_r̂z(-1.0)   # OK
validate_r̂z(2.0)    # Throws ArgumentError
validate_r̂z(-2.0)   # Throws ArgumentError
```
"""
function validate_r̂z(r̂z::Real)
	if r̂z < -1.0 || r̂z > 1.0
		throw(ArgumentError("r̂z must be in the range [-1, 1] (got r̂z = $r̂z)"))
	end
	return nothing
end

end
