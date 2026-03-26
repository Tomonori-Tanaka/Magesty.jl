"""
	module MySphericalHarmonics

This module provides functions to compute spherical harmonics ( Y_{l,m} ) and related derivatives, following the formalism described in Drautz (Phys. Rev. B 102, 024104, 2020). It includes:

- Normalized associated Legendre polynomials ( bar{P}_{l,m} ).
- Spherical harmonics ( Y_{l,m} ) and tesseral harmonics ( Z_{l,m} ).
- Partial derivatives with respect to Cartesian coordinates.

# Functions
- `P̄ₗₘ(l, m, r̂z)`: Compute the normalized associated Legendre polynomial.
- `Yₗₘ(l, m, uvec)`: Compute the spherical harmonic (validates inputs).
- `Yₗₘ_unsafe(l, m, uvec)`: Same as `Yₗₘ` without validation (for hot paths).
- `Zₗₘ(l, m, uvec)`: Compute the tesseral harmonic (validates inputs).
- `Zₗₘ_unsafe(l, m, uvec)`: Same as `Zₗₘ` without validation (for hot paths).
- `∂Yₗₘ_∂r̂x`, `∂Yₗₘ_∂r̂y`, `∂Yₗₘ_∂r̂z`, `yₗₘ`, `dP̄ₗₘ`, `∂Zₗₘ_∂r̂x`, `∂Zₗₘ_∂r̂y`, `∂Zₗₘ_∂r̂z`, `zzₗₘ`, `∂Zₗₘ_∂x`, `∂Zₗₘ_∂y`, `∂Zₗₘ_∂z`, `∂ᵢZlm`: validate then compute; each has a `…_unsafe` twin for hot paths.
- `d_Zlm` / `d_Zlm_unsafe`: length-3 lists of the Cartesian `∂Z` callbacks (safe vs unsafe).
"""
module MySphericalHarmonics

using LegendrePolynomials
using LinearAlgebra
using StaticArrays

# abstract type SphericalHarmonicsProduct end
export Zₗₘ, d_Zlm, ∂ᵢZlm
export dP̄ₗₘ_unsafe, Yₗₘ_unsafe, ∂Yₗₘ_∂r̂x_unsafe, ∂Yₗₘ_∂r̂y_unsafe, ∂Yₗₘ_∂r̂z_unsafe, yₗₘ_unsafe
export Zₗₘ_unsafe, ∂Zₗₘ_∂r̂x_unsafe, ∂Zₗₘ_∂r̂y_unsafe, ∂Zₗₘ_∂r̂z_unsafe, zzₗₘ_unsafe
export ∂Zₗₘ_∂x_unsafe, ∂Zₗₘ_∂y_unsafe, ∂Zₗₘ_∂z_unsafe, d_Zlm_unsafe, ∂ᵢZlm_unsafe

# Fast integer parity: (-1)^n without float exponentiation
@inline _parity(n::Integer) = isodd(n) ? -1 : 1

# Compute sqrt((2l+1)/(4π) * (l-m)!/(l+m)!) avoiding large factorial intermediates.
# Uses (l+m)!/(l-m)! = (l-m+1)*(l-m+2)*...*(l+m) to stay in Float64 throughout.
@inline function _plm_norm(l::Int, m::Int)::Float64
    acc = (2 * l + 1) / (4π)
    for i in (l - m + 1):(l + m)
        acc /= i
    end
    return sqrt(acc)
end

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
	am = abs(m)
	return _parity(m) * _plm_norm(l, am) * dnPl(r̂z, l, am)
end

function P̄ₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	return P̄ₗₘ(l, m, uvec[3])
end

"""
	dP̄ₗₘ_unsafe(l::Integer, m::Integer, r̂z::Real) -> Float64

Derivative ``\\mathrm{d}\\bar{P}_{\\ell m}/\\mathrm{d}\\hat{r}_z`` with **no** validation of
`l`, `m`, or `r̂z`. See also [`dP̄ₗₘ`](@ref).
"""
function dP̄ₗₘ_unsafe(l::Integer, m::Integer, r̂z::Real)::Float64
	am = abs(m)
	return _parity(m) * _plm_norm(l, am) * dnPl(r̂z, l, am + 1)
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
- For hot paths with valid `(l, m)` and `r̂z`, use [`dP̄ₗₘ_unsafe`](@ref).

# Examples
```julia
dP̄ₗₘ(2, 1, 0.5)              # Scalar input
dP̄ₗₘ(2, 1, [0.0, 0.0, 1.0])  # Vector input
```
"""
function dP̄ₗₘ(l::Integer, m::Integer, r̂z::Real)::Float64
	validate_lm(l, m)
	validate_r̂z(r̂z)
	return dP̄ₗₘ_unsafe(l, m, r̂z)
end

function dP̄ₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return dP̄ₗₘ_unsafe(l, m, uvec[3])
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


# Notes
- For repeated calls with inputs already checked, use [`Yₗₘ_unsafe`](@ref).

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
	return Yₗₘ_unsafe(l, m, uvec)
end

"""
	Yₗₘ_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex

Same as [`Yₗₘ`](@ref) but **does not** validate `l`, `m`, or `uvec`. Caller must ensure
`-l ≤ m ≤ l`, `l ≥ 0`, and `uvec` is a length-3 unit vector (within numerical tolerance
you care about). Violations may produce wrong results or errors from lower-level code.

Use after a single upfront `validate_lm` / `validate_uvec`, or when inputs come from
invariants elsewhere.
"""
function Yₗₘ_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex
	n = abs(m)
	plm = P̄ₗₘ(l, n, uvec[3])
	if m < 0
		return _parity(n) * ComplexF64(uvec[1], -uvec[2])^n * plm
	else
		return ComplexF64(uvec[1], uvec[2])^n * plm
	end
end

"""
	∂Yₗₘ_∂r̂x_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Same as [`∂Yₗₘ_∂r̂x`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Yₗₘ_∂r̂x_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	m == 0 && return Complex{Float64}(0)
	n = abs(m)
	plm = P̄ₗₘ(l, n, uvec[3])
	z_xy = if m < 0
		ComplexF64(uvec[1], -uvec[2])
	else
		ComplexF64(uvec[1], uvec[2])
	end
	return (m < 0 ? (-1)^n * n : m) * z_xy^(n - 1) * plm
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

# Notes
- For hot paths, use [`∂Yₗₘ_∂r̂x_unsafe`](@ref).
"""
function ∂Yₗₘ_∂r̂x(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Yₗₘ_∂r̂x_unsafe(l, m, uvec)
end

"""
	∂Yₗₘ_∂r̂y_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Same as [`∂Yₗₘ_∂r̂y`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Yₗₘ_∂r̂y_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	m == 0 && return Complex{Float64}(0)
	im_factor = m < 0 ? -im : im
	return im_factor * ∂Yₗₘ_∂r̂x_unsafe(l, m, uvec)
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

# Notes
- For hot paths, use [`∂Yₗₘ_∂r̂y_unsafe`](@ref).
"""
function ∂Yₗₘ_∂r̂y(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Yₗₘ_∂r̂y_unsafe(l, m, uvec)
end

"""
	∂Yₗₘ_∂r̂z_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Same as [`∂Yₗₘ_∂r̂z`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Yₗₘ_∂r̂z_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	n = abs(m)
	dplm = dP̄ₗₘ_unsafe(l, n, uvec[3])
	z_xy = if m < 0
		ComplexF64(uvec[1], -uvec[2])
	else
		ComplexF64(uvec[1], uvec[2])
	end
	phase = m < 0 ? (-1)^n : 1
	return phase * z_xy^n * dplm
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

The angular dependence on φ is carried by the complex exponential term (r̂x ± ir̂y)ᵐ,
while the θ dependence is in the derivative of P̄ₗₘ.

Reference: Equation (D20) in R. Drautz, Phys. Rev. B 102, 024104 (2020)

# Notes
- For hot paths, use [`∂Yₗₘ_∂r̂z_unsafe`](@ref).
"""
function ∂Yₗₘ_∂r̂z(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex{Float64}
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Yₗₘ_∂r̂z_unsafe(l, m, uvec)
end

"""
	yₗₘ_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex

Same as [`yₗₘ`](@ref) without validating `l`, `m`, or `uvec`.
"""
function yₗₘ_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex
	return uvec[1] * ∂Yₗₘ_∂r̂x_unsafe(l, m, uvec) +
		   uvec[2] * ∂Yₗₘ_∂r̂y_unsafe(l, m, uvec) +
		   uvec[3] * ∂Yₗₘ_∂r̂z_unsafe(l, m, uvec)
end

"""
	yₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Complex{Float64}

Compute the directional derivative (uvec ⋅ ∇)Yₗₘ.

# Mathematical Details
yₗₘ = x ∂Yₗₘ/∂x̂ + y ∂Yₗₘ/∂ŷ + z ∂Yₗₘ/∂ẑ

Reference: Equation (D22) in R. Drautz, Phys. Rev. B 102, 024104 (2020)

# Notes
- For hot paths, use [`yₗₘ_unsafe`](@ref).
"""
function yₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Complex
	validate_lm(l, m)
	validate_uvec(uvec)
	return yₗₘ_unsafe(l, m, uvec)
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

# Notes
- For repeated calls with inputs already checked, use [`Zₗₘ_unsafe`](@ref).
"""
function Zₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return Zₗₘ_unsafe(l, m, uvec)
end

"""
	Zₗₘ_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Same as [`Zₗₘ`](@ref) but **does not** validate `l`, `m`, or `uvec`. Caller must ensure
`-l ≤ m ≤ l`, `l ≥ 0`, and `uvec` is a length-3 unit vector. Violations may produce wrong
results or errors from lower-level code.
"""
function Zₗₘ_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	m == 0 && return P̄ₗₘ(l, 0, uvec[3])

	n = abs(m)
	plm = P̄ₗₘ(l, n, uvec[3])
	c = _parity(n) * √2 * plm
	z_pow = ComplexF64(uvec[1], uvec[2])^n

	return m > 0 ? c * real(z_pow) : c * imag(z_pow)
end

"""
	∂Zₗₘ_∂r̂x_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Same as [`∂Zₗₘ_∂r̂x`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Zₗₘ_∂r̂x_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	m == 0 && return 0.0
	n = abs(m)
	plm = P̄ₗₘ(l, n, uvec[3])
	c = _parity(n) * √2 * n * plm
	z_pow_n1 = ComplexF64(uvec[1], uvec[2])^(n - 1)
	return m > 0 ? c * real(z_pow_n1) : c * imag(z_pow_n1)
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

# Notes
- For hot paths, use [`∂Zₗₘ_∂r̂x_unsafe`](@ref).
"""
function ∂Zₗₘ_∂r̂x(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Zₗₘ_∂r̂x_unsafe(l, m, uvec)
end

"""
	∂Zₗₘ_∂r̂y_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Same as [`∂Zₗₘ_∂r̂y`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Zₗₘ_∂r̂y_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	m == 0 && return 0.0
	n = abs(m)
	plm = P̄ₗₘ(l, n, uvec[3])
	c = _parity(n) * √2 * n * plm
	z_pow_n1 = ComplexF64(uvec[1], uvec[2])^(n - 1)
	return m > 0 ? -c * imag(z_pow_n1) : c * real(z_pow_n1)
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
For hot paths, use [`∂Zₗₘ_∂r̂y_unsafe`](@ref).
"""
function ∂Zₗₘ_∂r̂y(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Zₗₘ_∂r̂y_unsafe(l, m, uvec)
end

"""
	∂Zₗₘ_∂r̂z_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Same as [`∂Zₗₘ_∂r̂z`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Zₗₘ_∂r̂z_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	n = abs(m)
	dplm = dP̄ₗₘ_unsafe(l, n, uvec[3])
	m == 0 && return dplm
	c = _parity(n) * √2 * dplm
	z_pow = ComplexF64(uvec[1], uvec[2])^n
	return m > 0 ? c * real(z_pow) : c * imag(z_pow)
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

# Notes
- For hot paths, use [`∂Zₗₘ_∂r̂z_unsafe`](@ref).
"""
function ∂Zₗₘ_∂r̂z(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Zₗₘ_∂r̂z_unsafe(l, m, uvec)
end

"""
	zzₗₘ_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Same as [`zzₗₘ`](@ref) (directional derivative ``\\hat{r}\\cdot\\nabla Z_{\\ell m}`` in the
tangent frame) without validating `l`, `m`, or `uvec`.
"""
function zzₗₘ_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	return uvec[1] * ∂Zₗₘ_∂r̂x_unsafe(l, m, uvec) +
		   uvec[2] * ∂Zₗₘ_∂r̂y_unsafe(l, m, uvec) +
		   uvec[3] * ∂Zₗₘ_∂r̂z_unsafe(l, m, uvec)
end

# zz imply "small z"
function zzₗₘ(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return zzₗₘ_unsafe(l, m, uvec)
end

"""
	∂Zₗₘ_∂x_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Same as [`∂Zₗₘ_∂x`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Zₗₘ_∂x_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	return ∂Zₗₘ_∂r̂x_unsafe(l, m, uvec) - uvec[1] * zzₗₘ_unsafe(l, m, uvec)
end

function ∂Zₗₘ_∂x(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Zₗₘ_∂x_unsafe(l, m, uvec)
end

"""
	∂Zₗₘ_∂y_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Same as [`∂Zₗₘ_∂y`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Zₗₘ_∂y_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	return ∂Zₗₘ_∂r̂y_unsafe(l, m, uvec) - uvec[2] * zzₗₘ_unsafe(l, m, uvec)
end

function ∂Zₗₘ_∂y(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Zₗₘ_∂y_unsafe(l, m, uvec)
end

"""
	∂Zₗₘ_∂z_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Float64

Same as [`∂Zₗₘ_∂z`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂Zₗₘ_∂z_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	return ∂Zₗₘ_∂r̂z_unsafe(l, m, uvec) - uvec[3] * zzₗₘ_unsafe(l, m, uvec)
end

function ∂Zₗₘ_∂z(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Float64
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂Zₗₘ_∂z_unsafe(l, m, uvec)
end

d_Zlm = [∂Zₗₘ_∂x, ∂Zₗₘ_∂y, ∂Zₗₘ_∂z]
d_Zlm_unsafe = [∂Zₗₘ_∂x_unsafe, ∂Zₗₘ_∂y_unsafe, ∂Zₗₘ_∂z_unsafe]

"""
	∂ᵢZlm_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Vector{Float64}

Same as [`∂ᵢZlm`](@ref) without validating `l`, `m`, or `uvec`.
"""
function ∂ᵢZlm_unsafe(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::SVector{3,Float64}
	x, y, z = uvec[1], uvec[2], uvec[3]
	n = abs(m)
	plm  = P̄ₗₘ(l, n, z)
	dplm = dP̄ₗₘ_unsafe(l, n, z)

	if m == 0
		# ∂r̂x = ∂r̂y = 0, ∂r̂z = dplm  →  zzₗₘ = z * dplm
		zz = z * dplm
		return SVector{3,Float64}(-x * zz, -y * zz, dplm - z * zz)
	end

	c = _parity(n) * √2
	z_xy     = ComplexF64(x, y)
	z_pow_n  = z_xy^n
	z_pow_n1 = z_xy^(n - 1)
	rn  = real(z_pow_n);  in_  = imag(z_pow_n)
	rn1 = real(z_pow_n1); in1  = imag(z_pow_n1)

	# ∂r̂ components of Zₗₘ
	dZx, dZy, dZz = if m > 0
		(c * n * plm * rn1, -c * n * plm * in1, c * dplm * rn)
	else
		(c * n * plm * in1,  c * n * plm * rn1, c * dplm * in_)
	end

	# zzₗₘ = r̂ ⋅ ∂r̂Z  (computed once, used three times)
	zz = x * dZx + y * dZy + z * dZz
	return SVector{3,Float64}(dZx - x * zz, dZy - y * zz, dZz - z * zz)
end

"""
	∂ᵢZlm(l::Integer, m::Integer, uvec::AbstractVector{<:Real}) -> Vector{Float64}

Cartesian gradient ``(\\partial_x Z_{\\ell m}, \\partial_y Z_{\\ell m}, \\partial_z Z_{\\ell m})``.

# Notes
- For hot paths, use [`∂ᵢZlm_unsafe`](@ref).
"""
function ∂ᵢZlm(l::Integer, m::Integer, uvec::AbstractVector{<:Real})::Vector{Float64}
	validate_lm(l, m)
	validate_uvec(uvec)
	return ∂ᵢZlm_unsafe(l, m, uvec)
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
