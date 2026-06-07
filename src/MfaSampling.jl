"""
	MfaSampling

Code-agnostic spin-configuration sampler based on the Mean-Field Approximation
(MFA). Given an initial spin matrix (`3 × n_atoms`, columns are atoms), it
draws thermally conditioned configurations whose per-atom directions follow a
von Mises-Fisher (vMF) distribution on the unit sphere and whose magnitudes are
preserved.

The sampler knows nothing about file formats; DFT-code adapters (e.g. VASP
INCAR) read the initial spin matrix, call [`mfa_sweep`](@ref), and write the
resulting configurations back. This keeps the physics reusable across codes.

## Conventions

- The control variable is either `"tau"` (scaled temperature ``τ = T/T_c``) or
  `"m"` (magnetization). Both are mapped onto the same self-consistent
  magnetization ``m`` and concentration ``κ = 3m/τ``.
- MFA self-consistency: ``m = \\coth(3m/τ) - τ/3m``, clamped to
  ``τ ∈ (MIN\\_TEMP, MAX\\_TEMP)``.
- Boundary behavior matches the reference sampler: ``τ < MIN\\_TEMP`` returns
  the input unchanged (fully ordered), ``τ > MAX\\_TEMP`` uses a vanishing
  concentration (near-uniform draw).
"""
module MfaSampling

import Roots
using LinearAlgebra
using Random
using StaticArrays

# Control variables accepted by the sampler.
const VALID_VARIABLES = ("tau", "m")

# Numerical guards on the scaled temperature τ = T/Tc.
const MIN_TEMP = 1.0e-5
const MAX_TEMP = 0.99999
# Concentration used in the near-uniform (τ > MAX_TEMP) limit; mirrors the
# reference sampler so the high-temperature distribution is reproduced exactly.
const KAPPA_UNIFORM = 1.0e-6
# Below this concentration the inverse-CDF draw degenerates, so fall back to a
# fully isotropic direction.
const KAPPA_MIN = 1.0e-12
# Norm below which a spin is treated as zero and left untouched.
const ZERO_NORM_ATOL = 1.0e-10

"""
	thermal_averaged_m(τ::Real)::Float64

Solve the MFA self-consistency equation ``m = \\coth(3m/τ) - τ/3m`` for the
thermally averaged magnetization ``m`` at scaled temperature `τ` (`T/Tc`).

Returns `1.0` for ``τ < MIN\\_TEMP`` (fully ordered) and `0.0` for
``τ > MAX\\_TEMP`` (fully disordered).
"""
function thermal_averaged_m(τ::Real)::Float64
	if τ < MIN_TEMP
		return 1.0
	elseif τ > MAX_TEMP
		return 0.0
	end
	f(m) = m - coth(3m / τ) + τ / 3m
	return Roots.find_zero(f, (MIN_TEMP, MAX_TEMP), Roots.Brent())
end

"""
	tau_from_magnetization(m::Real)::Float64

Invert the MFA self-consistency equation: return the scaled temperature
``τ = T/Tc`` whose thermally averaged magnetization is `m`.

Returns `1.0` for ``m ≤ 0`` (disordered) and `0.0` for ``m ≥ 1`` (ordered).
Magnetizations whose temperature falls outside ``[MIN\\_TEMP, MAX\\_TEMP]`` are
clamped to the ordered (`0.0`) or disordered (`1.0`) limit, which the sampler
handles directly.
"""
function tau_from_magnetization(m::Real)::Float64
	m <= 0.0 && return 1.0
	m >= 1.0 && return 0.0
	g(τ) = m - coth(3m / τ) + τ / 3m
	g_lo = g(MIN_TEMP)
	g_hi = g(MAX_TEMP)
	# Root outside the bracket (very ordered or very disordered m): clamp to the
	# ordered / disordered limit. g is monotonic in τ, so the boundary with the
	# smaller |g| is the one the root lies beyond.
	if g_lo * g_hi > 0
		return abs(g_lo) <= abs(g_hi) ? 0.0 : 1.0
	end
	return Roots.find_zero(g, (MIN_TEMP, MAX_TEMP), Roots.Brent())
end

"""
	sample_vmf_direction(mean_dir::SVector{3, Float64}, κ::Real)::SVector{3, Float64}

Draw a unit vector on ``S^2`` from the von Mises-Fisher distribution with mean
direction `mean_dir` (assumed unit) and concentration `κ ≥ 0`.

Uses the exact closed-form inverse-CDF construction for ``p = 3`` (Ulrich 1984,
Wood 1994): the cosine to `mean_dir` is ``w = 1 + κ^{-1}\\log(u + (1-u)e^{-2κ})``
with ``u ∼ U(0,1)``, combined with a uniform azimuth in the tangent plane. For
``κ < KAPPA\\_MIN`` the draw is isotropic.
"""
function sample_vmf_direction(mean_dir::SVector{3, Float64}, κ::Real)::SVector{3, Float64}
	if κ < KAPPA_MIN
		return _random_unit_vector()
	end
	# Cosine to the mean direction via the inverse CDF. exp(-2κ) underflows
	# cleanly to 0 for large κ, leaving w = 1 + log(u)/κ.
	u = rand()
	w = 1.0 + log(u + (1.0 - u) * exp(-2κ)) / κ
	w = clamp(w, -1.0, 1.0)
	s = sqrt(max(0.0, 1.0 - w * w))
	# Uniform azimuth in the plane orthogonal to mean_dir.
	φ = 2π * rand()
	e1, e2 = _orthonormal_basis(mean_dir)
	tangent = cos(φ) * e1 + sin(φ) * e2
	return w * mean_dir + s * tangent
end

"""
	mfa_sample(spin_matrix::AbstractMatrix{<:Real}, variable::AbstractString,
	           value::Real)::Matrix{Float64}

Draw one sampled configuration from `spin_matrix` (`3 × n_atoms`) at the given
control `value`. `variable` is `"tau"` (scaled temperature) or `"m"`
(magnetization). Per-atom magnitudes are preserved and zero-norm columns are
left untouched. For ``τ < MIN\\_TEMP`` the input is returned unchanged.
"""
function mfa_sample(
	spin_matrix::AbstractMatrix{<:Real},
	variable::AbstractString,
	value::Real,
)::Matrix{Float64}
	τ = _tau_from_value(variable, value)
	# Fully ordered limit: return the input unchanged.
	τ < MIN_TEMP && return Matrix{Float64}(spin_matrix)
	return _sample_vmf_columns(spin_matrix, _kappa_from_tau(τ))
end

"""
	mfa_sweep(spin_matrix::AbstractMatrix{<:Real}; variable, start, stop,
	          num_points, num_samples=1, randomize=false,
	          fixed_indices=Int[], uniform_indices=Int[]) -> Vector{Matrix{Float64}}

Sample configurations from `spin_matrix` (`3 × n_atoms`) over an evenly spaced
sweep of the control variable.

# Keyword arguments
- `variable::AbstractString`: `"tau"` (scaled temperature) or `"m"`
  (magnetization).
- `start::Real`, `stop::Real`, `num_points::Integer`: the sweep values are
  `range(start, stop; length = num_points)`.
- `num_samples::Integer = 1`: configurations drawn per sweep value.
- `randomize::Bool = false`: apply a single random global rotation
  (quantization-axis randomization) to each drawn configuration.
- `fixed_indices::AbstractVector{<:Integer} = Int[]`: 1-based atom indices kept
  at their input directions (rotated by the same global rotation when
  `randomize`), i.e. not sampled.
- `uniform_indices::AbstractVector{<:Integer} = Int[]`: 1-based atom indices
  whose direction is redrawn uniformly on the sphere (the disordered `κ → 0`
  limit) instead of from the vMF distribution, independently of the sweep value.
  The isotropic draw is applied after sampling and overrides it; indices also in
  `fixed_indices` are subsequently reset, so `fixed_indices` takes precedence.

# Returns
- `Vector{Matrix{Float64}}`: sampled configurations in `(value, sample)` order
  (value outer), each `3 × n_atoms` with per-atom magnitudes preserved.
"""
function mfa_sweep(
	spin_matrix::AbstractMatrix{<:Real};
	variable::AbstractString,
	start::Real,
	stop::Real,
	num_points::Integer,
	num_samples::Integer = 1,
	randomize::Bool = false,
	fixed_indices::AbstractVector{<:Integer} = Int[],
	uniform_indices::AbstractVector{<:Integer} = Int[],
)::Vector{Matrix{Float64}}
	_check_variable(variable)
	start <= stop || throw(ArgumentError("start ($start) must be <= stop ($stop)"))
	num_points >= 1 || throw(ArgumentError("num_points must be >= 1; got $num_points"))
	num_samples >= 1 || throw(ArgumentError("num_samples must be >= 1; got $num_samples"))
	size(spin_matrix, 1) == 3 ||
		throw(ArgumentError("spin_matrix must be 3 × n_atoms; got $(size(spin_matrix))"))

	n_atoms = size(spin_matrix, 2)
	_check_indices(fixed_indices, n_atoms, "fixed_indices")
	_check_indices(uniform_indices, n_atoms, "uniform_indices")

	values = range(start, stop; length = num_points)
	configs = Vector{Matrix{Float64}}(undef, num_points * num_samples)
	k = 0
	for value in values
		# τ and κ are constant across the inner sample loop; resolve once.
		τ = _tau_from_value(variable, value)
		ordered = τ < MIN_TEMP
		κ = ordered ? 0.0 : _kappa_from_tau(τ)
		for _ = 1:num_samples
			out = ordered ? Matrix{Float64}(spin_matrix) :
				_sample_vmf_columns(spin_matrix, κ)
			_overwrite_uniform!(out, spin_matrix, uniform_indices)
			# Optional global quantization-axis randomization. Split on
			# `randomize` so the rotation matrix stays concretely typed.
			if randomize
				R = _random_rotation_matrix()
				out = Matrix(R * out)
				for idx in fixed_indices
					out[:, idx] = R * _col_svec(spin_matrix, idx)
				end
			else
				for idx in fixed_indices
					out[:, idx] = _col_svec(spin_matrix, idx)
				end
			end
			k += 1
			configs[k] = out
		end
	end
	return configs
end

"""
	parse_atom_index_spec(spec::AbstractString; max_index::Integer)::Vector{Int}

Parse a 1-based atom-index specification such as `"1-10"`, `"1-10,12,20-22"`,
or `"1 2 5-8"` (commas and/or whitespace separate tokens; `a-b` and `a:b` are
inclusive ranges). Returns sorted, unique indices. An empty spec yields `Int[]`.
Errors on malformed tokens or indices outside `1:max_index`.
"""
function parse_atom_index_spec(spec::AbstractString; max_index::Integer)::Vector{Int}
	s = strip(spec)
	isempty(s) && return Int[]

	tokens = filter(!isempty, split(s, r"[,\s]+"))
	indices = Int[]
	for t in tokens
		if occursin(r"^\d+$", t)
			push!(indices, parse(Int, t))
			continue
		end
		m = match(r"^(\d+)[-:](\d+)$", t)
		if m !== nothing
			a = parse(Int, something(m.captures[1]))
			b = parse(Int, something(m.captures[2]))
			append!(indices, min(a, b):max(a, b))
			continue
		end
		throw(ArgumentError("invalid atom-index token: \"$t\". Example: \"1-10,12,20-22\""))
	end

	unique!(sort!(indices))
	if any(i -> i < 1 || i > max_index, indices)
		throw(ArgumentError("atom indices must be within 1:$max_index; got $(join(indices, ","))"))
	end
	return indices
end

# --- internal helpers ------------------------------------------------------

# Map a control (variable, value) pair to the scaled temperature τ.
function _tau_from_value(variable::AbstractString, value::Real)::Float64
	_check_variable(variable)
	return variable == "tau" ? Float64(value) : tau_from_magnetization(value)
end

# Concentration κ = 3m/τ for MIN_TEMP <= τ; the near-uniform limit above
# MAX_TEMP mirrors the reference sampler.
function _kappa_from_tau(τ::Real)::Float64
	return τ > MAX_TEMP ? KAPPA_UNIFORM : 3 * thermal_averaged_m(τ) / τ
end

# Draw one vMF-sampled configuration at fixed concentration κ. Magnitudes are
# preserved; zero-norm columns are left untouched.
function _sample_vmf_columns(spin_matrix::AbstractMatrix{<:Real}, κ::Real)::Matrix{Float64}
	n_atoms = size(spin_matrix, 2)
	out = zeros(Float64, 3, n_atoms)
	for i = 1:n_atoms
		spin = _col_svec(spin_matrix, i)
		mag = norm(spin)
		isapprox(mag, 0.0; atol = ZERO_NORM_ATOL) && continue
		out[:, i] = sample_vmf_direction(spin / mag, κ) * mag
	end
	return out
end

# Overwrite the given columns with isotropic directions, preserving magnitude.
function _overwrite_uniform!(
	out::AbstractMatrix{Float64},
	spin_matrix::AbstractMatrix{<:Real},
	uniform_indices::AbstractVector{<:Integer},
)
	for idx in uniform_indices
		mag = norm(_col_svec(spin_matrix, idx))
		isapprox(mag, 0.0; atol = ZERO_NORM_ATOL) && continue
		out[:, idx] = _random_unit_vector() * mag
	end
	return nothing
end

# Random rotation taking +z to a uniformly random axis on the sphere.
function _random_rotation_matrix()::SMatrix{3, 3, Float64}
	return _rotation_matrix_from_vectors(SVector(0.0, 0.0, 1.0), _random_unit_vector())
end

# Rotation matrix R with R*v1 = v2 (both normalized) via Rodrigues' formula;
# returns ±I for (anti)parallel inputs.
function _rotation_matrix_from_vectors(
	v1::SVector{3, Float64},
	v2::SVector{3, Float64},
)::SMatrix{3, 3, Float64}
	a = v1 / norm(v1)
	b = v2 / norm(v2)
	axis = cross(a, b)
	sin_θ = norm(axis)
	cos_θ = dot(a, b)
	if isapprox(sin_θ, 0.0; atol = ZERO_NORM_ATOL)
		return cos_θ > 0 ? SMatrix{3, 3, Float64}(I) : SMatrix{3, 3, Float64}(-1.0 * I)
	end
	v = axis / sin_θ
	K = @SMatrix [
		0.0   -v[3]  v[2]
		v[3]   0.0  -v[1]
		-v[2]  v[1]   0.0
	]
	return SMatrix{3, 3, Float64}(I) + sin_θ * K + (1 - cos_θ) * (K * K)
end

# Column i of a 3 × n_atoms matrix as a stack-resident SVector.
@inline function _col_svec(M::AbstractMatrix{<:Real}, i::Integer)::SVector{3, Float64}
	return SVector{3, Float64}(M[1, i], M[2, i], M[3, i])
end

# Uniformly random unit vector on S^2 (Gaussian-normalize; isotropic).
function _random_unit_vector()::SVector{3, Float64}
	v = SVector{3, Float64}(randn(), randn(), randn())
	return v / norm(v)
end

# An orthonormal pair spanning the plane orthogonal to the unit vector `n`.
function _orthonormal_basis(
	n::SVector{3, Float64},
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
	# Pick the global axis least aligned with n to avoid a degenerate cross.
	a = abs(n[1]) <= abs(n[2]) && abs(n[1]) <= abs(n[3]) ? SVector(1.0, 0.0, 0.0) :
		abs(n[2]) <= abs(n[3]) ? SVector(0.0, 1.0, 0.0) : SVector(0.0, 0.0, 1.0)
	e1 = cross(n, a)
	e1 = e1 / norm(e1)
	e2 = cross(n, e1)
	return e1, e2
end

function _check_variable(variable::AbstractString)
	variable in VALID_VARIABLES ||
		throw(ArgumentError("variable must be \"tau\" or \"m\"; got \"$variable\""))
	return nothing
end

function _check_indices(
	indices::AbstractVector{<:Integer},
	max_index::Integer,
	name::AbstractString,
)
	for idx in indices
		(idx < 1 || idx > max_index) &&
			throw(ArgumentError("$name entry $idx is outside 1:$max_index"))
	end
	return nothing
end

end # module MfaSampling
