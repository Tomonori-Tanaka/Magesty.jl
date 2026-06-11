"""
	Fitting.jl

This module contains functions for optimizing the SCE coefficients.
"""
module Fitting

using Base.Threads
using LinearAlgebra
using Printf
using GLMNet
using Statistics
using StaticArrays
using ..TesseralHarmonics
# Hot-path `*_unsafe` variants are not exported by TesseralHarmonics; import the
# ones this module calls explicitly.
using ..TesseralHarmonics: Zₗₘ_unsafe, ∂ᵢZlm_unsafe, Zₗₘ_grad_unsafe
using ..Structures
using ..Symmetries
using ..Clusters
using ..CoupledBases
using ..SALCBases
using ..SpinConfigs
using ..ProgressReporting: with_progress, tick!

export AbstractEstimator, OLS, Ridge, ElasticNet, Lasso, AdaptiveLasso, PrecomputedPilot,
	AdaptiveRidge

"""
	_cluster_scaling(n_sites::Integer) -> Float64

Return the SCE basis normalization factor `(4π)^(n_sites/2)` for an
`n_sites`-body cluster contribution. The factor cancels the
`(4π)^(-n_sites/2)` carried by the per-site tesseral harmonics. Being
dimensionless it does not change the unit of the fitted coefficients
`Jφ` (always the input energy unit, typically eV); it fixes their
normalization so their magnitudes map onto conventional spin-model
parameters.

See the Magesty.jl technical notes for the derivation. This helper is
internal; callers in this module invoke it directly as
`_cluster_scaling(n_C)` to keep the scaling convention in one place.
"""
@inline _cluster_scaling(n_sites::Integer)::Float64 = (4π)^(n_sites / 2)

"""
	_format_bytes(n::Integer) -> String

Format a byte count as a short human-readable string ("12.3 MB", "4.7 GB").
Uses base-1024 units up to TB and one decimal of precision; falls back to
plain bytes below 1 KB. Internal helper for the design-matrix memory log.
"""
function _format_bytes(n::Integer)::String
	# Negative input means upstream integer overflow on the byte-count
	# multiplication; assert rather than silently print a negative size.
	n >= 0 || throw(ArgumentError("_format_bytes expects a non-negative count; got $n"))
	x = Float64(n)
	x < 1024.0 && return @sprintf("%d B", n)
	# Climb units until the value fits in the current unit, or we hit the
	# TB cap (a value that exceeds 1024 TB is rendered in TB without an
	# extra unit).
	x /= 1024.0
	x < 1024.0 && return @sprintf("%.1f KB", x)
	x /= 1024.0
	x < 1024.0 && return @sprintf("%.1f MB", x)
	x /= 1024.0
	x < 1024.0 && return @sprintf("%.1f GB", x)
	x /= 1024.0
	return @sprintf("%.1f TB", x)
end

"""
	_log_design_matrix_memory_estimate(label, rows, num_salcs)

Print the pre-construction memory estimate for a design matrix of size
`rows × num_salcs` plus the `num_salcs × num_salcs` Cholesky Gram matrix.
`label` distinguishes the energy and torque builders in the output.
"""
function _log_design_matrix_memory_estimate(
	label::AbstractString,
	rows::Integer,
	num_salcs::Integer,
)
	bytes      = widen(rows) * num_salcs * sizeof(Float64)
	gram_bytes = widen(num_salcs) * num_salcs * sizeof(Float64)
	println(@sprintf(
		"Memory estimate: %s design matrix %s (%d x %d Float64), Cholesky Gram %s (%d x %d). Total ~%s.",
		label,
		_format_bytes(bytes), rows, num_salcs,
		_format_bytes(gram_bytes), num_salcs, num_salcs,
		_format_bytes(bytes + gram_bytes),
	))
	return nothing
end

"""
	_log_design_matrix_built(label, matrix)

Print the post-construction actual size of a freshly built design matrix
via `Base.summarysize`. Paired with `_log_design_matrix_memory_estimate`.
"""
function _log_design_matrix_built(label::AbstractString, matrix::AbstractMatrix)
	println(@sprintf(
		"%s design matrix built: %s actual.",
		label,
		_format_bytes(Base.summarysize(matrix)),
	))
	return nothing
end

"""
	_normal_equations(X, y) -> (XtX::Matrix{Float64}, Xty::Vector{Float64})

Build the symmetric Gram matrix `X'X` and the right-hand side `X'y` used
by the Cholesky-based normal-equation solvers. Promotes `X` and `y` to
`Float64` without copying when they already are.
"""
@inline function _normal_equations(
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
)
	Xf = X isa Matrix{Float64} ? X : Matrix{Float64}(X)
	yf = y isa Vector{Float64} ? y : Vector{Float64}(y)
	return Xf' * Xf, Xf' * yf
end

"""
	AbstractEstimator

Abstract type for SCE coefficient estimation methods. Concrete subtypes
are passed to `fit(SCEFit, dataset, estimator; torque_weight)` and carry
the estimator's hyperparameters in their fields; estimators without
hyperparameters are zero-field singleton types (`OLS`), and those with
hyperparameters expose them through keyword constructors (`Ridge`).
"""
abstract type AbstractEstimator end

"""
	OLS()

Ordinary least-squares estimator (no regularization). `OLS` is a
zero-field singleton because the OLS solver has no hyperparameter to
tune — there is no `lambda`-like knob, by design. Construct as `OLS()`;
contrast with `Ridge(lambda = ...)`.

The solve goes through Cholesky on the normal equations:
`cholesky(Symmetric(X'X)) \\ (X'y)`. Cholesky is the fastest stable
factorization for symmetric positive-definite systems and keeps memory
to the `num_salcs × num_salcs` Gram matrix.

If the design matrix is rank-deficient or numerically near-collinear,
`OLS` throws `ArgumentError` whose message explains the cause and
recommends `Ridge(lambda = ε)`. An unregularized fit on such data is
physically meaningless (the SCE coefficients are not identifiable), so
the right answer is an explicit error rather than a silent fallback.

# Examples
```julia
est = OLS()
f   = fit(SCEFit, dataset, OLS(); torque_weight = 0.3)
```
"""
struct OLS <: AbstractEstimator end

"""
	Ridge(; lambda::Real = 0.0)

L2-regularized least-squares (ridge) estimator. The penalty applies
uniformly to every SCE coefficient; the bias term `j0` does not need to
be excluded explicitly because it is eliminated analytically before the
solve (see `assemble_weighted_problem` / `extract_j0_jphi`).

The solve is `cholesky(Symmetric(X'X + lambda * I)) \\ (X'y)` — the same
Cholesky-on-normal-equations route as `OLS`, with the `lambda * I` shift
guaranteeing strict positive-definiteness for `lambda > 0` (so Cholesky
cannot fail here). When `lambda ≈ 0`, the call delegates to the OLS
solver and inherits its `PosDefException` → `ArgumentError` behavior.

Unlike `OLS`, `Ridge` carries one hyperparameter, `lambda`, and is
therefore a regular struct rather than a singleton. Constructing
`Ridge(lambda = 0.0)` is exactly equivalent to `OLS` numerically — the
extra type exists so that estimator sweeps can iterate over `OLS()` and
`Ridge(lambda = λ)` without special-casing.

# Fields
- `lambda::Float64`: Regularization strength `λ ∈ [0, ∞)`. `λ = 0`
  reduces to OLS; larger values shrink `jphi` toward zero.

# Examples
```julia
# Default keyword form (lambda = 0.0 -> equivalent to OLS).
est = Ridge()

# Typical regularized fit.
est = Ridge(lambda = 1e-4)
```
"""
struct Ridge <: AbstractEstimator
	lambda::Float64
end

function Ridge(; lambda::Real = 0.0)
	(isfinite(lambda) && lambda >= 0.0) ||
		throw(ArgumentError("Ridge lambda must be finite and non-negative; got $lambda"))
	return Ridge(Float64(lambda))
end

"""
	ElasticNet(; alpha::Real, lambda::Real, standardize::Bool = true)

Elastic-Net estimator backed by GLMNet.jl. Covers Lasso (`alpha = 1`),
GLMNet-style L2 (`alpha = 0`), and honest Elastic Net (mixed norm).
`standardize = true` divides each column of `X` by its empirical
standard deviation before the solve. The per-cluster `(4π)^(N/2)`
basis normalization already puts the columns on an equal footing in
the sphere-averaged (population) sense — that is what makes the fitted
coefficients map onto conventional spin-model parameters — so
standardization corrects only the *residual* per-column scale that
finite, non-uniform sampling leaves behind. L1 and mixed-norm
selection are sensitive to that residual scale (it shifts which
clusters enter the active set), so the default is `true`.

The bias term `j0` is eliminated analytically inside
`assemble_weighted_problem` (the energy block of `X` is mean-centered
before the solve) and re-fit afterward by `extract_j0_jphi` from the
un-scaled energy residual. This estimator therefore calls GLMNet with
`intercept = false`; adding GLMNet's own intercept on top of the
already-centered system would re-introduce a uniform offset across
both energy and torque rows and bias `jphi`.

# Fields
- `alpha::Float64`: Mixing parameter, `0 ≤ alpha ≤ 1`. `alpha = 1` is
  pure Lasso, `alpha = 0` is pure L2 (GLMNet's coordinate-descent
  variant), in between is Elastic Net.
- `lambda::Float64`: Penalty strength, `λ ≥ 0`. `λ = 0` reduces to OLS
  up to GLMNet's coordinate-descent precision.
- `standardize::Bool`: Forwarded to GLMNet, which divides each column
  by its empirical standard deviation. Default `true`: it equalizes
  the residual per-column scale left by finite sampling, to which
  L1 / mixed-norm selection is sensitive. (The `(4π)^(N/2)` basis
  factor already handles the population-level normalization.)

# Examples
```julia
est = ElasticNet(alpha = 0.5, lambda = 1e-3)

# Lasso(λ = ...) is a convenience function returning ElasticNet(alpha = 1, ...).
est = Lasso(lambda = 1e-3)
```
"""
struct ElasticNet <: AbstractEstimator
	alpha::Float64
	lambda::Float64
	standardize::Bool
end

function ElasticNet(; alpha::Real, lambda::Real, standardize::Bool = true)
	0.0 <= alpha <= 1.0 ||
		throw(ArgumentError("ElasticNet alpha must satisfy 0 <= alpha <= 1; got $alpha"))
	lambda >= 0.0 ||
		throw(ArgumentError("ElasticNet lambda must be non-negative; got $lambda"))
	return ElasticNet(Float64(alpha), Float64(lambda), standardize)
end

"""
	Lasso(; lambda::Real, standardize::Bool = true) -> ElasticNet

Convenience function (not a type) returning
`ElasticNet(alpha = 1.0, lambda = lambda, standardize = standardize)`.

There is no separate `Lasso` struct: `Lasso(lambda = ...) isa ElasticNet`
is `true`, and code that needs to detect the Lasso case should test
`e isa ElasticNet && e.alpha == 1.0`. A function (rather than a `const`
alias) prevents `Lasso(alpha = 0.5, ...)` from parsing, which would
contradict the intended meaning "α = 1 only".

# Examples
```julia
est = Lasso(lambda = 1e-3)
est isa ElasticNet  # true
est.alpha           # 1.0
```
"""
Lasso(; lambda::Real, standardize::Bool = true) =
	ElasticNet(alpha = 1.0, lambda = lambda, standardize = standardize)

"""
	AdaptiveLasso(; pilot::AbstractEstimator = OLS(),
	                lambda::Real,
	                gamma::Real = 1.0,
	                epsilon::Real = eps(Float64),
	                standardize::Bool = true)

One-shot Adaptive Lasso (Zou 2006). Runs `pilot` on `(X, y)` to obtain
`beta_pilot`, then solves the weighted-L1 Lasso

    min_b ||y - X * b||^2 / (2 n) + lambda * sum_j w_j * |b_j|

with `w_j = 1 / max(|beta_pilot[j]|, epsilon)^gamma`. Backed by GLMNet.jl
with `intercept = false`; the energy block of `X` is already mean-centered
upstream by `assemble_weighted_problem` and `j0` is recovered downstream by
`extract_j0_jphi`, the same post-processing `OLS`, `Ridge`, and `ElasticNet`
use.

Defaults match the ALAMODE adaptive-LASSO recipe (OLS pilot, `gamma = 1`).
`gamma = 0` reduces to plain Lasso, which the test suite exploits as a
correctness anchor.

Note: GLMNet internally rescales `penalty_factor` so the supplied weights
sum to `nvars`. The user-supplied `lambda` therefore interacts with the
*rescaled* weights, and matched-`lambda` comparison against plain Lasso is
not apples-to-apples once `gamma > 0`. Per-side `lambda` tuning is the
recommended pattern.

# Fields
- `pilot::AbstractEstimator`: First-stage estimator producing `beta_pilot`.
  Default `OLS()` (Zou 2006 verbatim; matches ALAMODE). For SCE designs
  that are rank-deficient or near-collinear (e.g. `num_salcs >=
  num_spinconfigs` at `torque_weight = 0`), use `pilot = Ridge(lambda =
  small)` instead -- the OLS minimum-norm solution populates null-space
  directions with ~1e-10 noise that is not clipped by `eps(Float64)` and
  miscalibrates the adaptive weights.
- `lambda::Float64`: Final L1 penalty strength, `lambda >= 0`.
- `gamma::Float64`: Weight exponent, `gamma >= 0`. Default `1.0`.
  `gamma = 0` reduces to plain Lasso.
- `epsilon::Float64`: Floor on `|beta_pilot[j]|` before reciprocation,
  `epsilon > 0`. Default `eps(Float64)`. Prevents `penalty_factor = Inf`
  when a pilot coefficient is numerically zero (e.g. when `pilot` is
  itself a Lasso).
- `standardize::Bool`: Forwarded to GLMNet. Default `true` for consistency
  with `ElasticNet` / `Lasso`; on SCE designs this is partially redundant
  with the adaptive reweighting but does not hurt.

# Examples
```julia
# Default: OLS pilot, gamma = 1.0 (ALAMODE-style).
est = AdaptiveLasso(lambda = 1e-3)

# Recommended for rank-deficient designs.
est = AdaptiveLasso(pilot = Ridge(lambda = 1e-4), lambda = 1e-3)

# Sanity check: gamma = 0 reduces to plain Lasso.
est = AdaptiveLasso(lambda = 1e-3, gamma = 0.0)

# Reuse a previously fitted SCEFit / SCEModel as the pilot, skipping
# the pilot regression. See `PrecomputedPilot` for the underlying
# adapter and `AdaptiveLasso(::SCEFit; ...)` /
# `AdaptiveLasso(::SCEModel; ...)` for the convenience constructors.
est = AdaptiveLasso(pilot = PrecomputedPilot(coef(prior_fit)), lambda = 1e-3)
```
"""
struct AdaptiveLasso <: AbstractEstimator
	pilot::AbstractEstimator
	lambda::Float64
	gamma::Float64
	epsilon::Float64
	standardize::Bool
end

function AdaptiveLasso(;
	pilot::AbstractEstimator = OLS(),
	lambda::Real,
	gamma::Real = 1.0,
	epsilon::Real = eps(Float64),
	standardize::Bool = true,
)
	lambda >= 0.0 ||
		throw(ArgumentError("AdaptiveLasso lambda must be non-negative; got $lambda"))
	gamma >= 0.0 ||
		throw(ArgumentError("AdaptiveLasso gamma must be non-negative; got $gamma"))
	epsilon > 0.0 ||
		throw(ArgumentError("AdaptiveLasso epsilon must be strictly positive; got $epsilon"))
	return AdaptiveLasso(pilot, Float64(lambda), Float64(gamma), Float64(epsilon), standardize)
end

"""
	PrecomputedPilot(beta::AbstractVector{<:Real})

Estimator adapter that returns a fixed coefficient vector from
`solve_coefficients`, ignoring the supplied `(X, y)` except for a length
check against `size(X, 2)`. Designed as an `AdaptiveLasso.pilot` choice:
lets the adaptive call reuse coefficients from a previous fit instead of
running a fresh pilot regression.

The input vector is copied at construction (enforced by the inner
constructor), so later mutation of the caller's storage does not leak
into `PrecomputedPilot.beta`.

# Fields
- `beta::Vector{Float64}`: Pilot coefficient vector. Named `beta` rather
  than `coef` to avoid visual collision with the `StatsAPI.coef`
  function that Magesty extends.

# Examples
```julia
# Reuse an existing fit's coefficients as the AdaptiveLasso pilot.
est = AdaptiveLasso(pilot = PrecomputedPilot(coef(fit)), lambda = 1e-3)
```
"""
struct PrecomputedPilot <: AbstractEstimator
	beta::Vector{Float64}
	PrecomputedPilot(b::AbstractVector{<:Real}) = new(Vector{Float64}(b))
end

# Compact display: the default struct printer dumps the whole `beta` vector
# (hundreds--thousands of entries), which floods the fit summary. Summarize
# by length instead. `AdaptiveLasso`'s default printer calls `show` on each
# field, so a labeled form keeps the nested pilot readable too.
Base.show(io::IO, p::PrecomputedPilot) =
	print(io, "PrecomputedPilot(", length(p.beta), " coefficients)")

Base.show(io::IO, e::AdaptiveLasso) =
	print(io, "AdaptiveLasso(pilot=", e.pilot, ", lambda=", e.lambda,
		", gamma=", e.gamma, ", epsilon=", e.epsilon,
		", standardize=", e.standardize, ")")

"""
	AdaptiveRidge(; lambda::Real,
	                epsilon::Real = 1e-8,
	                max_iter::Integer = 50,
	                tol::Real = 1e-6)

Iterative Adaptive Ridge estimator (Frommlet & Nuel 2016). Approximates
an L0-penalized fit by repeatedly refitting a per-coefficient weighted
ridge problem

    min_b  ||y - X b||^2 + lambda * sum_j w_j * b_j^2

and updating the weights `w_j = 1 / (b_j^2 + epsilon)` between
iterations. Iteration zero is a plain ridge solve (uniform weights);
each subsequent step rebuilds the weights from the current coefficients.
Large coefficients receive a light penalty and small ones a heavy
penalty, so iterating drives the small coefficients toward zero -- an L0
approximation.

Each weighted ridge subproblem is solved analytically via the closed
form `b = (X'X + lambda * Diagonal(w)) \\ (X'y)`, the same analytic
family as `Ridge`. Unlike `ElasticNet` / `AdaptiveLasso` no GLMNet call
is involved, so there is no `standardize` keyword: the penalty acts
directly on the coefficients, where the per-cluster `(4π)^(N/2)` basis
normalization already places them on the conventional spin-model scale,
so a single `epsilon` is a roughly uniform magnitude floor across
clusters. Because the reweighting approximates L0 (it selects), residual
per-column scale from finite, non-uniform sampling can still influence
which coefficients survive; rescale the design upstream if that matters
under strongly non-uniform sampling.

The iteration stops when the relative infinity-norm change in the
coefficient vector drops below `tol`, or after `max_iter` reweighting
steps, whichever comes first.

# Fields
- `lambda::Float64`: Ridge penalty strength, `lambda >= 0`. `lambda = 0`
  reduces to OLS -- the penalty term vanishes and the iteration is a
  no-op.
- `epsilon::Float64`: Floor added to `b_j^2` before reciprocation,
  `epsilon > 0`. Default `1e-8`. Keeps the weights finite when a
  coefficient is numerically zero, and sets the scale below which a
  coefficient is treated as negligible.
- `max_iter::Int`: Maximum number of reweighting iterations,
  `max_iter >= 1`. Default `50`.
- `tol::Float64`: Convergence threshold on the relative infinity-norm
  coefficient change, `tol > 0`. Default `1e-6`.

# Examples
```julia
# Default iterative Adaptive Ridge.
est = AdaptiveRidge(lambda = 1e-3)

# Tighter convergence, more iterations.
est = AdaptiveRidge(lambda = 1e-3, tol = 1e-8, max_iter = 200)
```
"""
struct AdaptiveRidge <: AbstractEstimator
	lambda::Float64
	epsilon::Float64
	max_iter::Int
	tol::Float64
end

function AdaptiveRidge(;
	lambda::Real,
	epsilon::Real = 1e-8,
	max_iter::Integer = 50,
	tol::Real = 1e-6,
)
	lambda >= 0.0 ||
		throw(ArgumentError("AdaptiveRidge lambda must be non-negative; got $lambda"))
	epsilon > 0.0 ||
		throw(ArgumentError("AdaptiveRidge epsilon must be strictly positive; got $epsilon"))
	max_iter >= 1 ||
		throw(ArgumentError("AdaptiveRidge max_iter must be at least 1; got $max_iter"))
	tol > 0.0 ||
		throw(ArgumentError("AdaptiveRidge tol must be strictly positive; got $tol"))
	return AdaptiveRidge(Float64(lambda), Float64(epsilon), Int(max_iter), Float64(tol))
end

Base.show(io::IO, e::AdaptiveRidge) =
	print(io, "AdaptiveRidge(lambda=", e.lambda, ", epsilon=", e.epsilon,
		", max_iter=", e.max_iter, ", tol=", e.tol, ")")

"""
	SHCache

Per-spinconfig precomputed tesseral spherical-harmonic values (and, for the
torque path, their direction gradients) for every atom in the structure.
Built once per spinconfig and reused across every SALC, coupled basis, and
cluster image during design-matrix construction.

# Layout

Values are addressed by `(lm_idx, atom)` with the standard flat tesseral
ordering

```
lm_idx(l, m) = l^2 + l + m + 1,    0 ≤ l ≤ l_max,   -l ≤ m ≤ l,
```

so the range `[l^2 + 1, (l + 1)^2]` covers all `2l + 1` m-values for a
given `l`. Column-major storage keeps a single atom's values contiguous,
matching the access pattern (inner loop varies `lm_idx`).

# Fields

- `Z::Matrix{Float64}` — shape `((l_max + 1)^2, num_atoms)`. Entry
  `Z[lm_idx(l, m), atom]` holds `Zₗₘ_unsafe(l, m, spin_directions[:, atom])`.
- `∂Z::Matrix{SVector{3, Float64}}` — shape `((l_max + 1)^2, num_atoms)`
  for the torque path, or `0×0` for the energy path. When populated, entry
  `∂Z[lm_idx(l, m), atom]` holds `∂ᵢZlm_unsafe(l, m, spin_directions[:, atom])`.
- `l_max::Int` — the maximum `l` value the cache covers.

Construct with [`build_sh_cache_energy`](@ref) (Z only) or
[`build_sh_cache_torque`](@ref) (Z + ∂Z). Not exported.
"""
struct SHCache
	Z::Matrix{Float64}
	∂Z::Matrix{SVector{3, Float64}}
	l_max::Int
end

# Largest `l` value appearing in any `cbc.ls` across the SALC list. The hot
# path indexes the cache by `(l^2 + l + m + 1, atom)`, so the cache must
# cover `0:l_max`.
function _compute_l_max(
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
)::Int
	l_max = 0
	for key_group in salc_list
		for cbc in key_group
			for l in cbc.ls
				if l > l_max
					l_max = l
				end
			end
		end
	end
	return l_max
end

"""
	build_sh_cache_energy(spin_directions, l_max) -> SHCache

Precompute the tesseral spherical harmonics `Zₗₘ(spin_directions[:, atom])`
for every `(l, m, atom)` with `0 ≤ l ≤ l_max` and `-l ≤ m ≤ l`. The
returned cache has `∂Z` left empty (`0×0`).

Allocates one `legendre_buf` of length `l_max + 1`, the upper bound on
`Zₗₘ_unsafe`'s `(l - |m| + 1)` requirement, and reuses it across every
`(l, m, atom)` call.

# Arguments
- `spin_directions::AbstractMatrix{<:Real}`: `3 × num_atoms` matrix of
  unit-vector spin directions (rows = x, y, z; columns = atoms).
- `l_max::Int`: largest `l` value the cache must cover. Typically the
  result of `_compute_l_max(salc_list)`.

# Returns
- `SHCache`: cache with populated `Z::Matrix{Float64}` of shape
  `((l_max+1)^2, num_atoms)` and empty `∂Z`.

Not exported.
"""
function build_sh_cache_energy(
	spin_directions::AbstractMatrix{<:Real},
	l_max::Int,
)::SHCache
	num_atoms = size(spin_directions, 2)
	nlm = (l_max + 1)^2
	Z = Matrix{Float64}(undef, nlm, num_atoms)
	legendre_buf = Vector{Float64}(undef, l_max + 1)
	@inbounds for atom = 1:num_atoms
		dir = SVector{3, Float64}(
			spin_directions[1, atom],
			spin_directions[2, atom],
			spin_directions[3, atom],
		)
		for l = 0:l_max
			base = l * l
			for m_idx = 1:(2*l+1)
				m = m_idx - l - 1
				Z[base + m_idx, atom] = Zₗₘ_unsafe(l, m, dir, legendre_buf)
			end
		end
	end
	return SHCache(Z, Matrix{SVector{3, Float64}}(undef, 0, 0), l_max)
end

"""
	build_sh_cache_torque(spin_directions, l_max) -> SHCache

Precompute both `Zₗₘ` values and their `∂ᵢZₗₘ` direction gradients for
every `(l, m, atom)` with `0 ≤ l ≤ l_max` and `-l ≤ m ≤ l`. The two
buffered overloads share the same `legendre_buf` (write-only scratch:
the Legendre cache is fully recomputed on every call from `(l, m, z)`,
so the prior contents are not consulted).

# Arguments
- `spin_directions::AbstractMatrix{<:Real}`: `3 × num_atoms` matrix of
  unit-vector spin directions (rows = x, y, z; columns = atoms).
- `l_max::Int`: largest `l` value the cache must cover. Typically the
  result of `_compute_l_max(salc_list)`.

# Returns
- `SHCache`: cache with populated `Z::Matrix{Float64}` and
  `∂Z::Matrix{SVector{3, Float64}}`, both of shape
  `((l_max+1)^2, num_atoms)`.

Not exported.
"""
function build_sh_cache_torque(
	spin_directions::AbstractMatrix{<:Real},
	l_max::Int,
)::SHCache
	num_atoms = size(spin_directions, 2)
	nlm = (l_max + 1)^2
	Z = Matrix{Float64}(undef, nlm, num_atoms)
	∂Z = Matrix{SVector{3, Float64}}(undef, nlm, num_atoms)
	legendre_buf = Vector{Float64}(undef, l_max + 1)
	@inbounds for atom = 1:num_atoms
		dir = SVector{3, Float64}(
			spin_directions[1, atom],
			spin_directions[2, atom],
			spin_directions[3, atom],
		)
		for l = 0:l_max
			base = l * l
			for m_idx = 1:(2*l+1)
				m = m_idx - l - 1
				z_val, grad_val = Zₗₘ_grad_unsafe(l, m, dir, legendre_buf)
				Z[base + m_idx, atom] = z_val
				∂Z[base + m_idx, atom] = grad_val
			end
		end
	end
	return SHCache(Z, ∂Z, l_max)
end

"""
	build_design_matrix_energy(salc_list, spinconfig_list, symmetry;
	                           verbosity = true) -> Matrix{Float64}

Build the energy design matrix `X_E` used for regression. One row per
spin configuration, one column per SALC key group; the bias term `j0` is
not represented as a column — it is recovered analytically after the
solve by `extract_j0_jphi`.

# Arguments
- `salc_list`: Vector of SALC key groups
  (`Vector{Vector{CoupledBases.CoupledBasis_with_coefficient}}`).
- `spinconfig_list`: Vector of spin configurations.
- `symmetry`: Symmetry information (accepted for API stability; orbit
  images are read from `cbc.clusters` and SH values from a per-config
  cache, so this argument is currently unused inside the kernel).

# Keyword arguments
- `verbosity::Bool = true`: When `true`, prints a pre-construction
  memory estimate (design matrix + Cholesky Gram), a thread-count line,
  a progress bar over SALC columns (live bar on TTY, a single
  "Done (X.XX sec)." line on non-TTY), and a post-construction actual
  size via `Base.summarysize`.

# Returns
- `Matrix{Float64}`: Energy design matrix of shape
  `(num_spinconfigs, num_salcs)`.
"""
function build_design_matrix_energy(
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
	spinconfig_list::AbstractVector{SpinConfig},
	symmetry::Symmetry;
	verbosity::Bool = true,
)::Matrix{Float64}
	num_salcs = length(salc_list)  # Number of key groups
	num_spinconfigs = length(spinconfig_list)

	if verbosity
		_log_design_matrix_memory_estimate("energy", num_spinconfigs, num_salcs)
	end

	# Construct design matrix A in Ax = b. One column per SALC; the bias
	# term `j0` is not represented as a column here — it is recovered
	# analytically after the solve by `extract_j0_jphi`.
	design_matrix = zeros(Float64, num_spinconfigs, num_salcs)

	_ = symmetry  # accepted for API stability; orbit images live on cbc.clusters, SH values come from sh_caches

	# Precompute the per-spinconfig SH cache up front, in parallel. The
	# energy kernel below threads over SALCs and reads every spinconfig
	# inside that loop, so building one cache per `j` shared read-only
	# across the SALC threads is strictly cheaper than rebuilding inside
	# the inner loop. Memory is small: `num_spinconfigs * (l_max+1)^2 *
	# num_atoms * 8 B` (FeGe 2×2×2 light: ~1.8 MB at l_max=5).
	l_max = _compute_l_max(salc_list)
	sh_caches = Vector{SHCache}(undef, num_spinconfigs)
	@threads for j = 1:num_spinconfigs
		sh_caches[j] = build_sh_cache_energy(spinconfig_list[j].spin_directions, l_max)
	end

	# Parallel over key-group columns; each thread writes to a disjoint column.
	# Rank-erasing annotations on `key_group` / `cbc` are intentionally absent so
	# that Julia specializes `design_matrix_energy_element` on each element's
	# concrete `CoupledBasis_with_coefficient{R}` type via call-site dispatch.
	if verbosity
		println(@sprintf(
			"Threading: %d SALC columns across %d thread(s).",
			num_salcs, nthreads(),
		))
	end
	with_progress(num_salcs, "Building energy design matrix"; verbosity = verbosity) do prog
		@threads for i = 1:num_salcs
			key_group = salc_list[i]
			n_C = length(key_group[1].atoms)  # Number of sites in the cluster
			scaling_factor = _cluster_scaling(n_C)
			@inbounds for j = 1:num_spinconfigs
				sh_cache = sh_caches[j]
				# Sum contributions from all CoupledBasis_with_coefficient in this key group
				group_value = 0.0
				for cbc in key_group
					group_value += design_matrix_energy_element(cbc, sh_cache)
				end
				design_matrix[j, i] = group_value * scaling_factor
			end
			tick!(prog)
		end
	end

	if verbosity
		_log_design_matrix_built("Energy", design_matrix)
	end

	return design_matrix
end


"""
	design_matrix_energy_element(cbc, sh_cache) -> Float64

Compute one energy-design feature for a given `CoupledBasis_with_coefficient`
and a precomputed [`SHCache`](@ref).

# Description
- Contracts the precomputed folded tensor (Mf axis already collapsed) with
  tesseral spherical-harmonic values read from `sh_cache`, summed over all
  symmetry-equivalent clusters in `cbc.clusters` and weighted by
  `cbc.multiplicity`.
- Equivalent to one column entry (excluding bias) in the energy design
  matrix, before the cluster scaling factor is applied by the caller.

# Arguments
- `cbc::CoupledBases.CoupledBasis_with_coefficient`: the coupled basis with
  its SALC coefficient and pre-enumerated cluster orbit.
- `sh_cache::SHCache`: per-spinconfig SH cache produced by
  [`build_sh_cache_energy`](@ref) (or `build_sh_cache_torque`). Only the
  `Z` field is consulted here; `∂Z` may be empty.

# Returns
- `Float64`: Feature value for the `CoupledBasis_with_coefficient`.
"""
function design_matrix_energy_element(
	cbc::CoupledBases.CoupledBasis_with_coefficient{R, N},
	sh_cache::SHCache,
)::Float64 where {R, N}
	# N = number of sites (also the folded_tensor rank); R = N + 1 is the
	# coeff_tensor rank. Translated cluster images are read from
	# `cbc.clusters` (precomputed by `CoupledBases.enumerate_orbit_clusters`
	# at SALC-build time); tesseral SH values come from the per-spinconfig
	# `sh_cache`.
	result::Float64 = 0.0
	dims_t = ntuple(i -> 2 * cbc.ls[i] + 1, Val(R - 1))
	other_dims_t = ntuple(i -> dims_t[i], Val(R - 2))
	other_site_indices = CartesianIndices(other_dims_t)
	idx_buf = MVector{R - 1, Int}(undef)

	# `l_offsets[site_idx] = ls[site_idx]^2` is the flat-index base for the
	# `Z[lm_idx, atom]` lookup at site `site_idx`. Depends only on `cbc.ls`,
	# so hoist out of the cluster loop.
	l_offsets = ntuple(i -> cbc.ls[i]^2, Val(R - 1))
	base_last = l_offsets[R - 1]
	dim_last = dims_t[R - 1]

	Z = sh_cache.Z

	@inbounds for cluster_atoms in cbc.clusters
		atom_last = cluster_atoms[R - 1]

		# Contract folded_tensor with cached spherical harmonics.
		# folded_tensor[m₁,…,m_N] = Σ_Mf coefficient[Mf] · coeff_tensor[…, Mf]
		# is precomputed at SALC-build time, so the Mf axis is already
		# collapsed; we only iterate the site m-indices here. Reuse the
		# product over the first N-1 sites for each last-site index;
		# `idx_buf::MVector{R-1, Int}` makes the splat indexing
		# `folded_tensor[idx_buf...]` statically resolvable.
		tensor_result = 0.0
		for other_tuple in other_site_indices
			product_other = 1.0
			for site_idx = 1:(R - 2)
				m_idx_other = other_tuple.I[site_idx]
				idx_buf[site_idx] = m_idx_other
				product_other *=
					Z[l_offsets[site_idx] + m_idx_other, cluster_atoms[site_idx]]
			end
			for m_idx_last = 1:dim_last
				idx_buf[R - 1] = m_idx_last
				tensor_result +=
					cbc.folded_tensor[idx_buf...] *
					(product_other * Z[base_last + m_idx_last, atom_last])
			end
		end

		result += tensor_result * cbc.multiplicity
	end

	return result
end

"""
	build_design_matrix_torque(salc_list, spinconfig_list, num_atoms, symmetry) -> Matrix{Float64}

Build the torque design matrix used for regression.

# Description
- For each spin configuration, constructs a block of size (3·num_atoms × num_salcs)
  whose rows are per-atom XYZ components of `cross(spin_dir, ∇ₑu)`.
- Blocks are vertically concatenated across configurations.

# Arguments
- `salc_list`: List of CoupledBasis_with_coefficient
- `spinconfig_list`: Vector of spin configurations
- `num_atoms`: Number of atoms in the structure
- `symmetry`: Symmetry information

# Keyword arguments
- `verbosity::Bool = true`: When `true`, prints a pre-construction
  memory estimate (design matrix + Cholesky Gram), a thread-count line,
  a progress bar over spin configurations (live bar on TTY, a single
  "Done (X.XX sec)." line on non-TTY), and a post-construction actual
  size via `Base.summarysize`.

# Returns
- `Matrix{Float64}`: Torque design matrix
"""
function build_design_matrix_torque(
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
	spinconfig_list::AbstractVector{SpinConfig},
	num_atoms::Integer,
	symmetry::Symmetry;
	verbosity::Bool = true,
)::Matrix{Float64}
	num_salcs = length(salc_list)  # Number of key groups
	num_spinconfigs = length(spinconfig_list)
	scaling_factors = Vector{Float64}(undef, num_salcs)
	for (salc_idx, key_group) in enumerate(salc_list)
		n_C = length(key_group[1].atoms)  # Number of sites in the cluster
		scaling_factors[salc_idx] = _cluster_scaling(n_C)
	end

	_ = symmetry  # accepted for API stability; orbit images live on cbc.clusters, SH and gradients come from sh_cache

	# Threading is over spinconfigs already, so each thread builds and
	# owns one `SHCache` for the spinconfig it is processing and reuses
	# it across every (salc_idx, cbc, cluster) accumulation.
	l_max = _compute_l_max(salc_list)

	# Preallocate the full design matrix and let each thread write into its
	# disjoint row block (sc_idx → rows [block_size*(sc_idx-1)+1 : block_size*sc_idx]).
	# Avoids per-thread block allocation and the final vcat copy.
	block_size = 3 * num_atoms
	if verbosity
		_log_design_matrix_memory_estimate("torque", num_spinconfigs * block_size, num_salcs)
	end
	design_matrix = Matrix{Float64}(undef, num_spinconfigs * block_size, num_salcs)

	# Thread-local accumulation buffers, one per possible threadid, pre-
	# allocated outside the `@threads` loop. Each spinconfig iteration
	# zeros its buffer with `fill!` instead of going back to the
	# allocator, which turns `num_spinconfigs` heap allocations of
	# ~1.9 MB each into `maxthreadid()` upfront allocations and removes
	# the corresponding GC pressure from the hot path. Sized by
	# `Threads.maxthreadid()` because Julia 1.10+ exposes thread pools
	# (default / interactive) whose ids can exceed `Threads.nthreads()`.
	# Layout `(3, num_atoms, num_salcs)` puts the `xyz` axis innermost
	# so SVector loads come off a single cache line.
	grad_bufs = [zeros(Float64, 3, num_atoms, num_salcs) for _ in 1:Threads.maxthreadid()]

	if verbosity
		println(@sprintf(
			"Threading: %d spin configurations across %d thread(s).",
			num_spinconfigs, nthreads(),
		))
	end
	with_progress(num_spinconfigs, "Building torque design matrix"; verbosity = verbosity) do prog
		@threads for sc_idx = 1:num_spinconfigs
			spinconfig = spinconfig_list[sc_idx]
			row_offset = block_size * (sc_idx - 1)
			sh_cache = build_sh_cache_torque(spinconfig.spin_directions, l_max)

			grad_buf = grad_bufs[Threads.threadid()]
			fill!(grad_buf, 0.0)

			# Cluster-major sweep: for each (salc, cbc, cluster), one pass
			# over the folded-tensor multi-index distributes the gradient
			# to all N sites of the cluster simultaneously, so the
			# per-(cbc, cluster) work scales with the orbit size rather
			# than with `num_atoms`.
			@inbounds for (salc_idx, key_group) in enumerate(salc_list)
				for cbc in key_group
					_accumulate_grad_torque_cluster!(grad_buf, cbc, salc_idx, sh_cache)
				end
			end

			# Reduction: for each (salc_idx, iatom) cell, take one
			# `cross(spin[iatom], grad_buf[:, iatom, salc_idx])`, apply
			# the per-cluster scaling, and write three contiguous design-
			# matrix rows. Looping salc_idx outermost matches the
			# column-major design-matrix layout (one column = one salc).
			@inbounds for salc_idx = 1:num_salcs
				scaling = scaling_factors[salc_idx]
				for iatom = 1:num_atoms
					dir_iatom_svec = SVector{3, Float64}(
						spinconfig.spin_directions[1, iatom],
						spinconfig.spin_directions[2, iatom],
						spinconfig.spin_directions[3, iatom],
					)
					grad_vec = SVector{3, Float64}(
						grad_buf[1, iatom, salc_idx],
						grad_buf[2, iatom, salc_idx],
						grad_buf[3, iatom, salc_idx],
					)
					torque_vec = cross(dir_iatom_svec, grad_vec) * scaling
					row_base = row_offset + 3 * (iatom - 1)
					design_matrix[row_base + 1, salc_idx] = torque_vec[1]
					design_matrix[row_base + 2, salc_idx] = torque_vec[2]
					design_matrix[row_base + 3, salc_idx] = torque_vec[3]
				end
			end
			tick!(prog)
		end
	end

	if verbosity
		_log_design_matrix_built("Torque", design_matrix)
	end

	return design_matrix
end


"""
	_accumulate_grad_torque_cluster!(grad_buf, cbc, salc_idx, sh_cache)

Cluster-major reverse-mode gradient accumulator for the torque design
matrix. For one `CoupledBasis_with_coefficient` and one SALC column,
sweep every cluster image in `cbc.clusters` and every multi-index of the
folded tensor, distributing the gradient contribution
``\\partial \\Phi_\\nu / \\partial \\hat{\\boldsymbol{e}}_{a_j}`` to all
``N`` sites ``a_j`` of the cluster in one pass. The per-(cbc, cluster)
work therefore scales with the orbit size, not with `num_atoms`: each
site `j` receives its gradient contribution directly via
`grad_buf[:, a_j, salc_idx]`.

# Arguments
- `grad_buf::Array{Float64, 3}`: Accumulation buffer of shape
  `(3, num_atoms, num_salcs)`. The `xyz` axis is innermost so reading a
  single 3-vector costs one cache line. Caller zeros the buffer at the
  start of each spinconfig.
- `cbc::CoupledBases.CoupledBasis_with_coefficient`: the coupled basis,
  with its folded tensor, SALC coefficient, multiplicity, and
  pre-enumerated cluster orbit.
- `salc_idx::Integer`: Column index into `grad_buf` (= SALC key-group index).
- `sh_cache::SHCache`: per-spinconfig SH cache produced by
  [`build_sh_cache_torque`](@ref). Both `Z` and `∂Z` are read.

# Returns
- `nothing` (mutates `grad_buf`).

Not exported.
"""
@inline function _accumulate_grad_torque_cluster!(
	grad_buf::Array{Float64, 3},
	cbc::CoupledBases.CoupledBasis_with_coefficient{R, N},
	salc_idx::Integer,
	sh_cache::SHCache,
) where {R, N}
	# N = number of sites (folded_tensor rank); R = N + 1 = coeff_tensor rank.
	dims_t = ntuple(i -> 2 * cbc.ls[i] + 1, Val(R - 1))
	site_indices = CartesianIndices(dims_t)
	idx_buf = MVector{R - 1, Int}(undef)

	# `l_offsets[k] = ls[k]^2` is the flat-index base for the
	# `Z[lm_idx, atom]` and `∂Z[lm_idx, atom]` lookups at site `k`. Depends
	# only on `cbc.ls`, so hoist out of both the cluster and multi-index
	# loops.
	l_offsets = ntuple(i -> cbc.ls[i]^2, Val(R - 1))

	Z = sh_cache.Z
	∂Z = sh_cache.∂Z
	mult::Float64 = cbc.multiplicity

	@inbounds for cluster_atoms in cbc.clusters
		# Cluster atoms are pairwise distinct: `Clusters.jl` filters every
		# seed cluster by `length(atom_list) == length(unique(atom_list))`,
		# and pure translations are bijections of the supercell, so the
		# distinctness propagates to every orbit image. The `+=`
		# accumulation here is robust to repeated atoms too (their per-site
		# contributions would sum to the correct gradient), but the
		# upstream filter is the contract we assert. Guarded by
		# `@boundscheck` so that production calls inside an `@inbounds`
		# context (and `@inline`'d into this kernel) elide the scan.
		@boundscheck @assert allunique(cluster_atoms) "cluster atoms must be distinct"

		for multi_idx in site_indices
			# Populate the folded-tensor indexer from the CartesianIndex.
			# `idx_buf::MVector{R-1, Int}` makes `folded_tensor[idx_buf...]`
			# statically resolvable.
			for k = 1:(R - 1)
				idx_buf[k] = multi_idx.I[k]
			end
			folded_x_mult = cbc.folded_tensor[idx_buf...] * mult

			# For each site `j` in the cluster, the leave-one-out product
			# `Π_{k ≠ j} Z[l_k, m_k][a_k]` weighted by the folded tensor
			# is the chain-rule coefficient on `∂Z[l_j, m_j][a_j]`.
			# Recomputing the product per `j` is O(N²) at fixed
			# multi-index but `N ≤ 3` for the current basis, so the
			# constant-factor cost is dominated by the indirect loads
			# rather than the multiplications.
			for j = 1:(R - 1)
				m_j = multi_idx.I[j]
				a_j = cluster_atoms[j]
				product_other = 1.0
				for k = 1:(R - 1)
					if k != j
						product_other *=
							Z[l_offsets[k] + multi_idx.I[k], cluster_atoms[k]]
					end
				end
				grad_jm = ∂Z[l_offsets[j] + m_j, a_j]
				coeff = folded_x_mult * product_other
				grad_buf[1, a_j, salc_idx] += coeff * grad_jm[1]
				grad_buf[2, a_j, salc_idx] += coeff * grad_jm[2]
				grad_buf[3, a_j, salc_idx] += coeff * grad_jm[3]
			end
		end
	end
	return nothing
end


"""
	_predict_energy(j0, jphi, salc_list, symmetry, spin_directions) -> Float64

Internal SCE energy evaluation from raw pieces. The `SCEModel`-typed
`predict_energy` method (in `Magesty`) unpacks a model and delegates
here. Returns `j0 + Σ jphi · ϕ(spin)` in eV.
"""
function _predict_energy(
	j0::Float64,
	jphi::AbstractVector{<:Real},
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
	symmetry::Symmetry,
	spin_directions::AbstractMatrix{<:Real},
)::Float64
	_ = symmetry  # accepted for API parity with `build_design_matrix_energy`
	num_salcs = length(salc_list)
	design_vector = Vector{Float64}(undef, num_salcs)
	# Single-threaded inference path; build one cache and reuse across SALCs.
	l_max = _compute_l_max(salc_list)
	sh_cache = build_sh_cache_energy(spin_directions, l_max)

	for i = 1:num_salcs
		key_group = salc_list[i]
		n_C = length(key_group[1].atoms)
		scaling_factor = _cluster_scaling(n_C)
		group_value = 0.0
		for cbc in key_group
			group_value += design_matrix_energy_element(cbc, sh_cache)
		end
		design_vector[i] = group_value * scaling_factor
	end

	return dot(design_vector, jphi) + j0
end

"""
	_predict_torque(jphi, salc_list, symmetry, spin_directions) -> Matrix{Float64}

Internal SCE per-atom torque evaluation from raw pieces. The
`SCEModel`-typed `predict_torque` method (in `Magesty`) unpacks a model
and delegates here. Returns a `3×N` torque matrix (eV); column `i` is
the torque on atom `i`. The sign convention matches
`build_design_matrix_torque` (`cross(spin, ∇ₑu)`).

Throws `ArgumentError` if `spin_directions` is not a `3×N` matrix.
"""
function _predict_torque(
	jphi::AbstractVector{<:Real},
	salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
	symmetry::Symmetry,
	spin_directions::AbstractMatrix{<:Real},
)::Matrix{Float64}
	_ = symmetry  # accepted for API parity with `build_design_matrix_torque`
	if size(spin_directions, 1) != 3
		throw(ArgumentError(
			"spin_directions must be a 3×N matrix (got $(size(spin_directions, 1)) rows)"
		))
	end
	num_atoms = size(spin_directions, 2)
	# Single-threaded inference path; build one cache and reuse across
	# all (salc, cbc, cluster) accumulations.
	l_max = _compute_l_max(salc_list)
	sh_cache = build_sh_cache_torque(spin_directions, l_max)

	# Accumulate the jphi-weighted angular gradient per atom in a single
	# pass over the cluster orbit, then take the cross product with the
	# atom's spin direction once at the end. Folding `scaling * jphi[salc]`
	# into the per-cbc coefficient avoids carrying a `(3, num_atoms,
	# num_salcs)` buffer just to weight by jphi afterwards.
	grad_per_atom = zeros(Float64, 3, num_atoms)
	@inbounds for (salc_idx, key_group) in enumerate(salc_list)
		n_C = length(key_group[1].atoms)
		coeff::Float64 = _cluster_scaling(n_C) * jphi[salc_idx]
		for cbc in key_group
			_accumulate_grad_torque_scaled!(grad_per_atom, cbc, coeff, sh_cache)
		end
	end

	torque = Matrix{Float64}(undef, 3, num_atoms)
	@inbounds for iatom = 1:num_atoms
		dir_iatom = SVector{3, Float64}(@view spin_directions[:, iatom])
		grad_vec = SVector{3, Float64}(
			grad_per_atom[1, iatom],
			grad_per_atom[2, iatom],
			grad_per_atom[3, iatom],
		)
		t = cross(dir_iatom, grad_vec)
		torque[1, iatom] = t[1]
		torque[2, iatom] = t[2]
		torque[3, iatom] = t[3]
	end

	return torque
end

# Cluster-major gradient accumulator with a scalar coefficient (used by
# inference / `_predict_torque`). Shares the algorithm of
# `_accumulate_grad_torque_cluster!` but writes into a `(3, num_atoms)`
# buffer with `coeff` (typically `scaling * jphi[salc]`) baked into the
# per-(cbc, multi-index) factor, so the prediction path does not need a
# `(3, num_atoms, num_salcs)` buffer.
@inline function _accumulate_grad_torque_scaled!(
	grad_per_atom::Matrix{Float64},
	cbc::CoupledBases.CoupledBasis_with_coefficient{R, N},
	coeff::Float64,
	sh_cache::SHCache,
) where {R, N}
	dims_t = ntuple(i -> 2 * cbc.ls[i] + 1, Val(R - 1))
	site_indices = CartesianIndices(dims_t)
	idx_buf = MVector{R - 1, Int}(undef)
	l_offsets = ntuple(i -> cbc.ls[i]^2, Val(R - 1))
	Z = sh_cache.Z
	∂Z = sh_cache.∂Z
	mult_x_coeff::Float64 = cbc.multiplicity * coeff

	@inbounds for cluster_atoms in cbc.clusters
		# See `_accumulate_grad_torque_cluster!` for the upstream
		# cluster-distinctness contract; same `@boundscheck` discipline.
		@boundscheck @assert allunique(cluster_atoms) "cluster atoms must be distinct"
		for multi_idx in site_indices
			for k = 1:(R - 1)
				idx_buf[k] = multi_idx.I[k]
			end
			folded_x_factor = cbc.folded_tensor[idx_buf...] * mult_x_coeff
			for j = 1:(R - 1)
				m_j = multi_idx.I[j]
				a_j = cluster_atoms[j]
				product_other = 1.0
				for k = 1:(R - 1)
					if k != j
						product_other *=
							Z[l_offsets[k] + multi_idx.I[k], cluster_atoms[k]]
					end
				end
				grad_jm = ∂Z[l_offsets[j] + m_j, a_j]
				c = folded_x_factor * product_other
				grad_per_atom[1, a_j] += c * grad_jm[1]
				grad_per_atom[2, a_j] += c * grad_jm[2]
				grad_per_atom[3, a_j] += c * grad_jm[3]
			end
		end
	end
	return nothing
end


"""
	assemble_weighted_problem(design_matrix_energy, design_matrix_torque,
	                          observed_energy_list, observed_torque, weight)
	    -> (X::Matrix{Float64}, y::Vector{Float64})

Build the augmented `(X, y)` linear system that all estimators solve.

`observed_torque` may be passed either as a vector of `3×n_atoms`
matrices (one per config) or as the already-flattened
`3 * n_atoms * n_config` vector; the matrix form is flattened internally
and both produce identical output.

The bias term `j0` is *not* represented as a column of `X`. Instead, the
energy block is mean-centered (column-wise on `design_matrix_energy` and
on `observed_energy_list`) before row scaling, which is the closed-form
result of eliminating `j0` analytically from the mixed energy / torque
objective via `∂L/∂j0 = 0`. The solver returns only `jphi`; the bias is
recovered afterward by `extract_j0_jphi` from the unscaled, un-centered
inputs as `j0 = mean(y_E - X_E * jphi)`.

Energy rows are scaled by `√((1 - weight) / n_E)` and torque rows
(after flattening) by `√(weight / n_T)`, where `n_E` is the number of
energy samples (configs) and `n_T = 3 * n_atoms * n_E` the number of
torque components. This per-sample normalization makes the augmented
least-squares objective equal to `(1 - weight) * MSE_energy + weight *
MSE_torque`, so `weight` selects a convex combination of the two mean
squared errors independent of how many rows each block contributes
(mirroring the energy/force loss convention used by ML interatomic
potentials).

# Arguments
- `design_matrix_energy`: Energy design matrix (no bias column).
- `design_matrix_torque`: Torque design matrix (no bias column).
- `observed_energy_list`: Observed energies.
- `observed_torque`: Observed torques, either as `3×n_atoms` matrices per
  config or as the flattened `3 * n_atoms * n_config` vector.
- `weight`: Trade-off; energy weight = `1 - weight`, torque weight = `weight`.

# Returns
- `X`: Augmented design matrix (centered, weighted energy rows stacked
  above flattened weighted torque rows).
- `y`: Augmented observation vector (centered, weighted energy entries
  stacked above weighted torque entries).
"""
function assemble_weighted_problem(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_torque::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_torque_list::AbstractVector{<:AbstractMatrix{<:Real}},
	weight::Real,
)
	return assemble_weighted_problem(
		design_matrix_energy,
		design_matrix_torque,
		observed_energy_list,
		convert(Vector{Float64}, vcat(vec.(observed_torque_list)...)),
		weight,
	)
end

function assemble_weighted_problem(
	design_matrix_energy::AbstractMatrix{<:Real},
	design_matrix_torque::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
	observed_torque_flattened::AbstractVector{<:Real},
	weight::Real,
)
	# weight parameters
	w_e = 1 - weight
	w_m = weight

	# Per-sample MSE normalization. The objective we want to minimize is
	#   L = (1 - weight) * MSE_energy + weight * MSE_torque
	#     = (1 - weight) / n_E * Σ ε²  +  weight / n_T * Σ τ²
	# but a linear least-squares solver minimizes only the plain residual
	# sum of squares ‖X β - y‖². Scaling a block's rows by a factor c
	# turns its residual sum of squares into c² * Σ(residual²), so to give
	# a block the weight w / n we scale its rows (and the matching entries
	# of y) by √(w / n). This is the standard weighted-least-squares
	# whitening transform W^(1/2) applied to (X, y); the √ appears because
	# the solver squares the residuals.
	# n_E = number of energy samples (configs),
	# n_T = 3 * n_atoms * n_E = number of torque components.
	n_E = length(observed_energy_list)
	n_T = length(observed_torque_flattened)
	scale_e = sqrt(w_e / n_E)
	scale_m = sqrt(w_m / n_T)

	# Eliminate `j0` analytically before row scaling. The closed form
	# `j0_star(jphi) = mean(y_E - X_E * jphi)` (from `∂L/∂j0 = 0`)
	# substituted back into the objective leaves an OLS-in-`jphi`
	# problem on column-centered inputs `(y_E - mean(y_E),
	# X_E - mean(X_E, dims = 1))`. The torque block is unchanged because
	# `j0` does not enter it.
	mean_X_E = mean(design_matrix_energy, dims = 1)
	mean_y_E = mean(observed_energy_list)
	centered_energy = design_matrix_energy .- mean_X_E
	centered_observed_energy = observed_energy_list .- mean_y_E

	# Row-whitening for the per-sample MSE normalization.
	normalized_design_matrix_energy = centered_energy .* scale_e
	normalized_design_matrix_torque = design_matrix_torque .* scale_m
	normalized_observed_energy_list = centered_observed_energy .* scale_e
	normalized_observed_torque_flattened = observed_torque_flattened .* scale_m

	X = vcat(
		normalized_design_matrix_energy,
		normalized_design_matrix_torque,
	)
	y = vcat(
		normalized_observed_energy_list,
		normalized_observed_torque_flattened,
	)

	return X, y
end

"""
	extract_j0_jphi(j_values, design_matrix_energy, observed_energy_list)
	    -> (j0::Float64, jphi::Vector{Float64})

Recover `(j0, jphi)` from the solver output. `j_values` is the
coefficient vector returned by `solve_coefficients` and has length
`num_salcs` — there is no bias slot to discard, because `j0` was
eliminated analytically inside `assemble_weighted_problem`. The bias is
re-estimated from the unscaled, un-centered energy data via
`mean(y_E - X_E * jphi)`, the closed form of `j0_star(jphi)` from
`∂L/∂j0 = 0`. Re-fitting on the original inputs keeps `j0` in the input
energy unit and independent of `torque_weight`.

# Arguments
- `j_values`: Coefficient vector returned by `solve_coefficients`
  (length = `num_salcs`).
- `design_matrix_energy`: Unscaled, un-centered energy design matrix
  (no bias column; shape `(n_E, num_salcs)`).
- `observed_energy_list`: Unscaled observed energies.

# Returns
- `(j0, jphi)`: Bias and SCE coefficients.
"""
function extract_j0_jphi(
	j_values::AbstractVector{<:Real},
	design_matrix_energy::AbstractMatrix{<:Real},
	observed_energy_list::AbstractVector{<:Real},
)
	jphi = collect(j_values)
	j0 = mean(observed_energy_list .- design_matrix_energy * jphi)
	return j0, jphi
end

"""
	solve_coefficients(estimator::AbstractEstimator, X, y) -> Vector{Float64}

Solve the augmented linear system `(X, y)` for the SCE coefficient
vector `jphi`. The bias term `j0` is *not* an entry of the returned
vector: it has been eliminated analytically before this point (the
energy block of `X` is mean-centered inside `assemble_weighted_problem`),
and is recovered afterward by `extract_j0_jphi`.

Each `AbstractEstimator` subtype defines exactly one method of this
function. Unknown estimator types raise `MethodError`.
"""
function solve_coefficients end

function solve_coefficients(
	::OLS,
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
)::Vector{Float64}
	XtX, Xty = _normal_equations(X, y)
	try
		return cholesky(Symmetric(XtX)) \ Xty
	catch e
		e isa PosDefException || rethrow()
		throw(ArgumentError(
			"OLS solve failed: the normal-equation matrix X'X is not " *
			"positive definite, which means the design matrix has " *
			"linearly dependent or numerically ill-conditioned columns. " *
			"Typical causes: (1) the basis size (number of SALCs) exceeds " *
			"the number of spin configurations when `torque_weight = 0` " *
			"(no torque rows to add rank); (2) the basis contains " *
			"linearly dependent SALCs — this can fire even with " *
			"`torque_weight > 0`; (3) the design is full-rank but " *
			"numerically ill-conditioned (κ(X) > 1e8 ish), which Cholesky " *
			"on X'X surfaces as a PosDef failure because the normal-" *
			"equation route squares the condition number. " *
			"Remedy: use `Ridge(lambda = 1e-6)` (or a similarly small " *
			"value) to regularize, or reduce the basis size."
		))
	end
end

"""
	solve_coefficients(estimator::Ridge, X, y) -> Vector{Float64}

L2-regularized solve `cholesky(Symmetric(X'X + lambda * I)) \\ (X'y)`.
The scalar `lambda` applies uniformly to every column — `X` contains
only SCE coefficient columns; the bias has already been removed by
centering in `assemble_weighted_problem`. For `lambda > 0` the shifted
Gram matrix is strictly positive definite, so Cholesky cannot fail.
When `estimator.lambda ≈ 0`, delegates to the OLS solver to inherit its
`PosDefException` → `ArgumentError` behavior.
"""
function solve_coefficients(
	e::Ridge,
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
)::Vector{Float64}
	if e.lambda ≈ 0.0
		return solve_coefficients(OLS(), X, y)
	end
	XtX, Xty = _normal_equations(X, y)
	# In-place diagonal shift `XtX[j,j] += lambda` avoids allocating a
	# second `p × p` matrix (cf. `XtX + lambda * I`). With `lambda > 0`
	# the shifted matrix is strictly SPD, so Cholesky cannot fail.
	@inbounds for j in axes(XtX, 1)
		XtX[j, j] += e.lambda
	end
	return cholesky(Symmetric(XtX)) \ Xty
end

"""
	solve_coefficients(estimator::ElasticNet, X, y) -> Vector{Float64}

Elastic-Net solve via GLMNet's coordinate descent. The augmented `(X, y)`
arriving here has already been row-scaled and energy-mean-centered inside
`assemble_weighted_problem`, so GLMNet sees the right problem directly with
`intercept = false`; `extract_j0_jphi` then recovers `j0` downstream from the
un-scaled energy data — the same post-processing OLS and Ridge already use.
"""
function solve_coefficients(
	e::ElasticNet,
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
)::Vector{Float64}
	return _glmnet_solve(
		X, y;
		alpha = e.alpha,
		lambda = e.lambda,
		standardize = e.standardize,
	)
end

"""
	solve_coefficients(estimator::AdaptiveLasso, X, y) -> Vector{Float64}

One-shot Adaptive Lasso: fit `e.pilot` on `(X, y)`, build per-column
weights `pf[j] = 1 / max(|beta_pilot[j]|, e.epsilon)^e.gamma`, then call
GLMNet with `alpha = 1.0`, `lambda = e.lambda`, and `penalty_factor = pf`.
The augmented `(X, y)` arriving here has already been row-scaled and
energy-mean-centered inside `assemble_weighted_problem`; `extract_j0_jphi`
recovers `j0` downstream from the un-scaled energy data -- the same
post-processing OLS / Ridge / ElasticNet already use.
"""
function solve_coefficients(
	e::AdaptiveLasso,
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
)::Vector{Float64}
	beta_pilot = solve_coefficients(e.pilot, X, y)
	pf = inv.(max.(abs.(beta_pilot), e.epsilon) .^ e.gamma)
	return _glmnet_solve(
		X, y;
		alpha = 1.0,
		lambda = e.lambda,
		standardize = e.standardize,
		penalty_factor = pf,
	)
end

"""
	solve_coefficients(estimator::PrecomputedPilot, X, y) -> Vector{Float64}

Return the stored coefficient vector unchanged, after a length check
against `size(X, 2)`. The supplied `y` is ignored. Throws
`DimensionMismatch` when the stored vector's length disagrees with the
design-matrix column count, which most commonly indicates that the
pilot fit used a different `SCEBasis`.
"""
function solve_coefficients(
	e::PrecomputedPilot,
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
)::Vector{Float64}
	length(e.beta) == size(X, 2) || throw(DimensionMismatch(
		"PrecomputedPilot coefficient length $(length(e.beta)) does not " *
		"match design-matrix column count $(size(X, 2)); the pilot was " *
		"likely fit on a different SCEBasis."))
	return e.beta
end

"""
	solve_coefficients(estimator::AdaptiveRidge, X, y) -> Vector{Float64}

Iterative Adaptive Ridge (Frommlet & Nuel 2016). A fixed-point iteration
that approximates an L0 penalty by repeatedly refitting a per-coefficient
weighted ridge problem. Iteration zero is a plain ridge solve; each
subsequent step rebuilds the weights `w[j] = 1 / (beta[j]^2 + epsilon)`
from the current coefficients and solves the analytic weighted ridge
`beta = (X'X + lambda * Diagonal(w)) \\ (X'y)`. The Gram matrix `X'X` and
`X'y` are formed once and reused; only the penalty diagonal changes
between iterations.

The augmented `(X, y)` arriving here has already been row-scaled and
energy-mean-centered inside `assemble_weighted_problem`; `extract_j0_jphi`
recovers `j0` downstream from the un-scaled energy data -- the same
post-processing OLS / Ridge / ElasticNet already use.

# Arguments
- `estimator::AdaptiveRidge`: Supplies `lambda`, `epsilon`, `max_iter`,
  `tol`.
- `X::AbstractMatrix{<:Real}`: Design matrix, already centered / scaled.
- `y::AbstractVector{<:Real}`: Observation vector, already centered /
  scaled.

# Returns
- `Vector{Float64}`: The converged coefficient vector, or the last
  iterate when `max_iter` reweighting steps are exhausted first.
"""
function solve_coefficients(
	e::AdaptiveRidge,
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
)::Vector{Float64}
	if e.lambda ≈ 0.0
		# The penalty term vanishes; the reweighting iteration is a no-op.
		return solve_coefficients(OLS(), X, y)
	end
	XtX, Xty = _normal_equations(X, y)
	p = size(XtX, 1)
	# Iteration zero: plain ridge with uniform weights, numerically the
	# same fit as `Ridge(e.lambda)`. The reweighting steps start here.
	# `lambda > 0` makes the Gram matrix strictly positive definite, so
	# Cholesky cannot fail. We later rebuild the diagonal in place via
	# `copyto!(A, XtX)`; the iteration-zero shift below uses a temporary
	# matrix because `XtX` itself must stay unmodified for the loop.
	A0 = copy(XtX)
	@inbounds for j in 1:p
		A0[j, j] += e.lambda
	end
	beta = cholesky(Symmetric(A0)) \ Xty
	w = Vector{Float64}(undef, p)
	# Penalty-augmented Gram matrix. Only its diagonal changes between
	# iterations, so it is rebuilt in place: `XtX` is copied back and the
	# weighted `e.lambda * w` is added to the diagonal each step.
	A = Matrix{Float64}(undef, p, p)
	iters = 0
	rel = Inf
	for outer iters in 1:e.max_iter
		@. w = 1.0 / (beta^2 + e.epsilon)
		copyto!(A, XtX)
		@inbounds for j in 1:p
			A[j, j] += e.lambda * w[j]
		end
		# X'X is positive semidefinite and lambda * diag(w) is strictly
		# positive (epsilon > 0 keeps every weight finite and positive),
		# so A is symmetric positive definite -- Cholesky is safe here.
		beta_new = cholesky(Symmetric(A)) \ Xty
		delta = mapreduce((a, b) -> abs(a - b), max, beta_new, beta)
		# The eps(Float64) floor only guards the 0/0 case where the whole
		# coefficient vector is zero; it is not a convergence scale.
		rel = delta / max(maximum(abs, beta_new), eps(Float64))
		beta = beta_new
		rel < e.tol && break
	end
	@debug "AdaptiveRidge solve" iterations = iters relative_change = rel
	return beta
end

# Single-lambda GLMNet wrapper. Returns the coefficient vector for the
# requested (alpha, lambda) point on the regularization path. `intercept`
# is fixed to `false` here because the energy block of `X` has already
# been centered upstream in `assemble_weighted_problem`.
#
# `penalty_factor` is forwarded to GLMNet's per-column weight when given,
# and omitted otherwise so the ElasticNet / Lasso call site stays
# byte-for-byte identical to the pre-AdaptiveLasso form. GLMNet rescales
# the supplied weights internally so that they sum to `nvars`, so the
# user-supplied `lambda` interacts with the rescaled vector.
function _glmnet_solve(
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real};
	alpha::Real,
	lambda::Real,
	standardize::Bool,
	penalty_factor::Union{Nothing, AbstractVector{<:Real}} = nothing,
)::Vector{Float64}
	Xf = Matrix{Float64}(X)
	yf = Vector{Float64}(y)
	path = if penalty_factor === nothing
		glmnet(
			Xf, yf;
			alpha = Float64(alpha),
			lambda = [Float64(lambda)],
			standardize = standardize,
			intercept = false,
		)
	else
		glmnet(
			Xf, yf;
			alpha = Float64(alpha),
			lambda = [Float64(lambda)],
			standardize = standardize,
			intercept = false,
			penalty_factor = Vector{Float64}(penalty_factor),
		)
	end
	return vec(Matrix(path.betas))
end

# --- GCV diagnostics: core numerics --------------------------------------
#
# Generalized cross-validation on the weighted, energy-mean-centered augmented
# system `(X, y)` that `assemble_weighted_problem` produces (energy rows stacked
# above torque rows). For a linear estimator the fitted values obey `ŷ = H y`
# with hat matrix `H`, and
#
#     GCV = (‖r‖² / N) / (1 − tr(H)/N)²,
#
# where `r = y − ŷ` is the augmented weighted residual, `tr(H)` is the effective
# degrees of freedom, and `N` is the number of *live* rows. A block whose
# whitening scale is zero (energy at `torque_weight == 1`, torque at
# `torque_weight == 0`) contributes only dead all-zero rows, which add nothing to
# `r` or to `tr(H)`; they must not inflate `N`. The eliminated bias `j0` costs
# one degree of freedom only when the energy block is live (`torque_weight < 1`).

"""
	_is_linear_estimator(estimator::AbstractEstimator) -> Bool

Whether `estimator` is a linear smoother, i.e. produces fitted values
`ŷ = H y` with a hat matrix `H` independent of `y` (`OLS`, `Ridge`,
`AdaptiveRidge`). GCV is defined only for these. `ElasticNet` (and its `Lasso`
alias, which constructs an `ElasticNet`) and `AdaptiveLasso` are non-linear and
return `false`. The default method returns `false`, so any future estimator is
treated as non-linear until it opts in explicitly.
"""
_is_linear_estimator(::AbstractEstimator)::Bool = false
_is_linear_estimator(::OLS)::Bool = true
_is_linear_estimator(::Ridge)::Bool = true
_is_linear_estimator(::AdaptiveRidge)::Bool = true

"""
	_require_linear_estimator(estimator::AbstractEstimator) -> Nothing

Throw an `ArgumentError` naming `estimator` unless it is a linear smoother
(see `_is_linear_estimator`). Used as a guard by the GCV entry points.
"""
function _require_linear_estimator(estimator::AbstractEstimator)::Nothing
	_is_linear_estimator(estimator) && return nothing
	throw(ArgumentError(
		"GCV requires a linear estimator (OLS, Ridge, or AdaptiveRidge); got " *
		"$(typeof(estimator)). Non-linear estimators -- ElasticNet (including " *
		"the Lasso alias) and AdaptiveLasso -- have no exact hat matrix, so " *
		"the generalized cross-validation score is undefined for them."))
end

"""
	_gcv_sample_count(n_energy, n_torque, torque_weight)
	    -> (n_eff::Int, intercept_dof::Int)

Live-row count `n_eff` and intercept degrees of freedom `intercept_dof` for the
weighted GCV system. A block whose whitening scale is zero contributes only dead
all-zero rows (`scale_e = 0` at `torque_weight == 1`; `scale_m = 0` at
`torque_weight == 0`), so it is excluded from `n_eff`. The eliminated bias `j0`
is identifiable, and costs one degree of freedom, only when the energy block is
live (`torque_weight < 1`).
"""
function _gcv_sample_count(
	n_energy::Integer,
	n_torque::Integer,
	torque_weight::Real,
)::Tuple{Int, Int}
	energy_live = torque_weight < 1.0
	torque_live = torque_weight > 0.0
	n_eff = (energy_live ? Int(n_energy) : 0) + (torque_live ? Int(n_torque) : 0)
	intercept_dof = energy_live ? 1 : 0
	return n_eff, intercept_dof
end

"""
	_gcv_value(rss::Real, dof::Real, n_eff::Integer) -> Float64

Assemble the GCV score `(rss / n_eff) / (1 - dof/n_eff)^2`. Returns `NaN` when
the denominator factor `1 - dof/n_eff` is non-positive — i.e. the effective
degrees of freedom meet or exceed the live-row count, so the model is
(numerically) saturated and GCV is undefined. Callers surface the `NaN` rather
than emit a spurious value.
"""
@inline function _gcv_value(rss::Real, dof::Real, n_eff::Integer)::Float64
	denom = 1.0 - dof / n_eff
	denom > 0.0 || return NaN
	return (rss / n_eff) / denom^2
end

"""
	_gcv_msy(y::AbstractVector{<:Real}, n_eff::Integer) -> Float64

Null-model mean square `‖y‖² / n_eff` of the weighted, energy-mean-centered
augmented target `y` — the prediction error of the `β = 0` baseline (energy
predicted at its mean, torque predicted as zero). Dead all-zero rows carry no
contribution to `‖y‖²`, so `n_eff` (live rows) is the correct divisor. It is the
denominator of the GCV-based predictive R² (see `_gcv_r2`).
"""
@inline function _gcv_msy(y::AbstractVector{<:Real}, n_eff::Integer)::Float64
	return sum(abs2, y) / n_eff
end

"""
	_gcv_r2(gcv::Real, msy::Real) -> Float64

GCV-based predictive R², `1 - gcv / msy`, normalizing the GCV score against the
null-model mean square `msy` (see `_gcv_msy`). Both quantities are in the same
weighted-objective unit, so the ratio is dimensionless: `1` is a perfect fit,
`0` matches the null (mean-energy / zero-torque) model, and negative values mean
the fit predicts worse than the null. Returns `NaN` when `gcv` is `NaN` (a
saturated model) or `msy` is non-positive (a degenerate all-zero target).
"""
@inline function _gcv_r2(gcv::Real, msy::Real)::Float64
	(isnan(gcv) || msy <= 0.0) && return NaN
	return 1.0 - gcv / msy
end

"""
	_effective_dof(X, estimator, beta, intercept_dof) -> Float64

Effective degrees of freedom `tr(H) = intercept_dof + tr(H_β)` of the linear
smoother on the weighted, energy-mean-centered design `X`. `intercept_dof` is
`1` when the eliminated bias `j0` is live (energy block present, i.e.
`torque_weight < 1`) and `0` otherwise. `beta` is the fitted coefficient vector
(`coef`), used only by `AdaptiveRidge`. Dead all-zero rows of `X` do not affect
`rank(X)` or its nonzero singular values, so they leave `tr(H_β)` unchanged.

- `OLS`: `intercept_dof + rank(X)`.
- `Ridge`: `intercept_dof + Σ_k σ_k² / (σ_k² + λ)` over the singular values
  `σ_k` of `X` (reduces to `OLS` when `λ ≈ 0`).
- `AdaptiveRidge`: conditional dof `intercept_dof + tr((X'X + λ D)^{-1} X'X)`
  with the converged diagonal reweighting
  `D = Diagonal(1 ./ (beta.^2 .+ epsilon))`. Recovering `D` from `beta` is exact
  at convergence (the iteration's fixed point sets the weights from the same
  coefficients); only treating `D` as fixed in `H` is the conditional
  approximation.
"""
function _effective_dof(
	X::AbstractMatrix{<:Real},
	::OLS,
	::AbstractVector{<:Real},
	intercept_dof::Integer,
)::Float64
	return intercept_dof + rank(X)
end

function _effective_dof(
	X::AbstractMatrix{<:Real},
	e::Ridge,
	::AbstractVector{<:Real},
	intercept_dof::Integer,
)::Float64
	e.lambda ≈ 0.0 && return intercept_dof + rank(X)
	s2 = abs2.(svdvals(X))
	return intercept_dof + sum(σ2 / (σ2 + e.lambda) for σ2 in s2)
end

function _effective_dof(
	X::AbstractMatrix{<:Real},
	e::AdaptiveRidge,
	beta::AbstractVector{<:Real},
	intercept_dof::Integer,
)::Float64
	e.lambda ≈ 0.0 && return intercept_dof + rank(X)
	length(beta) == size(X, 2) || throw(DimensionMismatch(
		"AdaptiveRidge dof: coefficient length $(length(beta)) does not match " *
		"design-matrix column count $(size(X, 2))."))
	w = inv.(abs2.(beta) .+ e.epsilon)
	Xf = X isa Matrix{Float64} ? X : Matrix{Float64}(X)
	XtX = Xf' * Xf
	A = XtX + e.lambda * Diagonal(w)
	# tr((X'X + λD)^{-1} X'X) via a single SPD solve; cyclic-trace identity
	# tr(A^{-1} XtX) avoids forming the hat matrix.
	return intercept_dof + tr(cholesky(Symmetric(A)) \ XtX)
end

"""
	_gcv_single(X, y, residuals, estimator, beta, n_eff, intercept_dof)
	    -> (gcv::Float64, dof::Float64)

Combined energy+torque GCV score and effective dof for one fit on the
weighted, energy-mean-centered augmented system `(X, y)`. `residuals` is the
augmented weighted residual `y - X*beta` (stored on the fit); dead rows
contribute zero to it. `n_eff` / `intercept_dof` come from `_gcv_sample_count`.
Returns `NaN` for `gcv` when the model is numerically saturated (see
`_gcv_value`).
"""
function _gcv_single(
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
	residuals::AbstractVector{<:Real},
	estimator::AbstractEstimator,
	beta::AbstractVector{<:Real},
	n_eff::Integer,
	intercept_dof::Integer,
)::Tuple{Float64, Float64}
	_require_linear_estimator(estimator)
	rss = sum(abs2, residuals)
	dof = _effective_dof(X, estimator, beta, intercept_dof)
	return _gcv_value(rss, dof, n_eff), dof
end

"""
	_gcv_lambda_path(X, y, lambdas, n_eff, intercept_dof)
	    -> (gcv::Vector{Float64}, dof::Vector{Float64})

Ridge GCV for every `lambda` from a single economy SVD of `X`. With
`X = U diag(σ) Vᵀ` and `a = Uᵀy`:

- `dof(λ) = intercept_dof + Σ_k σ_k² / (σ_k² + λ)`
- `RSS(λ) = ‖y_⊥‖² + Σ_k (λ / (σ_k² + λ))² a_k²`, where
  `‖y_⊥‖² = ‖y‖² − Σ_k a_k²` is the part of `y` outside `col(X)`.

`X` / `y` are the weighted, energy-mean-centered augmented system; dead rows
(zero whitening scale) carry zeros in `y` and `U`, so they drop out of `‖y‖²`
and `a`. `n_eff` / `intercept_dof` come from `_gcv_sample_count`. The caller is
responsible for validating that `lambdas` are non-negative.
"""
function _gcv_lambda_path(
	X::AbstractMatrix{<:Real},
	y::AbstractVector{<:Real},
	lambdas::AbstractVector{<:Real},
	n_eff::Integer,
	intercept_dof::Integer,
)::Tuple{Vector{Float64}, Vector{Float64}}
	svd_fac = svd(X)
	s2 = svd_fac.S .^ 2
	a = svd_fac.U' * y
	y_perp_sq = max(sum(abs2, y) - sum(abs2, a), 0.0)
	gcvs = Vector{Float64}(undef, length(lambdas))
	dofs = Vector{Float64}(undef, length(lambdas))
	@inbounds for (k, λ) in enumerate(lambdas)
		rss = y_perp_sq
		dof = Float64(intercept_dof)
		for j = 1:length(s2)
			denom = s2[j] + λ
			if denom == 0.0
				# σ = 0 and λ = 0: this direction is unfittable. It is fully
				# residual (hat value 0) and adds nothing to the dof.
				rss += a[j]^2
			else
				hat_k = s2[j] / denom            # σ²/(σ²+λ): hat value
				resid_factor_k = λ / denom       # λ/(σ²+λ): 1 − hat value
				rss += resid_factor_k^2 * a[j]^2
				dof += hat_k
			end
		end
		gcvs[k] = _gcv_value(rss, dof, n_eff)
		dofs[k] = dof
	end
	return gcvs, dofs
end

"""
	_calc_rmse(observed_list, predicted_list) -> Float64

Root-mean-square error between two equal-length vectors,
`sqrt(mean((observed_list .- predicted_list) .^ 2))`. Internal helper behind the
public `rmse_*` evaluation verbs.

Throws `ArgumentError` if the lengths differ.
"""
function _calc_rmse(
	observed_list::AbstractVector{<:Real},
	predicted_list::AbstractVector{<:Real},
)::Float64
	if length(observed_list) != length(predicted_list)
		throw(ArgumentError("The lengths of the two lists must be equal."))
	end
	return sqrt(mean((observed_list .- predicted_list) .^ 2))
end

"""
	_calc_r2score(observed_list, predicted_list) -> Float64

Coefficient of determination ``R^2 = 1 - SS_{res}/SS_{tot}`` between observed
and predicted values, with `SS_res = Σ(observed - predicted)²` and
`SS_tot = Σ(observed - mean(observed))²`. Internal helper behind the public
`r2_*` evaluation verbs.
"""
function _calc_r2score(
	observed_list::AbstractVector{<:Real},
	predicted_list::AbstractVector{<:Real},
)::Float64
	ss_res = sum((observed_list .- predicted_list) .^ 2)
	ss_tot = sum((observed_list .- mean(observed_list)) .^ 2)
	return 1 - ss_res / ss_tot
end

end

