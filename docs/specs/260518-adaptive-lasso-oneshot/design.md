# Design: Adaptive Lasso (oneshot) estimator

Status: draft (2026-05-18)
Spec: [`requirements.md`](requirements.md)

## Summary

Add `AdaptiveLasso <: AbstractEstimator` to `src/Fitting.jl`. The
estimator runs the user-supplied `pilot` to get `beta_pilot`, then
calls GLMNet with `alpha = 1.0`, the user-supplied `lambda`, and a
per-column `penalty_factor[j] = 1 / max(|beta_pilot[j]|, epsilon)^gamma`.
The Lasso / ElasticNet path's private helper `_glmnet_solve` gains an
optional `penalty_factor` keyword; the existing `solve_coefficients(::ElasticNet, ...)`
method passes `nothing`, so its behaviour is byte-for-byte unchanged.

`pilot::AbstractEstimator` (concrete instance, not a Symbol). This
keeps the AdaptiveLasso struct symmetric with the rest of the
dispatch hierarchy: any `AbstractEstimator` already in scope can be
plugged in. The default is `OLS()` (the Zou 2006 / ALAMODE choice).
The docstring recommends `pilot = Ridge(lambda = small)` for
rank-deficient SCE designs, since `OLS()` returns Julia's minimum-norm
solution there and the null-space noise miscalibrates the adaptive
weights.

`gamma` is kept exposed (default `1.0`, ALAMODE's fixed value). The
`gamma = 0` corner gives a uniform `penalty_factor = 1` and therefore
reduces to plain Lasso, which the test plan exploits to write a
direct agreement test against `Lasso(lambda = e.lambda)`.

## Module layout

| Target | Change |
|---|---|
| `Project.toml` | No change (`GLMNet` already on `[deps]` from spec 260517). |
| `src/Fitting.jl` | Add `AdaptiveLasso` struct + ctor (`0.0 <= gamma`, `0.0 <= lambda`, `0.0 < epsilon`). Add `solve_coefficients(::AdaptiveLasso, X, y)`. Extend `_glmnet_solve` with `penalty_factor::Union{Nothing, AbstractVector{<:Real}} = nothing`. Update `export`. |
| `src/Magesty.jl` | Re-export `AdaptiveLasso`. |
| `test/component/test_fitting_estimators.jl` | Append a new `@testset "AdaptiveLasso estimator" begin ... end` block with the five property tests. |
| `docs/src/api.md` | List `AdaptiveLasso` in the `## Estimators` `@docs` block. |
| `SPEC.md` | Update the `Fitting` module description, the `# Estimators` API line, and (no) the Primary external libraries table (no new library). |
| `CHANGELOG.md` | `[Unreleased]` `### Added` entry covering `AdaptiveLasso` and the `_glmnet_solve` `penalty_factor` extension. |

`assemble_weighted_problem` and `extract_j0_jphi` are **unchanged**.

## API

```julia
"""
    AdaptiveLasso(; pilot::AbstractEstimator = OLS(),
                    lambda::Real,
                    gamma::Real = 1.0,
                    epsilon::Real = eps(Float64),
                    standardize::Bool = true)

One-shot Adaptive Lasso (Zou 2006). Runs `pilot` on `(X, y)` to obtain
`beta_pilot`, then solves the weighted-L1 Lasso

    min_b ||y - X * b||^2 / (2 n) + lambda * sum_j w_j * |b_j|

with `w_j = 1 / max(|beta_pilot[j]|, epsilon)^gamma`. Backed by
GLMNet.jl with `intercept = false`; the energy block of `X` is already
mean-centered upstream by `assemble_weighted_problem` and `j0` is
recovered downstream by `extract_j0_jphi`, the same post-processing
`OLS`, `Ridge`, and `ElasticNet` use.

# Fields
- `pilot::AbstractEstimator`: First-stage estimator producing
  `beta_pilot`. Default `OLS()` (Zou 2006 verbatim; matches ALAMODE).
  For SCE designs that are rank-deficient or near-collinear
  (e.g. `num_salcs >= num_spinconfigs` at `torque_weight = 0`), use
  `pilot = Ridge(lambda = small)` instead -- the OLS minimum-norm
  solution populates null-space directions with noise that
  miscalibrates the adaptive weights.
- `lambda::Float64`: Final L1 penalty strength, `lambda >= 0`.
- `gamma::Float64`: Weight exponent, `gamma >= 0`. Default `1.0`.
  `gamma = 0` reduces to plain Lasso.
- `epsilon::Float64`: Floor on `|beta_pilot[j]|` before reciprocation,
  `epsilon > 0`. Default `eps(Float64)`. Prevents
  `penalty_factor = Inf` when a pilot coefficient is numerically
  zero (e.g. when `pilot` is itself a Lasso).
- `standardize::Bool`: Forwarded to GLMNet. Default `true` for
  consistency with `Lasso` / `ElasticNet`; on SCE designs this is
  partly redundant with the adaptive reweighting but does not hurt.

# Examples
```julia
# Default: OLS pilot, gamma = 1.0 (ALAMODE-style).
est = AdaptiveLasso(lambda = 1e-3)

# Recommended for rank-deficient designs.
est = AdaptiveLasso(pilot = Ridge(lambda = 1e-4), lambda = 1e-3)

# Sanity check: gamma = 0 reduces to plain Lasso.
est = AdaptiveLasso(lambda = 1e-3, gamma = 0.0)
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
```

### `solve_coefficients` method

```julia
function solve_coefficients(
    e::AdaptiveLasso,
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
)
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
```

### Extended `_glmnet_solve`

```julia
function _glmnet_solve(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real};
    alpha::Real,
    lambda::Real,
    standardize::Bool,
    penalty_factor::Union{Nothing, AbstractVector{<:Real}} = nothing,
)::Vector{Float64}
    kwargs = (
        alpha = Float64(alpha),
        lambda = [Float64(lambda)],
        standardize = standardize,
        intercept = false,
    )
    path = isnothing(penalty_factor) ?
        glmnet(Matrix{Float64}(X), Vector{Float64}(y); kwargs...) :
        glmnet(Matrix{Float64}(X), Vector{Float64}(y);
               kwargs..., penalty_factor = Vector{Float64}(penalty_factor))
    return vec(Matrix(path.betas))
end
```

The branch keeps the `ElasticNet` / `Lasso` path call-site
byte-for-byte identical (it doesn't pass `penalty_factor` and the
helper omits the argument from the GLMNet call).

## Types and conventions

- `AdaptiveLasso` is a concrete subtype of `AbstractEstimator`. No
  new abstract type. No new struct for the Lasso pilot case
  (`Lasso(lambda = ...)` already returns an `ElasticNet` instance, so
  `AdaptiveLasso(pilot = Lasso(lambda = ...), ...)` parses; it is a
  documented but non-default usage).
- `pilot` is held as `::AbstractEstimator` (abstract type). This
  costs a small amount of type-stability inside `solve_coefficients`,
  but `solve_coefficients` runs once per fit (cold path), so the
  cost is negligible. The trade-off is API symmetry: any
  `AbstractEstimator` already in scope can be plugged in without
  going through a type-parameter dance.
- The closed-form energy / torque convention is unchanged. Torque
  rows still do not see `j0`; the energy block arriving at
  `solve_coefficients` is already mean-centered upstream.

## Impact on linked sites

Sites named in CLAUDE.md "Linked sites":

- [ ] Spherical-harmonics convention (`TesseralHarmonics`):
      unaffected.
- [ ] SCE coefficient XML (`save` / `load`): unaffected -- only the
      fit step changes.
- [ ] `Fitting` <-> `SALCBasis`: unaffected -- column order,
      key-group order, and SALC indexing are unchanged.
      `AdaptiveLasso` consumes the same design matrix the rest of
      the pipeline produces.
- [ ] `.claude/agents/` references: the `code-reviewer` and
      `test-runner` agent prompts currently mention `OLS` / `Ridge` /
      `ElasticNet` / `Lasso`. Sweep to add `AdaptiveLasso` if
      any list is exhaustive.
- [ ] `SPEC.md` / `docs/src/api.md`: add `AdaptiveLasso` to the
      public estimator list.

## Test strategy

All tests are appended to `test/component/test_fitting_estimators.jl`
inside a new `@testset "AdaptiveLasso estimator" begin ... end`
block. Synthetic Gaussian fixtures with fixed RNG seeds, as in the
existing ElasticNet tests.

1. **`gamma = 0` reduces to plain Lasso** (strict).
   `AdaptiveLasso(pilot = OLS(), lambda = lambda, gamma = 0.0,
   standardize = false)` matches `Lasso(lambda = lambda, standardize
   = false)` on a centered `(X, y)` pair to within an `atol` chosen
   for GLMNet's coordinate-descent precision. Initial target
   `1e-7`; record the measured worst case in the test docstring and
   `CHANGELOG.md`.
2. **Weight construction is exact** (strict, `==` comparison).
   Manually compute `pf = inv.(max.(abs.(beta_pilot), eps) .^
   gamma)` from a known `beta_pilot` (obtained by calling
   `solve_coefficients(OLS(), X, y)` directly), then call
   `_glmnet_solve(X, y; alpha = 1, lambda, standardize = false,
   penalty_factor = pf)` and compare against
   `solve_coefficients(AdaptiveLasso(pilot = OLS(), lambda, gamma,
   epsilon, standardize = false), X, y)`. Both call sites go through
   the same `_glmnet_solve(... ; penalty_factor)` with identical
   inputs, so the two coefficient vectors must be element-wise
   `==` (or `isequal`), not merely `isapprox`. This test pins the
   exact construction of `pf` from `beta_pilot` -- any future edit
   to the weight formula (or to GLMNet's call shape) trips it.
3. **`epsilon` clip is active for a Lasso pilot** (strict).
   Build a synthetic `(X, y)` where `Lasso(lambda = large)` has
   several exactly-zero entries. Pass `pilot = Lasso(lambda = large)`
   into `AdaptiveLasso(pilot = ..., lambda = small, epsilon = 1e-12)`
   and confirm: (a) the call returns without erroring, (b) the
   resulting fit has all entries finite (`all(isfinite, b_fit)`),
   (c) the columns at which `beta_pilot == 0` are exactly zero in
   the adaptive fit (the very high `1/eps^gamma` penalty drove them
   to zero).
4. **`pilot = Ridge(...)` path works** (qualitative).
   On the synthetic problem from test 5 below, `AdaptiveLasso(pilot
   = Ridge(lambda = 1e-3), lambda = ..., gamma = 1)` returns a
   sparse fit that selects (most of) the true support. No exact
   numerical claim; just confirm the path runs and selects sensibly.
5. **Sparse recovery superior to plain Lasso, best-`lambda` per side**
   (qualitative).
   Because GLMNet rescales `penalty_factor` to sum to `nvars`, a
   matched-`lambda` comparison between plain Lasso and Adaptive Lasso
   does not represent an apples-to-apples penalty strength when
   `gamma > 0`. This test therefore tunes `lambda` separately for
   each estimator and compares the *best* support-recovery
   outcome each can achieve.
   Synthetic problem: `n = 100`, `p = 15`, sparse ground truth with
   only 3 truly-nonzero coefficients well-separated from zero
   (magnitudes ~1).
   - Sweep `lambda` over a small geometric grid (e.g.
     `10 .^ range(-3, 0, length = 13)`) for each estimator.
   - For each `lambda`, record the recovered support
     `Set(findall(!iszero, b_fit))`.
   - Pass condition: there exists some `lambda` at which
     `AdaptiveLasso(pilot = OLS(), gamma = 1, ...)` recovers
     *exactly* the true support `{j: beta_true[j] != 0}`, while
     plain `Lasso(...)` recovers the true support at *no* `lambda`
     in the same grid (it always either omits a true column or
     keeps a spurious one).
   This formulation is robust to the rescale issue and is the
   standard "Adaptive beats Lasso at support recovery" claim in
   Zou (2006).

(Tests 6 and 7 from the parent spec -- standardise effect,
integration smoke -- are not duplicated here; standardise is
already exercised by tests 1-3, and a Lasso integration path is
already smoke-tested by spec 260517. An `AdaptiveLasso` integration
smoke is added as a single line in one existing integration test
that already calls `fit(SCEFit, ...)`.)

## Risks and open items

- **GLMNet rescales `penalty_factor` internally.** The Fortran kernel
  (and so `GLMNet.jl`) normalises the supplied `penalty_factor` so
  that `sum(pf) == nvars`, then runs the path against the rescaled
  vector. Consequences this spec must own:
  - The user-supplied `lambda` interacts with the *rescaled* weights,
    so for `gamma > 0` the same numerical `lambda` value does not
    represent the same penalty strength as in plain
    `Lasso(lambda = ...)`. The `gamma = 0` case is unaffected (all
    weights are `1`, rescale is the identity to within rounding), so
    test 1 below remains a bit-level agreement check.
  - When `pilot` produces exact zeros (e.g. `pilot = Lasso(...)`),
    those columns receive `pf_j = 1 / epsilon^gamma` (a very large
    number). After rescale, the effective `lambda` on the
    non-clipped columns shrinks toward zero, so those columns
    behave closer to an OLS fit than to a uniform-penalty Lasso.
    Test 3 below specifically verifies the clipped columns end up
    at exactly zero; the non-clipped columns' values are
    intentionally not constrained.
  - When `pilot = OLS()` runs on a rank-deficient design, Julia's
    minimum-norm solution populates the null space with noise of
    magnitude ~ 1e-10..1e-12. Those entries are *not* clipped by
    `eps(Float64)` (~ 2.2e-16), so they produce `pf_j ~ 1e10` and
    drive the same rescale collapse described above. This is the
    concrete failure mode the docstring's "use `pilot = Ridge(...)`
    on rank-deficient designs" recommendation prevents.
  An implementation alternative -- compute the weights, drop the
  `pilot == 0` columns from `X` entirely before calling GLMNet, then
  re-insert zeros at those slots -- is not pursued in this spec; the
  rescale behaviour is documented and the test plan works with it.
- **Pilot-cost vs. SCE fit cost.** Running an OLS or Ridge pilot
  inside `solve_coefficients(::AdaptiveLasso, ...)` doubles the fit
  cost vs. plain Lasso. For SCE the dominant cost is the SALC
  construction (cached in `SCEBasis`), so this doubling is
  negligible in absolute terms. No bench entry needed.
- **Lasso-as-pilot recursion.** `AdaptiveLasso(pilot = Lasso(...),
  ...)` parses and works (test 3 covers it). Nesting deeper -- e.g.
  `AdaptiveLasso(pilot = AdaptiveLasso(...), ...)` -- also parses
  but is out of scope; we do not test it and do not document it as
  a supported pattern.
- **`standardize = true` choice.** The adaptive penalty already
  absorbs the per-column scale via `1/|beta_pilot|^gamma`, so on
  paper `standardize = true` is partially redundant. We default to
  `true` for consistency with `Lasso` / `ElasticNet` and document
  the redundancy. Users who want the literal Zou 2006 path can
  pass `standardize = false`.
- **`pilot::AbstractEstimator` field type.** Holding an abstract-
  type field forfeits some type stability inside
  `solve_coefficients(::AdaptiveLasso, ...)`. `solve_coefficients`
  is a cold path (called once per fit), so the cost is negligible.
  The alternative -- a type parameter `AdaptiveLasso{P<:AbstractEstimator}`
  -- adds construction noise (`AdaptiveLasso{OLS}` instead of plain
  `AdaptiveLasso`) for no measurable benefit on this code path.
