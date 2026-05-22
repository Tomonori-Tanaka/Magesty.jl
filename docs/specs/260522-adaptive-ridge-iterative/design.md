# Design: iterative Adaptive Ridge estimator

Status: draft (2026-05-22)

## Summary

Add `AdaptiveRidge <: AbstractEstimator`, an iterative L0-approximation
estimator following Frommlet & Nuel (2016). The fit is a fixed-point
iteration over a per-coefficient weighted ridge problem

    min_b  ||y - X b||^2 + lambda * sum_j w_j * b_j^2

with weights updated each iteration as `w_j = 1 / (b_j^2 + epsilon)`.
Large coefficients get a light penalty, small coefficients a heavy one;
iterating drives the small ones toward zero, approximating an L0
penalty.

Each weighted ridge subproblem is solved by the **analytic closed form**

    b = (XᵀX + lambda * Diagonal(w)) \ (Xᵀy)

This is the same analytic family as the existing `Ridge` estimator
(which calls `MultivariateStats.ridge`, not GLMNet).
`MultivariateStats.ridge` only accepts a scalar penalty, so the
weighted solve is written inline. `XᵀX` and `Xᵀy` are formed once and
reused across iterations; only the diagonal changes.

Alternatives considered and rejected:

- **GLMNet `alpha = 0` with `penalty_factor = w`.** Consistent in
  *backend* with `ElasticNet` / `AdaptiveLasso`, but GLMNet's ridge path
  carries ~1e-3 solver noise (an agreement test in spec `260517` only
  confirmed `rtol ~ 1e-3` against analytic `Ridge`), which floors a
  `tol = 1e-6` outer convergence test. GLMNet also rescales
  `penalty_factor` so the supplied weights sum to `nvars`, shifting the
  effective `lambda` — the same wrinkle documented in the `AdaptiveLasso`
  docstring. The analytic backend avoids both.
- **`PrecomputedPilot` inner-loop reuse.** The on-hold design note
  sketched feeding each iteration's coefficients back through
  `PrecomputedPilot` + the `AdaptiveLasso` one-shot kernel. That kernel
  is L1 (`alpha = 1`) with weights `1/max(|b|,eps)^gamma`; Adaptive Ridge
  needs L2 with weights `1/(b^2+eps)`. The kernels do not coincide, so
  `AdaptiveRidge` owns its iteration loop and `PrecomputedPilot` is not
  involved.

## Module layout

`AdaptiveRidge` has only scalar hyperparameters and no dependency on
`SCEFit` / `SCEModel`, so the whole feature lives inside
`src/Fitting.jl` next to the other estimators. No parent-module
constructors are needed (unlike `AdaptiveLasso`'s `SCEFit` / `SCEModel`
convenience ctors).

| Target | Change |
|---|---|
| `src/Fitting.jl` | Add `AdaptiveRidge` struct + keyword ctor (after `PrecomputedPilot`, ~line 295), labeled `Base.show`, and `solve_coefficients(::AdaptiveRidge, X, y)` (near the other `solve_coefficients` methods, ~line 1191). Add `AdaptiveRidge` to the **`Fitting` `export` line** (currently `src/Fitting.jl:23`). |
| `src/Magesty.jl` | Add `AdaptiveRidge` to the **`Magesty` public `export` line** (currently `src/Magesty.jl:82`) — a separate, distinct line from the `Fitting` export above. |
| `test/component/test_fitting_estimators.jl` | New `@testset "AdaptiveRidge"`. |
| `test/integration/febcc_2x2x2_pm/test.jl` | Add an `AdaptiveRidge` end-to-end smoke testset. |
| `docs/src/api.md` | Add `AdaptiveRidge` to the `## Estimators` `@docs` block. |
| `SPEC.md` | Update the `# Estimators` API line. |
| `CHANGELOG.md` | New `[Unreleased] ### Added` entry. |
| `DESIGN_NOTES.md` + `docs/design-notes/adaptive-ridge-iterative.md` | Point the index row at this spec; delete the design-note body on completion. |

## API

```julia
# src/Fitting.jl
struct AdaptiveRidge <: AbstractEstimator
    lambda::Float64
    epsilon::Float64
    max_iter::Int
    tol::Float64
end

function AdaptiveRidge(; lambda::Real,
                         epsilon::Real = 1e-8,
                         max_iter::Integer = 50,
                         tol::Real = 1e-6)
    lambda  >= 0.0 || throw(ArgumentError("AdaptiveRidge lambda must be non-negative; got $lambda"))
    epsilon  > 0.0 || throw(ArgumentError("AdaptiveRidge epsilon must be strictly positive; got $epsilon"))
    max_iter >= 1  || throw(ArgumentError("AdaptiveRidge max_iter must be at least 1; got $max_iter"))
    tol      > 0.0 || throw(ArgumentError("AdaptiveRidge tol must be strictly positive; got $tol"))
    return AdaptiveRidge(Float64(lambda), Float64(epsilon), Int(max_iter), Float64(tol))
end

function solve_coefficients(e::AdaptiveRidge,
                            X::AbstractMatrix{<:Real},
                            y::AbstractVector{<:Real})::Vector{Float64}
    # ... iteration body — see the algorithm section below.
end
```

Defaults follow the on-hold design note: `epsilon = 1e-8`,
`max_iter = 50`, `tol = 1e-6`. The `AdaptiveRidge` struct docstring uses
the standard Julia struct format (`# Fields`, `# Examples`) and
documents the per-cluster column-scale caveat shared with `Ridge` (no
`standardize`). The `solve_coefficients(::AdaptiveRidge, ...)` method
docstring follows the `# Arguments` / `# Returns` format, consistent
with the existing `solve_coefficients` docstrings for `Ridge` /
`ElasticNet` / `AdaptiveLasso` in `src/Fitting.jl`.

### `solve_coefficients` algorithm

1. If `e.lambda ≈ 0.0`, return `X \ y` — the penalty term vanishes and
   the iteration is a no-op (matches the `Ridge` short-circuit).
2. `Xf = Matrix{Float64}(X)`, `yf = Vector{Float64}(y)`,
   `XtX = Symmetric(Xf' * Xf)`, `Xty = Xf' * yf`. Formed once.
3. Iteration 0 (plain ridge initializer):
   `beta = (XtX + e.lambda * I) \ Xty`. Numerically equals
   `Ridge(e.lambda)`.
4. For `iter in 1:e.max_iter`:
   - `w_j = 1 / (beta_j^2 + e.epsilon)` — finite and strictly positive,
     so `XtX + e.lambda * Diagonal(w)` is symmetric positive definite
     whenever `e.lambda > 0` (PSD `XtX` plus a positive diagonal). No
     rank-deficiency handling is needed.
   - `beta_new = (XtX + e.lambda * Diagonal(w)) \ Xty`.
   - `rel = norm(beta_new - beta, Inf) / max(norm(beta_new, Inf), eps(Float64))`.
   - `beta = beta_new`; `break` when `rel < e.tol`.
5. `@debug` the iteration count and final `rel`; return `beta`.

## Types and conventions

- Fields are `Float64` / `Int`, converted in the keyword constructor —
  consistent with `Ridge` / `ElasticNet` / `AdaptiveLasso`.
- No physics-convention surface area is touched. `solve_coefficients`
  consumes the already-centered, already-weighted `(X, y)` from
  `assemble_weighted_problem`; `j0` is recovered downstream by
  `extract_j0_jphi`, exactly as for the other estimators.
- `lambda` carries the inverse-energy-squared scaling implied by the
  `(XᵀX + lambda·diag(w))` form, the same convention as the existing
  `Ridge` (`MultivariateStats.ridge(X, y, lambda)`).
- **New documented caveat** (not a new invariant): like `Ridge`,
  `AdaptiveRidge` does not standardize columns, so the per-cluster
  `(4π)^(N/2)` scale enters the weights. The docstring states this.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): unaffected.
- [ ] SCE coefficient XML (`save` / `load`): unaffected — `AdaptiveRidge`
      is an estimator recipe, never serialized.
- [ ] `Fitting` <-> `SALCBasis`: unaffected — the new code path runs
      strictly inside `solve_coefficients`; the design matrix and SALC
      ordering are untouched.
- [ ] `.claude/agents/` references: no module / Makefile change; sweep
      for stale estimator-list claims.
- [ ] `SPEC.md` / `docs/src/api.md` updates: add `AdaptiveRidge`.

## Test strategy

Component (`test/component/test_fitting_estimators.jl`,
`@testset "AdaptiveRidge"`):

1. **`lambda = 0` -> OLS.** `solve_coefficients(AdaptiveRidge(lambda = 0),
   X, y) == X \ y`.
2. **Analytic anchor (derived, not observed).** With
   `epsilon` chosen so `epsilon >> beta_j^2` for every `j`, all weights
   collapse to `≈ 1/epsilon`, so the converged fit equals plain
   `Ridge(lambda = lambda/epsilon)`:
   `(XᵀX + lambda·Diagonal(1/epsilon)) = (XᵀX + (lambda/epsilon)·I)`.
   Assert `AdaptiveRidge(lambda = L, epsilon = BIG)` agrees with
   `Ridge(lambda = L/BIG)` to a tight `rtol`.
3. **L0-approximation behavior.** On a fixture with a few large and
   several near-zero true coefficients, the converged `AdaptiveRidge`
   shrinks the small coefficients strictly more than plain `Ridge` at
   matched `lambda` (effective support no larger than plain ridge's).
4. **Convergence.** On a well-conditioned fixture the iteration
   converges before `max_iter`; all coefficients finite.
5. **`max_iter = 1` smoke.** One reweight yields a definite, finite
   result.
6. **Constructor validation.** Bad `lambda` / `epsilon` / `max_iter` /
   `tol` each throw `ArgumentError`.

Integration (`test/integration/febcc_2x2x2_pm/test.jl`): a new
`AdaptiveRidge` smoke testset fits a real `SCEFit` end-to-end and
asserts `coef` is finite, alongside the existing `AdaptiveLasso` smoke.

## Risks and open items

- **Numerical results.** This is a new estimator; no existing result
  changes. The analytic anchor (test 2) and the `Ridge` initializer
  (iteration 0) pin the new path to the existing `Ridge` numerics.
- **Non-convergence.** The fixed-point map is not guaranteed to be a
  contraction for all `(X, lambda, epsilon)`. On non-convergence within
  `max_iter` the method returns the last iterate and emits `@debug`;
  callers tune `max_iter` / `tol`. This matches the design-note intent
  (return only the final fit).
- **Column scaling.** Shared with `Ridge`: no standardization, so the
  per-cluster `(4π)^(N/2)` factor enters the weights. Documented in the
  docstring; not addressed here.
- **Cost.** Each iteration is a `p×p` solve (`p = num_salcs`) on a
  matrix that is rebuilt only in its diagonal; `XᵀX` / `Xᵀy` are formed
  once. Not a hot path (the basis-evaluation loop is untouched), so no
  `bench_log.md` entry is required.
