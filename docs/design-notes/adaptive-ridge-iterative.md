# Adaptive Ridge (iterative)

**Status**: on hold (2026-05-19)

Pre-spec sketch for an iterative Adaptive Ridge estimator (Frommlet &
Nuel 2016 L0 approximation) layered on the existing GLMNet-based
estimator hierarchy in `src/Fitting.jl`.

## Background

The L1 / mixed-norm / Adaptive families have already landed:

- `ElasticNet` (with a `Lasso` convenience constructor) via GLMNet.jl —
  spec [`260517-elasticnet-estimator`](../specs/260517-elasticnet-estimator/).
- `AdaptiveLasso` (oneshot, single-stage reweighting) using
  `_glmnet_solve`'s `penalty_factor` plumbing — spec
  [`260518-adaptive-lasso-oneshot`](../specs/260518-adaptive-lasso-oneshot/).
- A `PrecomputedPilot <: AbstractEstimator` adapter that returns a fixed
  coefficient vector from `solve_coefficients` — spec
  [`260519-adaptive-lasso-precomputed-pilot`](../specs/260519-adaptive-lasso-precomputed-pilot/).

`PrecomputedPilot` is the primitive the iterative variant calls inside
its inner loop: each iteration's coefficients become the next
`PrecomputedPilot.beta`, and the existing oneshot kernel handles the
GLMNet refit. The iterative Adaptive Ridge variant remains out of scope
for those specs because it adds an outer iteration loop with its own
convergence criteria.

## Estimator mapping

| Estimator | `alpha` | `penalty_factor[j]` | Note |
|---|---|---|---|
| AdaptiveRidge (`mode = :iterative`) | 0.0 | `1/(β_j² + ε)` updated iteratively | Frommlet & Nuel 2016. Refits GLMNet each iteration as an L0 approximation |

The behavior-descriptive symbol `:iterative` is preferred over an
author-name symbol, since the latter assumes the reader knows the paper.

## Open items for the follow-up

- `AdaptiveRidge(:iterative)` defaults: `epsilon = 1e-8`,
  `max_iter = 50`, `tol = 1e-6`.
- Convergence criterion: max coefficient change vs RSS change vs
  effective-support change (pick one, document the choice).
- Whether to expose the per-iteration coefficient trajectory as a
  debugging aid, or return only the final fit.
- Whether to default to `standardize = true` (consistent with
  AdaptiveLasso) — likely yes, by the same per-cluster `(4π)^(N/2)`
  column-scale argument.
- Whether to fold the existing analytic `Ridge` into GLMNet later
  (separate spec). The agreement test in spec `260517` already confirmed
  that GLMNet `α=0, λ_g = λ_a/n` agrees with analytic `Ridge(λ_a)` to
  `rtol ~ 1e-3`.

## Workflow

Mid-sized feature plus multiple design decisions — a spec is needed
before implementation. Implementation runs on a `refactor/<slug>` branch
([[feedback_branch_for_medium_refactor]]).
