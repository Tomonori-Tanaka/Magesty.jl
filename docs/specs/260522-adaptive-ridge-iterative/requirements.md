# Requirements: iterative Adaptive Ridge estimator

Status: draft (2026-05-22)

## Goal

Add an iterative Adaptive Ridge estimator (`AdaptiveRidge`) to the
`Fitting` estimator hierarchy. It approximates an L0-penalized fit
(Frommlet & Nuel 2016) by repeatedly refitting a per-coefficient
weighted ridge problem and updating the weights between iterations.

## Background

The L1 / mixed-norm / adaptive family already landed: `ElasticNet`
(with the `Lasso` convenience constructor), `AdaptiveLasso` (one-shot
reweighting), and the `PrecomputedPilot` adapter. The iterative
Adaptive Ridge variant was deferred because it adds an outer iteration
loop with its own convergence criteria. This spec implements it.

Design decisions settled with the user before drafting:

- **Solver backend: analytic weighted ridge.** Each refit solves the
  closed form `(XᵀX + lambda·Diagonal(w)) \ (Xᵀy)` directly. This is the
  same analytic family as the existing `Ridge`, which uses
  `MultivariateStats.ridge` (not GLMNet). `MultivariateStats.ridge`
  cannot take per-column weights, so the weighted solve is implemented
  inline. The analytic backend is exact and deterministic, has no solver
  noise floor, and avoids the GLMNet `penalty_factor` rescaling
  ambiguity that `AdaptiveLasso` carries.
- **Convergence: relative max coefficient change.** Stop when
  `norm(beta_new - beta, Inf) / max(norm(beta_new, Inf), eps(Float64))
  < tol`.
- **Trajectory: final solution only.** `solve_coefficients` keeps its
  `Vector{Float64}` return; iteration count / convergence surface via
  `@debug`.
- **No `mode` field.** Only one mode exists; a future one-shot variant
  can add `mode` then.
- **No `standardize` field.** The analytic backend has no GLMNet
  `standardize` knob, matching the existing `Ridge`.

## Scope

Includes:

- New `AdaptiveRidge <: AbstractEstimator` struct, keyword constructor,
  and `solve_coefficients(::AdaptiveRidge, X, y)` method in
  `src/Fitting.jl`.
- Export from `Fitting` and re-export from `Magesty`.
- Component + integration tests.
- `docs/src/api.md`, `SPEC.md`, `CHANGELOG.md` updates.

Excludes:

- Folding the existing analytic `Ridge` into a shared weighted-ridge
  kernel (separate spec, if ever).
- A one-shot Adaptive Ridge variant.
- Exposing the per-iteration coefficient trajectory.

## Invariants

- Existing XML round-trips byte-for-byte (`AdaptiveRidge` is an
  estimator recipe, never serialized).
- Spherical-harmonics convention unchanged.
- SALC ordering and the `Fitting` <-> `SALCBasis` design-matrix
  contract unchanged.
- The existing estimators (`OLS`, `Ridge`, `ElasticNet`, `Lasso`,
  `AdaptiveLasso`, `PrecomputedPilot`) and their numerical results are
  untouched.
- Public API is additive only.

## Completion criteria

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean.
- [ ] `AdaptiveRidge` reflected in `docs/src/api.md` and `SPEC.md`.
- [ ] Component test includes an analytically derived anchor
      (`AdaptiveRidge(lambda = L, epsilon = BIG)` collapses to
      `Ridge(lambda = L/BIG)`).
- [ ] Tier 2 four-axis review panel run and findings resolved.

## References

- Related specs: `260517-elasticnet-estimator`,
  `260518-adaptive-lasso-oneshot`,
  `260519-adaptive-lasso-precomputed-pilot`.
- Frommlet, F. & Nuel, G. (2016), "An adaptive ridge procedure for L0
  regularization", PLoS ONE 11(2): e0148620.
