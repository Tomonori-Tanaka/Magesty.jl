# Requirements: ElasticNet estimator (with Lasso convenience)

Status: draft (2026-05-17, rebased on top of 260518 energy-centered design matrix)

## Goal

Add an `ElasticNet` estimator to the `AbstractEstimator` dispatch
hierarchy in `Fitting.jl`, backed by GLMNet.jl, covering Lasso
(`alpha=1`), GLMNet-style L2 (`alpha=0`), and honest Elastic Net
(mixed norm). Provide a `Lasso(; lambda, standardize=true)`
convenience constructor that returns an `ElasticNet`. The existing
analytic-solver `Ridge` is kept; this spec also adds an agreement
test between `ElasticNet(alpha=0, lambda=λ)` and `Ridge(lambda=λ)`
to document the relationship.

`ElasticNet` must produce correct `jphi` for the **full mixed
energy/torque objective** (`0 ≤ weight ≤ 1`), not only the
energy-only case. Since spec 260518 (merged 2026-05-18), the
augmented design matrix from `assemble_weighted_problem` no longer
carries a bias column — `j0` is eliminated analytically by
mean-centering the energy block before the solve, and `extract_j0_jphi`
recovers `j0` afterward from the un-scaled energy residual. The
`ElasticNet` method therefore reduces to a direct GLMNet call with
`intercept=false`; no extra centering or column-stripping is needed
in the estimator's `solve_coefficients`.

`ElasticNet` defaults to `standardize=true` to neutralise the
per-cluster `(4π)^(N/2)` prefactor baked into the SCE design matrix,
which would otherwise bias L1 (or mixed) selection against high-N
clusters.

## Background

The current `Fitting.jl` supports `OLS` and `Ridge` only. Sparse
effective spin models require L1-style regularization. A prior design
note surveyed Lasso / Adaptive LASSO / Adaptive Ridge together
(selecting GLMNet.jl as the backend); this spec narrows the first
deliverable to plain `ElasticNet` (with a `Lasso` convenience
constructor) so the GLMNet plumbing and standardisation policy can be
settled and validated before the Adaptive variants are added.

Historical context: an earlier `ElasticNet(alpha, lambda)` struct
existed but silently ignored `alpha`, so it was renamed to `Ridge` in
[`260513-estimator-dispatch`](../260513-estimator-dispatch/), with the
explicit intent of freeing the `ElasticNet` name for a real mixed-norm
estimator later. This spec is that follow-up.

Why `standardize=true` by default (departs from `Ridge`'s implicit
"no standardisation"): SCE design-matrix columns carry a per-cluster
prefactor `(4π)^(N/2)` (≈ 3.5 for N=1 up to ≈ 158 for N=4). A uniform
L1 penalty on unscaled columns under-penalises high-N clusters and
biases selection. The analytic `Ridge` solver is far less sensitive to
column scale (continuous shrinkage), so its default stays unchanged.

Why `intercept=false` is the right call (post-260518): the energy
block of `X` arrives at `solve_coefficients` already centered by
`assemble_weighted_problem`. Asking GLMNet to *also* fit an
intercept would re-introduce a uniform offset across both energy
**and** torque rows, contradicting the physics "torques don't see
`j0`" and biasing `jphi`. The correct factorisation is the one the
codebase already uses for OLS / Ridge: fit `jphi` against the
centered, weighted system, and re-fit `j0` afterward from the
un-scaled energy data via `extract_j0_jphi`.

## Scope

Includes:

- New estimator `ElasticNet <: AbstractEstimator` with fields
  `alpha::Float64`, `lambda::Float64`, `standardize::Bool`, validated
  by a keyword constructor that enforces `0 ≤ alpha ≤ 1`.
- Convenience constructor `Lasso(; lambda, standardize=true)` that
  returns an `ElasticNet(alpha=1.0, ...)`. No new struct.
- `solve_coefficients(::ElasticNet, X, y)` method that calls a private
  `_glmnet_solve` helper with GLMNet's `intercept=false` and returns
  the resulting coefficient vector unchanged. No column stripping, no
  in-method centering: the centered, no-bias-column `X` produced by
  `assemble_weighted_problem` is exactly what GLMNet should see.
- `using GLMNet` in `src/Fitting.jl`; `GLMNet` added to `Project.toml`
  `[deps]` and `[compat]`.
- Re-export of `ElasticNet` and `Lasso` from `src/Magesty.jl`.
- Component tests covering: OLS limit (`λ → 0`), GLMNet-`α=0` vs
  analytic-Ridge agreement, Lasso sparsity monotonicity,
  standardisation effect, and **`weight > 0` correctness**
  (the mixed energy/torque case).
- At least one integration smoke test exercising the public
  `fit(SCEFit, dataset, Lasso(lambda=...))` pipeline.
- `docs/src/api.md` and `CHANGELOG.md` `[Unreleased]` updates.

Excludes:

- Adaptive LASSO, Adaptive Ridge `:oneshot` / `:iterative` (deferred
  to a follow-up spec; pre-spec sketch retained in the design note).
- Removing or replacing the existing analytic `Ridge`. This spec
  documents the GLMNet `α=0` vs analytic-Ridge relationship via a
  test; any deprecation decision is a separate spec.
- Lambda-path / cross-validated `λ` selection.
- Folding the existing `Ridge` into GLMNet.
- Changes to `assemble_weighted_problem` or `extract_j0_jphi`. Both
  are already in their post-260518 form (no bias column, energy
  block mean-centered).

## Invariants

- `OLS` and `Ridge` numerical behaviour are byte-for-byte unchanged
  (no edits to their `solve_coefficients` methods).
- `assemble_weighted_problem` returns the same `(X, y)` as today (no
  signature change; `ElasticNet` consumes the same pair as OLS /
  Ridge).
- `extract_j0_jphi` is unchanged. `ElasticNet` returns a `jphi`
  vector of length `num_salcs`; `extract_j0_jphi` then computes
  `j0 = mean(y_e - X_e · jphi)` against the un-scaled energy data —
  identical post-processing to OLS / Ridge.
- Existing XML round-trips for `SCEBasis` / `SCEModel` continue
  byte-for-byte; `ElasticNet` only affects the fit step.
- `Jφ` coefficients returned by `ElasticNet` are in the input energy
  unit (typically eV); GLMNet's internal standardisation is fully
  reversed before the coefficients leave `solve_coefficients`.
- Torque rows do not see `j0`: this is enforced by the
  `assemble_weighted_problem` centering (set up in spec 260518) plus
  GLMNet's `intercept=false`. See `design.md` for the derivation.
- Spin direction stays unit-vector with `3 × n_atoms` layout
  (untouched).
- Tesseral spherical-harmonics conventions (untouched).
- SALC ordering and the `Fitting` ↔ `SALCBasis` ↔ XML link (see
  CLAUDE.md "Linked sites") are untouched.

## Completion criteria

- [ ] `make test-all` passes.
- [ ] `make test-jet`, `make test-aqua` clean (no new warnings).
- [ ] New component tests in `test/component/test_fitting_estimators.jl`
      (or appended) cover the five `ElasticNet` properties listed in
      the design's "Test strategy" section.
- [ ] GLMNet `α=0` vs analytic `Ridge` agreement test passes with the
      `rtol` value decided during implementation, and that value is
      recorded in the test docstring and in `CHANGELOG.md`.
- [ ] The `weight > 0` correctness test passes: `Lasso(λ→0)` with
      `weight=0.5` matches `OLS()` with `weight=0.5` (within the
      same `atol` as the energy-only OLS-limit test), confirming
      the centered-block + `intercept=false` path is correct.
- [ ] One integration test path exercises `Lasso(...)` end-to-end via
      the public `fit` API.
- [ ] `ElasticNet` and `Lasso` appear in `docs/src/api.md` alongside
      `OLS` / `Ridge`.
- [ ] `CHANGELOG.md` `[Unreleased]` records the new estimator, the
      `GLMNet` dependency, and the Ridge-agreement test outcome.
- [ ] `Status:` line in this spec's `tasklist.md` and the matching
      row in `docs/specs/README.md` updated to `complete` together.

## References

- Related spec (analytic `j0` elimination + bias-column removal):
  [`260518-energy-centered-design-matrix/`](../260518-energy-centered-design-matrix/)
- Related spec (prior estimator-dispatch refactor):
  [`260513-estimator-dispatch/`](../260513-estimator-dispatch/)
- GLMNet.jl: https://github.com/JuliaStats/GLMNet.jl
- Friedman, Hastie, Tibshirani (2010), "Regularization Paths for
  Generalized Linear Models via Coordinate Descent", J. Stat. Softw.
