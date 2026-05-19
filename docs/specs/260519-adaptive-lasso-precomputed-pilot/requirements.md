# Requirements: AdaptiveLasso precomputed-pilot path

Status: draft (2026-05-19)

## Goal

Let `AdaptiveLasso` accept a pre-computed pilot coefficient vector (and,
by extension, an already-fitted `SCEFit` / `SCEModel`) instead of always
running a fresh pilot estimator inside `solve_coefficients`. This both
saves a redundant pilot fit when the user already has a usable model and
gives the future iterative Adaptive variant a natural inner-loop
primitive (each iteration's `coef_{k-1}` becomes the next pilot).

## Background

Currently `solve_coefficients(::AdaptiveLasso, X, y)` always calls
`solve_coefficients(e.pilot, X, y)` to obtain `beta_pilot`. The oneshot
spec (`260518-adaptive-lasso-oneshot`) intentionally kept this minimal,
but two follow-up needs surface immediately:

1. **Reuse a fitted model**: A user who has already fit `OLS` or `Ridge`
   for diagnostics shouldn't pay the pilot cost a second time when
   switching to `AdaptiveLasso`. On the production SCE pipeline the
   pilot is the dominant cost of the adaptive call.
2. **Iterative Adaptive variant**: The Frommlet & Nuel iterative L0
   approximation reuses the previous iteration's coefficients as the
   next pilot. A "precomputed pilot" hook in the oneshot path lets the
   iterative variant be a thin outer loop over the same kernel.

## Scope

Includes:

- New `PrecomputedPilot <: AbstractEstimator` carrying a fixed
  coefficient vector; `solve_coefficients(::PrecomputedPilot, X, y)`
  returns the stored vector after a length check against `size(X, 2)`.
- Convenience constructors `AdaptiveLasso(model::SCEModel; ...)` and
  `AdaptiveLasso(fit::SCEFit; ...)` that wrap `coef(model)` /
  `coef(fit)` in `PrecomputedPilot` and forward the remaining kwargs.
- Unit tests covering: pilot reuse parity (oneshot with
  `PrecomputedPilot(coef(prefit))` equals oneshot with the original
  estimator pilot, bit-for-bit), length-mismatch error path, and the
  convenience constructors.
- Doc updates: `AdaptiveLasso` docstring, `docs/src/api.md`,
  `CHANGELOG.md`, and the iterative-variant sketch in
  `docs/design-notes/`.

Excludes:

- The iterative Adaptive variant itself. That remains a separate spec;
  this spec only provides the primitive it will need.
- Any change to the default `AdaptiveLasso` behavior or to existing
  estimators. Numerics must be unchanged when no precomputed pilot is
  passed.
- Persisting `PrecomputedPilot` to XML or any other on-disk format.

## Invariants

- Default `AdaptiveLasso(lambda = ...)` (and all existing call sites)
  produces numerically identical coefficients before and after this
  change.
- SCE XML round-trip is unchanged (`PrecomputedPilot` is an estimator
  recipe, not a fitted artifact, and is never serialized).
- Physics conventions, spin-direction layout, and SALC ordering are
  untouched.
- `AbstractEstimator` continues to be the only type the
  `AdaptiveLasso.pilot` field accepts; no `Union` widening.

## Completion criteria

- [ ] `make test-all` passes.
- [ ] New component tests cover the pilot-reuse parity case, the
      convenience constructors, and the length-mismatch error.
- [ ] `AdaptiveLasso` docstring, `docs/src/api.md`, and `SPEC.md` list
      `PrecomputedPilot` and the new constructors.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] Iterative-variant design note updated to point at
      `PrecomputedPilot` as the reusable kernel.

## References

- Related specs: [`260518-adaptive-lasso-oneshot`](../260518-adaptive-lasso-oneshot/)
- Related design notes: `docs/design-notes/` (Lasso-family adaptive
  estimators sketch)
