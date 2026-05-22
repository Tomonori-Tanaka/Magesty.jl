# Requirements: `refit` — post-selection refit on the selected basis support

Status: draft (2026-05-22)

## Goal

Add a `refit` function that takes a fitted `SCEFit`, reads its selected
basis support (the SCE coefficients judged non-negligible), and re-solves
the fit restricted to that support with a chosen estimator. The dropped
bases keep coefficient `0.0`; the `SCEBasis` length and ordering are
unchanged.

## Background

The sparse estimators (`Lasso`, `AdaptiveLasso`, `AdaptiveRidge`) select a
basis subset, but the surviving coefficients are biased: L1 shrinkage
pulls them toward zero, and `AdaptiveRidge` applies a residual L2 penalty.
A standard remedy is post-selection refit (a.k.a. post-Lasso / debiasing):
fix the support chosen by the sparse fit, then re-estimate on that support
with an unpenalized (or lightly penalized) estimator.

A truncated basis is not removed from the instance — it persists in the
coefficient vector `jphi` as `0.0`. So `refit` selects the support by
thresholding the coefficient magnitudes of the input fit, and the dropped
bases simply stay at `0.0` in the output.

## Scope

Includes:

- New exported function `refit(fit::SCEFit, estimator; threshold)` in the
  `Fitting` pipeline, returning a fresh `SCEFit`.
- Selection criterion: scaled coefficient magnitude
  `abs(jphi[j]) * norm(X[:, j])`, where `X` is the assembled weighted
  design matrix.
- Documentation (`docs/src/api.md`, theory page) and `CHANGELOG.md`.

Excludes:

- A `Refit` estimator type / one-shot `fit(SCEFit, dataset, Refit(...))`
  path. `refit` is post-hoc on an existing `SCEFit` only.
- A reusable column-restriction estimator (`OnSupport`) or a public
  `support` helper. (Considered and deferred; see design.md.)
- `refit` on a bare `SCEModel`: it carries no dataset, so it cannot be
  re-solved.
- Cross-validated / automatic threshold selection.

## Invariants

- `SCEBasis` and the `jphi` length / ordering are unchanged; dropped
  bases stay at coefficient `0.0`. The `Fitting` <-> `SALCBasis`
  key-group order is untouched.
- `j0` is still recovered by `extract_j0_jphi` from the unscaled,
  un-centered energy data; it stays in the input energy unit and
  independent of `torque_weight`.
- No change to existing estimators or to `fit(SCEFit, ...)` — `refit`
  is purely additive.
- XML `save` / `load` round-trips are unaffected (a refit result is an
  ordinary `SCEFit` / `SCEModel`).

## Completion criteria

- [ ] `refit(fit, estimator; threshold)` implemented and exported.
- [ ] `make test-all` passes; new component tests cover the L1 exact-zero
      case, the `AdaptiveRidge` thresholded case, the empty-support case,
      and OLS idempotence (`refit` of an OLS fit at `threshold = 0`
      reproduces it).
- [ ] `make test-aqua` / `make test-jet` clean.
- [ ] New API reflected in `docs/src/api.md` and the fitting theory page.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.

## References

- Related specs: `260513-estimator-dispatch` (the `solve_coefficients`
  dispatch reused here), `260517-elasticnet-estimator`,
  `260518-adaptive-lasso-oneshot`, `260522-adaptive-ridge-iterative`.
