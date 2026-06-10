# Tasklist: GCV diagnostics (hat-matrix generalized cross-validation)

Status: implemented, pending commit (2026-06-10)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Core numerics (`Fitting.jl`)

- [x] `_is_linear_estimator` predicate; `ArgumentError` for non-linear.
- [x] `_gcv_single` (single score) and `_gcv_lambda_path` (SVD closed form).
- [x] `_effective_dof` (OLS / Ridge / AdaptiveRidge conditional dof).
- [x] `_gcv_sample_count` (live-row count + conditional intercept dof; review
      fix for the zeroed-block case at `torque_weight` 0 / 1).

### M2 — Public API (`Magesty.jl`)

- [x] `GCVLambdaPath` (`gcv_scores` field) / `GCVSizeCurve` structs.
- [x] `gcv(f)`, `gcv_lambda(dataset, lambdas; ...)`,
      `gcv_learning_curve(dataset, estimator; ...)` with docstrings + exports.
- [x] Default `sizes` grid (6 points, `max(p+2,10)`..`length(dataset)`);
      `torque_weight ∈ [0, 1]` validation on both sweep entry points.

  A bad draw (rank-deficient OLS, or `N − dof ≤ 0`) records `NaN` for that
  repeat and warns, rather than emitting a bogus value.

### M3 — Output + docs

- [x] `write_gcv_lambda` / `write_gcv_learning_curve` in `FitCheckIO.jl`.
- [x] `tools/FitCheck_gcv_lambda.py` / `tools/FitCheck_gcv_learning_curve.py`.
- [x] Example script (`examples/04_gcv_diagnostics.jl`); `docs/src` narrative
      diagnostics section + `api.md` `@docs` + `examples.md`.

### M4 — Tests + review

- [x] `test/component/test_gcv.jl` (analytic / cross-check suite, 108 tests).
- [x] `make test-all` / `test-jet` / `test-aqua` green.
- [x] Tier 2 review panel (numerical / maintainability / performance / API);
      findings resolved (live-row `N` blocker, two-block test, renames,
      `torque_weight` validation, error-message attribution, style fixes).

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes (23,619).
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] If results changed: regression or validation test added. (Only new GCV
      numbers, validated against closed-form / dense recomputation.)
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [x] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`. (Not a fit hot path; baseline recorded anyway.)
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [x] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ (None changed.)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync. (Done on merge.)
- [ ] Implementation commit hash appended below. (Added at commit time.)
