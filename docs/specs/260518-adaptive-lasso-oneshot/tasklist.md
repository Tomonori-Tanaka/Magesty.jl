# Tasklist: Adaptive Lasso (oneshot) estimator

Status: complete (2026-05-19)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 -- Struct, helper extension, and dispatch (one commit)

The `AdaptiveLasso` struct, the `_glmnet_solve` `penalty_factor`
extension, and the `solve_coefficients(::AdaptiveLasso, ...)` method
are tightly coupled. Commit them together.

- [x] Add `struct AdaptiveLasso` with fields `pilot::AbstractEstimator`,
      `lambda::Float64`, `gamma::Float64`, `epsilon::Float64`,
      `standardize::Bool`. Keyword ctor validates `lambda >= 0`,
      `gamma >= 0`, `epsilon > 0`. Defaults: `pilot = OLS()`,
      `gamma = 1.0`, `epsilon = eps(Float64)`, `standardize = true`.
- [x] Extend `_glmnet_solve` in `src/Fitting.jl` with
      `penalty_factor::Union{Nothing, AbstractVector{<:Real}} = nothing`.
      When `nothing`, omit the `penalty_factor` argument from the
      `glmnet` call so the existing `ElasticNet` / `Lasso` path is
      byte-for-byte unchanged.
- [x] Implement `solve_coefficients(::AdaptiveLasso, X, y)` per the
      design (pilot fit -> weights `inv.(max.(abs.(beta_pilot),
      epsilon) .^ gamma)` -> `_glmnet_solve(...; penalty_factor)`
      with `alpha = 1.0`, `intercept = false`).
- [x] Export `AdaptiveLasso` from `Fitting`; re-export from `Magesty`.
- [x] `make test-unit` green: existing OLS / Ridge / ElasticNet /
      Lasso paths remain identical numerically (no edits to their
      methods; the helper's `penalty_factor = nothing` default keeps
      them on the same GLMNet call path).

### M2 -- Component tests

- [x] Tests 1-5 from `design.md` "Test strategy" added inside a new
      `@testset "AdaptiveLasso estimator" begin ... end` block in
      `test/component/test_fitting_estimators.jl`. Tests 1 (gamma=0
      = plain Lasso agreement) and 2 (exact weight construction) are
      the load-bearing strict tests; tests 3-5 cover the
      Lasso-pilot epsilon clip, the Ridge-pilot path, and the
      sparse-recovery superiority claim. Test 3 was weakened during
      implementation: GLMNet's `sum(pf) = nvars` rescale puts
      pilot-zero columns at moderate effective `lambda` rather than
      infinite, so "pilot-zero implies adaptive-zero" is data-
      dependent; the test now pins finite-output and alive-pilot-
      column survival instead. `design.md` was updated to match.
- [x] Decide and record the `gamma = 0` agreement tolerance in the
      test docstring. Recorded as `atol = 1e-6` (empirical worst
      ~3e-8, ~30x headroom).
- [x] `make test-unit` green (21738 tests).

### M3 -- Integration smoke and docs

- [x] Add a one-line `AdaptiveLasso(...)` smoke path inside
      `test/integration/febcc_2x2x2_pm/test.jl` (the multi-cluster,
      well-conditioned fixture; chosen because dimer / chain are
      too sparse to exercise sparsity behaviour and FeGe is
      rank-deficient). Confirm `fit(SCEFit, dataset,
      AdaptiveLasso(...))` runs and `Magesty.save` / `Magesty.load`
      round-trip the resulting model. No existing integration test
      currently exercises `Lasso` or `ElasticNet`; the smoke is
      additive.
- [x] Update `docs/src/api.md` to list `AdaptiveLasso` alongside
      `ElasticNet` / `Lasso`.
- [x] Update `SPEC.md`: extend the `Fitting` module description and
      the `# Estimators` API line. No new entry in the Primary
      external libraries table (GLMNet is already there).
- [x] Append a `[Unreleased]` entry in `CHANGELOG.md` covering the
      new estimator and the `_glmnet_solve` extension, with the
      `gamma = 0` agreement tolerance recorded.
- [x] `make test-integration` green (817 tests).

### M4 -- Wrap-up

- [x] `make test-jet` clean.
- [x] Update `Status:` line in this file -> `complete (2026-05-19)`
      and the matching row in `docs/specs/README.md` in the same
      commit.
- [x] Remove the now-superseded Adaptive Lasso (oneshot) section
      from `docs/design-notes/lasso-adaptive-estimators.md` -- its
      content is fully absorbed into this spec, and lingering
      basename references in long-lived docs are discouraged.
      Retain the `:iterative` Adaptive Ridge sketch only, and
      update the note's status line to reflect "oneshot complete,
      iterative next-up".

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes (21738 unit + 817 integration).
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] `gamma = 0` agreement test (Adaptive vs. plain Lasso) passes;
      tolerance recorded in test docstring and `CHANGELOG.md`
      (`atol = 1e-6`).
- [x] Sparse-recovery test passes: on a correlated-feature synthetic
      problem (`0.85 * target + 0.15 * noise`), `AdaptiveLasso(pilot
      = OLS(), gamma = 1)` recovers the true support exactly at some
      `lambda` in the geometric grid while plain `Lasso` misses it
      everywhere on the same grid.
- [x] If results changed: regression or validation test added.
      `AdaptiveLasso` is additive -- existing OLS / Ridge /
      ElasticNet / Lasso paths are untouched, so no regression test
      for them; new tests cover `AdaptiveLasso` directly.
- [x] Public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [x] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ Not a hot path; fit cost roughly
      doubles vs. plain Lasso (one extra pilot fit per call) but
      `solve_coefficients` runs once per `fit(SCEFit, ...)` and the
      SCE-side cost dwarfs both. No bench entry.
- [x] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ No module or Makefile
      change; `.claude/agents/` lists are descriptive ("OLS / Ridge
      fitting" in `test-runner.md`) rather than exhaustive, so no
      update was forced.
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [x] Implementation commit hash appended below.

## Implementation commits

- `cf1a6e2`  docs(specs): add Adaptive Lasso (oneshot) estimator spec
- `5a11206`  feat(Fitting): add AdaptiveLasso oneshot estimator
- `3149596`  test(Fitting): add AdaptiveLasso component tests
- `1c66a1d`  docs: document AdaptiveLasso estimator and add integration smoke
