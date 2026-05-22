# Tasklist: iterative Adaptive Ridge estimator

Status: complete (2026-05-22)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — `AdaptiveRidge` struct + `solve_coefficients`

- [x] Add the `AdaptiveRidge` struct, keyword constructor (with
      `ArgumentError` validation), labeled `Base.show`, and standard
      docstring in `src/Fitting.jl` (after `PrecomputedPilot`).
- [x] Add `solve_coefficients(::AdaptiveRidge, X, y)` with the analytic
      weighted-ridge iteration.
- [x] Export `AdaptiveRidge` from `Fitting` and re-export from
      `Magesty`.
- [x] `make test-all` (existing tests stay green).

### M2 — Tests

- [x] Component `@testset "AdaptiveRidge"` in
      `test/component/test_fitting_estimators.jl`: `lambda = 0` -> OLS,
      analytic `Ridge(L/BIG)` anchor, L0-shrinkage vs plain `Ridge`,
      convergence, `max_iter = 1` smoke, constructor validation.
- [x] Integration `AdaptiveRidge` smoke testset in
      `test/integration/febcc_2x2x2_pm/test.jl`.
- [x] `make test-all`.

### M3 — Docs and wrap-up

- [x] `docs/src/api.md` — add `AdaptiveRidge` to the Estimators block.
- [x] `docs/src/theory/design_matrix_and_fitting.md` — add `AdaptiveRidge`
      to the estimator table, the objective-function list, the parameter
      prose, the cross-estimator `lambda` note, and the references
      (Frommlet & Nuel 2016).
- [x] `docs/src/input_keys.md` — add `AdaptiveRidge` to the estimator
      enumeration.
- [x] `SPEC.md` — update the `# Estimators` API line.
- [x] `CHANGELOG.md` `[Unreleased] ### Added` entry.
- [x] `AdaptiveRidge` docstring cross-references the estimator family.
- [x] Delete `docs/design-notes/adaptive-ridge-iterative.md` and remove
      its row from `DESIGN_NOTES.md` (spec is now canonical).
- [x] Update `Status:` here and the row in `docs/specs/README.md`.

### M4 — Tier 2 review panel

- [x] Launch the four-axis review panel (numerical / maintainability /
      performance / API) in parallel on the final diff.
- [x] Apply all numerical-reviewer findings; resolve blockers / major;
      escalate any material perf-vs-maintainability tradeoff to the
      user. Panel result: 0 blockers; numerical axis clean (0 major);
      maintainability 1 major (over-length `@debug` line), API 2 major
      (`(4π)` consistency, sibling return annotations), performance 1
      major (per-iteration `p×p` allocation, tagged
      `[contention: maintainability]`). All applied except the API
      finding to retrofit `::Vector{Float64}` return annotations onto
      the pre-existing `OLS` / `Ridge` / `ElasticNet` / `AdaptiveLasso`
      `solve_coefficients` methods — held as out of scope (a separate
      cleanup, would muddy this diff). The performance finding was not
      escalated: `solve_coefficients` is not a hot path and the in-place
      fix carries negligible readability cost, so the tradeoff is not
      material.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes (22848 tests).
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] If results changed: regression or validation test added.
      (New estimator; analytic anchor + `Ridge` initializer pin it.)
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [x] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ (Not a hot path.)
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [x] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ (No module / Makefile change.)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [x] Implementation commit hash appended below.

## Commit log

- `231207f` docs(specs): add iterative AdaptiveRidge spec
- `b7469ee` feat(Fitting): add AdaptiveRidge estimator
- `446089a` test(Fitting): add AdaptiveRidge component and integration tests
- `7585307` docs: document AdaptiveRidge across api, theory, and changelog
