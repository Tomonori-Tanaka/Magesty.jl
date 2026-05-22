# Tasklist: `refit` — post-selection refit on the selected basis support

Status: implementing (2026-05-23)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — `refit` implementation

- [x] Add `refit(fit::SCEFit, estimator = OLS(); threshold, verbosity)`
      in `src/Magesty.jl`, reusing `assemble_weighted_problem` /
      `solve_coefficients` / `extract_j0_jphi`.
- [x] Scaled-magnitude support selection (`abs(jphi[j]) * norm(X[:,j])`),
      empty-support `@warn` path, `threshold >= 0` `ArgumentError`.
- [x] Export `refit`; docstring in standard Julia format.
- [x] `PrecomputedPilot`-backed estimator rejection with a clear
      `ArgumentError` (covers `AdaptiveLasso(::SCEFit)` too).

### M2 — Tests

- [x] `test/component/test_refit.jl`: OLS idempotence, masked-support
      OLS solve, L1 exact-zero support, empty support warn,
      `torque_weight` reuse, `PrecomputedPilot` rejection, negative
      `threshold` rejection.
- [x] Register the new test file in the component test runner.

### M3 — Docs

- [x] `docs/src/api.md` `@docs` block includes `refit`.
- [x] `SPEC.md`: `refit` listed in the public API section.
- [x] `docs/src/theory/design_matrix_and_fitting.md`: post-selection
      refit section + the scaled criterion.
- [x] `CHANGELOG.md` `[Unreleased]` entry.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] If results changed: regression or validation test added.
      (`refit` is additive; covered by new `test_refit.jl`.)
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ (No hot-path change.)
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [ ] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ (No such change.)
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
