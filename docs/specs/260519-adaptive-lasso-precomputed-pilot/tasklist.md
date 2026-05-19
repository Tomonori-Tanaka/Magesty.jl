# Tasklist: AdaptiveLasso precomputed-pilot path

Status: draft (2026-05-19)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 -- `PrecomputedPilot` + `solve_coefficients`

- [ ] Add the `PrecomputedPilot` struct (with copying inner ctor) and
      `solve_coefficients(::PrecomputedPilot, X, y)` method in
      `src/Fitting.jl`.
- [ ] Export `PrecomputedPilot` from `Fitting` and re-export from
      `Magesty`.
- [ ] Run `make test-all` (existing tests must stay green; no new
      tests yet).

### M2 -- Convenience constructors

- [ ] Add `AdaptiveLasso(model::SCEModel; kwargs...)` and
      `AdaptiveLasso(fit::SCEFit; kwargs...)` in `src/Magesty.jl`,
      placed after the `coef(::SCEFit)` / `coef(::SCEModel)`
      definitions (Fitting cannot host these: SCEFit / SCEModel are
      defined in the parent module *after* `include("Fitting.jl")`).
- [ ] Run `make test-all`.

### M3 -- Tests

- [ ] Component tests in `test/component/test_fitting_estimators.jl`:
      pilot-reuse parity (bit-exact), length-mismatch error,
      convenience-ctor smoke.
- [ ] Integration test extension in
      `test/integration/febcc_2x2x2_pm/test.jl`: reuse a fitted
      `SCEFit` via `AdaptiveLasso(fitted; lambda = ...)`.
- [ ] `make test-all`.

### M4 -- Docs and wrap-up

- [ ] `AdaptiveLasso` docstring updated to mention the convenience
      constructors and `PrecomputedPilot`.
- [ ] `docs/src/api.md` -- add `PrecomputedPilot` to the Estimators
      block.
- [ ] `SPEC.md` -- update the `# Estimators` API line.
- [ ] `CHANGELOG.md` `[Unreleased] ### Added` entry.
- [ ] Iterative-variant design note: point at `PrecomputedPilot` as
      the reusable kernel.
- [ ] Update `Status:` in this file and the row in
      `docs/specs/README.md`.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] ~~If results changed: regression or validation test added.~~
      (no numerical change; pilot-reuse parity test added anyway as a
      regression anchor.)
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ (no hot-path change.)
- [ ] If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.

## Commit log

<!-- Appended as milestones land. -->
