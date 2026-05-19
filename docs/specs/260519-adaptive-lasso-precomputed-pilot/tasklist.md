# Tasklist: AdaptiveLasso precomputed-pilot path

Status: complete (2026-05-19)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 -- `PrecomputedPilot` + `solve_coefficients`

- [x] Add the `PrecomputedPilot` struct (with copying inner ctor) and
      `solve_coefficients(::PrecomputedPilot, X, y)` method in
      `src/Fitting.jl`.
- [x] Export `PrecomputedPilot` from `Fitting` and re-export from
      `Magesty`.
- [x] Run `make test-all` (existing tests must stay green; no new
      tests yet).

### M2 -- Convenience constructors

- [x] Add `AdaptiveLasso(model::SCEModel; kwargs...)` and
      `AdaptiveLasso(fit::SCEFit; kwargs...)` in `src/Magesty.jl`,
      placed after the `coef(::SCEFit)` / `coef(::SCEModel)`
      definitions (Fitting cannot host these: SCEFit / SCEModel are
      defined in the parent module *after* `include("Fitting.jl")`).
- [x] Run `make test-all`.

### M3 -- Tests

- [x] Component tests in `test/component/test_fitting_estimators.jl`:
      pilot-reuse parity (bit-exact), defensive-copy at construction,
      length-mismatch `DimensionMismatch` error. Convenience-ctor smoke
      lives in the integration test (real `SCEFit` / `SCEModel` are
      not constructible from a component test fixture).
- [x] Integration test extension in
      `test/integration/febcc_2x2x2_pm/test.jl`: reuse a fitted
      `SCEFit` via `AdaptiveLasso(fitted; lambda = ...)` end-to-end,
      plus duplicate-`pilot`-kwarg guard.
- [x] `make test-all`.
- [x] Discovered during M3 testing: Julia's kwarg-splat semantics
      *silently* override the precomputed pilot when a `pilot`
      keyword is passed in addition to the positional fit/model
      (the spec had assumed Julia would raise duplicate-keyword
      error). Folded a `haskey(kwargs, :pilot) && throw(ArgumentError(...))`
      guard into M3 (rather than amending M2) so the test that
      demonstrates the issue ships in the same commit.

### M4 -- Docs and wrap-up

- [x] `AdaptiveLasso` docstring updated to mention the convenience
      constructors and `PrecomputedPilot`.
- [x] `docs/src/api.md` -- add `PrecomputedPilot` to the Estimators
      block.
- [x] `SPEC.md` -- update the `# Estimators` API line and list the
      two convenience constructors.
- [x] `CHANGELOG.md` `[Unreleased] ### Added` entry.
- [x] Iterative-variant design note: point at `PrecomputedPilot` as
      the reusable kernel.
- [x] Update `Status:` in this file and the row in
      `docs/specs/README.md`.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes.
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] ~~If results changed: regression or validation test added.~~
      (no numerical change; pilot-reuse parity test added anyway as a
      regression anchor.)
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [x] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ (no hot-path change.)
- [x] If module names or Makefile targets changed:
      `.claude/agents/` swept and updated. (Sweep confirmed no stale
      AdaptiveLasso / pilot references in `.claude/agents/`.)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [x] Implementation commit hash appended below.

## Commit log

- `aa724df` docs(specs): add AdaptiveLasso precomputed-pilot spec
- `6d54189` feat(Fitting): add PrecomputedPilot estimator adapter  (M1)
- `7aa9074` feat(Magesty): add AdaptiveLasso(::SCEFit) / AdaptiveLasso(::SCEModel) convenience ctors  (M2)
- `c574296` test(Fitting): add PrecomputedPilot + AdaptiveLasso reuse tests  (M3, incl. silent-override guard)
- (M4 docs commit hash appended after this commit lands.)
