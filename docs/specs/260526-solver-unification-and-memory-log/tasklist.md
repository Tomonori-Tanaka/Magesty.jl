# Tasklist: Solver unification (Cholesky) and design-matrix memory logging

Status: complete (2026-05-26)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Solver switch (OLS / Ridge / AdaptiveRidge)

- [ ] Replace `solve_coefficients(::OLS, X, y)` body with
      `cholesky(Symmetric(X'X)) \ (X'y)` plus `PosDefException` rethrow
      as `ArgumentError` pointing to `Ridge`.
- [ ] Replace `solve_coefficients(::Ridge, X, y)` body with explicit
      `cholesky(Symmetric(X'X + λI))`; `λ ≈ 0` path delegates to the OLS
      Cholesky path. Drop the `MultivariateStats.ridge` call.
- [ ] Replace the two `Symmetric \` calls in
      `solve_coefficients(::AdaptiveRidge, ...)` with `cholesky`. Confirm
      no `PosDefException` catch is needed (both call sites are SPD by
      construction: `λ > 0` is checked upfront, `w > 0` always). Route
      `λ ≈ 0` through the OLS Cholesky path for consistency.
- [ ] Drop `using MultivariateStats` from `Fitting.jl` and remove
      `MultivariateStats` from `Project.toml` `[deps]` and `[compat]`.
- [ ] Refresh OLS / Ridge docstrings to describe Cholesky and the new
      `ArgumentError` path.

### M2 — Memory logging

- [ ] Add `_format_bytes(::Integer) -> String` private helper.
- [ ] Add pre-construction + post-construction memory log to
      `build_design_matrix_energy` under `verbosity = true`.
- [ ] Same for `build_design_matrix_torque`.

### M3 — Tests

- [ ] Cholesky-OLS vs prior-QR-OLS cross-check on a well-conditioned
      `(X, y)`.
- [ ] `PosDefException` → `ArgumentError` rethrow test on a duplicated
      column.
- [ ] Cholesky-Ridge vs `MultivariateStats.ridge` cross-check.
- [ ] Smoke test confirming the memory log line is emitted under
      `verbosity = true`.

### M4 — Docs and bench

- [ ] Refresh the OLS / Ridge paragraphs in `SPEC.md` and
      `docs/src/api.md`.
- [ ] Add a `[Unreleased]` entry in `CHANGELOG.md` noting the OLS
      behavioral change.
- [ ] Add a before/after `solve_coefficients` timing to
      `.claude/bench_log.md`.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] If results changed: regression or validation test added.
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
      (API surface unchanged but docstrings updated; refresh the
      narrative paragraphs too.)
- [x] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved. The numerical axis is primary here:
      condition-number squaring and `PosDefException` semantics.
- [ ] If module names or Makefile targets changed:
      `.claude/agents/` swept and updated. (Not applicable.)
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
