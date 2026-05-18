# Tasklist: ElasticNet estimator (with Lasso convenience)

Status: **complete (2026-05-18)** -- branch `refactor/lasso-estimator`,
commits `1634fe8..ab3eee0`.

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Plumb GLMNet into the package

- [ ] Add `GLMNet` to `Project.toml` `[deps]` with `[compat]`
      starting at `"0.7"`; adjust after the API check below.
- [ ] Verify the GLMNet API actually used by the design
      (`glmnet(X, y; alpha, lambda=[…], standardize, intercept)`,
      `path.betas`). At a REPL with `using GLMNet`,
      check `methods(glmnet)` and `?glmnet`; record any keyword
      naming or return-shape surprises and adjust `_glmnet_solve`
      in M2 accordingly.
- [ ] `using GLMNet` in `src/Fitting.jl`; confirm `Magesty`
      precompiles.
- [ ] Confirm `make test-aqua` is still clean (Aqua does not flag
      stray bindings).

### M2 — Add the `ElasticNet` struct and dispatch method

- [ ] Add `struct ElasticNet`, keyword constructor with
      `0 ≤ alpha ≤ 1` and `lambda ≥ 0` validation, and `Lasso(...)`
      convenience constructor.
- [ ] Add private helper `_glmnet_solve`.
- [ ] Implement `solve_coefficients(::ElasticNet, X, y)` per the
      design doc — a thin wrapper around `_glmnet_solve` with
      `intercept=false`. No in-method centering, no column
      stripping (both are already handled by
      `assemble_weighted_problem` since spec 260518).
- [ ] Export `ElasticNet` and `Lasso` from `Fitting`; re-export
      from `Magesty`.
- [ ] `make test-unit` green: existing `OLS` / `Ridge` paths are
      untouched, so they should remain identical numerically.

### M3 — Component tests

- [ ] Tests 1–6 from `design.md` "Test strategy" added (or
      extended) in `test/component/test_fitting_estimators.jl`.
      Tests 1 and 2 (energy-only and mixed-weight OLS limit) are
      the load-bearing regression tests for the
      `intercept=false` + upstream-centering wiring; the rest
      cover behaviour.
- [ ] Decide and record the GLMNet-`α=0` vs analytic-Ridge
      tolerance (`rtol`) in the test docstring.
- [ ] `make test-unit` green.

### M4 — Integration smoke and docs

- [ ] Add a `Lasso(lambda=…)` smoke path to one existing
      `test/integration/*` test or alongside it.
- [ ] Update `docs/src/api.md` to list `ElasticNet` and `Lasso`.
- [ ] Append a `[Unreleased]` entry in `CHANGELOG.md` covering the
      new estimator, the GLMNet dependency, and the Ridge-agreement
      tolerance number.
- [ ] `make test-integration` green.

### M5 — Wrap-up

- [ ] `make test-jet` clean.
- [ ] Update `Status:` line in this file → `complete (YYYY-MM-DD)`
      and the matching row in `docs/specs/README.md` in the same
      commit.
- [ ] Mark the on-hold Adaptive sketch in the related design note
      as the next-up follow-up (no other content change).

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] GLMNet-`α=0` vs analytic-Ridge agreement test passes; the
      tolerance value is recorded in the test docstring and in
      `CHANGELOG.md` `[Unreleased]`.
- [ ] Mixed-weight (`weight = 0.5`) OLS-limit test for `ElasticNet`
      passes within the same `atol` as the energy-only OLS-limit
      test — confirms the upstream energy-only centering plus
      `intercept=false` is correctly wired.
- [ ] If results changed: regression or validation test added.
      `ElasticNet` is additive — existing `OLS` / `Ridge` paths are
      untouched, so no regression test for them; the new tests
      cover `ElasticNet` directly.
- [ ] Public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ Not a hot path. If
      `make bench-fitting` exists, a courtesy entry is appended.
- [ ] If module names or Makefile targets changed:
      `.claude/agents/` swept and updated. (`code-reviewer` and
      `test-runner` agent prompts: confirm no estimator-name list
      needs the `ElasticNet` / `Lasso` addition.)
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.

## Implementation commits

- `1634fe8`  chore(deps): add GLMNet for ElasticNet estimator
- `391934f`  feat(Fitting): add ElasticNet estimator backed by GLMNet
- `a094b32`  test(Fitting): add ElasticNet component tests
- `45b97a9`  feat(Magesty): show nonzero coef count in fit summary
- `ab3eee0`  docs: document ElasticNet/Lasso estimators and nonzero coef line
