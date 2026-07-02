# Tasklist: SALC projection time-reversal dedup

Status: draft (2026-07-02)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Baseline capture

- [ ] Capture pre-change snapshots on the development machine
      (scratch, not committed): `projection_mat` outputs and full SALC
      coefficients for the dimer and `fege_2x2x2` fixtures.
- [ ] Land the committed regression probe in
      `test/component/test_SALC_projection.jl` (projection matrix vs
      committed baseline, `isapprox(atol = 1e-8)`;
      `salc_fingerprint` with `==`) against the *current*
      implementation — test must pass before any loop change.
- [ ] Record the "before" benchmark
      (`bench/benchmark_salcbasis_hotspots.jl`) numbers.

### M2 — Loop restructure

- [ ] Implement the single-pass dual-buffer loop in
      `_projection_matrix_coupled_basis` per `design.md`.
- [ ] One-time bit-identity gate: re-capture on the same machine and
      compare `==` against the M1 snapshots. Any difference: stop and
      consult.
- [ ] Regression probe green; `make test-all` / `test-jet` /
      `test-aqua` green.

### M3 — Bench + review

- [ ] Record the "after" benchmark; append the before/after entry
      (including the bit-identity gate outcome) to
      `.claude/bench_log.md`.
- [ ] Numerical-reviewer sign-off on the accumulation-order argument.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] If results changed: regression or validation test added.
      (Results must NOT change — the bit-identity test is the guard.)
- [ ] ~~If public API changed: `SPEC.md` and `docs/src/api.md`
      updated.~~ (No API change.)
- [ ] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.
- [ ] Tier 2 review panel run (numerical / maintainability /
      performance / API axes) and findings resolved.
- [ ] ~~If module names or Makefile targets changed: `.claude/agents/`
      swept and updated.~~ (None changed.)
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
