# Tasklist: SALC projection time-reversal dedup

Status: implemented (2026-07-02)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Baseline capture

- [x] Capture pre-change snapshots on the development machine
      (scratch, not committed): `projection_mat` outputs and full SALC
      coefficients for the dimer and `fege_2x2x2` fixtures.
- [x] Land the committed regression probe in
      `test/component/test_SALC_projection.jl` (projection matrix vs
      committed baseline, `isapprox(atol = 1e-8)`) against the
      *current* implementation — passed before the loop change.
      Baselines: one fege group per distinct dimension (4, 12, 20, 36,
      60); the dimer fixture has a single trivial 1×1 group and is
      covered by the snapshot gate instead. `salc_fingerprint`
      equality is exercised by the snapshot gate and the existing
      fresh-build baseline test in `test_save_load.jl`.
- [x] Record the "before" benchmark
      (`bench/benchmark_salcbasis_hotspots.jl`) numbers.

### M2 — Loop restructure

- [x] Implement the single-pass dual-buffer loop in
      `_projection_matrix_coupled_basis` per `design.md`.
- [x] One-time bit-identity gate: re-captured on the same machine and
      compared `==` against the M1 snapshots — PASS (54 fege
      projection matrices, 81,192 SALC coefficient values,
      fingerprints; dimer + fege).
- [x] Regression probe green; `make test-all` (23,752) / `test-jet` /
      `test-aqua` green.

### M3 — Bench + review

- [x] Record the "after" benchmark; append the before/after entry
      (including the bit-identity gate outcome) to
      `.claude/bench_log.md`. Projection stage 1.95× / −47% allocs;
      constructor 1.19× / −29% allocs (adjacent-session A/B).
- [x] Numerical-reviewer sign-off on the accumulation-order argument —
      clean (all four invariants verified against the diff; the
      committed probe's `isapprox(atol = 1e-8)` tolerance confirmed
      numerically sound with ~5 orders of margin over cross-platform
      trig round-off).

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes (23,752).
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] If results changed: regression or validation test added.
      (Results did NOT change — bit-identity gate PASS; the committed
      projection-value probe guards against future drift.)
- [x] ~~If public API changed: `SPEC.md` and `docs/src/api.md`
      updated.~~ (No API change.)
- [x] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.
- [x] Tier 2 review panel run (numerical / maintainability /
      performance / API axes) and findings resolved. Numerical: clean
      sign-off (all four invariants verified). Maintainability: 1
      major (baseline regeneration recipe) + 3 minor (buffer naming,
      haskey guard, per-group testsets) — all applied. Performance: 2
      minor, both documented no-action (deferred `reorder_atoms`
      follow-up; second scratch buffer acceptable). API: 2 minor
      (CHANGELOG heading moved to `### Internal`; `sum(ls)`
      permutation-invariance confirmed by the numerical axis).
- [x] ~~If module names or Makefile targets changed: `.claude/agents/`
      swept and updated.~~ (None changed.)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [x] Implementation commit hash appended below.

Implementation commit: b5cd2bc
