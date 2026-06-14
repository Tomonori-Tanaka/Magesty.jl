# Tasklist: Magesty.jl main-module file split

Status: complete (2026-06-14)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Extract `Evaluation.jl`

- [x] Create `src/Evaluation.jl` with a leading comment stating its
      required-in-scope dependencies (core types, `SpinConfig`, `Fitting`).
- [x] Move verbatim: `SCEEvalData`, `_check_basis`, `_eval_energy`,
      `_eval_torque`, `predict_energy` / `predict_torque` (all overloads),
      `r2_*`, `rss_*`, `residuals_*`, `rmse_*`.
- [x] Replace the removed block in `Magesty.jl` with `include("Evaluation.jl")`.
- [x] `make test-all` + `make test-jet` green.

### M2 — Extract `GCV.jl`

(Independent of M1 — the GCV block makes no calls into the Evaluation
block, so order is a matter of preference.)

- [x] Create `src/GCV.jl` with a leading dependency comment.
- [x] Move verbatim: `GCVLambdaPath`, `GCVSizeCurve` (+ `show`),
      `_check_gcv_torque_weight`, `_gcv_core`, `gcv`, `gcv_r2`,
      `gcv_lambda`, `gcv_learning_curve`, `_gcv_subset`,
      `_default_gcv_sizes`, `_argmin_ignore_nan`.
- [x] Replace the removed block in `Magesty.jl` with `include("GCV.jl")`.
- [x] `make test-all` + `make test-jet` green.

### M3 — Verify + housekeeping

- [x] `git diff -w` audited: `Magesty.jl` shows 845 deletions + the 11-line
      include snippet; the bodies were `sed`-sliced verbatim into the new
      files (no content edits).
- [x] `make test-aqua` green (10/10).
- [x] Check `SPEC.md` for an `src/` file-layout list: the module table lists
      submodules only (`FitCheckIO.jl` and the other include fragments are
      not listed), so `Evaluation.jl` / `GCV.jl` need no row. No update.
- [x] Sweep `.claude/agents/` for old single-file-layout references: none
      found (agent files are generic, no line/file-layout citations).
- [x] Update `docs/design-notes/investigations/package-review-backlog.md`
      (marked "main module size" done).

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes (23677 passed).
- [x] `make test-aqua` / `make test-jet` clean (aqua 10/10, jet 1 passed).
- [x] ~~If results changed: regression or validation test added.~~ Verbatim
      move — no result change, no new test (confirmed via `git diff -w`).
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
      (API unchanged; `SPEC.md` module table is submodule-only — no update;
      `api.md` is binding-based — no update.)
- [x] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ No hot path touched; AST equivalent.
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved. Design-phase panel consulted (all
      four axes); a Tier 1 `code-reviewer` pass on the landed verbatim diff
      stands in for a full rerun, since no logic changed.
- [x] If module names or Makefile targets changed:
      `.claude/agents/` swept and updated. (No module names change; swept
      for stale file-layout references — none found.)
- [x] `CHANGELOG.md` `[Unreleased]` updated (`### Internal` entry added).
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
