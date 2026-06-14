# Requirements: Magesty.jl main-module file split

Status: complete (2026-06-14)

## Goal

Split the oversized main module file `src/Magesty.jl` (~2088 lines) into
smaller `include`d files that splice into the same `Magesty` module
namespace, so the file is readable and navigable without changing any
public API or numerical result.

## Background

A prior whole-package review flagged the main module size as a
maintainability item: `Magesty.jl` mixes the core types, fitting,
prediction, ~24 evaluation metrics, GCV diagnostics, and I/O in one file.
A four-axis design panel (numerical / maintainability / performance / API)
was consulted on the approach and converged on:

- Mechanism: plain `include`d files in the same `Magesty` namespace (the
  idiom already used by the trailing includes `FitCheckIO.jl`,
  `VaspConvert.jl`, `VaspSampling.jl`, `SunnyExport.jl`), **not** Julia
  submodules. Submodules would force `Fitting` internal-export churn and
  dependency-ordering pain for no benefit.
- Scope: extract an `Evaluation.jl` (prediction + metrics) and a `GCV.jl`.
  Keep `save`/`load`, `show` methods, and `_print_fit_summary` in the main
  file.

This is a behavior-preserving, AST-equivalent textual move. The panel
confirmed zero numerical risk, zero public-API impact, and no benchmark
re-baseline needed.

## Scope

Includes:

- Extract `src/Evaluation.jl`: `SCEEvalData`, `_check_basis`,
  `_eval_energy`, `_eval_torque`, `predict_energy`, `predict_torque`, and
  the eight metric verbs `r2_energy/torque`, `rss_energy/torque`,
  `residuals_energy/torque`, `rmse_energy/torque`.
- Extract `src/GCV.jl`: `GCVLambdaPath`, `GCVSizeCurve` (with their `show`
  methods), `gcv`, `gcv_r2`, `gcv_lambda`, `gcv_learning_curve`, and the
  private GCV helpers (`_gcv_core`, `_gcv_subset`, `_default_gcv_sizes`,
  `_argmin_ignore_nan`, `_check_gcv_torque_weight`).
- Wire the two new files via `include(...)` in `Magesty.jl` at positions
  that respect the existing load order.

Excludes:

- `save`/`load` (thin XMLIO delegators, ~115 lines) — stay in `Magesty.jl`.
- Main-type `show` methods (`SCEBasis`, `SCEDataset`, `SCEFit`, `SCEModel`)
  and `_print_fit_summary` — stay in `Magesty.jl`. (The GCV-type `show`
  methods for `GCVLambdaPath` / `GCVSizeCurve` move with their types into
  `GCV.jl`.)
- Core types, `fit`/`refit`, StatsAPI verbs — stay in `Magesty.jl`.
- Any logic edit, rename, signature change, or docstring rewrite. This is
  a verbatim move only.
- The hot paths (`Fitting.jl`, `SALCBases.jl`, `TesseralHarmonics.jl`) —
  untouched.

## Invariants

- No public name moves out of module `Magesty`; all exports
  (`Magesty.jl:92-109`) and the non-exported `Magesty.save` /
  `Magesty.load` keep identical signatures, docstrings, and behavior.
- No numerical result changes: summation order, signs, units,
  normalization all byte-for-byte identical (verbatim move).
- Existing XML round-trips byte-for-byte (save/load not moved; XMLIO
  untouched).
- SALC key-group / design-matrix column ordering unchanged.
- Spin direction stays unit-vector with `3 × n_atoms` layout (not touched).
- `git diff -w` on `src/` shows only line moves plus the new `include`
  lines (auditable verbatim-move claim).

## Completion criteria

- [ ] `make test-all` passes (unit + integration).
- [ ] `make test-jet` clean (catches any include-order / scope breakage).
- [ ] `make test-aqua` clean.
- [ ] `git diff -w` confirms the move is verbatim (no content deltas beyond
      `include` lines).
- [ ] `Magesty.jl` reduced to ~1250 lines; `Evaluation.jl` and `GCV.jl`
      created.
- [ ] No new test needed (existing suite already covers predict / metrics /
      GCV); none added unless a gap is found.
- [ ] No `.claude/bench_log.md` entry needed (no hot path touched, AST
      equivalent).
- [ ] `docs/specs/README.md` table + this `tasklist.md` `Status:` updated.

## References

- Related specs / design notes: whole-package review backlog
  (`docs/design-notes/investigations/package-review-backlog.md`,
  "Maintainability — main module size").
- Prior layout specs: 260516-src-layout-refactor, 260601-src-refactor.
