# Design: Magesty.jl main-module file split

Status: complete (2026-06-14)

## Summary

Move two cohesive blocks out of `src/Magesty.jl` into new `include`d files
that splice into the same `Magesty` module namespace. This is the existing
trailing-include idiom (`FitCheckIO.jl` etc.) applied to the body of the
module. Because `include` is a textual splice evaluated in the enclosing
module, the result is AST-equivalent: same dispatch, same codegen, same
numbers, same exported bindings. No submodules.

Two files are created:

- `src/Evaluation.jl` — "use the model on data": the prediction verbs and
  the accuracy metrics that compare observed vs predicted. The shared
  helpers `SCEEvalData`, `_check_basis`, `_eval_energy`, `_eval_torque`
  live here too, since they bridge prediction and metrics.
- `src/GCV.jl` — the generalized-cross-validation diagnostics, already the
  most self-contained block (own types, own private helpers, one-directional
  calls into `Fitting`).

Rejected alternatives:

- **Submodules.** The moved code makes 15+ calls into `Fitting` internals
  (`Fitting._calc_r2score`, `Fitting._gcv_*`, etc.) and references the core
  types (`SCEModel`/`SCEFit`/`SCEDataset`). A submodule would require
  exporting those internals from `Fitting` or duplicating `using`/`import`
  per submodule — churn with no benefit.
- **Separate `Prediction.jl` + `Metrics.jl`.** The seam between prediction
  and metrics is artificial: `_eval_*` (the metrics' shared front-end)
  calls `predict_*` directly, so the two are one cohesive unit. Keeping
  them together avoids a cross-file private dependency. (Considered and
  declined in favor of a single `Evaluation.jl`.)
- **Extracting `save`/`load`.** Only ~115 lines of thin XMLIO delegation;
  a separate file saves little reading effort for another file to track.

## Module layout

| Target | Change |
|---|---|
| `src/Magesty.jl` | Remove the evaluation block (current ~1035–1470) and the GCV block (current ~1473–1872); add `include("Evaluation.jl")` and `include("GCV.jl")` at the correct load positions. Keep types, `fit`/`refit`, StatsAPI verbs, `SCEModel(::SCEFit)`, `save`/`load`, `show`, `_print_fit_summary`. ~2088 → ~1250 lines. |
| `src/Evaluation.jl` (new) | `SCEEvalData` const; `predict_energy` / `predict_torque` (all overloads); `_check_basis`; `_eval_energy` / `_eval_torque`; `r2_energy/torque`, `rss_energy/torque`, `residuals_energy/torque`, `rmse_energy/torque`. |
| `src/GCV.jl` (new) | `GCVLambdaPath`, `GCVSizeCurve` + their `Base.show`; `_check_gcv_torque_weight`, `_gcv_core`; `gcv`, `gcv_r2`, `gcv_lambda`, `gcv_learning_curve`; `_gcv_subset`, `_default_gcv_sizes`, `_argmin_ignore_nan`. |

## Load-order constraints

`include` injects top-level expressions in place. The hard constraint is
only on parse-time-evaluated top level (`const`, `struct`, and return-type
annotations on method definitions); method *bodies* resolve at call time
and are order-tolerant.

- `Evaluation.jl` must be included **after** the four core types
  (`SCEBasis`/`SCEModel`/`SCEDataset`/`SCEFit`) and `SpinConfig` are in
  scope (`const SCEEvalData = Union{SCEDataset, AbstractVector{SpinConfig},
  AbstractString}` is parse-time), and **before** the trailing domain
  includes (`FitCheckIO.jl`/`VaspConvert.jl`/`VaspSampling.jl`/
  `SunnyExport.jl`) that consume `predict_*`.
- `GCV.jl` must be included **after** `Fitting` (`using .Fitting`), the
  core types, and `MersenneTwister`/`randperm`/`mean`/`std` are in scope
  (`GCVSizeCurve.estimator::AbstractEstimator` field + the return-type
  annotations are parse-time).
- `_print_fit_summary` stays in `Magesty.jl` and calls `_eval_*`; since it
  is only *called* at runtime (from `fit`/`refit`), it does not need
  `_eval_*` defined before its own definition — but the cleanest placement
  is to `include("Evaluation.jl")` before `_print_fit_summary`'s definition
  anyway, so the dependency is visibly satisfied.

Chosen placement: keep each `include` at the position the moved block
currently occupies in `Magesty.jl` (evaluation block → `include
("Evaluation.jl")`; GCV block → `include("GCV.jl")`). This preserves the
current top-to-bottom order exactly, so no relative ordering changes. Each
new file opens with a one-line comment stating its required-in-scope
dependencies.

## API

No API additions, changes, or deletions. Every public binding keeps its
current signature and docstring; only the source file holding it changes.
Docstrings attach to bindings, not files, so `?fn` and Documenter `@ref`
links are unaffected.

## Types and conventions

No new types, no convention changes, no unit/normalization/sign touch.
`GCVLambdaPath` / `GCVSizeCurve` definitions move verbatim (with their
`show` methods, improving locality).

## Impact on linked sites

- [x] Spherical-harmonics convention (`TesseralHarmonics`): **not touched.**
      Moved code reaches harmonics only transitively via the unmoved
      `Fitting._predict_*`.
- [x] SCE coefficient XML (`save` / `load`): **not touched.** `save`/`load`
      stay in `Magesty.jl`; the format/order contract lives in `XMLIO.jl`
      (unmoved).
- [x] `Fitting` <-> `SALCBasis`: **not touched.** Moved code consumes the
      prebuilt `salc_list` / `X_E` / `X_T` ordering; it does not define it.
- [ ] `.claude/agents/` references: **deferred to M3** — sweep for any
      agent file naming the removed line ranges or asserting "all
      metrics/GCV live in Magesty.jl" (exit-checklist item).
- [ ] `SPEC.md` / `docs/src/api.md`: `api.md` is binding-based (no file
      paths) — no update needed. The `SPEC.md` file-layout table lists only
      submodule/utility files and does not enumerate the trailing includes
      (`FitCheckIO.jl` etc.), so it needs no row for `Evaluation.jl` /
      `GCV.jl` either — **no update expected**. Confirm during M3.

## Test strategy

No new tests. The existing suite already exercises `predict_energy`/
`predict_torque`, the metric verbs, the GCV verbs, and the XML round-trip.
A verbatim include-move cannot change a number; if include order or symbol
scope is wrong, the module fails to load or a test fails — exactly the
signal wanted.

Evidence gate:
1. `make test-all` green.
2. `make test-jet` green (primary detector of scope/include-order breakage).
3. `make test-aqua` green.
4. `git diff -w` audited to confirm verbatim move (only `include` lines as
   real content deltas).

## Risks and open items

- **Include-order mistake** is the only realistic failure mode (a
  parse-time `const`/`struct`/return-annotation referenced before its
  definition). Mitigated by keeping includes at the blocks' current
  positions and by `make test-jet`.
- **`SPEC.md` file-layout drift** — verify whether `SPEC.md` enumerates
  `src/` files; update if so.
- **`.claude/agents/` drift** — sweep for any reference to the old
  single-file layout (exit checklist item).
- No deferred alternatives carry numerical implications.
