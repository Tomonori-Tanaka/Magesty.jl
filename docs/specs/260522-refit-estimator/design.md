# Design: `refit` — post-selection refit on the selected basis support

Status: draft (2026-05-22)

## Summary

`refit` is a post-hoc function on a fitted `SCEFit`. It reads the input
fit's coefficient vector, decides which bases form the *support* (the
non-negligible set), re-assembles the same weighted least-squares problem
the original `fit` used, and re-solves it restricted to the support
columns with a caller-chosen estimator. The result is a fresh `SCEFit`
whose `jphi` is full length: support columns carry the re-estimated
coefficients, dropped columns stay at `0.0`.

The selection criterion is the **scaled coefficient magnitude**

    keep basis j  iff  abs(jphi[j]) * norm(X[:, j]) > threshold

where `X` is the assembled weighted design matrix. `abs(jphi[j]) *
norm(X[:, j])` is the magnitude of basis `j`'s contribution to the
weighted objective the fit minimizes — the column scale `norm(X[:, j])`
neutralizes the per-cluster `(4π)^(N/2)` factor that otherwise makes raw
`jphi` magnitudes incomparable across cluster sizes.

A single criterion (no mode flag) covers both selector families:

- **L1 selectors** (`Lasso`, `AdaptiveLasso`): coefficients are exactly
  `0.0`, so `abs(jphi[j]) * norm(X[:, j]) = 0`. The default
  `threshold = 0.0` keeps exactly the non-zero support — no tuning
  needed.
- **`AdaptiveRidge`**: coefficients shrink but never reach exact zero, so
  `threshold = 0.0` keeps everything (`refit` becomes a no-op). The
  caller passes an explicit positive `threshold` to drop the
  near-negligible bases.

Alternatives considered and rejected:

- A `Refit <: AbstractEstimator` wrapper with a one-shot
  `fit(SCEFit, dataset, Refit(selector, refit_est))` path. Rejected: the
  use case is post-hoc inspection of an existing fit; a wrapper would
  re-run the (expensive) selection. The post-hoc function composes the
  same way (`fit` then `refit`) without that cost.
- A public `OnSupport(estimator, columns)` column-restriction estimator
  plus a `support(fit; ...)` helper. Rejected for the first version to
  keep the API surface minimal; `refit` is the only new public name.
- A `relative` / `:absolute` mode flag. Rejected: the scaled criterion
  alone serves both families, and `threshold = 0.0` makes it invisible
  to L1 users.

## Module layout

| Target | Change |
|---|---|
| `src/Magesty.jl` | New `refit` function next to `fit(::Type{SCEFit}, ...)`; add `refit` to the `export` list. |
| `src/Fitting.jl` | No new kernel functions — `refit` reuses `assemble_weighted_problem`, `solve_coefficients`, `extract_j0_jphi`. A small internal helper for column-norm support selection may live here. |
| `docs/src/api.md` | Add `refit` to the fitting `@docs` block. |
| `docs/src/theory/design_matrix_and_fitting.md` | Describe post-selection refit and the scaled criterion. |
| `CHANGELOG.md` | `[Unreleased]` entry. |
| `test/component/test_refit.jl` | New component test file. |

## API

```julia
"""
    refit(fit::SCEFit, estimator::AbstractEstimator = OLS();
          threshold::Real = 0.0, verbosity::Bool = true) -> SCEFit

Post-selection refit on the basis support of `fit`.

# Arguments
- `fit::SCEFit`: The fitted model whose coefficient support is reused.
- `estimator::AbstractEstimator = OLS()`: Estimator for the support
  re-solve. Must adapt to the column count it is given (see the
  `PrecomputedPilot` restriction below).
- `threshold::Real = 0.0`: Scaled-magnitude cutoff; bases with
  `abs(jphi[j]) * norm(X[:,j]) <= threshold` are dropped.
- `verbosity::Bool = true`: Print the fit summary, as `fit` does.

# Returns
- `SCEFit`: A fresh fit; dropped bases carry coefficient `0.0`.

# Examples
(L1 selection + OLS debiasing; AdaptiveRidge + explicit threshold.)
"""
function refit(
    fit::SCEFit,
    estimator::AbstractEstimator = OLS();
    threshold::Real = 0.0,
    verbosity::Bool = true,
)::SCEFit
```

Algorithm:

1. Validate: `threshold >= 0` else `ArgumentError`; reject a
   `PrecomputedPilot`-backed `estimator` (see edge cases) with
   `ArgumentError`.
2. `X, y = assemble_weighted_problem(fit.dataset.X_E, fit.dataset.X_T,
   fit.dataset.y_E, fit.dataset.y_T, fit.torque_weight)` — the same call
   `fit` makes, reusing `fit.dataset` and `fit.torque_weight`.
3. `scale[j] = norm(@view X[:, j])` for every column.
4. `support = findall(j -> abs(coef(fit)[j]) * scale[j] > threshold,
   eachindex(coef(fit)))`.
5. **If `isempty(support)`**: short-circuit — emit `@warn`, set
   `j_values = zeros(num_salcs)`, skip the solve. Otherwise
   `j_sub = solve_coefficients(estimator, X[:, support], y)` and scatter
   into `j_values = zeros(num_salcs); j_values[support] = j_sub`.
6. `j0, jphi = extract_j0_jphi(j_values, fit.dataset.X_E,
   fit.dataset.y_E)`.
7. `residuals = y .- X * j_values`.
8. Return `SCEFit(fit.dataset, j0, jphi, estimator, fit.torque_weight,
   residuals)`.

`estimator` defaults to `OLS()` — classic post-selection debiasing. Pass
`Ridge(lambda = small)` when the selected support is still
rank-deficient or near-collinear.

Edge cases / errors:

- **Empty support** (every basis below `threshold`): step 5 short-circuits
  *before* `solve_coefficients`, so no estimator ever sees a zero-column
  matrix (a GLMNet-backed estimator would error on one). `@warn`, return
  an all-zero `jphi`; `extract_j0_jphi` then yields
  `j0 = mean(fit.dataset.y_E)`. Not an error — a valid, if degenerate,
  fit.
- **`PrecomputedPilot`-backed `estimator` is rejected.** `refit` passes
  the sub-matrix `X[:, support]` (only `|support|` columns) to
  `solve_coefficients`. `PrecomputedPilot` stores a fixed vector of
  length `num_salcs` and length-checks against the column count, so it —
  and any `AdaptiveLasso` whose `pilot` is a `PrecomputedPilot` — would
  throw `DimensionMismatch` deep in the solve. `refit` detects this
  upfront and raises a clear `ArgumentError` instead. (A precomputed
  pilot is meaningless for a refit anyway: the support has already been
  chosen.)
- The `threshold >= 0` precondition is checked with `ArgumentError`.
- `verbosity` mirrors `fit`: prints the standard fit summary via
  `_print_fit_summary`.

The returned `SCEFit.estimator` is the *refit* estimator; `refit`
intentionally does not record which selection produced the support
(documented in the docstring — provenance lives in the caller's code).

## Types and conventions

- No physics-convention change. `jphi` length and SALC key-group order
  are untouched; dropped bases are `0.0`, exactly as a truncated basis
  already appears.
- `j0` recovery is unchanged (`extract_j0_jphi`), so `j0` stays in the
  input energy unit and independent of `torque_weight`.
- New invariant: `refit` reuses the input fit's `torque_weight`, so the
  refit minimizes the *same* weighted objective as the original `fit` —
  the support selection and the refit are mutually consistent.
- Selection criterion `abs(jphi[j]) * norm(X[:, j])` uses the assembled
  (centered, weighted) `X`. An all-zero column (a basis with no signal)
  has `norm = 0` and is dropped regardless of `jphi[j]` — correct, since
  it contributes nothing.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): none.
- [ ] SCE coefficient XML (`save` / `load`): none — a refit result is an
      ordinary `SCEFit`; `SCEModel(refit_result)` saves unchanged.
- [x] `Fitting` <-> `SALCBasis`: `refit` must preserve `jphi` length and
      key-group order. It only zeroes a subset and re-solves the rest;
      no reordering. Verified by test.
- [ ] `.claude/agents/` references: none (no module/Makefile changes).
- [x] `SPEC.md` / `docs/src/api.md` updates: `api.md` `@docs` block plus
      the fitting theory page.

## Test strategy

New `test/component/test_refit.jl`:

- **L1 exact-zero support**: fit with `Lasso`, `refit` with `OLS` at the
  default `threshold = 0.0`. Assert the refit support equals
  `findall(!iszero, coef(lasso_fit))`, and the refit coefficients equal
  an independent `OLS` solve restricted to those columns.
- **`AdaptiveRidge` thresholded**: fit with `AdaptiveRidge`, `refit` with
  a positive `threshold`; assert the support shrinks and dropped
  coefficients are exactly `0.0`.
- **OLS idempotence**: `refit` of an `OLS` fit at `threshold = 0` with
  `OLS` reproduces `coef` / `intercept` to round-off (the support is the
  full set when no coefficient is exactly zero).
- **Empty support**: a `threshold` above every scaled magnitude warns and
  returns an all-zero `jphi` with `j0 == mean(y_E)`.
- **`torque_weight` consistency**: `refit` of a fit made with
  `torque_weight = 0.5` reuses `0.5`.
- Expected values are derived analytically (independent `\` solves on
  the column subset), not read back from the implementation.

No hot-path change; a `bench_log.md` entry is not required, but the
refit cost (one extra `assemble_weighted_problem` + solve) is noted in
the docstring.

## Risks and open items

- The result's `estimator` field loses selection provenance. Accepted;
  documented. A future `Refit`-record wrapper could carry it if needed.
- `AdaptiveRidge` threshold choice is left to the caller; there is no
  automatic scale. The docstring gives guidance (start from the scaled
  magnitude distribution of `coef(fit)`).
- Deferred: `OnSupport` estimator + `support` helper for a reusable
  column-restriction path and a one-shot selection+refit. Recorded here
  in case a later spec wants it.
