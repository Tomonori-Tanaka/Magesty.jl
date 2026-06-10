# Design: GCV diagnostics (hat-matrix generalized cross-validation)

Status: draft (2026-06-10)

## Summary

GCV is computed on the **weighted, energy-mean-centered augmented system**
`(X, y)` that the fit already minimizes — `X = vcat(scale_e·centered(X_E),
scale_m·X_T)`, `y = vcat(scale_e·centered(y_E), scale_m·y_T)` — produced by
`Fitting.assemble_weighted_problem`. This makes the combined energy+torque GCV
the exactly-defined quantity for linear estimators, with no change to the fit.

```
GCV = (‖r‖² / N) / (1 − tr(H)/N)²
```

`r` is the augmented weighted residual (already stored on `SCEFit.residuals`,
or recomputed for a subset). `tr(H)` is the effective degrees of freedom. `N` is
the number of **live** rows: a block whose whitening scale is zero
(`scale_e = 0` at `torque_weight = 1`; `scale_m = 0` at `torque_weight = 0`)
contributes only dead all-zero rows that add nothing to `r` or `tr(H)`, so it
must not inflate `N`. Concretely `N = n_E + n_T` for `0 < w < 1`, `N = n_T` for
`w = 1`, `N = n_E` for `w = 0`. The eliminated `j0` costs one degree of freedom
only when the energy block is live (`w < 1`), i.e. the `+1` is conditional.
`_gcv_sample_count(n_E, n_T, w)` returns `(N, intercept_dof)`.

Two sweep axes share the single-score primitive:

- **`gcv_lambda`** (penalty selection): one economy SVD `X = U·diag(σ)·Vᵀ`
  yields every `lambda`:
  - `dof(λ) = intercept_dof + Σ_k σ_k²/(σ_k²+λ)`
  - `RSS(λ) = ‖y_⊥‖² + Σ_k (λ/(σ_k²+λ))² aₖ²`, `aₖ = (Uᵀy)_k`,
    `‖y_⊥‖² = ‖y‖² − Σ aₖ²` (dead rows carry zeros in `y` and `U`, so they drop
    out of `‖y‖²` and `a`)
  - `GCV(λ) = (RSS(λ)/N)/(1 − dof(λ)/N)²`; `lambda_best = argmin`.
- **`gcv_learning_curve`** (data sufficiency): for each size `n` draw `repeats` random
  config subsets (seeded RNG), build the sub-dataset via existing
  `dataset[idx]`, `fit` it, take `gcv(f)`; report mean ± std.

Rejected alternatives: putting the logic in an external script (it must reach
into `assemble_weighted_problem` anyway, so core is cleaner); k-fold CV
(heavier, not requested); per-block (energy-only / torque-only) GCV (the fit
minimizes the combined objective, so the combined GCV is the consistent target).

## Module layout

| Target | Change |
|---|---|
| `src/Fitting.jl` | Add internal numerics: `_is_linear_estimator`, `_gcv_from_svd`, `_gcv_lambda_path(X, y, lambdas)`, `_effective_dof_adaptive_ridge`. They consume the assembled `(X, y)`. |
| `src/Magesty.jl` | Add user-facing `gcv` / `gcv_lambda` / `gcv_learning_curve`, structs `GCVLambdaPath` / `GCVSizeCurve`, and `export` lines (after the `rmse_*` block). `gcv_learning_curve` uses `dataset[idx]` + `fit` + `gcv`. |
| `src/FitCheckIO.jl` | Add `write_gcv_lambda` / `write_gcv_learning_curve` mirroring `_write_energy_file`. |
| `tools/FitCheck_gcv_lambda.py`, `tools/FitCheck_gcv_learning_curve.py` | New plot scripts (GCV-vs-λ log-x; GCV-vs-n with std error bars). |
| `examples/` | One end-to-end example script. |
| `docs/src/` | Narrative diagnostics page + `api.md` `@docs` entries. |

## API

```julia
# Single combined GCV score for an existing fit (linear estimators only).
gcv(f::SCEFit)::Float64

# Ridge penalty path — one SVD, all lambdas.
struct GCVLambdaPath
    lambdas::Vector{Float64}
    gcv_scores::Vector{Float64}    # named to avoid shadowing the `gcv` function
    dof::Vector{Float64}
    lambda_best::Float64
    torque_weight::Float64
end
gcv_lambda(dataset::SCEDataset, lambdas::AbstractVector{<:Real};
           torque_weight::Real = 1.0)::GCVLambdaPath

# Data-sufficiency sweep — random subsets, mean ± std.
struct GCVSizeCurve
    sizes::Vector{Int}
    gcv_mean::Vector{Float64}
    gcv_std::Vector{Float64}
    repeats::Int
    seed::Int
    estimator::AbstractEstimator
    torque_weight::Float64
end
gcv_learning_curve(dataset::SCEDataset, estimator::AbstractEstimator = OLS();
             sizes::AbstractVector{<:Integer} = <auto grid>,
             repeats::Integer = 5, seed::Integer = 0,
             torque_weight::Real = 1.0)::GCVSizeCurve

# Table writers (FitCheck convention).
write_gcv_lambda(path::GCVLambdaPath, filename::AbstractString)::Nothing
write_gcv_learning_curve(curve::GCVSizeCurve, filename::AbstractString)::Nothing
```

All public functions get standard-format docstrings
(`# Arguments` / `# Returns` / `# Examples`) and explicit type annotations.

## Types and conventions

- GCV value is in the **weighted-objective unit** (dimensionless mixture of eV²
  energy and eV² torque scaled by `(1−w)/n_E` and `w/n_T`), **not eV²**. This is
  documented; the diagnostic is for comparing shapes/plateaus, not absolute
  physical error. The docstrings state this explicitly.
- Estimator support: `OLS` / `Ridge` exact; `AdaptiveRidge` **conditional**
  (treats its converged per-column reweighting `D` as fixed in the hat matrix:
  `dof = intercept_dof + tr(X(XᵀX+λD)⁻¹Xᵀ)`). `ElasticNet` (incl. the `Lasso`
  alias) and `AdaptiveLasso` → `ArgumentError` naming the unsupported type.
- `torque_weight` is validated to `[0, 1]` in `gcv_lambda` /
  `gcv_learning_curve` (matching `fit`'s documented convention), with the error
  attributed to the public entry point.
- **`_is_linear_estimator` dispatch (avoid the `Lasso` trap).** `Lasso` is *not*
  a struct — `Lasso(lambda=...)` constructs an `ElasticNet(alpha=1.0, …)`
  instance (see `src/Fitting.jl`). The predicate must branch on the concrete
  types `OLS` / `Ridge` / `AdaptiveRidge` returning `true` and otherwise
  `false`, so that `ElasticNet` (and therefore any `Lasso`-built instance) and
  `AdaptiveLasso` are rejected. Never branch on a non-existent `Lasso` type.
- **`AdaptiveRidge` weight recovery (no `SCEFit` change).** The converged
  weights are `D = Diagonal(w)`, `w[j] = 1/(coef(f)[j]^2 + f.estimator.epsilon)`
  — `coef(f)` *is* the converged `beta` of the same weighted, centered system
  the iteration solved, so this recovers `D` **exactly** at convergence (it is
  not an approximation of the weights; only the GCV treatment of `D` as fixed is
  the "conditional" part). `lambda` and `epsilon` come from `f.estimator`. No
  field is added to `SCEFit` and the solver return value is unchanged.
- `gcv_lambda` is a scalar-L2 (ridge) path by construction.
- Default `sizes` grid: 6 integer points (`unique`, ascending) linearly spaced
  from `max(p + 2, 10)` to `length(dataset)` inclusive (rounded to nearest
  integer), where `p = number of SALCs`. The lower bound keeps `N − dof > 0`
  (torque rows add `3·n_atoms` per config, so even the smallest size is
  comfortably overdetermined). If `length(dataset) < max(p + 2, 10)` the grid
  collapses to the available range and a warning is emitted.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): none.
- [ ] SCE coefficient XML (`save` / `load`): none (no format change).
- [x] `Fitting` <-> `SALCBasis`: read-only reuse of `assemble_weighted_problem`
      and `SCEDataset` slicing; no reordering, no behavior change.
- [ ] `.claude/agents/` references: none (no module/Makefile-target rename).
- [x] `SPEC.md` / `docs/src/api.md` updates: add the new public API.

## Test strategy

`test/component/test_gcv.jl` (analytic / cross-check; never tailored to output):

- OLS GCV on a tiny synthetic design equals `(RSS/N)/(1−p/N)²` by hand.
- SVD-based ridge dof equals dense `1 + tr(X(XᵀX+λI)⁻¹Xᵀ)` (several λ, ~1e-10).
- λ-path closed-form `RSS(λ)` equals re-solving ridge + computing residuals.
- `gcv_lambda` at a λ equals `gcv(fit(..., Ridge(lambda=λ)))` on full data.
- `gcv_learning_curve` reproducible under fixed `seed`; mean/std correct over repeats.
- Non-linear estimators raise `ArgumentError`.
- Writer round-trip: file columns parse back to the struct fields.

## Risks and open items

- No numerical-result change is intended; the only "new numbers" are the GCV
  diagnostics themselves, validated against closed-form / dense recomputation.
- `gcv_learning_curve` at very small `n` can hit a rank-deficient OLS solve or a
  non-positive `N − dof` denominator; the implementation records `NaN` + warns
  rather than emitting a bogus value, and the default grid avoids that region.
- `AdaptiveRidge` GCV is conditional (an approximation); documented as such.
