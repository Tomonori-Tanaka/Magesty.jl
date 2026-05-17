# LASSO / Adaptive LASSO / Adaptive Ridge estimators

**Status**: not started (2026-05-16)

Design note for adopting GLMNet.jl and adding the L1 / Adaptive families
to the estimator-dispatch hierarchy in `Fitting.jl`. Pre-spec draft.

## Background

The SCE coefficient estimator in `src/Fitting.jl` currently supports
only `OLS` (`X \ y`) and `Ridge` (`MultivariateStats.ridge` with a
per-column penalty that excludes the bias column). To obtain sparse
effective spin models we want L1-style regularization (LASSO, Adaptive
LASSO), and as an L0 approximation, Adaptive Ridge (the iterative
reweighting scheme of Frommlet & Nuel 2016).

## Library selection

`MultivariateStats.jl` does not provide LASSO. Comparison of candidates:

| Aspect | Lasso.jl | GLMNet.jl |
|---|---|---|
| Implementation | Pure Julia (coordinate descent) | Wrapper around the Friedman/Hastie/Tibshirani Fortran `glmnet` (reference implementation) |
| Reliability | JuliaStats-official, actively maintained | Bit-identical to R `glmnet` and scikit-learn; 30-year track record |
| Julia dependencies | Pulls GLM, StatsModels, StatsBase, Distributions, etc. | Distributions, SparseArrays, plus `GLMNet_jll` |
| JLL dependency | None | Yes (we already use Spglib_jll, so this is no new burden) |
| `penalty_factor` (bias exclusion / Adaptive weights) | Supported | Supported |
| Unified `alpha` for LASSO / Ridge / Elastic Net | Partial | Complete |

**Choice: GLMNet.jl.** Reasons: (a) numerical reproducibility as the
reference implementation (aligns with CLAUDE.md "Numerical correctness
and reproducibility first"), (b) minimal Julia dependencies,
(c) structurally simple: `alpha` plus `penalty_factor` covers LASSO,
Adaptive LASSO, Adaptive Ridge, and Elastic Net with one implementation.

The existing `Ridge` keeps using MultivariateStats (regression-risk
avoidance). Only new estimators go through GLMNet.

## Conventions (match the existing Ridge)

| Item | Policy |
|---|---|
| Bias column | `penalty_factor[bias_col] = 0` (excluded from penalty) |
| `j0` extraction | Reuse the existing `extract_j0_jphi` |
| Standardize | OFF by default (SCE scale factors give columns physically meaningful units already); enable via `standardize::Bool` kwarg |
| `λ` selection | User-specified single `λ` (matches Ridge). Automatic CV is a separate spec |
| Augmented system | Use the existing `assemble_weighted_problem` output `(X, y, bias_col)` as-is |

## Estimator mapping

| Estimator | `alpha` | `penalty_factor[j != bias]` | Note |
|---|---|---|---|
| Lasso | 1.0 | 1 | — |
| AdaptiveLasso | 1.0 | `1/|β̂_init,j|^γ` | `β̂_init` from OLS or Ridge (Zou 2006) |
| AdaptiveRidge (`mode = :oneshot`) | 0.0 | `1/|β̂_init,j|^γ` | Single-stage reweighting |
| AdaptiveRidge (`mode = :iterative`) | 0.0 | `1/(β_j² + ε)` updated iteratively | Frommlet & Nuel 2016. Refits GLMNet each iteration as an L0 approximation |

We use behavior-descriptive symbols `:oneshot` / `:iterative` rather
than author-name symbols, since the latter assume the reader knows the
paper.

## Draft API

```julia
struct Lasso <: AbstractEstimator
    lambda::Float64
    standardize::Bool
end

struct AdaptiveLasso <: AbstractEstimator
    lambda::Float64
    gamma::Float64           # typically 1.0
    init::Symbol             # :ols | :ridge
    init_lambda::Float64     # only used when init = :ridge
    standardize::Bool
end

struct AdaptiveRidge <: AbstractEstimator
    lambda::Float64
    mode::Symbol             # :oneshot | :iterative
    gamma::Float64           # :oneshot only
    init::Symbol             # :oneshot only (:ols | :ridge)
    init_lambda::Float64
    epsilon::Float64         # :iterative only (added to β_j²)
    max_iter::Int            # :iterative only
    tol::Float64             # :iterative only (L-infinity change in coefficients)
    standardize::Bool
end
```

Each `solve_coefficients` method calls a shared helper
`_glmnet_solve(X, y, alpha, penalty_factor; ...)`. Internal helpers
`_initial_estimate` / `_adaptive_penalty_factor` /
`_iterative_penalty_factor`.

Export `Lasso` / `AdaptiveLasso` / `AdaptiveRidge`; re-export from
`Magesty.jl`.

## Test plan

1. **Agrees with OLS as `λ -> 0` (requirement)**: on small,
   well-conditioned `X, y`, `Lasso(λ=1e-10)` / `AdaptiveLasso(λ=1e-10)` /
   `AdaptiveRidge(λ=1e-10, mode=:oneshot)` /
   `AdaptiveRidge(λ=1e-10, mode=:iterative)` match `OLS()` within
   `atol`. Includes the bias coefficient (`j0`).
2. **Bias column excluded from penalty**: even with an extreme penalty,
   the bias_col coefficient equals `mean(y)` (in the limit where the
   others collapse to zero).
3. **Adaptive LASSO oracle property**: on synthetic sparse-truth data,
   zero-coefficients actually go to zero.
4. **AdaptiveRidge `:iterative` convergence**: meets `tol` within
   `max_iter`.
5. **GLMNet vs MultivariateStats.ridge agreement**: with `alpha=0,
   penalty_factor=unit`, the GLMNet solution agrees with the existing
   Ridge (numerical reproducibility).
6. `make test-all`, `make test-jet`, `make test-aqua` all pass.

## Files to change

- `Project.toml` — add `GLMNet` to `[deps]` and `[compat]`.
- `src/Fitting.jl` — `using GLMNet`; three types, methods, helpers,
  exports.
- `src/Magesty.jl` — re-exports.
- `test/component/test_fitting_estimators.jl` (new file or appended to
  an existing one).

## Workflow

Mid-sized feature + multiple design decisions + new external dependency
= spec needed. Once agreed, create
`docs/specs/260516-lasso-adaptive-estimators/`
(requirements / design / tasklist) and implement on a
`refactor/lasso-adaptive-estimators` branch
([[feedback_branch_for_medium_refactor]]).

## Open items

- `standardize` defaults to OFF (consistent with Ridge); confirm.
- `AdaptiveLasso` / `AdaptiveRidge(:oneshot)` defaults: `init = :ols`,
  `γ = 1.0`.
- `AdaptiveRidge(:iterative)` defaults: `epsilon=1e-8`, `max_iter=50`,
  `tol=1e-6`.
- Whether to fold the existing `Ridge` into GLMNet later (separate
  spec).
