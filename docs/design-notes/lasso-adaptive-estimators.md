# ElasticNet / LASSO / Adaptive estimators

**Status**: ElasticNet (incl. Lasso convenience constructor) → spec
started (2026-05-17); Adaptive variants → on hold

Design note for adopting GLMNet.jl and adding the L1 / mixed-norm /
Adaptive families to the estimator-dispatch hierarchy in `Fitting.jl`.
The first concrete deliverable is `ElasticNet` (with a `Lasso`
convenience constructor) — see the active spec for that work. Adaptive
LASSO / Adaptive Ridge remain as pre-spec sketches at the bottom of
this note.

## Background

The SCE coefficient estimator in `src/Fitting.jl` currently supports
only `OLS` (`X \ y`) and `Ridge` (`MultivariateStats.ridge` with a
per-column penalty that excludes the bias column). To obtain sparse
effective spin models we want L1-style regularization (LASSO) and the
continuous interpolation to L2 (Elastic Net). Adaptive LASSO and
Adaptive Ridge are natural follow-ups but are out of scope for the
first spec.

Historical note: an earlier `ElasticNet(alpha, lambda)` struct existed
but silently ignored `alpha`, so it was renamed to `Ridge` in
[`docs/specs/260513-estimator-dispatch/`](../specs/260513-estimator-dispatch/),
with the explicit intent of freeing the `ElasticNet` name for a real
mixed-norm estimator later. This spec is that follow-up.

## Library selection

`MultivariateStats.jl` does not provide LASSO or honest Elastic Net.
Comparison of candidates:

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

The existing `Ridge` keeps using `MultivariateStats.ridge` — it is the
analytic solver (`(XᵀX + λI)β = Xᵀy`), already validated by the test
suite, and worth preserving as a numerically clean reference for the L2
problem. `ElasticNet(alpha=0, lambda=λ)` should be numerically close to
`Ridge(lambda=λ)` but **not necessarily byte-identical** — GLMNet uses
coordinate descent with its own internal `λ`-scaling convention. The
spec includes an explicit agreement test (see below); regardless of the
outcome, both estimators stay in the public API.

## Conventions (current scope: ElasticNet + Lasso convenience)

| Item | Policy |
|---|---|
| Bias column | Excluded from the augmented `X` passed to GLMNet; recovered via GLMNet's own `intercept=true` (option A). See "Bias-column handling" below. |
| `j0` extraction | Read directly from `GLMNetPath.a0` (no `extract_j0_jphi` call for ElasticNet) |
| Standardize | **ON by default** for `ElasticNet`. Rationale: the per-cluster prefactor `(4π)^(N/2)` baked into each column of the SCE design matrix scales columns by ~3.5×–158× across N=1..4 body terms; a uniform L1 (or mixed) penalty on unscaled columns under-penalises high-N clusters and biases selection. Standardisation removes that scale artefact. Exposed as `standardize::Bool` kwarg so the choice is auditable per fit. The existing `Ridge` keeps its current behaviour (no standardisation) — its analytic solver is unchanged. |
| `α` range | `0 ≤ alpha ≤ 1`. Constructor validates. `alpha=1` is Lasso, `alpha=0` is the GLMNet Ridge analog, intermediate values are honest Elastic Net. |
| `λ` selection | User-specified single `λ` (matches Ridge). Lambda-path / CV are explicitly out of scope; revisit in a separate spec. |
| Augmented system | The existing `assemble_weighted_problem` output `(X, y, bias_col)` is consumed by **stripping the bias column** before handing the rest to GLMNet (option A). `bias_col` is still threaded through `solve_coefficients` so OLS / Ridge methods stay identical. |

### Bias-column handling — why option A (GLMNet `intercept=true`)

The augmented `X` from `assemble_weighted_problem` has a constant-1 bias
column at `bias_col`. With `standardize=true`, that column has zero
standard deviation, which GLMNet would either reject or special-case
unpredictably. Two viable approaches:

- **(A) Strip the bias column and use `intercept=true`.** GLMNet centers
  `y` and the remaining columns, fits Elastic Net, then restores both
  the intercept (`a0`) and the slope coefficients in original units.
  This is the canonical glmnet usage pattern and is what R/Python users
  get by default. `j0 = a0`; `jphi = beta`. `extract_j0_jphi` is
  bypassed for the ElasticNet path.

- **(B) Keep the augmented `X` and write a custom standardiser** that
  centers/scales all columns except `bias_col`, calls GLMNet with
  `intercept=false, standardize=false`, and inverts the scaling on the
  returned coefficients.

Option A is chosen: it is shorter, matches the reference implementation
exactly, and avoids reinventing well-tested glmnet bookkeeping. Option
B is recorded here so the trade-off is on file.

A side-effect of option A is that the energy-row weighting `√(1 - w)`
is already baked into the rows of `X` and `y` from
`assemble_weighted_problem`, so what GLMNet sees is the same weighted
least-squares problem the other estimators see. The reset
`Xe[:, 1] .= 1.0` step in `assemble_weighted_problem` becomes a no-op
for the ElasticNet path (we drop that column), but is harmless.

## Draft API (current scope)

```julia
struct ElasticNet <: AbstractEstimator
    alpha::Float64       # 0 ≤ alpha ≤ 1; 0=Ridge, 1=Lasso
    lambda::Float64
    standardize::Bool
end

ElasticNet(; alpha::Real, lambda::Real, standardize::Bool = true) =
    ElasticNet(Float64(alpha), Float64(lambda), standardize)

# Convenience constructor — returns an ElasticNet, not a new struct.
Lasso(; lambda::Real, standardize::Bool = true) =
    ElasticNet(alpha = 1.0, lambda = lambda, standardize = standardize)
```

`solve_coefficients(::ElasticNet, X, y; bias_col)` strips
`X[:, bias_col]`, calls a thin internal helper
`_glmnet_solve(X_nobias, y; alpha, lambda, standardize)`, and
reassembles the full-length coefficient vector with `a0` placed at
`bias_col`. The `_glmnet_solve` helper is private and intentionally
narrow so Adaptive variants can reuse it later.

Export `ElasticNet` and `Lasso`; re-export from `Magesty.jl`.

## Test plan (ElasticNet)

1. **Agrees with OLS as `λ → 0` (requirement)**: on small,
   well-conditioned `X, y`, `ElasticNet(alpha=1, lambda=1e-10, standardize=false)`
   and `ElasticNet(alpha=0, lambda=1e-10, standardize=false)` both match
   `OLS()` within `atol`. Includes the bias coefficient (`j0`). A
   separate looser-tolerance variant confirms the `standardize=true`
   path also converges to OLS up to GLMNet's standardisation/back-
   transform round-trip precision.
2. **GLMNet `α=0` vs `Ridge` agreement**: on a synthetic dataset,
   compare `ElasticNet(alpha=0, lambda=λ, standardize=false)` against
   `Ridge(lambda=λ)` for a handful of `λ` values. Document the
   observed disagreement (expected ~ GLMNet `tol`, but may include a
   small constant offset from GLMNet's `λ`-scaling convention). Goal:
   confirm "close enough" qualitatively and pin down the discrepancy
   model. Pass condition: max relative difference of `jphi` below
   some agreed `rtol` (e.g. `1e-4`) — exact threshold to be set after
   first measurement. If they cannot be reconciled within a sensible
   tolerance, the spec records this and both estimators stay
   independent (no Ridge deprecation either way; that's out of
   scope).
3. **Bias column not penalised**: with a synthetic dataset where `j0`
   is far from zero, even a large `λ` leaves `j0` close to `mean(y)`
   (the un-penalised intercept), while non-bias coefficients shrink
   toward zero.
4. **Sparsity sanity check (Lasso path)**: on synthetic data with a
   known sparse ground truth, increasing `λ` at `alpha=1` monotonically
   increases the number of exactly-zero coefficients.
5. **Standardisation effect**: on a fabricated SCE-like design matrix
   with deliberately large between-column scale spread, compare
   selected supports for `standardize=true` vs `false` at the same
   `λ` and document the qualitative difference. Not a strict numerical
   assertion — this test is for human inspection / regression
   detection of API behaviour.
6. **Integration smoke**: at least one existing `test/examples/*` is
   exercised through `fit(SCEFit, dataset, Lasso(lambda=...))` to
   confirm the public API plumbs through correctly.
7. `make test-all`, `make test-jet`, `make test-aqua` all pass.

## Files to change (ElasticNet)

- `Project.toml` — add `GLMNet` to `[deps]` and `[compat]`.
- `src/Fitting.jl` — `using GLMNet`; `ElasticNet` struct, `Lasso`
  convenience constructor, `solve_coefficients(::ElasticNet, ...)`,
  private `_glmnet_solve` helper, exports.
- `src/Magesty.jl` — re-export `ElasticNet` and `Lasso`.
- `test/component/test_fitting_estimators.jl` (or appended to existing
  fitting tests).
- `docs/src/api.md` — list `ElasticNet` and `Lasso` alongside `OLS` /
  `Ridge`.
- `CHANGELOG.md` — `[Unreleased]` entry recording the new estimator,
  the GLMNet dependency, and the Ridge-agreement test outcome.

## Workflow

Mid-sized feature + multiple design decisions + new external
dependency = spec needed. Active spec for ElasticNet lives at
[`docs/specs/260517-elasticnet-estimator/`](../specs/260517-elasticnet-estimator/);
implementation runs on a `refactor/lasso-estimator` branch
([[feedback_branch_for_medium_refactor]]). Branch name is kept as
`refactor/lasso-estimator` rather than renaming — the spec is the
authoritative scope record.

---

## Adaptive variants (on hold, pre-spec sketch)

The original design note also covered Adaptive LASSO and Adaptive
Ridge. They are deferred until the plain ElasticNet path ships and we
have practical experience with GLMNet on real SCE data. Sketch
retained below for the follow-up spec.

### Estimator mapping

| Estimator | `alpha` | `penalty_factor[j != bias]` | Note |
|---|---|---|---|
| ElasticNet (incl. Lasso at α=1) | 0..1 | 1 | shipping in current spec |
| AdaptiveLasso | 1.0 | `1/|β̂_init,j|^γ` | `β̂_init` from OLS or Ridge (Zou 2006) |
| AdaptiveRidge (`mode = :oneshot`) | 0.0 | `1/|β̂_init,j|^γ` | Single-stage reweighting |
| AdaptiveRidge (`mode = :iterative`) | 0.0 | `1/(β_j² + ε)` updated iteratively | Frommlet & Nuel 2016. Refits GLMNet each iteration as an L0 approximation |

Behavior-descriptive symbols `:oneshot` / `:iterative` are preferred
over author-name symbols, since the latter assume the reader knows the
paper.

### Open items for the follow-up

- Whether Adaptive variants should also default to `standardize=true`
  (likely yes, by the same column-scale argument; confirm).
- `AdaptiveLasso` / `AdaptiveRidge(:oneshot)` defaults: `init = :ols`,
  `γ = 1.0`.
- `AdaptiveRidge(:iterative)` defaults: `epsilon=1e-8`, `max_iter=50`,
  `tol=1e-6`.
- Whether to fold the existing `Ridge` into GLMNet later (separate
  spec). Outcome depends on the agreement test in this spec.
