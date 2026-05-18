# Design: ElasticNet estimator (with Lasso convenience)

Status: draft (2026-05-17, rebased on top of 260518 energy-centered design matrix)
Spec: [`requirements.md`](requirements.md)

## Summary

Add `ElasticNet <: AbstractEstimator` to `src/Fitting.jl`, backed by
GLMNet.jl. The struct exposes `alpha`, `lambda`, and `standardize`,
covering Lasso (`alpha=1`), GLMNet-style Ridge (`alpha=0`), and honest
Elastic Net (mixed norm). A `Lasso(; lambda, standardize=true)`
convenience constructor returns an `ElasticNet(alpha=1.0, ...)` to
keep call sites readable.

Since spec 260518 (merged 2026-05-18), `assemble_weighted_problem`
already mean-centers the energy block and produces an `X` with no
bias column, with `j0` recovered analytically by `extract_j0_jphi`
after the solve. `ElasticNet` joins this convention: its
`solve_coefficients` method is a thin wrapper around `glmnet(...,
intercept=false)`. No in-method centering, no column stripping, no
`bias_col` plumbing.

`OLS` and the analytic `Ridge` are unchanged. An agreement test
between `ElasticNet(alpha=0, lambda=λ, standardize=false)` and
`Ridge(lambda=λ)` is included so the relationship between the GLMNet
coordinate-descent L2 solver and the analytic L2 solver is documented
and trackable.

## Why `intercept=false`

`assemble_weighted_problem` builds the augmented system
`[√((1-w)/n_E) · (X_e - μ_X_e); √(w/n_T) · X_t]` and
`[√((1-w)/n_E) · (y_e - μ_y_e); √(w/n_T) · y_t]`, with the energy
block centered by its own (energy-only) mean before scaling and the
torque block left alone. This is the closed-form substitution of
`j0(jphi) = mean_e(y_e - X_e · jphi)` into the weighted
least-squares + Elastic-Net loss

```
min_{j0, jphi} (1 - w)/n_E · ||y_e - j0 - X_e · jphi||²
             +       w/n_T · ||y_t        - X_t · jphi||²
             + λ · P_α(jphi)
```

after using `∂L/∂j0 = 0`. The torque block does not depend on `j0`
(no intercept enters the torque equations of motion), so the
substitution touches only the energy block; the `λ · P_α(jphi)`
term is untouched (`j0` is not penalised).

Asking GLMNet to *also* fit an intercept on top of this would
re-introduce a uniform offset across both energy and torque rows,
contradicting the structure of the augmented system. `intercept=false`
is therefore required, and the post-solve `extract_j0_jphi` step
recovers `j0` from the un-scaled energy residual exactly as it does
for OLS / Ridge.

(For the full derivation of the energy-only centering rule, see
[`260518-energy-centered-design-matrix/design.md`](../260518-energy-centered-design-matrix/design.md).)

## Module layout

| Target | Change |
|---|---|
| `Project.toml` | Add `GLMNet` to `[deps]`. In `[compat]`, pin `GLMNet = "0.7"` (current stable; bump after verifying the API used here is stable across the chosen range during M1). |
| `src/Fitting.jl` | `using GLMNet`. Add `ElasticNet` struct + ctor, `Lasso` convenience ctor, `solve_coefficients(::ElasticNet, X, y)`, and a private `_glmnet_solve` helper. Update `export`. |
| `src/Magesty.jl` | Re-export `ElasticNet` and `Lasso`. No call-site changes (the `solve_coefficients(estimator, X, y)` signature is unchanged). |
| `test/component/test_fitting_estimators.jl` | New or extended; tests listed in "Test strategy". |
| `test/integration/<one chosen example>/test.jl` | Add a `Lasso(...)` smoke path or a new integration test file alongside it. |
| `docs/src/api.md` | List `ElasticNet` / `Lasso` next to `OLS` / `Ridge`. |
| `CHANGELOG.md` | `[Unreleased]` entry. |

`assemble_weighted_problem` and `extract_j0_jphi` are **unchanged**.

## API

```julia
"""
    ElasticNet(; alpha::Real, lambda::Real, standardize::Bool = true)

Elastic Net estimator backed by GLMNet.jl. Covers Lasso (`alpha=1`),
GLMNet-style L2 (`alpha=0`), and honest Elastic Net (mixed norm).
`standardize=true` neutralises the per-cluster `(4π)^(N/2)` column
scale baked into the SCE design matrix.

`j0` is eliminated analytically inside `assemble_weighted_problem`
(energy-only centering) and re-fit from the un-scaled energy residual
by `extract_j0_jphi` downstream — the same post-processing OLS and
Ridge already use. This estimator therefore calls GLMNet with
`intercept=false`.

# Fields
- `alpha::Float64`  : `0 ≤ alpha ≤ 1`.
- `lambda::Float64` : `λ ≥ 0`. Penalty strength.
- `standardize::Bool` : passed to GLMNet.
"""
struct ElasticNet <: AbstractEstimator
    alpha::Float64
    lambda::Float64
    standardize::Bool
end

function ElasticNet(; alpha::Real, lambda::Real, standardize::Bool = true)
    0.0 <= alpha <= 1.0 ||
        throw(ArgumentError("ElasticNet alpha must satisfy 0 ≤ alpha ≤ 1; got $alpha"))
    lambda >= 0.0 ||
        throw(ArgumentError("ElasticNet lambda must be non-negative; got $lambda"))
    return ElasticNet(Float64(alpha), Float64(lambda), standardize)
end

"""
    Lasso(; lambda::Real, standardize::Bool = true) -> ElasticNet

Convenience **function** (not a type) that returns
`ElasticNet(alpha = 1.0, lambda, standardize)`. There is no separate
`Lasso` struct, so:

- `Lasso(lambda=…) isa ElasticNet` is `true`.
- `Lasso(lambda=…) isa Lasso` is **not valid** — `Lasso` is a
  function, not a type. Code that needs to detect the Lasso case
  should check `e isa ElasticNet && e.alpha == 1.0`.

Reason for keeping `Lasso` as a function rather than a `const`
type alias: an alias would make `Lasso(alpha = 0.5, lambda = …)`
parse, contradicting the meaning "α = 1 only". A standalone
function fixes `alpha = 1.0`.
"""
Lasso(; lambda::Real, standardize::Bool = true) =
    ElasticNet(alpha = 1.0, lambda = lambda, standardize = standardize)
```

### `solve_coefficients` method

The abstract contract is `solve_coefficients(estimator, X, y) ->
Vector{Float64}`, unchanged. `ElasticNet` adds its method:

```julia
function solve_coefficients(
    e::ElasticNet,
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
)
    return _glmnet_solve(
        X, y;
        alpha       = e.alpha,
        lambda      = e.lambda,
        standardize = e.standardize,
    )
end
```

The `X` arriving here has already been centered (energy block) and
row-scaled by `assemble_weighted_problem`, with no bias column.
GLMNet sees the right problem directly.

### Helper

```julia
# Private. Returns beta of length size(X, 2).
function _glmnet_solve(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real};
    alpha::Real,
    lambda::Real,
    standardize::Bool,
)
    path = glmnet(
        Matrix{Float64}(X), Vector{Float64}(y);
        alpha       = Float64(alpha),
        lambda      = [Float64(lambda)],
        standardize = standardize,
        intercept   = false,
    )
    return vec(Matrix(path.betas))   # (p, 1) -> (p,)
end
```

`intercept=false` is the critical flag: without it, GLMNet would
re-introduce a uniform intercept across all rows, undoing the
energy-only centering done upstream in
`assemble_weighted_problem`.

## Types and conventions

- `ElasticNet` is a concrete subtype of `AbstractEstimator`. No new
  abstract type.
- `Lasso` is a function, not a struct. Code that wants to test
  "is this Lasso?" should check
  `e isa ElasticNet && e.alpha == 1.0`.
- The convention "energy unit in, energy unit out" is preserved by
  GLMNet's `standardize=true` round-trip; the returned `beta` is in
  the original column units.

## Impact on linked sites

Sites named in CLAUDE.md "Linked sites":

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): unaffected.
- [ ] SCE coefficient XML (`save` / `load`): unaffected — `ElasticNet`
      only changes how `jphi` is fit, not how it is serialized.
- [ ] `Fitting` ↔ `SALCBasis`: unaffected — column order, key-group
      order, and SALC indexing are unchanged. `ElasticNet` consumes
      whatever design matrix the existing pipeline produces.
- [ ] `.claude/agents/` references: review the `code-reviewer` and
      `test-runner` agent prompts for any place that enumerates
      estimator names. Update to include `ElasticNet` / `Lasso` if so.
- [ ] `SPEC.md` / `docs/src/api.md`: add `ElasticNet` and `Lasso` to
      the public estimator list.

## Test strategy

New tests live in `test/component/test_fitting_estimators.jl` (extend
the existing file; if it does not yet exist, create it alongside the
other component tests). All tests use small synthetic data so they
run in well under a second.

1. **`λ → 0` limit agrees with OLS, energy-only (`weight = 0`)** (strict).
   `ElasticNet(alpha=1, lambda=1e-10, standardize=false)` and
   `ElasticNet(alpha=0, lambda=1e-10, standardize=false)` match `OLS()`
   on the same `(X, y)` produced by `assemble_weighted_problem`,
   within an `atol` chosen for GLMNet's coordinate-descent precision
   (start at `1e-7`; record the measured worst case in the test
   docstring). Cleanest sanity check that the wrapper is wired
   correctly.
2. **`λ → 0` limit, mixed objective (`weight = 0.5`)** (strict).
   Same as test 1 but with `weight = 0.5` so the augmented system
   contains both energy and torque rows. Confirms `intercept=false`
   is the right call given that the energy block has been centered
   by `assemble_weighted_problem`. Tolerance matches test 1.
3. **`λ → 0` limit with `standardize=true`** (loose).
   Same as test 1 but with `standardize=true`, looser `atol`
   (`1e-5`), documenting the standardisation/back-transform round-trip
   precision. One sub-case with `weight = 0` and one with `weight = 0.5`.
4. **GLMNet `α=0` vs analytic `Ridge` agreement.**
   For several `λ` values (e.g. `[1e-3, 1e-2, 1e-1, 1.0]`), compare
   `solve_coefficients(ElasticNet(alpha=0, lambda=λ, standardize=false), X, y)`
   against `solve_coefficients(Ridge(lambda=λ), X, y)`. Pass
   condition: max relative L∞ difference on `jphi` below
   `rtol = 1e-4` (initial guess — adjust in the implementation pass
   if GLMNet's internal `λ`-scaling shifts the effective penalty).
   Whatever tolerance ends up being used, it is recorded in
   (a) the test's docstring, (b) `CHANGELOG.md`'s `[Unreleased]`
   entry, (c) `tasklist.md` exit checklist.
5. **Lasso sparsity monotonicity.**
   On synthetic data with a known sparse ground truth, increasing
   `λ` at `alpha=1, standardize=true` monotonically (weakly)
   increases the number of exactly-zero coefficients.
6. **Standardisation qualitative effect.**
   On a fabricated SCE-like design matrix with columns scaled by
   factors of `[1, 12, 45, 158]` (mimicking `(4π)^(N/2)` for N=1..4)
   and a sparse ground truth that lives mostly on the large-scale
   columns, fit at `alpha=1` once with `standardize=true` and once
   with `standardize=false` at the same `λ`. Assert that the
   `standardize=true` fit retains the high-N columns at the chosen
   `λ` while the `standardize=false` fit drops them earlier. Use
   counts (`count(!iszero, beta)`) and whether specific known-true
   columns are non-zero.
7. **Integration smoke test.**
   Pick one existing `test/integration/*` whose fit currently uses
   OLS or Ridge with a clean L2 result. Add a second pass that calls
   `fit(SCEFit, dataset, Lasso(lambda=…))` and confirms (a) it runs,
   (b) the returned `SCEModel` has the expected field shapes, (c)
   `Magesty.save` / `Magesty.load` round-trip the resulting model
   byte-for-byte. Numerical fit quality is not asserted (Lasso
   selection is a separate property tested above).
8. **Static analysis.** `make test-jet`, `make test-aqua` remain
   clean (no new warnings introduced by `using GLMNet`).

(The "bias not penalised" test from the pre-260518 draft is dropped:
there is no longer a bias column to penalise or exempt. The
analogous physics — `j0` recovered as `mean(y_e - X_e · jphi)` — is
already covered by spec 260518's regression test.)

## Risks and open items

- **GLMNet `λ` convention.** GLMNet internally normalises the
  penalty by `sum(weights) / n` (or similar) before applying it. The
  agreement test (test 4) may fail at the initial `rtol = 1e-4` and
  need either tolerance relaxation or a small constant correction
  factor. If the latter is required, the spec records the factor
  rather than hiding it: callers see `λ` as GLMNet sees it, and the
  test docstring documents the offset.
  **Guard condition**: if the empirically required `rtol` exceeds
  `1e-2`, treat this as a sign that the GLMNet and analytic Ridge
  models disagree on more than rounding noise (different `λ`
  semantics, different scaling). In that case pause and revisit
  the spec — do not simply relax the tolerance further to make the
  test pass. The expected outcome is `rtol ∈ [1e-7, 1e-3]`; outside
  that window the design assumption "GLMNet `α=0` is the same
  problem as analytic Ridge" should be re-examined.
- **`weight = 1` edge case.** Strictly speaking this means "torque
  loss only, ignore energies", but the augmented system still
  contains the (centered) energy rows; their σ-scaled contribution
  is just zero. What `weight = 1` actually breaks is
  `extract_j0_jphi`: with `jphi` fit from torques alone, the
  un-scaled energy residual `mean(y_e - X_e · jphi)` is computed
  against potentially arbitrary energies, which is well-defined but
  physically meaningless. This is already the situation for OLS /
  Ridge today and is not in scope to fix here.
- **JET / Aqua warnings on GLMNet.** GLMNet may surface latent type
  instabilities or unused-binding warnings. If so, isolate them
  inside `_glmnet_solve` and add a targeted `@nowarn` (with a
  comment) rather than restructuring GLMNet's API in our code.

## Notes

- **Branch name.** The branch is `refactor/lasso-estimator` from
  the earlier scoping conversation; the rename to
  `refactor/elasticnet-estimator` is deliberately not done — the
  spec is the authoritative scope record and the branch name is a
  minor inconvenience.
