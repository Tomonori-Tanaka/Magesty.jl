# Requirements: Adaptive Lasso (oneshot) estimator

Status: draft (2026-05-18)

## Goal

Add an `AdaptiveLasso` estimator to the `AbstractEstimator` dispatch
hierarchy in `Fitting.jl`. The estimator performs the **oneshot
adaptive Lasso of Zou (2006)**: a pilot regression produces an
initial coefficient vector `beta_pilot`, then a weighted-L1 Lasso is
fit with per-column penalty weights `w_j = 1 / |beta_pilot[j]|^gamma`.
Backed by GLMNet.jl through the existing `_glmnet_solve` helper, which
is extended with an optional `penalty_factor` keyword.

The implementation matches the formulation used by ALAMODE
(<https://alamode.readthedocs.io/en/latest/almdir/formalism_alm.html>),
which fixes `gamma = 1` and uses OLS as the pilot. This spec keeps
the same defaults but exposes `gamma` and `pilot` as overrides so
unit tests can exercise `gamma = 0` (= plain Lasso) for agreement
checks, and so users can swap in `pilot = Ridge(lambda = ...)` when
the SCE design matrix is rank-deficient.

Adaptive Ridge (the Frommlet--Nuel 2016 iterative L0 approximation)
is a separate, follow-up spec.

## Background

The plain `Lasso` / `ElasticNet` paths shipped in spec 260517
(`refactor/lasso-estimator`, merged as PR #20). Two known weaknesses
of plain Lasso are addressed by the oneshot adaptive variant:

1. **Bias on large coefficients.** Plain L1 shrinks every active
   coefficient toward zero by a constant amount, leaving a residual
   bias on the truly-nonzero entries. The adaptive weights
   `1/|beta_pilot|^gamma` send the per-column effective penalty to
   zero on large pilot coefficients, recovering oracle unbiasedness
   on truly-nonzero entries at sufficiently large sample size.
2. **Inconsistent variable selection.** Plain Lasso is not
   selection-consistent in general; adaptive Lasso is, under mild
   conditions on the pilot (Zou 2006, Theorem 2). For SCE this means
   weak / spurious cluster contributions are pushed to exact zero
   more reliably than under plain Lasso.

ALAMODE uses this exact formulation for force-constant fitting:
`Phi_adalasso = argmin (1/(2 N_d)) ||F_DFT - A Phi||_2^2 + alpha
sum_i w_i |Phi_i|` with `w_i = 1 / |Phi_OLS,i|`. We mirror the recipe
(with `lambda` for `alpha` and an exposed `gamma` for the exponent),
and plug it into the existing `solve_coefficients` dispatch so the
post-processing (energy-mean-centering upstream, `j0` recovery
downstream) is identical to `OLS` / `Ridge` / `ElasticNet`.

The pilot default is `OLS()` (Zou 2006 verbatim, ALAMODE default).
The docstring will recommend `pilot = Ridge(lambda = small)` when
the SCE design is rank-deficient or near-collinear (e.g.
`num_salcs >= num_spinconfigs` at `torque_weight = 0`, the
situation that required dropping the FeGe (`fege_2x2x2` integration
fixture) Ridge cases in spec 260518). In that regime an OLS pilot
returns Julia's minimum-norm solution, which is mathematically
valid but populates null-space directions with noise of magnitude
~1e-10 -- not small enough to be clipped by the default
`epsilon = eps(Float64)`. That noise then drives wildly
miscalibrated `penalty_factor` entries, which interact poorly with
GLMNet's internal `penalty_factor` rescaling (see `design.md`
"Risks and open items").

## Scope

Includes:

- New estimator `AdaptiveLasso <: AbstractEstimator` with fields
  `pilot::AbstractEstimator`, `lambda::Float64`, `gamma::Float64`,
  `epsilon::Float64`, `standardize::Bool`. Keyword constructor
  validates `lambda >= 0`, `gamma >= 0`, `epsilon > 0`.
- Defaults: `pilot = OLS()`, `gamma = 1.0`, `epsilon = eps(Float64)`,
  `standardize = true`.
- `solve_coefficients(::AdaptiveLasso, X, y)` method that:
  1. Calls `solve_coefficients(e.pilot, X, y)` to get `beta_pilot`.
  2. Builds `penalty_factor[j] = 1 / max(|beta_pilot[j]|, epsilon)^gamma`.
  3. Calls the extended `_glmnet_solve(...; penalty_factor)` with
     `alpha = 1.0`, `lambda = e.lambda`, `standardize = e.standardize`,
     `intercept = false`.
- Extend `_glmnet_solve` (in `src/Fitting.jl`) with an optional
  keyword `penalty_factor::Union{Nothing, AbstractVector{<:Real}} =
  nothing`. When `nothing`, the call to GLMNet omits the argument
  (preserving current `ElasticNet` / `Lasso` behaviour
  byte-for-byte). When a vector is given, it is forwarded to
  GLMNet's `penalty_factor` argument.
- Re-export `AdaptiveLasso` from `src/Magesty.jl`.
- Component tests (extending `test/component/test_fitting_estimators.jl`)
  covering: `gamma = 0` reduces to plain Lasso; weight construction
  is exactly `1 / max(|beta_pilot|, eps)^gamma`; `epsilon` clip is
  active for a Lasso-as-pilot synthetic case; both `pilot = OLS()`
  and `pilot = Ridge(lambda = ...)` paths produce sensible fits;
  sparse-recovery superior to plain Lasso at matched `lambda`.
- Documentation: `docs/src/api.md` lists `AdaptiveLasso`; `SPEC.md`
  estimator list extends; `CHANGELOG.md` `[Unreleased]` entry.

Excludes:

- Adaptive Ridge (`:oneshot` and `:iterative` variants). Separate spec.
- Cross-validated `lambda` / `gamma` selection. Out of scope; users
  pass `lambda` directly. A user-side CV loop can wrap this
  estimator.
- Adaptive ElasticNet (mixed-norm adaptive). Strictly speaking the
  same machinery generalizes (`alpha != 1.0` with custom
  `penalty_factor`); deferred until a concrete use case appears.
- Modifying `assemble_weighted_problem` or `extract_j0_jphi`. Both
  are reused as-is (post-260518 form: energy-mean-centered, no bias
  column).

## Invariants

- `OLS`, `Ridge`, and `ElasticNet` / `Lasso` numerical behaviour are
  byte-for-byte unchanged. Extending `_glmnet_solve` with an optional
  `penalty_factor = nothing` default preserves the existing
  call-through to GLMNet (the `penalty_factor` argument is omitted
  when `nothing`).
- `assemble_weighted_problem` and `extract_j0_jphi` are unchanged.
  `AdaptiveLasso` consumes the same `(X, y)` shape as `OLS` / `Ridge`
  / `ElasticNet` and returns a `jphi`-only coefficient vector.
- `Jphi` coefficients returned by `AdaptiveLasso` are in the input
  energy unit (typically eV); GLMNet's `standardize` round-trip
  reverses the internal column scaling before returning, as for
  `ElasticNet`.
- `j0` recovery still happens downstream in `extract_j0_jphi` via
  `mean(y_E - X_E * jphi)`. Torque rows still do not see `j0`.
- Spin direction stays unit-vector with `3 x n_atoms` layout
  (untouched).
- Tesseral spherical-harmonics conventions (untouched).
- SALC ordering and the `Fitting` <-> `SALCBasis` <-> XML link (see
  CLAUDE.md "Linked sites") are untouched.

## Completion criteria

- [ ] `make test-all` passes.
- [ ] `make test-jet`, `make test-aqua` clean (no new warnings).
- [ ] New component tests verify the five properties listed in
      `design.md` "Test strategy".
- [ ] `gamma = 0` agreement test: at `gamma = 0`, `AdaptiveLasso(...)`
      matches `Lasso(lambda = e.lambda, standardize = e.standardize)`
      to within an `atol` chosen for GLMNet's coordinate-descent
      precision. Target ~`1e-7`; record the measured worst case in
      the test docstring and in `CHANGELOG.md`.
- [ ] Sparse-recovery test: on a synthetic problem where plain Lasso
      retains spurious columns at a chosen `lambda`, `AdaptiveLasso`
      (`pilot = OLS()`, `gamma = 1`, same `lambda`) selects the
      correct support.
- [ ] `AdaptiveLasso` appears in `docs/src/api.md` alongside `OLS` /
      `Ridge` / `ElasticNet` / `Lasso`.
- [ ] `CHANGELOG.md` `[Unreleased]` records the new estimator and the
      `_glmnet_solve` extension.
- [ ] `Status:` line in this spec's `tasklist.md` and the matching
      row in `docs/specs/README.md` updated to `complete` together.

## References

- Zou, H. (2006). "The Adaptive Lasso and Its Oracle Properties".
  *J. Am. Stat. Assoc.* 101, 1418-1429.
- ALAMODE documentation, "Adaptive LASSO" section:
  <https://alamode.readthedocs.io/en/latest/almdir/formalism_alm.html>
- Related spec (parent estimator path):
  [`260517-elasticnet-estimator/`](../260517-elasticnet-estimator/)
- Related spec (energy-centered design matrix that this estimator
  consumes via `assemble_weighted_problem`):
  [`260518-energy-centered-design-matrix/`](../260518-energy-centered-design-matrix/)
- Design note (Adaptive sketch retained as pre-spec context;
  superseded by this spec, and slated for removal in M4 wrap-up so
  the basename reference does not outlive the spec):
  [`../../design-notes/lasso-adaptive-estimators.md`](../../design-notes/lasso-adaptive-estimators.md)
