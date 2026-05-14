# Requirements: Optimize.jl estimator dispatch refactor

Status: draft (2026-05-13)
Owner: T. Tanaka
Branch: `refactor/estimator-dispatch`

## Goal

Restructure the regression dispatch in `src/Optimize.jl` so that adding a new
estimator (Ridge, Lasso, Bayesian, NNLS, ...) requires only:

1. A new `<: AbstractEstimator` struct carrying its hyperparameters.
2. A single new method of `solve_coefficients(::NewEstimator, X, y; bias_col)`.

The current `isa` chain in `_fit_sce_model_internal` and the
"assemble + solve + extract" monolith inside `elastic_net_regression` are
replaced by three orthogonal layers connected by multiple dispatch.

## Scope

In scope:

- `src/Optimize.jl`: `_fit_sce_model_internal`, `fit_sce_model_ols`,
  `fit_sce_model_elastic_net`, `elastic_net_regression`, the `Optimizer`
  outer constructor signature, and the `AbstractEstimator` type hierarchy.
- Internal helper introductions: `assemble_weighted_problem`,
  `solve_coefficients`, `extract_j0_jphi`.
- **Rename `ElasticNet` → `Ridge` (breaking).** The current
  `ElasticNet` struct is functionally a pure L2 ridge (its `alpha`
  field is silently ignored). The struct is renamed to `Ridge` with a
  `lambda::Float64` field only. The name `ElasticNet` is freed for a
  future spec that adds a real `α ∈ (0,1)` mixed-norm estimator.
- TOML schema cleanup tied to the rename: `[regression].alpha` becomes
  optional and emits a deprecation warning when present (cannot be
  fully removed in this spec without breaking every existing
  `input.toml`; full removal is a follow-up).
- All in-repo call sites updated: `test/examples/*/test.jl`,
  `test/benchmark_optimize.jl`, `tools/check_convergence_embset.jl`,
  `SPEC.md`, `docs/src/{api,tutorial,examples}.md`.

Out of scope (deferred to follow-up specs):

- Adding genuinely new estimators (Lasso / true ElasticNet / Bayesian /
  NNLS). This spec only rebuilds the existing OLS + (renamed) Ridge
  behavior on the new structure.
- StatsAPI / GLM.jl-style user-facing API (`fit(SCEModel, ...)`,
  `coef`, `intercept`, `predict`). Covered by the user-facing API
  proposal in `DESIGN_NOTES.md`.
- Full removal of `[regression].alpha` from the TOML schema (kept as
  a deprecated-but-accepted key in this spec; remove in a follow-up
  once downstream users have migrated).
- XML I/O format changes.
- `SALCBasis` / `Symmetry` / `Structure` / SALC layout — untouched.

## Invariants (must hold before vs. after)

1. **Numerical results unchanged.** For the same `(structure, basisset,
   spinconfig_list, alpha, lambda, weight)` inputs, the returned
   `(j0, jphi)` must match the current implementation to within
   floating-point noise (`≈ atol=1e-12 rtol=1e-10` on existing example
   datasets).
2. **Bias-column exclusion from regularization preserved.** The first
   column of the augmented design matrix is the bias term; its
   regularization weight stays zero regardless of estimator.
3. **Energy / torque combined regression with √weight scaling preserved.**
   `w_e = 1 - weight`, `w_m = weight`; energy rows scaled by `√w_e`,
   torque rows by `√w_m`; energy bias column reset to 1.0 after scaling.
4. **`(l, m, site)` basis ordering preserved.** `solve_coefficients` does
   not permute columns; downstream `extract_j0_jphi` and design-matrix
   contract with `SALCBases.jl` are unchanged.
5. **Public API surface, except the `ElasticNet` rename, is preserved.**
   `Optimizer`, `SCEModel`, `fit_sce_model`, `predict_energy`,
   `AbstractEstimator`, `OLS` remain exported with their current
   signatures. `ElasticNet` is removed from `export` and replaced by
   `Ridge` (single breaking change, explicitly approved as part of
   this spec). `Optimizer(..., alpha, lambda, weight, spinconfig_list;
   estimator=...)` keeps working — the positional `alpha`/`lambda`
   are honored for backward compatibility with the TOML path (a
   non-zero `alpha` emits a deprecation warning and is then ignored,
   matching today's silent behavior).
6. **Physics conventions untouched.** Real tesseral spherical harmonics,
   `3 × n_atoms` spin-direction layout, SCE coefficient units (eV)
   remain as-is.

## Non-goals

- No performance regression target beyond "no measurable slowdown" on
  `julia --project test/benchmark_optimize.jl --with-fit --samples 20`. Speedups are welcome but not required.
- No change in error-handling philosophy beyond replacing the explicit
  `throw(ArgumentError(...))` for unknown estimator with the implicit
  `MethodError` from missing `solve_coefficients` dispatch.

## Completion criteria

- [ ] `_fit_sce_model_internal` is a three-line orchestration:
      `assemble_weighted_problem` → `solve_coefficients` →
      `extract_j0_jphi`. No `isa` branches.
- [ ] `solve_coefficients(::OLS, X, y; bias_col)` and
      `solve_coefficients(::ElasticNet, X, y; bias_col)` exist and
      cover the current OLS / ElasticNet paths.
- [ ] `fit_sce_model_ols` and `fit_sce_model_elastic_net` either
      (a) become thin shims that build the estimator and call
      `_fit_sce_model_internal`, or (b) are deprecated in favor of
      passing an estimator. Decision recorded in `design.md`.
- [ ] `Optimizer` outer constructor: keep backward-compat positional
      `(alpha, lambda, weight)` signature, but document that
      `estimator` is the canonical knob. A non-zero `alpha` emits a
      deprecation warning and is then ignored.
- [ ] `ElasticNet` removed from source. `Ridge <: AbstractEstimator`
      added with a single `lambda::Float64` field. All in-repo
      call sites (`test/`, `tools/`, `docs/`, `SPEC.md`) updated.
- [ ] `Config4Optimize` keeps reading `[regression].alpha` for
      backward compatibility but emits a deprecation warning when
      the value is non-zero.
- [ ] `make test-all` passes. `make test-jet` and `make test-aqua`
      stay green.
- [ ] Numerical regression test added: a tiny synthetic fixture
      (energy + torque) where `(j0, jphi)` from the new code matches a
      golden value captured from `main` at the start of this branch.
- [ ] `julia --project test/benchmark_optimize.jl --with-fit --samples 20` shows no statistically significant
      regression (record before/after in `.claude/bench_log.md`).
- [ ] `DESIGN_NOTES.md` "設計案: estimator dispatch リファクタリング"
      section gets a one-line pointer to the implementing PR /
      "完了" marker.

## Linked sections (must update or check if touched)

- `src/Magesty.jl` — re-export list: replace `ElasticNet` with `Ridge`.
- `src/utils/ConfigParser.jl` — `Config4Optimize` field/default for
  `alpha` (kept, depwarn) and the constructor that builds the
  estimator (build `Ridge(lambda)` instead of
  `ElasticNet(alpha, lambda)`).
- `test/examples/*/test.jl` — they call both `fit_sce_model` and
  `Optimizer(...)`. Smoke-test all of them after the refactor.
  `fept_tetragonal_2x2x2/test.jl` constructs `ElasticNet(...)` and
  must be updated to `Ridge(lambda=...)`.
- `test/benchmark_optimize.jl` — update to `Ridge(lambda=...)`.
- `tools/check_convergence_embset.jl` — update to
  `Ridge(lambda=cfg_opt.lambda)`.
- `SPEC.md` / `docs/src/api.md` / `docs/src/tutorial.md` /
  `docs/src/examples.md` — replace `ElasticNet` with `Ridge` in
  prose and code samples.
- `DESIGN_NOTES.md` §"設計案: Optimize.jl の estimator dispatch
  リファクタリング" — keep as-is until merged; then mark "完了"
  with PR link.

## Risk register

| Risk | Mitigation |
|------|------------|
| Silent numerical drift from refactor | Golden-value regression test on small fixture; spot-check `test/examples/dimer` `(j0, jphi)`. |
| `bias_col` contract ambiguity for future estimators | Document the contract in `solve_coefficients` docstring: `bias_col::Int` is the column index whose regularization weight must be zero; estimators that intercept-fit internally may ignore it. |
| External callers depending on `elastic_net_regression` directly | Grep before removal. If any survive outside `src/Optimize.jl`, keep it as a deprecated shim. |
