# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `AdaptiveLasso` estimator implementing the Zou 2006 / ALAMODE-style
  one-shot reweighted L1. Runs a user-supplied
  `pilot::AbstractEstimator` (default `OLS()`) to obtain a pilot
  coefficient vector, then solves a weighted-L1 Lasso with
  `pf[j] = 1 / max(|beta_pilot[j]|, epsilon)^gamma`. Exposed kwargs:
  `pilot`, `lambda` (`>= 0`), `gamma` (`>= 0`, default `1.0`),
  `epsilon` (`> 0`, default `eps(Float64)`), `standardize`
  (default `true`). The docstring recommends `pilot = Ridge(lambda =
  small)` on rank-deficient SCE designs (where OLS's minimum-norm
  null-space noise miscalibrates the adaptive weights). The private
  `_glmnet_solve` helper gained an optional
  `penalty_factor::Union{Nothing, AbstractVector{<:Real}} = nothing`
  keyword; the existing `OLS` / `Ridge` / `ElasticNet` / `Lasso`
  paths are byte-for-byte unchanged (the `nothing` default routes
  through the original glmnet call form). A regression test confirms
  `AdaptiveLasso(gamma = 0)` reduces to plain Lasso to `atol = 1e-6`
  (empirical worst ~3e-8 on the test fixture, ~30x headroom). GLMNet
  internally rescales `penalty_factor` so the supplied weights sum to
  `nvars`, so the user-supplied `lambda` interacts with the rescaled
  vector -- matched-`lambda` comparison against plain Lasso is not
  apples-to-apples once `gamma > 0`; per-side `lambda` tuning is the
  recommended pattern.
- `ElasticNet` estimator and `Lasso` convenience function (returning
  `ElasticNet(alpha = 1.0, ...)`) for sparse / mixed-norm SCE fits.
  `ElasticNet` is backed by GLMNet.jl and exposes `alpha::Float64`
  (`0 <= alpha <= 1`), `lambda::Float64` (`lambda >= 0`), and
  `standardize::Bool` (default `true`, neutralising the per-cluster
  `(4pi)^(N/2)` column scale). The estimator hooks into the existing
  `solve_coefficients(estimator, X, y)` dispatch and calls GLMNet
  with `intercept = false`: the energy block of `X` is already
  mean-centered upstream by `assemble_weighted_problem`, and `j0` is
  recovered downstream by `extract_j0_jphi` -- the same post-processing
  `OLS` and `Ridge` already use. A regression test confirms
  `ElasticNet(alpha = 0, lambda = lambda_a/n, standardize = false)`
  reproduces analytic `Ridge(lambda = lambda_a)` to `rtol ~ 1e-3` over
  `lambda_a in {1e-3, 1e-2, 1e-1}` (GLMNet's coordinate-descent
  objective is `(1/2n)||y - X*beta||^2 + (lambda_g/2)||beta||^2`,
  equivalent to analytic ridge under `lambda_g = lambda_a / n`).
- The `fit(SCEFit, ...)` summary now prints an estimator-independent
  `nonzero coefs : k / num_coefs` line, counting `count(!iszero, jphi)`.
  For `OLS` / `Ridge` this is the full coefficient count (no exact
  zeros); for `Lasso` / `ElasticNet` it reports the active support at
  the chosen `lambda`.
- `SpinConfig` is now exported from `Magesty` as a public type. The
  docstring fully describes its fields, including the
  derived-at-construction `local_magfield_vertical` and `torques`
  observables. Build one directly with `SpinConfig(energy, magmom_size,
  spin_directions, local_magfield)` or read a sequence with
  `read_embset`.
- `predict_energy` and `predict_torque` now accept
  `AbstractVector{<:AbstractMatrix{<:Real}}` (list of spin-direction
  matrices) and `AbstractVector{SpinConfig}` in addition to the
  existing single-config / `SCEDataset` inputs. Both vector forms
  return `Vector{Float64}` / `Vector{Matrix{Float64}}` in input order,
  letting callers skip building an `SCEDataset` for ad-hoc batches.

### Changed

- **Breaking:** `save` and `load` are no longer exported from `Magesty`.
  Call them as `Magesty.save(...)` / `Magesty.load(...)`. This avoids the
  name clash with the generic `save` / `load` exported by JLD2, FileIO,
  CSV.jl, and other packages.
- **Breaking:** The default `torque_weight` of
  `fit(SCEFit, dataset, estimator)` changed from `0.5` to `1.0`
  (torque-only fit). The SCE coefficients `jphi` are best determined by
  torque residuals, and `j0` is recovered in closed form from the
  energy block, so a torque-only default produces better-conditioned
  fits on typical DFT datasets. Set `torque_weight = 0.5` explicitly to
  recover the previous behavior.

### Removed

- **Breaking:** Removed the exported `VERSION` constant and the
  `Magesty.versioninfo()` function. They duplicated functionality that
  Julia already provides — use `pkgversion(Magesty)` for the package
  version and `Base.versioninfo()` for the runtime environment dump.
  Internally, the `Version` submodule was folded into a private helper
  inside `XMLIO.jl` (still records the package version in the XML
  `<version>` element).
- **Breaking:** Removed `Magesty.install_tools()`. The CLI scripts
  under `tools/` are unchanged; invoke them directly with
  `julia --project=@vMAJOR.MINOR /path/to/script.jl` (or a shell
  alias of your own) instead of relying on a packaged installer.

### Internal

- Removed the bias column from `SCEDataset.X_E` / the energy design
  matrix returned by `Fitting.build_design_matrix_energy`. The
  reference energy `j0` is now eliminated analytically before the
  solver runs: `Fitting.assemble_weighted_problem` mean-centers the
  energy block (closed form from `∂L/∂j0 = 0`,
  `j0_star(jphi) = mean(y_E - X_E * jphi)`), the solver returns only
  `jphi`, and `Fitting.extract_j0_jphi` recovers `j0` afterward from
  the unscaled inputs via the same closed form. The public API
  (`fit`, `predict_energy`, `predict_torque`, `coef`, `intercept`,
  the estimator types, `Magesty.save` / `Magesty.load`) is
  signature-unchanged. For `OLS` and `Ridge` the formulation is
  mathematically equivalent to the previous bias-in-`X` setup;
  numerical output agrees up to floating-point rounding noise from
  the different operation sequences, with the regression test pinned
  at `rtol = 1e-9` (measured worst case `2.3e-11`).

## [0.1.0] - 2026-05-17

First public release.

### Added

- Spin-cluster expansion (SCE) pipeline: `SCEBasis`, `SCEDataset`, `SCEFit`,
  `SCEModel`.
- Four input paths for `SCEBasis`: TOML file, parsed dict, AtomsBase system,
  keyword arguments.
- Fitting estimators: `OLS`, `Ridge`.
- Prediction API: `predict_energy`, `predict_torque`.
- StatsAPI-compatible evaluation: `coef`, `intercept`, `nobs`, `dof`,
  `r2_energy`, `r2_torque`, `rmse_energy`, `rmse_torque`, `rss_*`,
  `residuals_*`.
- XML persistence: `save` / `load` for `SCEBasis` and `SCEModel`.
- EMBSET reader: `read_embset`.
- Tesseral real spherical harmonics with safe / unsafe variants.
- Symmetry-adapted linear combination (SALC) basis construction with
  isotropy option.
- Examples: `examples/01_basic_flow.jl`, `02_cif_input.jl`,
  `03_save_load.jl`.
- Documentation hosted at
  <https://Tomonori-Tanaka.github.io/Magesty.jl/>.

[Unreleased]: https://github.com/Tomonori-Tanaka/Magesty.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Tomonori-Tanaka/Magesty.jl/releases/tag/v0.1.0
