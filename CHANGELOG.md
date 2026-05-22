# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `refit(fit::SCEFit, estimator = OLS(); threshold = 0.0, verbosity = true)`: post-selection
  refit on the basis support of an existing fit. A basis is kept iff
  `abs(coef(f)[j]) * norm(X[:, j]) > threshold`, where `X` is the
  weighted, energy-centered design matrix that `f` was fitted against
  (`f.dataset` and `f.torque_weight` are reused). The column-norm factor
  cancels the per-cluster `(4π)^(N/2)` scaling so the scaled criterion
  is comparable across cluster sizes, and the default `threshold = 0.0`
  keeps an L1 fit's exact-zero support intact while `AdaptiveRidge`
  users pass a positive threshold. The support is re-solved with the
  caller-chosen `estimator` (default `OLS()` for classic post-selection
  debiasing); dropped bases stay at `0.0`, so `jphi` length and SALC
  ordering are preserved. `PrecomputedPilot`-backed estimators (and
  `AdaptiveLasso` whose pilot is one, including the
  `AdaptiveLasso(::SCEFit)` form) are rejected upfront because their
  fixed pilot vector has the original column count, not the refit
  support length. Purely additive — no impact on `fit` or any existing
  estimator.
- `AdaptiveRidge`: iterative Adaptive Ridge estimator (Frommlet & Nuel
  2016) that approximates an L0-penalized fit. It refits a per-coefficient
  weighted ridge problem and updates the weights
  `w_j = 1 / (beta_j^2 + epsilon)` between iterations, driving small
  coefficients toward zero. Each weighted ridge subproblem is solved by
  the analytic closed form `(X'X + lambda * Diagonal(w)) \ (X'y)` — the
  same analytic family as `Ridge`, no GLMNet call — so there is no
  `standardize` keyword. Keyword constructor:
  `AdaptiveRidge(; lambda, epsilon = 1e-8, max_iter = 50, tol = 1e-6)`;
  the iteration stops on a relative infinity-norm coefficient change
  below `tol`. Purely additive, no impact on existing estimators.
- `SCEDataset(model::SCEModel, ...)` / `SCEDataset(fit::SCEFit, ...)`
  constructors that accept a fitted `SCEModel` or `SCEFit` in place of an
  `SCEBasis`, reusing the basis embedded in it (`model.basis` /
  `fit.dataset.basis`) without rebuilding the SALCs. The second argument
  is the usual `Vector{SpinConfig}` or EMBSET file path; the methods
  forward to the basis-based constructors, so design matrices are built
  exactly as before (purely additive, no numerical-convention impact).
- `magesty` command-line interface, provided by the `MagestyCLI` package
  in the `cli/` subdirectory of the repository and built with Comonicon.jl.
  Adding `MagestyCLI` and running `Pkg.build("MagestyCLI")` writes the
  launcher into `~/.julia/bin`. Keeping the CLI in its own package leaves
  the core `Magesty` package free of the Comonicon dependency, so static
  analysis (JET) covers the whole core again. Subcommands:
  `magesty vasp extxyz` (convert a VASP run to extended XYZ),
  `magesty vasp toml` (convert a VASP POSCAR to a Magesty input TOML
  configuration), `magesty vasp embset` (convert VASP OUTCAR files to the
  EMBSET training-data format), and `magesty version`.
- `vasp_to_extxyz`: exported function that converts a VASP run
  (`vasprun.xml`, optionally `OSZICAR`) to extended XYZ and returns the
  extxyz text; the `magesty vasp extxyz` subcommand is a thin wrapper over
  it. The VASP parsing and extxyz writing previously living under `tools/`
  are now package submodules (`VaspIO`, `ExtXYZ`).
- `poscar_to_toml`: exported function that converts a VASP POSCAR
  structure file to a Magesty input TOML configuration and returns the
  TOML text; the `magesty vasp toml` subcommand is a thin wrapper over it.
  The generated configuration is a starting point with placeholder
  interaction settings (`lmax = 0`, `cutoff = -1`).
- `outcar_to_embset`: exported function that extracts the energy, per-atom
  magnetic moments, and per-atom constraining field from one or more VASP
  `OUTCAR` files and returns the result as EMBSET training-data text; the
  `magesty vasp embset` subcommand is a thin wrapper over it. Magnetic
  moments and fields are rotated by the `saxis` quantization-axis rotation.
- `write_energies` / `write_torques`: exported plain-text writers that dump
  observed (DFT) versus predicted (SCE) energies and per-atom torques to
  whitespace-separated files, consumed by the `FitCheck_energy.py` /
  `FitCheck_torque.py` visualization scripts under `tools/`. Each writer has
  a self-contained `write_*(f::SCEFit)` form (evaluates the fit's own
  training dataset) and a general `write_*(predictor, data)` form where
  `predictor` is `SCEModel` or `SCEFit` and `data` is an `SCEDataset`, a
  `Vector{SpinConfig}`, or an EMBSET file path -- so a held-out validation
  or test set can be checked too. The output format matches the legacy
  writers; predictions go through `predict_energy` / `predict_torque`
  unchanged (pure I/O, no numerical-convention impact).
- `PrecomputedPilot <: AbstractEstimator` adapter that returns a fixed
  coefficient vector from `solve_coefficients` (after a length check
  against `size(X, 2)`), letting `AdaptiveLasso.pilot` reuse a prior
  fit's coefficients instead of running a fresh pilot regression. The
  input vector is copied at construction (enforced by the inner
  constructor) so caller-side mutation cannot leak into the stored
  pilot. Designed as the primitive an iterative Adaptive variant will
  call inside its inner loop.
- `AdaptiveLasso(f::SCEFit; kwargs...)` and
  `AdaptiveLasso(model::SCEModel; kwargs...)` convenience constructors
  that wrap `coef(f)` / `coef(model)` in `PrecomputedPilot` and
  forward the remaining kwargs (`lambda`, `gamma`, `epsilon`,
  `standardize`) to the keyword constructor. The prior fit must have
  been produced on the same `SCEBasis`; only a length check is
  enforced (same-length, different-SALC-ordering mismatch is not
  detectable). Passing a `pilot` keyword alongside the positional
  `SCEFit` / `SCEModel` raises `ArgumentError` -- otherwise Julia's
  kwarg-splat semantics would let the splatted `pilot` silently
  override the precomputed one, defeating the constructor's purpose.
  Placed in `src/Magesty.jl` after `coef(::SCEFit)` / `coef(::SCEModel)`
  because `SCEFit` / `SCEModel` are not visible inside the `Fitting`
  submodule (defined in the parent module after `include("Fitting.jl")`);
  an explicit `import .Fitting: AdaptiveLasso` silences the Julia 1.12
  "constructor extended without explicit qualification" warning.
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
- `SpinConfig` inner constructor now validates that every column of
  `spin_directions` has unit norm, enforcing the documented unit-vector
  invariant that `Zₗₘ` semantics rely on. Tolerance is controlled by
  the new keyword `atol_unit_norm::Float64` (default `1e-6`). Columns
  containing `NaN` are rejected explicitly. Existing in-memory
  fixtures that already supplied unit-norm directions are unaffected.
- In-memory basis-identity check: `SCEBasis` gains a fifth field
  `salc_fingerprint::UInt64` derived from the structural identifiers
  of every SALC (`ls`, `Lf`, `Lseq`, `atoms`, `multiplicity`) in
  key-group order. `_check_basis` now short-circuits on either `===`
  identity or matching fingerprint, so a model reloaded from XML can
  be applied to a dataset built before save without a spurious
  `ArgumentError`. The XML schema is unchanged; the fingerprint is
  recomputed from the reconstructed SALC ordering on every load.
- `Magesty.save(f::SCEFit, path)` convenience that delegates to
  `Magesty.save(SCEModel(f), path)`. The fit-time dataset and
  estimator are not persisted -- saving an `SCEFit` and saving its
  `SCEModel(f)` produce byte-identical files.

### Changed

- `vcat` of `SCEDataset` objects now accepts any datasets whose
  `SCEBasis` objects share a `salc_fingerprint`, not only datasets built
  from the *same* `SCEBasis` instance. Two bases constructed from the
  same input — or reloaded from the same XML — are distinct objects with
  an identical fingerprint, so their datasets can now be concatenated.
  This matches the prediction path (`predict_energy` / `predict_torque`),
  which already accepted fingerprint-equal bases. Datasets from
  structurally different bases (mismatched fingerprint) still raise
  `ArgumentError`. Non-breaking: every `vcat` call that previously
  succeeded still succeeds.
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
- **Breaking:** Removed `Magesty.install_tools()`, superseded by the
  Comonicon-based `magesty` command-line interface (see Added).
- The standalone `tools/vasp/vasp2extxyz.jl` script — its single-directory
  VASP-to-extxyz conversion is now the `magesty vasp extxyz` command (and
  the `vasp_to_extxyz` API function).
- The standalone `tools/vasp/pos2toml.jl` script — its POSCAR-to-TOML
  conversion is now the `magesty vasp toml` command (and the
  `poscar_to_toml` API function).
- The standalone `tools/extract.jl` script — its OUTCAR-to-EMBSET
  extraction is now the `magesty vasp embset` command (and the
  `outcar_to_embset` API function). The old script's `--randomize` flag
  is not carried forward; shuffle the input OUTCAR list before invoking
  the command if a randomized configuration order is needed.
- The `tools/vasp/vasp2extxyz_recursive.jl` script (batch directory-tree
  conversion), removed as unused.
- The `make test-tools` target: the converter tests it ran are now package
  component tests covered by `make test-unit` / `make test-all`.

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
