# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- `r2_energy` / `r2_torque` now return `NaN` on a degenerate evaluation set
  (all observed values equal, e.g. a single configuration), matching the GCV
  convention, instead of an inconsistent `±Inf`/`NaN`.
- The raw-matrix `predict_energy` / `predict_torque` entry points validate
  that every spin-direction column is a unit vector (tolerance `1e-6`, same
  as the `SpinConfig` constructor) and throw an `ArgumentError` on non-unit
  or non-finite columns instead of silently returning wrong physics.

### Internal

- `SALCBasis` construction: the projection-matrix stage no longer recomputes
  the time-reversal-independent work (atom shift, translation fold, atom
  reorder, tensor inner products) twice per symmetry operation. Projection
  stage ~1.9× faster with −47% allocations on the FeGe 2×2×2 fixture
  (constructor overall ~1.2× faster, −29% allocations); output verified
  bit-identical on the development machine.

## [0.2.0] - 2026-06-25

### Added

- `sce_to_sunny` gains a `scaling` keyword (`:auto` / `:moment` / `:coupling`; CLI
  `--scaling`) that lets the exporter handle **itinerant / non-half-integer
  moments** (e.g. Fe `2.2 μB` ⇒ `S_eff = 1.1`). The `:coupling` route keeps Sunny's
  `Moment` at a placeholder `s₀ = 1` and folds the physical `S_eff` into the
  couplings (`J = M/(s₀·√(S_i S_j))`, single-ion `1/(s₀ S_i)`), so the magnon
  dispersion is physical for any positive real `S_eff`. Only the dispersion is
  preserved — the static `energy(sys)` is then no longer the SCE energy. Exact for a
  uniform `S_eff`; non-uniform `S_eff` keeps the off-diagonal exchange exact and
  warns that the on-site (Larmor) term is approximate. `:auto` uses `:moment` for
  half-integer spins and `:coupling` otherwise.

### Changed

- **BREAKING CHANGE:** `sce_to_sunny` now requires a `spin` keyword — the physical
  effective spin length `S_eff = m/(g μ_B)` of each magnetic species (a scalar, or
  a `species => value` Dict) — and accepts optional `g`, `mode`, and `scaling`. The
  SCE couplings absorb the spin magnitude (`J_SCE = J_phys·S²`), so the LSWT magnon
  dispersion needs the physical spin; the previous fixed `s = 1` inflated magnon
  frequencies by a factor `~S` (e.g. ~2.5× for MnTe, `S = 5/2`). With the default
  `scaling = :moment` (selected for half-integer spins), bilinear bonds are rescaled
  by `1/(s_i s_j)` and single-ion terms by a mode-dependent factor at emission, so
  `energy(sys)` is unchanged while the dispersion scales physically; non-half-integer
  spins use the new `:coupling` route (see *Added*). The `magesty sunny script` CLI
  gains required `--spin` and optional `--g` / `--mode` / `--scaling`.

## [0.1.1] - 2026-06-19

### Documentation

- Installation instructions now lead with `Pkg.add("Magesty")` from the
  General registry, and document installing the `magesty` command-line tool
  without a repository checkout via
  `Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl", subdir="cli")`.
- Citation metadata (`CITATION.cff`) now references the published article,
  Phys. Rev. Research 8, 023300 (2026), rather than the arXiv preprint.

## [0.1.0] - 2026-06-15

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
- **GCV diagnostics (spec 260610-gcv-diagnostics).** New exported
  generalized cross-validation diagnostics for SCE fitting, computed on the
  combined energy+torque weighted objective the fit minimizes via the hat
  matrix `H` (`GCV = (‖r‖²/N) / (1 − tr(H)/N)²`). `gcv(f::SCEFit)` returns the
  score for a fit and `gcv_r2(f::SCEFit)` returns the companion predictive R²
  (`1 − GCV/MSY`, with `MSY = ‖y‖²/N` the null-model mean square), which reads on
  a fixed scale (`1` perfect, `0` matches the null model) where the raw GCV score
  does not; `gcv_lambda(dataset, lambdas; torque_weight)` sweeps the
  ridge penalty from a single SVD and reports the GCV minimizer
  (`GCVLambdaPath`); `gcv_learning_curve(dataset, estimator; sizes, repeats, seed,
  torque_weight)` sweeps the training-set size with random subsets to check
  data sufficiency (`GCVSizeCurve`). Defined for the linear estimators `OLS`,
  `Ridge`, and `AdaptiveRidge` (the last with a conditional effective-dof);
  non-linear estimators raise `ArgumentError`. The predictive R² is also carried
  per-point on the sweep results (`GCVLambdaPath.gcv_r2`,
  `GCVSizeCurve.gcv_r2_mean` / `gcv_r2_std`). `write_gcv_lambda` /
  `write_gcv_learning_curve` write the sweeps to text (including the R² columns)
  for the new `tools/FitCheck_gcv_lambda.py` /
  `tools/FitCheck_gcv_learning_curve.py` plotters (`--r2` plots the predictive
  R²). No change to existing fit results.
- **MFA spin sampling (spec 260605-mfa-sampling-cli).** New exported
  `sample_mfa_incar(incar_path; variable, start, stop, num_points, num_samples,
  randomize, fix, uniform_atoms, outdir, prefix)` draws thermally conditioned
  spin configurations from a VASP INCAR with the Mean-Field Approximation and
  writes one INCAR per configuration, also available from the command line as
  `magesty vasp mfa`. Per-atom directions are sampled from a von Mises-Fisher
  distribution whose concentration `κ = 3m/τ` follows from the MFA
  self-consistency equation `m = coth(3m/τ) − τ/3m`; magnitudes are preserved
  and both `MAGMOM` and `M_CONSTR` are written. The control variable is `tau`
  (reduced temperature `T/Tc`) or `m` (magnetization); the sweep is specified by
  point count (`num_points`, evenly spaced) rather than step width. The
  code-agnostic sampler (`Magesty.MfaSampling`) and VASP INCAR reader/writer
  (`Magesty.IncarIO`) are new internal modules; the von Mises-Fisher draw uses an
  exact closed-form p=3 sampler, so the only new dependency is the lightweight
  `Roots`. This promotes the former `tools/sampling_mfa.jl` script.
- **Sunny.jl export (spec 260529-sce-sunny-export).** New exported
  `sce_to_sunny(model; output, placement, symprec)` turns a fitted
  `SCEModel` into a runnable [Sunny.jl](https://github.com/SunnySuite/Sunny.jl)
  linear-spin-wave-theory script, also available from the command line as
  `magesty sunny script`. Two-site `l₁ = l₂ = 1` SALCs map to 3×3 bilinear
  exchange (Heisenberg + Dzyaloshinskii–Moriya + anisotropic symmetric) and
  single-site `l = 2` SALCs to single-ion anisotropy; higher-order terms are
  skipped with a warning. `placement = :auto` emits the chemical primitive cell
  (unfolded dispersion) when the model is cleanly unfoldable (interaction range
  below half the supercell) and otherwise falls back to the exact training
  supercell (folded dispersion). Spins use the reduced `s = 1`, `g = 2`
  convention. Magesty gains no Sunny dependency (text generation only); a
  Sunny-backed round-trip energy check lives in `test/sunny/` and runs via the
  new `make test-sunny` target.

### Changed

- **Internal `src/` cleanup (spec 260601-src-refactor).** Behavior-preserving
  refactor and hot-path tuning surfaced by a review of `src/`. The
  angular-momentum coupling cache is now reached through a
  `CoupledBases.cached_coupling_results` accessor instead of `SALCBases`
  indexing the cache dict directly. The torque spherical-harmonics cache builds
  `Z` and `∇Z` from a single Legendre recursion (new `Zₗₘ_grad_unsafe`,
  ~2.4× faster per-atom fill), and `projection_matrix_coupled_basis` plus the
  SALC translational-equivalence check shed avoidable allocations. Results are
  bit-identical to the previous implementation.
- **Solver unification (spec 260526-solver-unification-and-memory-log).**
  `OLS`, `Ridge`, and `AdaptiveRidge` now all route through
  `cholesky(Symmetric(X'X + Λ)) \ (X'y)`, replacing the previous
  non-pivoted-QR (OLS) and `Symmetric \ ` Bunch-Kaufman (Ridge,
  AdaptiveRidge) paths. About **3× faster** on representative sizes
  (e.g. 19300 × 146: 32 ms → 10 ms; 1470 × 31: 260 μs → 68 μs) with
  numerical agreement to ~1e-16. **Behavioral change for OLS on
  rank-deficient designs**: the previous QR + pivoted-QR fallback
  silently returned a min-norm solution; the new Cholesky path catches
  `PosDefException` and rethrows it as an `ArgumentError` whose
  message directs the user to `Ridge(lambda = ε)`. This is by design:
  an OLS fit on a rank-deficient design is not physically identifiable.
  The `MultivariateStats` dependency is dropped from `Project.toml`
  (it was only used for `MultivariateStats.ridge`, now replaced by
  the explicit Cholesky path).
- **Design-matrix memory log.** `build_design_matrix_energy` and
  `build_design_matrix_torque` print a pre-construction memory
  estimate (matrix bytes + Cholesky Gram bytes) and a post-construction
  `Base.summarysize` line under the existing `verbosity = true` gate.
  Helps size the problem before waiting through the build.
- Design-matrix construction restructured for performance (spec
  260524-design-matrix-restructuring). `build_design_matrix_energy`
  and `build_design_matrix_torque` are about **12× faster** on the
  FeGe B20 2×2×2 light fixture (energy 1.218 s → 0.101 s, torque
  9.833 s → 0.833 s; 4 threads, 5-trial median). The improvements
  come from four staged changes: (A) folding the SALC `M_f`
  coefficient into the coupled tensor at SALC-build time
  (`CoupledBasis_with_coefficient` gains a precomputed
  `folded_tensor` field); (C) pre-enumerating each orbit's
  translation images at SALC-build time (new `clusters` field);
  (B) caching per-spinconfig tesseral spherical harmonics in an
  internal `SHCache` so the hot path reads cached values instead
  of recomputing per call; (D) rewriting the torque kernel as a
  cluster-major reverse-mode Jacobian that distributes the per-
  cluster gradient to all `N` sites in one pass. No public-API
  signature changes. Internal `calc_∇ₑu!` / `calc_∇ₑu` and the
  `EnergyWorkspace` / `GradWorkspace` scratch types are removed
  (they are replaced by `SHCache` and the cluster-major
  accumulator). XML round-trip preserves bit-identity: the new
  `folded_tensor` / `clusters` fields are recomputed
  deterministically on load. See the new theory pages
  `docs/src/theory/{folded_tensor,orbit_clusters,sh_cache,cluster_major_torque}.md`
  for the derivations.

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
  configuration), `magesty vasp embset` (convert VASP OSZICAR files to the
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
- `oszicar_to_embset`: exported function that extracts the energy, per-atom
  magnetic moments, and per-atom constraining field from one or more VASP
  `OSZICAR` files and returns the result as EMBSET training-data text; the
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
- The standalone `tools/extract.jl` script — its OSZICAR-to-EMBSET
  extraction is now the `magesty vasp embset` command (and the
  `oszicar_to_embset` API function). The old script's `--randomize` flag
  is not carried forward; shuffle the input OSZICAR list before invoking
  the command if a randomized configuration order is needed.
- The `tools/vasp/vasp2extxyz_recursive.jl` script (batch directory-tree
  conversion), removed as unused.
- The `make test-tools` target: the converter tests it ran are now package
  component tests covered by `make test-unit` / `make test-all`.

### Internal

- Split the oversized main module file `src/Magesty.jl` into two new files
  that are `include`d into the same `Magesty` module namespace (the idiom
  already used by the trailing includes such as `FitCheckIO.jl`), not Julia
  submodules. `src/Evaluation.jl` now holds the prediction verbs
  (`predict_energy` / `predict_torque`) and the accuracy metrics
  (`r2_*` / `rss_*` / `residuals_*` / `rmse_*`) with their shared
  `_eval_*` helpers; `src/GCV.jl` holds the generalized-cross-validation
  diagnostics (`gcv` / `gcv_r2` / `gcv_lambda` / `gcv_learning_curve` and the
  `GCVLambdaPath` / `GCVSizeCurve` result types). This is a verbatim code
  move: every public binding keeps its module, signature, and docstring, so
  `using Magesty`, `Magesty.save` / `Magesty.load`, and all `@ref` docstring
  links are unaffected. The full test suite (including JET and Aqua) passes
  unchanged; no numerical output changes.
- `Cluster` construction (`src/Clusters.jl`) is significantly faster on
  three-body systems. `irreducible_clusters` switches from an
  O(N_clusters^2) linear scan against accepted representatives to an
  O(N_clusters) lookup against a `Set` keyed by the lex-minimum
  translation image of each cluster (`_translation_canonical_form`). In
  the same change, `set_mindist_pairs` is computed once in the `Cluster`
  constructor and threaded into `generate_clusters` instead of being
  recomputed inside it. On three-body benchmark fixtures the combined
  effect lifts `make bench-cluster` runtime from minutes/hours down to
  sub-second: `fege_2x2x2_3body_fefe_open` drops from 0.50 s to 0.065 s
  (7.7x), and `fege_2x2x2_3body_all_open` — a configuration that used to
  take roughly 4000 s to build — now completes in 0.34 s. Output data
  structures (`cluster_dict`, `irreducible_cluster_dict`,
  `cluster_orbits_dict`, `min_distance_pairs`) are bit-identical; the
  full test suite including XML save/load round-trip and all integration
  fixtures passes unchanged.
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

[Unreleased]: https://github.com/Tomonori-Tanaka/Magesty.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Tomonori-Tanaka/Magesty.jl/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/Tomonori-Tanaka/Magesty.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/Tomonori-Tanaka/Magesty.jl/releases/tag/v0.1.0
