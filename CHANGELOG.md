# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
