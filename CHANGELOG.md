# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
