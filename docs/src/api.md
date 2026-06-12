# API Reference

```@meta
CurrentModule = Magesty
```

The user-facing API is built around four types: `SCEBasis`, `SCEDataset`,
`SCEFit`, and `SCEModel`. Names listed below are stable: removals and
signature-incompatible changes follow a deprecation cycle. For lower-level
building blocks used internally by the SCE pipeline (and not covered by
this stability guarantee), see the [Internal API](@ref).

## Main types

```@docs
SCEBasis
SCEDataset
SCEFit
SCEModel
```

## Fitting

```@docs
fit
refit
coef
intercept
nobs
dof
```

## Prediction

```@docs
predict_energy
predict_torque
```

## Evaluation

```@docs
r2_energy
r2_torque
rmse_energy
rmse_torque
rss_energy
rss_torque
residuals_energy
residuals_torque
```

## Fit-quality output

`write_energies` and `write_torques` dump observed (DFT) versus predicted
(SCE) values to whitespace-separated text files. The files are consumed by
the `FitCheck_energy.py` / `FitCheck_torque.py` visualization scripts under
`tools/`.

```@docs
write_energies
write_torques
```

## Cross-validation diagnostics

Generalized cross-validation (GCV) estimates out-of-sample prediction error
from a single fit, on the same combined energy+torque weighted objective that
`fit` minimizes. `gcv` returns the score for a fit and `gcv_r2` returns the
companion predictive R² (`1 − GCV/MSY`), which reads on a fixed scale (`1`
perfect, `0` matches the null model) where the raw GCV score does not;
`gcv_lambda` sweeps the
ridge penalty (one SVD serves the whole path) and reports the GCV minimizer;
`gcv_learning_curve` sweeps the training-set size with random subsets to check data
sufficiency. GCV is defined only for the linear estimators `OLS`, `Ridge`, and
`AdaptiveRidge`. `write_gcv_lambda` / `write_gcv_learning_curve` write the sweep
results to text for the `FitCheck_gcv_lambda.py` / `FitCheck_gcv_learning_curve.py`
scripts under `tools/`. See [Cross-validation diagnostics](@ref) in the theory
notes for the formula and conventions.

```@docs
gcv
gcv_r2
gcv_lambda
gcv_learning_curve
GCVLambdaPath
GCVSizeCurve
write_gcv_lambda
write_gcv_learning_curve
```

## VASP conversion

The `vasp_to_extxyz`, `poscar_to_toml`, and `oszicar_to_embset` functions
convert VASP output to, respectively, extended XYZ, a Magesty input TOML
configuration, and the EMBSET training-data format. Each is also available
from the command line under `magesty vasp` (see [Installation](@ref) and
[Tools](tools.md)).

```@docs
vasp_to_extxyz
poscar_to_toml
oszicar_to_embset
```

## Spin sampling

`sample_mfa_incar` draws thermally conditioned spin configurations from a VASP
INCAR with the Mean-Field Approximation (von Mises-Fisher direction sampling)
and writes one INCAR per configuration. It is also available from the command
line as `magesty vasp mfa` (see [Tools](tools.md)).

```@docs
sample_mfa_incar
```

## Sunny.jl export

`sce_to_sunny` turns a fitted `SCEModel` into a runnable
[Sunny.jl](https://github.com/SunnySuite/Sunny.jl) script that computes a
linear spin-wave-theory magnon dispersion. It is also available from the
command line as `magesty sunny script` (see [Tools](tools.md)). Magesty itself
gains no Sunny dependency — the function only emits text.

```@docs
sce_to_sunny
```

## Persistence

`save` and `load` are not exported — call them as `Magesty.save` /
`Magesty.load` to avoid clashing with the generic `save` / `load` exported
by JLD2, FileIO, CSV.jl, and others.

```@docs
Magesty.save
Magesty.load
```

## Estimators

```@docs
AbstractEstimator
OLS
Ridge
ElasticNet
Lasso
AdaptiveLasso
PrecomputedPilot
AdaptiveRidge
```

## Data loading

```@docs
read_embset
SpinConfig
```

## Version information

Use the standard Julia idioms:

- `pkgversion(Magesty)` returns the package version as a `VersionNumber`.
- `Base.versioninfo()` dumps the active Julia / platform / threading
  context.

Magesty does not provide its own `VERSION` constant or `versioninfo`
function — the standard library already covers both needs.
