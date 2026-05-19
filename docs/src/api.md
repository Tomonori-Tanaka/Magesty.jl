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
