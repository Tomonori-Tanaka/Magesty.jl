# API Reference

```@meta
CurrentModule = Magesty
```

The user-facing API is built around four types: `SCEBasis`, `SCEDataset`,
`SCEFit`, and `SCEModel`.

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
```

## Data loading

```@docs
read_embset
SpinConfig
```

## CLI tools

```@docs
install_tools
```

## Submodules

### Structures
```@docs
Structures.Structure
```

### Symmetries
```@docs
Symmetries.Symmetry
Symmetries.SymmetryOperation
Symmetries.Maps
```

### Clusters
```@docs
Clusters.Cluster
```

### SALCBases
```@docs
SALCBases.SALCBasis
```

## Utility types

### Spherical harmonics transforms
```@docs
SphericalHarmonicsTransforms.c2r_sph_harm_matrix
SphericalHarmonicsTransforms.r2c_sph_harm_matrix
```

### Atom cells
```@docs
AtomCells.AtomCell
```

### Input specs
```@docs
InputSpecs.SystemSpec
InputSpecs.InteractionSpec
InputSpecs.SymmetryOptions
InputSpecs.parse_toml_inputs
```

## Utility functions

### Spherical harmonics
```@docs
TesseralHarmonics.Zₗₘ
TesseralHarmonics.Zₗₘ_unsafe
TesseralHarmonics.∂ᵢZlm
TesseralHarmonics.∂ᵢZlm_unsafe
```

### Rotation matrices
```@docs
RotationMatrix.rotmat2euler
RotationMatrix.Δl
```

### Version information

Use the standard Julia idioms:

- `pkgversion(Magesty)` returns the package version as a `VersionNumber`.
- `Base.versioninfo()` dumps the active Julia / platform / threading
  context.

Magesty does not provide its own `VERSION` constant or `versioninfo`
function — the standard library already covers both needs.
