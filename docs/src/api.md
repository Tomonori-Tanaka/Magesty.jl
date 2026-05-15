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

```@docs
save
load
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

### Spin configurations
```@docs
SpinConfigs.SpinConfig
SpinConfigs.read_embset
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

### Configuration parsers
```@docs
ConfigParser.Config4System
ConfigParser.Config4Optimize
```

## Utility functions

### Spherical harmonics
```@docs
MySphericalHarmonics.Zₗₘ
MySphericalHarmonics.Zₗₘ_unsafe
MySphericalHarmonics.∂ᵢZlm
MySphericalHarmonics.∂ᵢZlm_unsafe
```

### Rotation matrices
```@docs
RotationMatrix.rotmat2euler
RotationMatrix.Δl
```

### Version information
```@docs
Version.version_string
```
