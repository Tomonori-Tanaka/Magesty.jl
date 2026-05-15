# API Reference

```@meta
CurrentModule = Magesty
```

The user-facing API is built around four types: `SCEBasis`, `SCEDataset`,
`SCEFit`, and `SCEModel`. The legacy `System` / `SpinCluster` /
`build_sce_basis` / `fit_sce_model` / `write_xml` entry points are still
exported for backwards compatibility and will be removed in a future
release.

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

## Legacy API (to be removed)

The legacy entry points remain exported during the public-API refactor.
New code should prefer the four-type API above.

### Legacy types

```@docs
System
SpinCluster
```

### Legacy construction

```@docs
build_sce_basis
build_sce_basis_from_xml
SpinCluster(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
SpinCluster(system::System, input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
SpinCluster(system::System, input_dict::AbstractDict{<:AbstractString, <:Any}, spinconfig_list::AbstractVector{SpinConfig}; verbosity::Bool = true)
```

### Legacy fitting and output

```@docs
fit_sce_model
write_xml
```

### Legacy non-exported helpers

After `using Magesty`, access these with the `Magesty.` prefix:

```julia
Magesty.calc_energy(sc, spin_config)
Magesty.calc_torque(sc, spin_config)
Magesty.get_j0(sc)
Magesty.get_jphi(sc)
Magesty.get_j0_jphi(sc)
Magesty.write_energies(sc, filename)
Magesty.write_torques(sc, filename)
```

```@docs
Magesty.calc_energy
Magesty.calc_torque
Magesty.get_j0
Magesty.get_jphi
Magesty.get_j0_jphi
Magesty.write_energies
Magesty.write_torques
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
Clusters.cluster_orbits
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

### Energy calculations
```@docs
EnergyTorque.calc_energy
```

### Version information
```@docs
Version.version_string
```
