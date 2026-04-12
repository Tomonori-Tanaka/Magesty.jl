# API Reference

```@meta
CurrentModule = Magesty
```

## Main Types

### System
```@docs
System
```

### SpinCluster
```@docs
SpinCluster
```

## Main Functions

### System Building
```@docs
build_sce_basis
build_sce_basis_from_xml
```

### SpinCluster Creation
```@docs
SpinCluster(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
SpinCluster(system::System, input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
SpinCluster(system::System, input_dict::AbstractDict{<:AbstractString, <:Any}, spinconfig_list::AbstractVector{SpinConfig}; verbosity::Bool = true)
```

### Model Fitting
```@docs
fit_sce_model
```

### Estimators
```@docs
AbstractEstimator
OLS
ElasticNet
```

### Output Functions
```@docs
write_xml
```

### Data Loading
```@docs
read_embset
```

## Non-Exported Convenience Functions

The following functions are defined in the `Magesty` module but are **not exported**.
After `using Magesty`, access them with the `Magesty.` prefix:

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

### BasisSets
```@docs
BasisSets.BasisSet
```

### Spin Configurations
```@docs
SpinConfigs.SpinConfig
SpinConfigs.read_embset
```

## Utility Types

### Spherical Harmonics Transforms
```@docs
SphericalHarmonicsTransforms.c2r_sph_harm_matrix
SphericalHarmonicsTransforms.r2c_sph_harm_matrix
```

### Atom Cells
```@docs
AtomCells.AtomCell
```

### Configuration Parsers
```@docs
ConfigParser.Config4System
ConfigParser.Config4Optimize
```

## Utility Functions

### Spherical Harmonics
```@docs
MySphericalHarmonics.Zₗₘ
MySphericalHarmonics.Zₗₘ_unsafe
MySphericalHarmonics.∂ᵢZlm
MySphericalHarmonics.∂ᵢZlm_unsafe
```

### Rotation Matrices
```@docs
RotationMatrix.rotmat2euler
RotationMatrix.Δl
```

### Energy Calculations
```@docs
EnergyTorque.calc_energy
```

### Version Information
```@docs
Version.version_string
```

## Container Types

### Sorted Containers
```@docs
SortedContainer.SortedVector
SortedContainer.SortedUniqueVector
SortedContainer.SortedCountingUniqueVector
```

### Counting Containers
```@docs
CountingContainer.CountingUniqueVector
CountingContainer.getcounts
```

