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

### System Creation
```@docs
System(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
System(toml_file::AbstractString; verbosity::Bool = true)
```

### SpinCluster Creation
```@docs
SpinCluster(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
SpinCluster(toml_file::AbstractString; verbosity::Bool = true)
SpinCluster(system::System, input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
SpinCluster(system::System, input_dict::AbstractDict{<:AbstractString, <:Any}, spinconfig_list::AbstractVector{SpinConfig}; verbosity::Bool = true)
```

### Energy and Torque Calculations
```@docs
calc_energy
calc_torque
```

### Results Access
```@docs
get_j0
get_jphi
get_j0_jphi
```

### Output Functions
```@docs
write_xml
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

### BasisSets
```@docs
BasisSets.BasisSet
```

### Optimization
```@docs
Optimize.Optimizer
```

### Spin Configurations
```@docs
SpinConfigs.SpinConfig
SpinConfigs.read_embset
```

## Utility Types

### Atomic Indices
```@docs
AtomicIndices.Indices
AtomicIndices.IndicesUniqueList
AtomicIndices.get_total_L
AtomicIndices.get_atom_l_list
```

### SALCs (Symmetry-Adapted Linear Combinations)
```@docs
SALCs.SALC
```

### Unitary Matrix
```@docs
UnitaryMatrixCl.UniMatCl
UnitaryMatrixCl.getindex_m
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
MySphericalHarmonics.Sₗₘ
MySphericalHarmonics.d_Slm
MySphericalHarmonics.∂ᵢSlm
```

### Rotation Matrices
```@docs
RotationMatrix.rotmat2euler
RotationMatrix.Δl
```

### Energy Calculations
```@docs
CalcEnergy.calc_energy
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
SortedContainer.getcount
SortedContainer.addcount!
```

### Counting Containers
```@docs
CountingContainer.CountingUniqueVector
CountingContainer.getcounts
```

## Constants

```@docs
VERSION
```