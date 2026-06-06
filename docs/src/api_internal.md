# Internal API

```@meta
CurrentModule = Magesty
```

!!! warning "Stability"
    The names on this page are **internal**: they are not exported and may
    be renamed, removed, or have their signatures changed in any release
    without a deprecation cycle. They are documented here for contributors
    and for advanced users who need access to lower-level building blocks
    (e.g. inspecting a `SALCBasis` after construction). For the stable,
    user-facing surface see the [API Reference](@ref).

---

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

---

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

---

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

### MFA spin sampling

Code-agnostic building blocks behind the exported [`sample_mfa_incar`](@ref);
see [Mean-Field Sampling](tips/mfa_sampling.md) for the theory.

```@docs
MfaSampling.thermal_averaged_m
MfaSampling.tau_from_magnetization
MfaSampling.sample_vmf_direction
MfaSampling.mfa_sample
MfaSampling.mfa_sweep
MfaSampling.parse_atom_index_spec
```
