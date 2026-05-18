# Magesty.jl — specification and architecture

Julia package for building effective spin models of magnetic materials.
Constructs general spin models from noncollinear spin DFT calculations
using the spin-cluster expansion (SCE).

## Main modules (`src/`)

All modules live directly under `src/` (flat layout). The `include` order
in `src/Magesty.jl` reflects the dependency order.

### SCE pipeline

| Module | File | Role |
|---|---|---|
| `Structures` | `src/Structures.jl` | Crystal-structure loading and supercell generation |
| `Symmetries` | `src/Symmetries.jl` | Space-group symmetry operations (Spglib wrapper) |
| `Clusters` | `src/Clusters.jl` | Cluster expansion and distance-matrix computation |
| `SALCBases` | `src/SALCBases.jl` | SALC basis construction (computationally expensive) |
| `Fitting` | `src/Fitting.jl` | Design-matrix construction and regression (OLS / Ridge) |
| `SpinConfigs` | `src/SpinConfigs.jl` | Spin-configuration loading and management |

### Utilities and data types

| Module / file | Role |
|---|---|
| `SphericalHarmonicsTransforms.jl` | Real spherical-harmonics transforms |
| `TesseralHarmonics.jl` | Tesseral spherical harmonics `Zₗₘ` (with unsafe variants) |
| `AngularMomentumCoupling.jl` | Wigner coefficients and angular-momentum coupling |
| `CoupledBases.jl` | Intermediate representation for angular-momentum-coupled bases |
| `InputSpecs.jl` | Typed values (`SystemSpec` / `InteractionSpec` / `SymmetryOptions`) and TOML / Dict parser (`parse_toml_inputs`, including wildcard species expansion) |
| `RotationMatrix.jl` | Rotation matrices |
| `XMLIO.jl` | XML I/O for bases and models |
| `AtomsBaseAdapter.jl` | AtomsBase / Unitful boundary adapter (`system_to_specs` / `kwargs_to_specs`) |
| `AtomCells.jl` | Lightweight type holding atom sites plus unit-cell info |
| `SortedCounters.jl` | Internal counter that iterates by sorted key |

## Main types (current API)

```
# Input typed-value gate (src/InputSpecs.jl). All four input paths
# (TOML file / Dict / AtomsBase / kwargs) assemble this 3-tuple and
# pass it to the SCEBasis inner constructor.

SystemSpec
├── name::String
├── num_atoms::Int
├── kd_name::Vector{String}           # Species names (order indexes kd_int_list)
├── kd_int_list::Vector{Int}          # Per-atom species index
├── lattice_vectors::Matrix{Float64}  # 3x3 (Å)
├── x_fractional::Matrix{Float64}     # 3 x num_atoms
└── is_periodic::Vector{Bool}         # length 3

InteractionSpec
├── nbody::Int
├── body1_lmax::Vector{Int}                  # Length = number of species
├── bodyn_lsum::OffsetArray{Int, 1}          # Axis 2:nbody
└── bodyn_cutoff::OffsetArray{Float64, 3}    # (2:nbody, 1:nkd, 1:nkd), symmetric

SymmetryOptions
├── tolerance_sym::Float64
└── isotropy::Bool

# Downstream SCE pipeline

SCEBasis
├── structure::Structure    # Crystal structure (unit cell + supercell)
├── symmetry::Symmetry      # Symmetry operations and translational symmetry
├── salcbasis::SALCBasis    # SALC basis list
└── isotropy::Bool          # Whether the basis was built with the isotropy restriction

SCEDataset
├── basis::SCEBasis
├── spinconfigs::Vector{SpinConfig}
├── X_E::Matrix{Float64}    # Energy design matrix (unweighted; column 1 is the bias)
├── X_T::Matrix{Float64}    # Torque design matrix (unweighted; no bias column)
├── y_E::Vector{Float64}    # Observed energies
└── y_T::Vector{Float64}    # Observed torques (flattened)

SCEFit <: StatsAPI.RegressionModel
├── dataset::SCEDataset
├── j0::Float64             # Estimated reference energy
├── jphi::Vector{Float64}   # Estimated SCE coefficients
├── estimator::AbstractEstimator
├── torque_weight::Float64
├── residuals::Vector{Float64}
└── metrics::Dict{Symbol, Any}  # In-sample RMSE / R²

SCEModel
├── basis::SCEBasis
├── j0::Float64
└── jphi::Vector{Float64}
```

`SCEModel(f::SCEFit)` produces a lightweight predictor from an `SCEFit`.

## Directory layout

```
src/               Package source
test/
  component/       Unit tests (per module)
  integration/     Integration tests (real computational examples)
tools/             Scripts independent of the package source
  vasp/            VASP I/O conversion utilities
  test/            Tests for tools (run via `make test-tools`)
  personal/        Personal scripts (not quality-assured)
docs/              Documenter.jl documentation
examples/          Usage examples (basic_flow / CIF input / save-load)
```

## Public API (exported)

```julia
# Construction
SCEBasis(toml_path)                # Build SCE basis from TOML (includes SALC)
SCEBasis(input_dict)               # From a parsed TOML dict
SCEBasis(system; interaction, ...) # From an AtomsBase.AbstractSystem
SCEBasis(; lattice, kd, kd_list, positions, periodicity, interaction, ...)

SCEDataset(basis, spinconfigs)     # From a SpinConfig vector
SCEDataset(basis, embset_path)     # From an EMBSET file
SCEDataset(system, spinconfigs; interaction, ...)  # AtomsBase shortcut
SCEDataset(toml_path, spinconfigs) # TOML shortcut

# Fitting
fit(SCEFit, dataset, estimator; torque_weight = 1.0) -> SCEFit
SCEModel(f::SCEFit) -> SCEModel    # Lightweight prediction-only conversion

# Prediction (data may be: AbstractMatrix, SpinConfig,
#                          AbstractVector{<:AbstractMatrix},
#                          AbstractVector{SpinConfig}, or SCEDataset.
#                          Single inputs return a scalar / 3 × num_atoms
#                          matrix; vector / dataset inputs return a
#                          Vector{Float64} / Vector{Matrix{Float64}}.)
predict_energy(model_or_fit, data)
predict_torque(model_or_fit, data)

# Evaluation (StatsAPI conventions + Magesty natives)
coef(f), intercept(f), nobs(f), dof(f)
r2_energy(f), r2_torque(f)
rmse_energy(f), rmse_torque(f)
rss_energy(f), rss_torque(f)
residuals_energy(f), residuals_torque(f)

# Persistence (XML) — not exported; call through the module to avoid
# clashing with `save` / `load` from JLD2, FileIO, CSV.jl, etc.
Magesty.save(basis_or_model, path) # path is .xml
Magesty.load(SCEBasis, path)       # Read SCEBasis (model XML also accepted)
Magesty.load(SCEModel, path)       # Read SCEModel

# Estimators
AbstractEstimator, OLS, Ridge

# Data reading
read_embset(path)                  # EMBSET -> Vector{SpinConfig}
SpinConfig(energy, magmom_size, spin_directions, local_magfield)
```

## Primary external libraries

| Library | Purpose |
|---|---|
| `AtomsBase` / `Unitful` | AtomsBase input path |
| `StatsAPI` | `fit` / `coef` / `nobs` / `dof` conventions |
| `Spglib` | Space-group analysis |
| `StaticArrays` | Stack-allocated arrays (performance) |
| `EzXML` | XML I/O |
| `WignerD`, `WignerSymbols`, `LegendrePolynomials` | Angular-momentum coupling |
| `LinearAlgebra`, `Statistics` | Linear algebra and statistics |
| `MultivariateStats` | Ridge regression (`ridge`) |
