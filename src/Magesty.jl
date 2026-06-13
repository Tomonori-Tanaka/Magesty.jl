"""
	Magesty

The main module of the Magesty package, providing the entry point for
spin-cluster expansion (SCE) analysis and fitting.

# Usage
```julia
using Magesty

basis   = SCEBasis("input.toml")
dataset = SCEDataset(basis, "EMBSET")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4))   # torque_weight = 1.0 by default
Magesty.save(SCEModel(f), "model.xml")                 # save / load are not exported; call via the module
```

# Main types
- `SCEBasis`: structure + symmetry + SALC basis
- `SCEDataset`: a basis paired with training data and design matrices
- `SCEFit <: StatsAPI.RegressionModel`: a fitted regression
- `SCEModel`: a basis plus fitted `(j0, jphi)` — the persistable predictor

# Submodules
- `Structures`: crystal structure processing
- `Symmetries`: symmetry operations processing
- `SALCBases`: SALC basis generation
- `Fitting`: design-matrix construction and regression internals
"""
module Magesty

using Printf
using TOML
using LinearAlgebra: norm, dot, det, cross, inv
using StaticArrays: SMatrix, SVector, MVector
using Statistics: mean, std
using Random: MersenneTwister, randperm
import AtomsBase
import DataStructures
import StatsAPI                          # for the StatsAPI.RegressionModel supertype
import StatsAPI: fit, coef, nobs, dof    # extended with SCEFit / SCEModel methods

include("SortedCounters.jl")
include("AtomCells.jl")
include("SpinConfigs.jl")
using .SpinConfigs: SpinConfig, read_embset

include("SphericalHarmonicsTransforms.jl")
include("AngularMomentumCoupling.jl")
include("CoupledBases.jl")
using .AngularMomentumCoupling
using .CoupledBases

include("InputSpecs.jl")
include("AtomsBaseAdapter.jl")
include("RotationMatrix.jl")
include("TesseralHarmonics.jl")
using .InputSpecs
using .AtomsBaseAdapter

include("Structures.jl")
include("Symmetries.jl")
include("Clusters.jl")
include("ProgressReporting.jl")
include("SALCBases.jl")
include("Fitting.jl")

using .Structures
using .Symmetries
using .Clusters
using .SALCBases
using .Fitting
# Explicit import so that the AdaptiveLasso convenience constructors
# added below (AdaptiveLasso(::SCEFit) / AdaptiveLasso(::SCEModel))
# extend Fitting.AdaptiveLasso without triggering the Julia 1.12
# "constructor extended without explicit qualification" warning.
import .Fitting: AdaptiveLasso

include("XMLIO.jl")
using .XMLIO

# VASP I/O and extxyz writer; reached as `Magesty.VaspIO` / `Magesty.ExtXYZ`
# and used by the `vasp_to_extxyz` conversion API.
include("ExtXYZ.jl")
include("VaspIO.jl")

# Code-agnostic MFA spin sampler and the VASP INCAR text reader/writer;
# reached as `Magesty.MfaSampling` / `Magesty.IncarIO` and used by the
# `sample_mfa_incar` API.
include("MfaSampling.jl")
include("IncarIO.jl")

export SCEBasis, SCEDataset, SCEFit, SCEModel
export AbstractEstimator, OLS, Ridge, ElasticNet, Lasso, AdaptiveLasso,
	PrecomputedPilot, AdaptiveRidge
export predict_energy, predict_torque
export write_energies, write_torques
export write_gcv_lambda, write_gcv_learning_curve
# `refit` is a Magesty-native verb (post-selection refit on the basis
# support); it is not imported from `StatsAPI`, which currently does not
# define `refit`. If `StatsAPI` ever adds it, revisit this export so the
# two definitions do not collide silently.
export fit, refit, coef, intercept, nobs, dof
export r2_energy, r2_torque, rss_energy, rss_torque
export residuals_energy, residuals_torque, rmse_energy, rmse_torque
export gcv, gcv_r2, gcv_lambda, gcv_learning_curve, GCVLambdaPath, GCVSizeCurve
export SpinConfig, read_embset
export vasp_to_extxyz, poscar_to_toml, oszicar_to_embset
export sample_mfa_incar
export sce_to_sunny

# Shared skeleton for the SCEBasis input-driven constructors.
# Returns the (structure, symmetry, cluster) triplet; callers append the
# SALCBasis — either computed via `SALCBasis(...)` or loaded from XML.
function _build_structure_skeleton(
	system::SystemSpec,
	options::SymmetryOptions,
	interaction::InteractionSpec;
	verbosity::Bool = true,
)
	structure::Structure = Structure(system, verbosity = verbosity)
	symmetry::Symmetry = Symmetry(structure, options, verbosity = verbosity)
	cluster::Cluster = Cluster(structure, symmetry, interaction, verbosity = verbosity)
	return structure, symmetry, cluster
end

"""
	SCEBasis

Material + basis: crystal structure, symmetry, and SALC basis. The heavy
part is the SALC construction; an `SCEBasis` can be persisted and
reused. `Cluster` is a construction step, not a stored field — it is
computed inside the constructor to build `salcbasis`, then discarded.

# Fields
- `structure::Structure`: Crystal structure information.
- `symmetry::Symmetry`: Symmetry operations.
- `salcbasis::SALCBasis`: SALC basis functions.
- `isotropy::Bool`: Whether the basis was built with the isotropic
  restriction (only `Lf = 0` terms). Provenance metadata: it is not
  derivable from `salcbasis` alone without inspecting every basis
  function.
- `salc_fingerprint::UInt64`: Structural fingerprint of `salcbasis`
  computed via `SALCBases.salc_fingerprint(salcbasis)`. Used by the
  basis-identity check between an `SCEModel` / `SCEFit` and an
  `SCEDataset` so that a basis reloaded from disk still matches the
  in-memory original. Stable across `Magesty.save` / `Magesty.load`
  because the recipe hashes only integer structural identifiers; see
  `SALCBases.salc_fingerprint` for the included / excluded fields.

# Examples
```julia
# Build from a TOML input file (most common path):
basis = SCEBasis("input.toml")

# The SALC construction is the expensive step; persist for reuse:
Magesty.save(basis, "basis.xml")
basis2 = Magesty.load(SCEBasis, "basis.xml")
```
"""
struct SCEBasis
	structure::Structure
	symmetry::Symmetry
	salcbasis::SALCBasis
	isotropy::Bool
	salc_fingerprint::UInt64
end

# Forwards four concrete components and derives `salc_fingerprint`
# from `salcbasis`. The default 5-arg constructor is also callable
# directly when a test needs to inject a deliberately mismatching
# fingerprint.
SCEBasis(
	structure::Structure,
	symmetry::Symmetry,
	salcbasis::SALCBasis,
	isotropy::Bool,
) = SCEBasis(structure, symmetry, salcbasis, isotropy,
             SALCBases.salc_fingerprint(salcbasis))

"""
	SCEBasis(input_dict::AbstractDict; verbosity = true) -> SCEBasis
	SCEBasis(toml_path::AbstractString; verbosity = true) -> SCEBasis
	SCEBasis(system::AtomsBase.AbstractSystem; interaction, name, tolerance_sym, isotropy, verbosity) -> SCEBasis
	SCEBasis(; lattice, kd, kd_list, positions, periodicity, interaction, name, tolerance_sym, isotropy, verbosity) -> SCEBasis
	SCEBasis(x::Union{SCEBasis, SCEModel, SCEDataset, SCEFit}) -> SCEBasis

Construct an `SCEBasis` from one of four input paths: a parsed TOML
dictionary, a TOML file, an `AtomsBase.AbstractSystem`, or raw Julia
keyword arguments. All paths funnel into the three typed value objects
(`SystemSpec`, `InteractionSpec`, `SymmetryOptions`) and then run the
same structure / symmetry / cluster / SALC construction.

The `interaction` argument is a nested `NamedTuple` keyed `body1`,
`body2`, ...; see the API documentation for the accepted format and
the per-path Unitful requirements.

The extractor form `SCEBasis(x)` returns the basis embedded in an
`SCEModel`, `SCEDataset`, or `SCEFit` (and is the identity on an
`SCEBasis`), so generic code can accept any of the four uniformly. The
returned basis shares storage with the input — `SCEBasis` is treated as
a value object and is not mutated by the public API.
"""
function SCEBasis(
	system_spec::SystemSpec,
	interaction::InteractionSpec,
	options::SymmetryOptions;
	verbosity::Bool = true,
)
	structure, symmetry, cluster =
		_build_structure_skeleton(system_spec, options, interaction;
		                          verbosity = verbosity)
	salcbasis::SALCBasis =
		SALCBasis(structure, symmetry, cluster, interaction, options;
		          verbosity = verbosity)
	# `cluster` is a construction step only; it is not stored in SCEBasis.
	return SCEBasis(structure, symmetry, salcbasis, options.isotropy)
end

function SCEBasis(
	input_dict::AbstractDict{<:AbstractString, <:Any};
	verbosity::Bool = true,
)
	specs = parse_toml_inputs(input_dict)
	return SCEBasis(specs...; verbosity = verbosity)
end

function SCEBasis(toml_path::AbstractString; verbosity::Bool = true)
	try
		open(toml_path) do io
			input_dict = TOML.parse(read(io, String))
			return SCEBasis(input_dict; verbosity = verbosity)
		end
	catch e
		if isa(e, SystemError)
			throw(SystemError("Failed to read file: $toml_path"))
		elseif isa(e, TOML.ParserError)
			throw(ErrorException("Failed to parse TOML file: $toml_path\n" *
				sprint(showerror, e)))
		else
			rethrow()
		end
	end
end

function SCEBasis(
	system::AtomsBase.AbstractSystem;
	interaction::NamedTuple,
	name::AbstractString = "system",
	tolerance_sym::Real = 1e-5,
	isotropy::Bool = false,
	verbosity::Bool = true,
)
	specs = AtomsBaseAdapter.system_to_specs(
		system,
		interaction;
		name = name,
		tolerance_sym = tolerance_sym,
		isotropy = isotropy,
	)
	return SCEBasis(specs...; verbosity = verbosity)
end

function SCEBasis(;
	lattice::AbstractMatrix,
	kd::AbstractVector{Symbol},
	kd_list::AbstractVector{<:Integer},
	positions::AbstractVector,
	periodicity = (true, true, true),
	interaction::NamedTuple,
	name::AbstractString = "system",
	tolerance_sym::Real = 1e-5,
	isotropy::Bool = false,
	verbosity::Bool = true,
)
	specs = AtomsBaseAdapter.kwargs_to_specs(
		lattice = lattice,
		kd = kd,
		kd_list = kd_list,
		positions = positions,
		periodicity = periodicity,
		interaction = interaction,
		name = name,
		tolerance_sym = tolerance_sym,
		isotropy = isotropy,
	)
	return SCEBasis(specs...; verbosity = verbosity)
end


"""
	SCEModel

A fitted SCE model: a `SCEBasis` together with the fitted reference
energy `j0` and SCE coefficients `jphi`. This is the lightweight,
persistable predictor — `predict_energy` / `predict_torque` and the
evaluation verbs accept it. Build one from a `SCEFit` via `SCEModel(f)`.

# Fields
- `basis::SCEBasis`: Structure, symmetry, and SALC basis the fit used.
- `j0::Float64`: Reference energy (bias term) in eV.
- `jphi::Vector{Float64}`: SCE coefficients.

The constructor checks that `length(jphi)` matches the number of SALCs
in `basis` — there is one coefficient per SALC.

# Examples
```julia
# From a fit:
model = SCEModel(f)

# Predict on any dataset (e.g., the held-out split):
ŷ_E = predict_energy(model, test)
ŷ_T = predict_torque(model, test)

# Persist and reload (round-trip is byte-stable):
Magesty.save(model, "model.xml")
model2 = Magesty.load(SCEModel, "model.xml")
```
"""
struct SCEModel
	basis::SCEBasis
	j0::Float64
	jphi::Vector{Float64}

	function SCEModel(basis::SCEBasis, j0::Real, jphi::AbstractVector{<:Real})
		n_salc = length(basis.salcbasis.salc_list)
		length(jphi) == n_salc || throw(ArgumentError(
			"length(jphi) ($(length(jphi))) must equal the number of SALCs ($n_salc)"))
		# `new` converts j0/jphi to the field types; a Vector{Float64} jphi
		# is kept as-is (convert is identity), preserving reference sharing.
		return new(basis, j0, jphi)
	end
end


"""
	SCEDataset

A `SCEBasis` paired with training data: the spin configurations and the
unweighted design matrices and observation vectors derived from them.

`X_E` and `X_T` are stored unweighted. The torque weight is applied later
at `fit` time, so a weight sweep reuses one `SCEDataset` without
rebuilding design matrices.

# Fields
- `basis::SCEBasis`: The SCE basis the design matrices were built from.
- `spinconfigs::Vector{SpinConfig}`: Training spin configurations.
- `X_E::Matrix{Float64}`: Energy design matrix, one row per
  configuration, one column per SALC. No bias column — the reference
  energy `j0` is recovered analytically at `fit` time and is not a
  column of `X_E`. Unweighted.
- `X_T::Matrix{Float64}`: Torque design matrix, no bias column, a
  `3 * num_atoms` block of rows per configuration. Unweighted.
- `y_E::Vector{Float64}`: Observed energies, length `n_configs`.
- `y_T::Vector{Float64}`: Observed torques, flattened, length
  `3 * num_atoms * n_configs`.

# Examples
```julia
# Pair a basis with training data read from an EMBSET file:
dataset = SCEDataset(basis, "EMBSET")

# Indexing yields a new SCEDataset — handy for train/test splits:
train = dataset[1:80]
test  = dataset[81:end]

# Concatenate datasets that share the same basis (e.g., for cross-validation):
combined = vcat(train, test)
```
"""
struct SCEDataset
	basis::SCEBasis
	spinconfigs::Vector{SpinConfig}
	X_E::Matrix{Float64}
	X_T::Matrix{Float64}
	y_E::Vector{Float64}
	y_T::Vector{Float64}
end

"""
	SCEDataset(basis::SCEBasis, spinconfigs::AbstractVector{SpinConfig}) -> SCEDataset
	SCEDataset(basis::SCEBasis, embset_path::AbstractString) -> SCEDataset
	SCEDataset(model::SCEModel, spinconfigs_or_embset_path) -> SCEDataset
	SCEDataset(f::SCEFit, spinconfigs_or_embset_path) -> SCEDataset
	SCEDataset(system::AtomsBase.AbstractSystem, spinconfigs; interaction, ...) -> SCEDataset
	SCEDataset(toml_path::AbstractString, spinconfigs::AbstractVector{SpinConfig}) -> SCEDataset

Build a `SCEDataset` from a basis and training data. The base method
takes an explicit `SCEBasis` and a vector of `SpinConfig`; it builds the
unweighted energy and torque design matrices once. The `embset_path`
method reads the configurations from an EMBSET file first. Passing a
fitted `SCEModel` or `SCEFit` in place of the basis reuses the basis
embedded in it (`model.basis` / `f.dataset.basis`) without rebuilding
it. The two sugar methods build a throwaway `SCEBasis` internally (from
an `AtomsBase.AbstractSystem` or a TOML file) and embed it in the
dataset; for workflows that share one basis across several datasets,
construct the `SCEBasis` explicitly and pass it to the base method.
"""
function SCEDataset(
	basis::SCEBasis,
	spinconfigs::AbstractVector{SpinConfig};
	verbosity::Bool = true,
)
	_check_active_atoms_have_moment(basis, spinconfigs)
	X_E = Fitting.build_design_matrix_energy(
		basis.salcbasis.salc_list,
		spinconfigs,
		basis.symmetry;
		verbosity = verbosity,
	)
	X_T = Fitting.build_design_matrix_torque(
		basis.salcbasis.salc_list,
		spinconfigs,
		basis.structure.supercell.num_atoms,
		basis.symmetry;
		verbosity = verbosity,
	)
	y_E::Vector{Float64} = [sc.energy for sc in spinconfigs]
	y_T::Vector{Float64} =
		isempty(spinconfigs) ? Float64[] :
		reduce(vcat, (vec(sc.torques) for sc in spinconfigs))
	return SCEDataset(basis, collect(SpinConfig, spinconfigs), X_E, X_T, y_E, y_T)
end

function SCEDataset(
	basis::SCEBasis,
	embset_path::AbstractString;
	verbosity::Bool = true,
)
	return SCEDataset(basis, SpinConfigs.read_embset(embset_path); verbosity = verbosity)
end

function SCEDataset(
	system::AtomsBase.AbstractSystem,
	spinconfigs::AbstractVector{SpinConfig};
	interaction::NamedTuple,
	name::AbstractString = "system",
	tolerance_sym::Real = 1e-5,
	isotropy::Bool = false,
	verbosity::Bool = true,
)
	basis = SCEBasis(
		system;
		interaction = interaction,
		name = name,
		tolerance_sym = tolerance_sym,
		isotropy = isotropy,
		verbosity = verbosity,
	)
	return SCEDataset(basis, spinconfigs; verbosity = verbosity)
end

function SCEDataset(
	toml_path::AbstractString,
	spinconfigs::AbstractVector{SpinConfig};
	verbosity::Bool = true,
)
	return SCEDataset(
		SCEBasis(toml_path; verbosity = verbosity),
		spinconfigs;
		verbosity = verbosity,
	)
end

# Atoms that appear in any SALC carry spin degrees of freedom in the model.
# A spin configuration that puts a (near-)zero moment on such an atom is
# ill-posed: the energy basis still reads the atom's now-meaningless direction,
# silently biasing the fit. Catching it here covers both a non-magnetic site
# wrongly pulled into the basis (fix the interaction cutoffs / per-species
# `lmax`) and a magnetic site that collapsed to zero moment in a configuration.
# Non-magnetic atoms that stay out of the basis are not referenced and are
# therefore allowed to carry zero moment.
function _check_active_atoms_have_moment(
	basis::SCEBasis,
	spinconfigs::AbstractVector{SpinConfig},
)::Nothing
	referenced = Set{Int}()
	for group in basis.salcbasis.salc_list, cbc in group
		for a in cbc.atoms
			push!(referenced, a)
		end
	end
	isempty(referenced) && return nothing

	num_atoms = basis.structure.supercell.num_atoms
	for (config_index, sc) in enumerate(spinconfigs)
		# A dimension mismatch is reported by the design-matrix build; skip
		# the moment check rather than raising a cryptic bounds error here.
		length(sc.magmom_size) == num_atoms || continue
		for a in referenced
			if sc.magmom_size[a] ≤ SpinConfigs.ZERO_MOMENT_ATOL
				throw(
					ArgumentError(
						"atom $a in spin configuration $config_index has a " *
						"(near-)zero magnetic moment (‖m‖ = $(sc.magmom_size[a])) " *
						"but is referenced by the SALC basis, so its undefined " *
						"spin direction would bias the fit. Exclude non-magnetic " *
						"sites from the interaction cutoffs / per-species lmax, or " *
						"drop configurations where a magnetic site collapses to " *
						"zero moment.",
					),
				)
			end
		end
	end
	return nothing
end

# --- SCEDataset slicing -------------------------------------------------

Base.length(d::SCEDataset)::Int = length(d.spinconfigs)
Base.lastindex(d::SCEDataset)::Int = length(d)
Base.firstindex(d::SCEDataset)::Int = 1

# Torque rows for config indices `idx`: each config owns a contiguous
# `block_size` block, mirroring `build_design_matrix_torque`'s layout.
function _torque_rows(d::SCEDataset, idx::AbstractVector{<:Integer})::Vector{Int}
	block_size = 3 * d.basis.structure.supercell.num_atoms
	rows = Vector{Int}(undef, block_size * length(idx))
	for (k, c) in enumerate(idx)
		dst = block_size * (k - 1)
		src = block_size * (c - 1)
		@inbounds for b = 1:block_size
			rows[dst+b] = src + b
		end
	end
	return rows
end

function Base.getindex(d::SCEDataset, r::AbstractVector{<:Integer})::SCEDataset
	trows = _torque_rows(d, r)
	return SCEDataset(
		d.basis,
		d.spinconfigs[r],
		d.X_E[r, :],
		d.X_T[trows, :],
		d.y_E[r],
		d.y_T[trows],
	)
end

Base.getindex(d::SCEDataset, i::Integer)::SCEDataset = d[i:i]
Base.getindex(d::SCEDataset, r::AbstractVector{Bool})::SCEDataset = d[findall(r)]

# Two SCEBasis objects are interchangeable for design-matrix purposes
# when they are the same object or share a structural fingerprint. The
# fingerprint hashes the SALC structure (see `SALCBases.salc_fingerprint`),
# so an equal fingerprint guarantees the same SALC count and the same
# design-matrix column ordering. Bases reloaded from the same XML are
# distinct objects but compare equal here.
_same_basis(b1::SCEBasis, b2::SCEBasis)::Bool =
	b1 === b2 || b1.salc_fingerprint == b2.salc_fingerprint

function Base.vcat(d1::SCEDataset, ds::SCEDataset...)::SCEDataset
	for d in ds
		_same_basis(d.basis, d1.basis) || throw(ArgumentError(
			"vcat requires all SCEDataset arguments to share the same SALC " *
			"basis: their SCEBasis fingerprints differ, so the design-matrix " *
			"columns follow different (l, m, site) orderings and cannot be " *
			"stacked. Build the datasets from the same SCEBasis, or from " *
			"bases reloaded from the same XML."))
	end
	all_d = (d1, ds...)
	return SCEDataset(
		d1.basis,
		reduce(vcat, (d.spinconfigs for d in all_d)),
		reduce(vcat, (d.X_E for d in all_d)),
		reduce(vcat, (d.X_T for d in all_d)),
		reduce(vcat, (d.y_E for d in all_d)),
		reduce(vcat, (d.y_T for d in all_d)),
	)
end


"""
	SCEFit

A fitted SCE regression. Holds the dataset it was fit on, the fitted
coefficients, the estimator and torque weight used, and the residuals
of the augmented (weighted) least-squares system.

`SCEFit <: StatsAPI.RegressionModel`. The response-block-independent
verbs `coef`, `intercept`, `nobs`, `dof` are defined for it.

# Fields
- `dataset::SCEDataset`: The training dataset.
- `j0::Float64`: Fitted reference energy (bias term).
- `jphi::Vector{Float64}`: Fitted SCE coefficients.
- `estimator::AbstractEstimator`: The estimator used (`OLS` / `Ridge`).
- `torque_weight::Float64`: Torque weight in `[0, 1]` used at fit time.
- `residuals::Vector{Float64}`: Residuals of the augmented *weighted*
  least-squares system (energy rows stacked above flattened torque rows).
  These carry the `torque_weight` scaling and are not in physical units;
  use `rmse_energy(f)` / `rmse_torque(f)` / `r2_energy(f)` /
  `r2_torque(f)` for interpretable in-sample errors.

# Examples
```julia
# Fit (see `fit(SCEFit, ...)` for full options):
f = fit(SCEFit, dataset, Ridge(lambda = 1e-4))

# Inspect:
coef(f)         # SCE coefficient vector
intercept(f)    # fitted reference energy j0 (eV)
rmse_torque(f)  # in-sample torque RMSE
r2_torque(f)    # in-sample torque R²

# Convert to the persistable predictor:
model = SCEModel(f)
```
"""
struct SCEFit <: StatsAPI.RegressionModel
	dataset::SCEDataset
	j0::Float64
	jphi::Vector{Float64}
	estimator::AbstractEstimator
	torque_weight::Float64
	residuals::Vector{Float64}
end

# `SCEDataset` constructors that reuse the basis embedded in a fitted
# `SCEModel` / `SCEFit`. Defined here because `SCEFit` is declared after
# `SCEDataset`. They forward to the basis-based base methods above.
SCEDataset(model::SCEModel, spinconfigs::AbstractVector{SpinConfig}) =
	SCEDataset(model.basis, spinconfigs)
SCEDataset(model::SCEModel, embset_path::AbstractString) =
	SCEDataset(model.basis, embset_path)
SCEDataset(f::SCEFit, spinconfigs::AbstractVector{SpinConfig}) =
	SCEDataset(f.dataset.basis, spinconfigs)
SCEDataset(f::SCEFit, embset_path::AbstractString) =
	SCEDataset(f.dataset.basis, embset_path)

# `SCEBasis` extractors for objects that embed one. Defined here because
# `SCEModel` / `SCEDataset` / `SCEFit` are declared after `SCEBasis`. The
# embedded basis is returned by reference — the same sharing semantics as
# accessing `.basis` directly. `SCEBasis` is treated as a value object
# (no public API mutates its contents), so this is safe in practice.
SCEBasis(basis::SCEBasis) = basis
SCEBasis(model::SCEModel) = model.basis
SCEBasis(dataset::SCEDataset) = dataset.basis
SCEBasis(f::SCEFit) = f.dataset.basis

"""
	fit(::Type{SCEFit}, dataset::SCEDataset, estimator::AbstractEstimator;
	    torque_weight::Real = 1.0, verbosity::Bool = true) -> SCEFit

Fit SCE coefficients on `dataset` with `estimator`, returning a `SCEFit`.

`torque_weight` in `[0, 1]` sets the convex combination of the per-sample
energy and torque mean squared errors that the augmented least-squares
problem minimizes:

```
loss(jphi, j0) = (1 - torque_weight) * MSE_energy + torque_weight * MSE_torque
```

with

```
MSE_energy = (1 / n_E) * Σ_{i=1..n_E}            (y_E[i] - ŷ_E[i])^2
MSE_torque = (1 / n_T) * Σ_{k=1..n_T}            (y_T[k] - ŷ_T[k])^2
```

where `n_E` is the number of configurations,
`n_T = 3 * num_atoms * n_E` is the total number of torque components,
`y_E` / `y_T` are the observed energies and flattened torques, and
`ŷ_E` / `ŷ_T` are the predictions of the model parameterised by
`(j0, jphi)`. Dividing each block by its sample count makes the two
terms commensurate, so the convex combination is meaningful regardless
of how many torque components a system has.

Limiting cases:

| `torque_weight` | Behavior                                                       |
|-----------------|----------------------------------------------------------------|
| `0`             | Energy-only fit; torques are ignored.                          |
| `1`             | Torque-only fit (default); energies enter through `j0` only.   |
| `0 < w < 1`     | Joint fit. `0.5` weighs both per-sample MSEs equally.          |

The default `1.0` is chosen on physical grounds: the SCE coefficients
`jphi` are best determined by torque residuals, which carry the per-atom
directional information that drives the response of the spin model.
Energies enter the fit through the closed-form reference-energy
recovery `j0 = mean(y_E - X_E * jphi)`, so a torque-only fit
still yields a usable `j0`. Set `torque_weight < 1` only when the
energies carry information that the torques do not — typically because
the dataset is energy-rich and torque-poor.

Both MSEs are computed on the raw (unscaled) energy and torque
residuals — the energy unit of the DFT input (typically eV) for the
energy term, and eV per unit spin direction for the torque term. The
two terms therefore live on the *same* physical scale (energy squared),
so the default and the limiting cases are meaningful out of the box.

The design matrices stored in `dataset` are unweighted, so a
`torque_weight` sweep reuses one `SCEDataset` without rebuilding them.

This is the StatsAPI `fit` verb; `using Magesty` re-exports it.

# Arguments
- `::Type{SCEFit}`: Target type — dispatch tag, written literally.
- `dataset::SCEDataset`: Training dataset (basis + design matrices +
  observations).
- `estimator::AbstractEstimator`: Regression estimator (`OLS()` or
  `Ridge(lambda=...)`).
- `torque_weight::Real = 1.0`: Convex weight described above.
- `verbosity::Bool = true`: Whether to print a summary of the fit
  (estimator, sizes, `j0`, in-sample RSS / RMSE / R² on energy and
  torque blocks, elapsed time) to stdout after solving.

# Returns
- `SCEFit`: The fitted model. Inspect with `coef`, `intercept`,
  `r2_energy`, `rmse_torque`, …; persist with
  `Magesty.save(SCEModel(f), path)`.

# Examples
```julia
# Default: torque-only fit with Ridge regularization.
f = fit(SCEFit, dataset, Ridge(lambda = 1e-4))

# Energy-only fit, ignoring torques entirely.
f0 = fit(SCEFit, dataset, OLS(); torque_weight = 0.0)

# Reuse one dataset across a `torque_weight` sweep without rebuilding
# the design matrices.
for w in (0.0, 0.25, 0.5, 0.75, 1.0)
    println(w, "  ", rmse_torque(fit(SCEFit, dataset, OLS();
                                     torque_weight = w)))
end
```
"""
function fit(
	::Type{SCEFit},
	dataset::SCEDataset,
	estimator::AbstractEstimator;
	torque_weight::Real = 1.0,
	verbosity::Bool = true,
)::SCEFit
	start_time = time_ns()
	X, y = Fitting.assemble_weighted_problem(
		dataset.X_E,
		dataset.X_T,
		dataset.y_E,
		dataset.y_T,
		torque_weight,
	)
	# `j0` is eliminated analytically before this point: the energy
	# block of `X` has been mean-centered inside
	# `assemble_weighted_problem`, so the solver returns only `jphi`.
	# `extract_j0_jphi` then recovers `j0` from the unscaled, un-centered
	# energy data via `mean(y_E - X_E * jphi)` -- the closed form of
	# `∂L/∂j0 = 0` -- which keeps `j0` in the input energy unit and
	# independent of `torque_weight` / estimator choice.
	j_values = Fitting.solve_coefficients(estimator, X, y)
	j0, jphi = Fitting.extract_j0_jphi(j_values, dataset.X_E, dataset.y_E)
	residuals::Vector{Float64} = y .- X * j_values
	f = SCEFit(
		dataset,
		j0,
		jphi,
		estimator,
		Float64(torque_weight),
		residuals,
	)
	if verbosity
		elapsed_time = (time_ns() - start_time) / 1e9
		_print_fit_summary(f, elapsed_time)
	end
	return f
end

"""
	refit(fit::SCEFit, estimator::AbstractEstimator = OLS();
	      threshold::Real = 0.0, verbosity::Bool = true) -> SCEFit

Post-selection refit on the basis support of `fit`.

A basis index `j` enters the *support* iff
`abs(coef(fit)[j]) * norm(X[:, j]) > threshold`, where `X` is the
weighted, energy-centered design matrix used by `fit` (reassembled from
`fit.dataset` with the input fit's `torque_weight`). The criterion is
scale-equivariant: multiplying column `j` by any positive constant
leaves the product invariant, because the OLS coefficient scales
inversely while the column norm scales linearly. This cancels both the
per-cluster `(4π)^(N/2)` factor and any other column-scale convention,
so a single criterion serves both selector families: L1 estimators
(`Lasso`, `AdaptiveLasso`) have exact-zero coefficients that fall out
at the default `threshold = 0.0`, while `AdaptiveRidge` users pass a
positive `threshold` to drop near-negligible bases.

The support columns are resolved on the *same* weighted problem with
`estimator`; dropped bases stay at `0.0`, so `jphi` keeps its full
length and SALC ordering — a refit `SCEFit` round-trips through
`SCEModel` and `Magesty.save` like any other fit. `j0` is recovered from
the uncentered energy data via `extract_j0_jphi`, identical to `fit`.

A `PrecomputedPilot`-backed `estimator` (including an `AdaptiveLasso`
whose `pilot` is a `PrecomputedPilot`, such as `AdaptiveLasso(::SCEFit)`)
is rejected upfront: the fixed pilot vector has the original column
count, not the refit support length. A precomputed pilot is also
meaningless once a support has been chosen.

# Arguments
- `fit::SCEFit`: The fitted model whose coefficient support is reused.
  `fit.dataset` and `fit.torque_weight` are reused verbatim, so the
  refit minimizes the same weighted objective as the original `fit`.
- `estimator::AbstractEstimator = OLS()`: Estimator for the support
  resolve. Defaults to `OLS()` — classic post-selection debiasing. Pass
  `Ridge(lambda = small)` when the selected support is still
  rank-deficient or near-collinear.
- `threshold::Real = 0.0`: Scaled-magnitude cutoff. Bases with
  `abs(coef(fit)[j]) * norm(X[:, j]) <= threshold` are dropped. Must be
  `>= 0`. The cutoff is in the scale of the row-whitened augmented
  residual block (energy rows scaled by `sqrt((1-torque_weight)/n_E)`,
  torque rows by `sqrt(torque_weight/n_T)`) rather than raw eV, so pick
  it relative to the largest scaled magnitude in the input fit, e.g.
  `maximum(abs.(coef(fit)) .* [norm(@view X[:, j]) for j in axes(X, 2)])`.
- `verbosity::Bool = true`: Whether to print the standard fit summary,
  exactly as `fit` does.

# Returns
- `SCEFit`: A fresh fit; dropped bases carry coefficient `0.0`. The
  returned `estimator` field is the *refit* estimator; the selection
  estimator is not recorded.

# Examples
Both snippets assume `dataset::SCEDataset` is already in scope (built
from an `SCEBasis` plus spin configurations).
```julia
# L1 selection + OLS debiasing. The default threshold = 0.0 keeps the
# Lasso's exact-zero support verbatim.
fit_lasso = fit(SCEFit, dataset, Lasso(lambda = 1e-3))
fit_db    = refit(fit_lasso, OLS())

# AdaptiveRidge fit: pass a positive threshold to drop near-zero bases,
# then refit with Ridge against residual collinearity.
fit_ar    = fit(SCEFit, dataset, AdaptiveRidge(lambda = 1e-2))
fit_db    = refit(fit_ar, Ridge(lambda = 1e-6); threshold = 1e-4)
```
"""
function refit(
	fit::SCEFit,
	estimator::AbstractEstimator = OLS();
	threshold::Real = 0.0,
	verbosity::Bool = true,
)::SCEFit
	threshold >= 0 || throw(ArgumentError(
		"refit threshold must be non-negative; got $threshold"))
	_reject_precomputed_pilot(estimator)

	start_time = time_ns()
	dataset = fit.dataset
	torque_weight = fit.torque_weight
	X, y = Fitting.assemble_weighted_problem(
		dataset.X_E,
		dataset.X_T,
		dataset.y_E,
		dataset.y_T,
		torque_weight,
	)

	jphi_in = coef(fit)
	num_salcs = length(jphi_in)
	support = findall(
		j -> abs(jphi_in[j]) * norm(@view X[:, j]) > threshold,
		eachindex(jphi_in),
	)

	j_values = zeros(Float64, num_salcs)
	if isempty(support)
		# Short-circuit before `solve_coefficients` -- a GLMNet-backed
		# estimator would error on a zero-column design matrix, and the
		# resulting all-zero `jphi` is well-defined: `extract_j0_jphi`
		# then recovers `j0 = mean(y_E)`. A degenerate but valid fit.
		@warn "refit: empty support — every basis is at or below the " *
			"scaled-magnitude threshold. Returning an all-zero `jphi`; " *
			"`j0` falls back to `mean(y_E)`." threshold
	else
		# `view` avoids a sub-matrix allocation on the OLS / Ridge paths.
		# GLMNet- and AdaptiveRidge-backed solvers materialise their input
		# to a dense `Matrix{Float64}` internally either way, so the view
		# costs nothing there.
		j_sub = Fitting.solve_coefficients(estimator, view(X, :, support), y)
		j_values[support] .= j_sub
	end
	j0, jphi = Fitting.extract_j0_jphi(j_values, dataset.X_E, dataset.y_E)
	residuals::Vector{Float64} = y .- X * j_values
	refit_fit = SCEFit(
		dataset,
		j0,
		jphi,
		estimator,
		torque_weight,
		residuals,
	)
	if verbosity
		elapsed_time = (time_ns() - start_time) / 1e9
		_print_fit_summary(refit_fit, elapsed_time)
	end
	return refit_fit
end

# A `PrecomputedPilot` returns a fixed coefficient vector whose length
# was set against the original design matrix. `refit` solves on a
# sub-matrix `X[:, support]`, so the column counts disagree and the
# pilot would throw `DimensionMismatch` deep in `solve_coefficients`.
# Reject upfront with a clear message — and also reject an
# `AdaptiveLasso` whose pilot is a `PrecomputedPilot`, the form returned
# by `AdaptiveLasso(::SCEFit)` / `AdaptiveLasso(::SCEModel)`.
function _reject_precomputed_pilot(e::AbstractEstimator)
	if e isa PrecomputedPilot
		throw(ArgumentError(
			"refit does not accept a PrecomputedPilot-backed estimator: " *
			"the pilot's fixed coefficient vector has the original column " *
			"count, not the refit support length. A precomputed pilot is " *
			"also meaningless once a support has been chosen."))
	elseif e isa AdaptiveLasso && e.pilot isa PrecomputedPilot
		throw(ArgumentError(
			"refit does not accept an AdaptiveLasso whose pilot is a " *
			"PrecomputedPilot (e.g. AdaptiveLasso(::SCEFit; ...)): the " *
			"pilot's fixed coefficient vector has the original column " *
			"count, not the refit support length. Pass a fresh pilot " *
			"such as OLS() instead."))
	end
	return nothing
end

# --- StatsAPI verbs (response-block independent) ------------------------

"""
	coef(f::SCEFit) -> Vector{Float64}
	coef(m::SCEModel) -> Vector{Float64}

The fitted SCE coefficients `jphi` — a single set shared by the energy
and torque models.
"""
coef(f::SCEFit)::Vector{Float64} = f.jphi
coef(m::SCEModel)::Vector{Float64} = m.jphi

"""
	AdaptiveLasso(f::SCEFit; kwargs...)
	AdaptiveLasso(model::SCEModel; kwargs...)

Build an `AdaptiveLasso` estimator that reuses `coef(f)` / `coef(model)`
as its pilot via `PrecomputedPilot`, skipping a fresh pilot regression
inside `solve_coefficients`. All remaining keyword arguments (`lambda`,
`gamma`, `epsilon`, `standardize`) are forwarded to the standard
`AdaptiveLasso(; ...)` keyword constructor.

The prior fit must have been produced on the same `SCEBasis`: the
adaptive call applies a length check against the new design-matrix
column count but cannot detect a same-length, different-SALC-ordering
mismatch. Passing a `pilot` keyword in addition to the positional
`SCEFit` / `SCEModel` raises `ArgumentError`; the precomputed pilot is
the whole point of these constructors, and Julia's kwarg-splat
semantics would otherwise silently override it.

# Examples
```julia
fit_ols = fit(SCEFit, dataset, OLS(); torque_weight = 0.5)
est = AdaptiveLasso(fit_ols; lambda = 1e-3)
fit_ada = fit(SCEFit, dataset, est; torque_weight = 0.5)
```
"""
function AdaptiveLasso(f::SCEFit; kwargs...)
	haskey(kwargs, :pilot) && throw(ArgumentError(
		"AdaptiveLasso(::SCEFit; ...) sets the pilot from `coef(f)`; " *
		"passing a `pilot` keyword as well would silently override it."))
	return AdaptiveLasso(; pilot = PrecomputedPilot(coef(f)), kwargs...)
end

function AdaptiveLasso(model::SCEModel; kwargs...)
	haskey(kwargs, :pilot) && throw(ArgumentError(
		"AdaptiveLasso(::SCEModel; ...) sets the pilot from `coef(model)`; " *
		"passing a `pilot` keyword as well would silently override it."))
	return AdaptiveLasso(; pilot = PrecomputedPilot(coef(model)), kwargs...)
end

"""
	intercept(f::SCEFit) -> Float64
	intercept(m::SCEModel) -> Float64

The fitted reference energy `j0` (bias term), in the energy unit of the
DFT input (typically eV). The returned value is exactly the `j0` field of
the underlying `SCEFit` / `SCEModel`, kept separate from the SCE
coefficients `jphi`; XML persistence stores it as the `j0` attribute of
the `<JPhi>` block. `intercept` is a Magesty-native verb adopted because
StatsAPI has no `intercept` concept; the two names refer to the same
quantity.
"""
function intercept end
intercept(f::SCEFit)::Float64 = f.j0
intercept(m::SCEModel)::Float64 = m.j0

"""
	nobs(f::SCEFit) -> Int

Number of observations (spin configurations) the fit was trained on, i.e.
the energy-block observation count `n_configs`.
"""
nobs(f::SCEFit)::Int = length(f.dataset)

"""
	dof(f::SCEFit) -> Int

Degrees of freedom consumed by the fit: `length(coef(f)) + 1` — the SCE
coefficients plus the intercept.
"""
dof(f::SCEFit)::Int = length(f.jphi) + 1


# --- SCEModel conversion ------------------------------------------------

"""
	SCEModel(f::SCEFit) -> SCEModel

Extract a lightweight `SCEModel` from a fitted `SCEFit`: the `SCEBasis`
together with the fitted `(j0, jphi)`. The training dataset and
residuals held by the `SCEFit` are dropped.
"""
SCEModel(f::SCEFit)::SCEModel = SCEModel(f.dataset.basis, f.j0, f.jphi)


# --- Prediction ---------------------------------------------------------

"""
	predict_energy(model::SCEModel, spin_directions::AbstractMatrix{<:Real}) -> Float64
	predict_energy(model::SCEModel, sc::SpinConfig) -> Float64
	predict_energy(model::SCEModel, sd_list::AbstractVector{<:AbstractMatrix{<:Real}}) -> Vector{Float64}
	predict_energy(model::SCEModel, configs::AbstractVector{SpinConfig}) -> Vector{Float64}
	predict_energy(model::SCEModel, dataset::SCEDataset) -> Vector{Float64}

`f::SCEFit` may be passed in place of `model`; the `SCEFit` overloads
delegate through `SCEModel(f)`.

Predict SCE energies for one or more spin configurations.

# Arguments
- `model::SCEModel` or `f::SCEFit`: Trained predictor. `SCEFit` inputs
  delegate through `SCEModel(f)`.
- `spin_directions::AbstractMatrix{<:Real}`: Spin direction matrix of size
  `3 × num_atoms` (rows = x, y, z; unit-vector columns).
- `sc::SpinConfig`: Single spin configuration; equivalent to passing
  `sc.spin_directions`.
- `sd_list::AbstractVector{<:AbstractMatrix{<:Real}}`: Sequence of spin
  direction matrices. The number of atoms must match across entries
  and the predictor's `SCEBasis`; the matrices themselves can be views
  or freshly allocated.
- `configs::AbstractVector{SpinConfig}`: Sequence of spin configurations
  (e.g. the output of `read_embset`); only `spin_directions` is read,
  the other fields are ignored.
- `dataset::SCEDataset`: Batch evaluation reusing the dataset's stored
  energy design matrix. Must share the predictor's `SCEBasis` (same
  `(l, m, site)` column ordering as the fitted coefficients).

# Returns

Return type depends on the second positional argument:

| Input form                                                | Return                                                                                       |
|-----------------------------------------------------------|----------------------------------------------------------------------------------------------|
| `spin_directions::AbstractMatrix`                         | `Float64` — energy in the unit of the training data (typically eV).                          |
| `sc::SpinConfig`                                          | `Float64` — same as above, evaluated at `sc.spin_directions`.                                |
| `sd_list::AbstractVector{<:AbstractMatrix}`               | `Vector{Float64}` of length `length(sd_list)`, in input order.                               |
| `configs::AbstractVector{SpinConfig}`                     | `Vector{Float64}` of length `length(configs)`, in input order.                               |
| `dataset::SCEDataset`                                     | `Vector{Float64}` of length `length(dataset)`, in dataset order.                             |
"""
predict_energy(model::SCEModel, spin_directions::AbstractMatrix{<:Real})::Float64 =
	Fitting._predict_energy(
		model.j0, model.jphi,
		model.basis.salcbasis.salc_list, model.basis.symmetry, spin_directions)
predict_energy(model::SCEModel, sc::SpinConfig)::Float64 =
	predict_energy(model, sc.spin_directions)
predict_energy(f::SCEFit, spin_directions::AbstractMatrix{<:Real})::Float64 =
	predict_energy(SCEModel(f), spin_directions)
predict_energy(f::SCEFit, sc::SpinConfig)::Float64 =
	predict_energy(SCEModel(f), sc)

predict_energy(
	model::SCEModel,
	sd_list::AbstractVector{<:AbstractMatrix{<:Real}},
)::Vector{Float64} =
	Float64[predict_energy(model, sd) for sd in sd_list]
predict_energy(
	model::SCEModel,
	configs::AbstractVector{SpinConfig},
)::Vector{Float64} =
	Float64[predict_energy(model, sc) for sc in configs]
predict_energy(
	f::SCEFit,
	sd_list::AbstractVector{<:AbstractMatrix{<:Real}},
)::Vector{Float64} =
	predict_energy(SCEModel(f), sd_list)
predict_energy(
	f::SCEFit,
	configs::AbstractVector{SpinConfig},
)::Vector{Float64} =
	predict_energy(SCEModel(f), configs)

function predict_energy(model::SCEModel, dataset::SCEDataset)::Vector{Float64}
	_check_basis(model, dataset)
	return dataset.X_E * model.jphi .+ model.j0
end
predict_energy(f::SCEFit, dataset::SCEDataset)::Vector{Float64} =
	predict_energy(SCEModel(f), dataset)

"""
	predict_torque(model::SCEModel, spin_directions::AbstractMatrix{<:Real}) -> Matrix{Float64}
	predict_torque(model::SCEModel, sc::SpinConfig) -> Matrix{Float64}
	predict_torque(model::SCEModel, sd_list::AbstractVector{<:AbstractMatrix{<:Real}}) -> Vector{Matrix{Float64}}
	predict_torque(model::SCEModel, configs::AbstractVector{SpinConfig}) -> Vector{Matrix{Float64}}
	predict_torque(model::SCEModel, dataset::SCEDataset) -> Vector{Matrix{Float64}}

`f::SCEFit` may be passed in place of `model`; the `SCEFit` overloads
delegate through `SCEModel(f)`.

Predict per-atom SCE torques for one or more spin configurations.

# Arguments
- `model::SCEModel` or `f::SCEFit`: Trained predictor. `SCEFit` inputs
  delegate through `SCEModel(f)`.
- `spin_directions::AbstractMatrix{<:Real}`: Spin direction matrix of size
  `3 × num_atoms` (rows = x, y, z; unit-vector columns).
- `sc::SpinConfig`: Single spin configuration; equivalent to passing
  `sc.spin_directions`.
- `sd_list::AbstractVector{<:AbstractMatrix{<:Real}}`: Sequence of spin
  direction matrices. The number of atoms must match across entries
  and the predictor's `SCEBasis`.
- `configs::AbstractVector{SpinConfig}`: Sequence of spin configurations
  (e.g. the output of `read_embset`); only `spin_directions` is read.
- `dataset::SCEDataset`: Batch evaluation reusing the dataset's stored
  torque design matrix. Must share the predictor's `SCEBasis` (same
  `(l, m, site)` column ordering as the fitted coefficients).

# Returns

Return type depends on the second positional argument:

| Input form                                                | Return                                                                                                          |
|-----------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| `spin_directions::AbstractMatrix`                         | `Matrix{Float64}` of size `3 × num_atoms` (rows = x, y, z).                                                     |
| `sc::SpinConfig`                                          | `Matrix{Float64}` of size `3 × num_atoms`, evaluated at `sc.spin_directions`.                                   |
| `sd_list::AbstractVector{<:AbstractMatrix}`               | `Vector{Matrix{Float64}}` of length `length(sd_list)`, in input order; each element is `3 × num_atoms`.          |
| `configs::AbstractVector{SpinConfig}`                     | `Vector{Matrix{Float64}}` of length `length(configs)`, in input order; each element is `3 × num_atoms`.          |
| `dataset::SCEDataset`                                     | `Vector{Matrix{Float64}}` of length `length(dataset)`, in dataset order; each element is `3 × num_atoms`.        |
"""
predict_torque(model::SCEModel, spin_directions::AbstractMatrix{<:Real})::Matrix{Float64} =
	Fitting._predict_torque(
		model.jphi, model.basis.salcbasis.salc_list, model.basis.symmetry, spin_directions)
predict_torque(model::SCEModel, sc::SpinConfig)::Matrix{Float64} =
	predict_torque(model, sc.spin_directions)
predict_torque(f::SCEFit, spin_directions::AbstractMatrix{<:Real})::Matrix{Float64} =
	predict_torque(SCEModel(f), spin_directions)
predict_torque(f::SCEFit, sc::SpinConfig)::Matrix{Float64} =
	predict_torque(SCEModel(f), sc)

predict_torque(
	model::SCEModel,
	sd_list::AbstractVector{<:AbstractMatrix{<:Real}},
)::Vector{Matrix{Float64}} =
	Matrix{Float64}[predict_torque(model, sd) for sd in sd_list]
predict_torque(
	model::SCEModel,
	configs::AbstractVector{SpinConfig},
)::Vector{Matrix{Float64}} =
	Matrix{Float64}[predict_torque(model, sc) for sc in configs]
predict_torque(
	f::SCEFit,
	sd_list::AbstractVector{<:AbstractMatrix{<:Real}},
)::Vector{Matrix{Float64}} =
	predict_torque(SCEModel(f), sd_list)
predict_torque(
	f::SCEFit,
	configs::AbstractVector{SpinConfig},
)::Vector{Matrix{Float64}} =
	predict_torque(SCEModel(f), configs)

function predict_torque(model::SCEModel, dataset::SCEDataset)::Vector{Matrix{Float64}}
	_check_basis(model, dataset)
	num_atoms = dataset.basis.structure.supercell.num_atoms
	block_size = 3 * num_atoms
	flat = dataset.X_T * model.jphi
	return [
		reshape(flat[((i-1)*block_size+1):(i*block_size)], 3, num_atoms)
		for i = 1:length(dataset)
	]
end
predict_torque(f::SCEFit, dataset::SCEDataset)::Vector{Matrix{Float64}} =
	predict_torque(SCEModel(f), dataset)


# --- Evaluation verbs ---------------------------------------------------

# Evaluation data accepted by the `(predictor, data)` verb family.
const SCEEvalData = Union{SCEDataset, AbstractVector{SpinConfig}, AbstractString}

# The evaluation data must come from the same SALC basis as the
# predictor: the (l, m, site) column ordering of the design matrices
# must match the fitted coefficients. Only checkable when the data is a
# SCEDataset (configs / EMBSET paths are evaluated through the
# predictor's own basis, so they are compatible by construction).
function _check_basis(model::SCEModel, dataset::SCEDataset)
	_same_basis(model.basis, dataset.basis) && return nothing
	throw(ArgumentError(
		"the evaluation SCEDataset was built from a different SCEBasis " *
		"than the predictor; design-matrix column ordering is set by " *
		"the SCEBasis, so combining objects from different bases would " *
		"silently mix incompatible `(l, m, site)` orderings. Build " *
		"the SCEDataset, SCEFit, and SCEModel from the *same* SCEBasis " *
		"instance, or, if reloading from disk, reuse a single " *
		"`Magesty.load(SCEBasis, path)` result everywhere."))
end
_check_basis(f::SCEFit, dataset::SCEDataset) = _check_basis(SCEModel(f), dataset)

# (observed, predicted) energy vectors for a predictor and evaluation data.
_eval_energy(predictor, dataset::SCEDataset) =
	(dataset.y_E, predict_energy(predictor, dataset))
_eval_energy(predictor, configs::AbstractVector{SpinConfig}) =
	(Float64[sc.energy for sc in configs],
	 Float64[predict_energy(predictor, sc) for sc in configs])
_eval_energy(predictor, embset_path::AbstractString) =
	_eval_energy(predictor, SpinConfigs.read_embset(embset_path))

# (observed, predicted) flattened torque vectors.
function _eval_torque(predictor, dataset::SCEDataset)
	predicted = isempty(dataset.spinconfigs) ? Float64[] :
		reduce(vcat, (vec(t) for t in predict_torque(predictor, dataset)))
	return (dataset.y_T, predicted)
end
function _eval_torque(predictor, configs::AbstractVector{SpinConfig})
	observed = isempty(configs) ? Float64[] :
		reduce(vcat, (vec(sc.torques) for sc in configs))
	predicted = isempty(configs) ? Float64[] :
		reduce(vcat, (vec(predict_torque(predictor, sc)) for sc in configs))
	return (observed, predicted)
end
_eval_torque(predictor, embset_path::AbstractString) =
	_eval_torque(predictor, SpinConfigs.read_embset(embset_path))

"""
	r2_energy(predictor, data) -> Float64
	r2_energy(f::SCEFit) -> Float64

Coefficient of determination (R²) of the SCE energy predictions.

Call as `r2_energy(f)` to evaluate in-sample on the training dataset
embedded in `f`; call as `r2_energy(predictor, data)` to evaluate on a
different dataset — typically a held-out test set.

# Arguments
- `predictor::Union{SCEModel, SCEFit}`: Trained predictor.
- `data`: Evaluation data — `SCEDataset`, `AbstractVector{SpinConfig}`,
  or an EMBSET file path (`AbstractString`).
- `f::SCEFit` (single-argument form): Evaluates `f` in-sample on its own
  training dataset.

# Returns
- `Float64`: R² over the observed and predicted energies (eV).
"""
function r2_energy(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_energy(predictor, data)
	return Fitting._calc_r2score(observed, predicted)
end
r2_energy(f::SCEFit)::Float64 = r2_energy(f, f.dataset)

"""
	r2_torque(predictor, data) -> Float64
	r2_torque(f::SCEFit) -> Float64

Coefficient of determination (R²) of the SCE torque predictions,
computed over the flattened torque components.

Call as `r2_torque(f)` to evaluate in-sample on the training dataset
embedded in `f`; call as `r2_torque(predictor, data)` to evaluate on a
different dataset — typically a held-out test set.

# Arguments
- `predictor::Union{SCEModel, SCEFit}`: Trained predictor.
- `data`: Evaluation data — `SCEDataset`, `AbstractVector{SpinConfig}`,
  or an EMBSET file path (`AbstractString`).
- `f::SCEFit` (single-argument form): Evaluates `f` in-sample on its own
  training dataset.

# Returns
- `Float64`: R² over the flattened observed and predicted torques (eV),
  length `3 * num_atoms * n_configs`.
"""
function r2_torque(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_torque(predictor, data)
	return Fitting._calc_r2score(observed, predicted)
end
r2_torque(f::SCEFit)::Float64 = r2_torque(f, f.dataset)

"""
	rss_energy(predictor, data) -> Float64
	rss_energy(f::SCEFit) -> Float64

Residual sum of squares of the SCE energy predictions.

Call as `rss_energy(f)` to evaluate in-sample on the training dataset
embedded in `f`; call as `rss_energy(predictor, data)` to evaluate on a
different dataset — typically a held-out test set.

# Arguments
- `predictor::Union{SCEModel, SCEFit}`: Trained predictor.
- `data`: Evaluation data — `SCEDataset`, `AbstractVector{SpinConfig}`,
  or an EMBSET file path (`AbstractString`).
- `f::SCEFit` (single-argument form): Evaluates `f` in-sample on its own
  training dataset.

# Returns
- `Float64`: `sum((observed - predicted).^2)` over the energies (eV²).
"""
function rss_energy(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_energy(predictor, data)
	return sum(abs2, observed .- predicted)
end
rss_energy(f::SCEFit)::Float64 = rss_energy(f, f.dataset)

"""
	rss_torque(predictor, data) -> Float64
	rss_torque(f::SCEFit) -> Float64

Residual sum of squares of the SCE torque predictions, over the
flattened torque components.

Call as `rss_torque(f)` to evaluate in-sample on the training dataset
embedded in `f`; call as `rss_torque(predictor, data)` to evaluate on a
different dataset — typically a held-out test set.

# Arguments
- `predictor::Union{SCEModel, SCEFit}`: Trained predictor.
- `data`: Evaluation data — `SCEDataset`, `AbstractVector{SpinConfig}`,
  or an EMBSET file path (`AbstractString`).
- `f::SCEFit` (single-argument form): Evaluates `f` in-sample on its own
  training dataset.

# Returns
- `Float64`: `sum((observed - predicted).^2)` over the flattened torques
  (eV²), length `3 * num_atoms * n_configs`.
"""
function rss_torque(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_torque(predictor, data)
	return sum(abs2, observed .- predicted)
end
rss_torque(f::SCEFit)::Float64 = rss_torque(f, f.dataset)

"""
	residuals_energy(predictor, data) -> Vector{Float64}
	residuals_energy(f::SCEFit) -> Vector{Float64}

Per-configuration energy residuals `observed - predicted` (eV).

Call as `residuals_energy(f)` to evaluate in-sample on the training
dataset embedded in `f`; call as `residuals_energy(predictor, data)` to
evaluate on a different dataset — typically a held-out test set.

# Arguments
- `predictor::Union{SCEModel, SCEFit}`: Trained predictor.
- `data`: Evaluation data — `SCEDataset`, `AbstractVector{SpinConfig}`,
  or an EMBSET file path (`AbstractString`).
- `f::SCEFit` (single-argument form): Evaluates `f` in-sample on its own
  training dataset.

# Returns
- `Vector{Float64}` of length `n_configs`, in input order.
"""
function residuals_energy(
	predictor::Union{SCEModel, SCEFit},
	data::SCEEvalData,
)::Vector{Float64}
	observed, predicted = _eval_energy(predictor, data)
	return observed .- predicted
end
residuals_energy(f::SCEFit)::Vector{Float64} = residuals_energy(f, f.dataset)

"""
	residuals_torque(predictor, data) -> Vector{Float64}
	residuals_torque(f::SCEFit) -> Vector{Float64}

Flattened torque residuals `observed - predicted` (eV).

Call as `residuals_torque(f)` to evaluate in-sample on the training
dataset embedded in `f`; call as `residuals_torque(predictor, data)` to
evaluate on a different dataset — typically a held-out test set.

# Arguments
- `predictor::Union{SCEModel, SCEFit}`: Trained predictor.
- `data`: Evaluation data — `SCEDataset`, `AbstractVector{SpinConfig}`,
  or an EMBSET file path (`AbstractString`).
- `f::SCEFit` (single-argument form): Evaluates `f` in-sample on its own
  training dataset.

# Returns
- `Vector{Float64}` of length `3 * num_atoms * n_configs`, flattened over
  (component, atom, configuration).
"""
function residuals_torque(
	predictor::Union{SCEModel, SCEFit},
	data::SCEEvalData,
)::Vector{Float64}
	observed, predicted = _eval_torque(predictor, data)
	return observed .- predicted
end
residuals_torque(f::SCEFit)::Vector{Float64} = residuals_torque(f, f.dataset)

"""
	rmse_energy(predictor, data) -> Float64
	rmse_energy(f::SCEFit) -> Float64

Root mean squared error of the SCE energy predictions (eV).

Call as `rmse_energy(f)` to evaluate in-sample on the training dataset
embedded in `f`; call as `rmse_energy(predictor, data)` to evaluate on a
different dataset — typically a held-out test set.

# Arguments
- `predictor::Union{SCEModel, SCEFit}`: Trained predictor.
- `data`: Evaluation data — `SCEDataset`, `AbstractVector{SpinConfig}`,
  or an EMBSET file path (`AbstractString`).
- `f::SCEFit` (single-argument form): Evaluates `f` in-sample on its own
  training dataset.

# Returns
- `Float64`: `sqrt(mean((observed - predicted).^2))` over the energies (eV).
"""
function rmse_energy(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_energy(predictor, data)
	return Fitting._calc_rmse(observed, predicted)
end
rmse_energy(f::SCEFit)::Float64 = rmse_energy(f, f.dataset)

"""
	rmse_torque(predictor, data) -> Float64
	rmse_torque(f::SCEFit) -> Float64

Root mean squared error of the SCE torque predictions (eV), over the
flattened torque components.

Call as `rmse_torque(f)` to evaluate in-sample on the training dataset
embedded in `f`; call as `rmse_torque(predictor, data)` to evaluate on a
different dataset — typically a held-out test set.

# Arguments
- `predictor::Union{SCEModel, SCEFit}`: Trained predictor.
- `data`: Evaluation data — `SCEDataset`, `AbstractVector{SpinConfig}`,
  or an EMBSET file path (`AbstractString`).
- `f::SCEFit` (single-argument form): Evaluates `f` in-sample on its own
  training dataset.

# Returns
- `Float64`: `sqrt(mean((observed - predicted).^2))` over the flattened
  torques (eV), length `3 * num_atoms * n_configs`.
"""
function rmse_torque(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_torque(predictor, data)
	return Fitting._calc_rmse(observed, predicted)
end
rmse_torque(f::SCEFit)::Float64 = rmse_torque(f, f.dataset)


# --- GCV diagnostics ----------------------------------------------------
#
# Generalized cross-validation on the combined energy+torque weighted
# objective the fit minimizes. See the `Fitting` GCV-core helpers for the
# `GCV = (‖r‖²/N) / (1 − tr(H)/N)²` formula and the effective-dof computation.
# All entry points are restricted to linear estimators (`OLS` / `Ridge` /
# `AdaptiveRidge`); non-linear estimators raise `ArgumentError`.

"""
	GCVLambdaPath

Result of a ridge GCV penalty sweep ([`gcv_lambda`](@ref)). The combined
energy+torque GCV score is evaluated at each penalty `lambda` from a single
SVD of the weighted design matrix.

# Fields
- `lambdas::Vector{Float64}`: The penalty values swept, in input order.
- `gcv_scores::Vector{Float64}`: GCV score at each `lambda` (`NaN` where the
  model is numerically saturated).
- `gcv_r2::Vector{Float64}`: GCV-based predictive R² at each `lambda`,
  `1 - gcv / msy` (`1` perfect, `0` matches the null model, `< 0` worse than
  null). Interpretable on a fixed scale, unlike the raw `gcv_scores`.
- `dof::Vector{Float64}`: Effective degrees of freedom `tr(H)` at each `lambda`.
- `lambda_best::Float64`: The `lambda` minimizing the GCV score.
- `torque_weight::Float64`: The torque weight the sweep used.
"""
struct GCVLambdaPath
	lambdas::Vector{Float64}
	gcv_scores::Vector{Float64}
	gcv_r2::Vector{Float64}
	dof::Vector{Float64}
	lambda_best::Float64
	torque_weight::Float64
end

function Base.show(io::IO, p::GCVLambdaPath)
	print(io, "GCVLambdaPath(", length(p.lambdas), " lambdas, lambda_best=",
		p.lambda_best, ", torque_weight=", p.torque_weight, ")")
end

"""
	GCVSizeCurve

Result of a data-sufficiency GCV sweep ([`gcv_learning_curve`](@ref)). At each
training-set size, `repeats` random config subsets are fit and scored; the
mean and standard deviation across draws are reported, so a flattening curve
signals that enough data is present.

# Fields
- `sizes::Vector{Int}`: Training-set sizes, ascending.
- `gcv_mean::Vector{Float64}`: Mean GCV over the random draws at each size.
- `gcv_std::Vector{Float64}`: Standard deviation over the draws at each size.
- `gcv_r2_mean::Vector{Float64}`: Mean GCV-based predictive R² over the draws at
  each size (`1 - gcv / msy`; `1` perfect, `0` matches the null model). Read on a
  fixed scale, unlike the raw `gcv_mean`.
- `gcv_r2_std::Vector{Float64}`: Standard deviation of the predictive R² over the
  draws at each size.
- `repeats::Int`: Random draws per size.
- `seed::Int`: RNG seed used for reproducibility.
- `estimator::AbstractEstimator`: The estimator fit on each subset.
- `torque_weight::Float64`: The torque weight the sweep used.
"""
struct GCVSizeCurve
	sizes::Vector{Int}
	gcv_mean::Vector{Float64}
	gcv_std::Vector{Float64}
	gcv_r2_mean::Vector{Float64}
	gcv_r2_std::Vector{Float64}
	repeats::Int
	seed::Int
	estimator::AbstractEstimator
	torque_weight::Float64
end

function Base.show(io::IO, c::GCVSizeCurve)
	print(io, "GCVSizeCurve(", length(c.sizes), " sizes ", first(c.sizes), "..",
		last(c.sizes), ", repeats=", c.repeats, ", ", c.estimator, ")")
end

# Validate the convex energy/torque weight, matching the [0, 1] convention `fit`
# documents. Named so the error attributes to the public GCV entry point.
function _check_gcv_torque_weight(caller::AbstractString, torque_weight::Real)
	0 <= torque_weight <= 1 || throw(ArgumentError(
		"$caller: torque_weight must be in [0, 1]; got $torque_weight"))
	return nothing
end

# Shared GCV setup. Assemble the weighted, energy-centered augmented problem
# `(X, y)` from `dataset` under `torque_weight`, and compute the effective live-row
# count `n_eff` plus the intercept degree of freedom `intercept_dof` (the
# eliminated `j0`, live only when the energy block survives). Every GCV entry
# point — single-fit scoring, the penalty path, and the learning-curve subsets —
# starts from exactly this state, so each caller is a thin wrapper over it.
function _gcv_core(
	dataset::SCEDataset,
	torque_weight::Real,
)::Tuple{Matrix{Float64}, Vector{Float64}, Int, Int}
	X, y = Fitting.assemble_weighted_problem(
		dataset.X_E, dataset.X_T, dataset.y_E, dataset.y_T, torque_weight)
	n_eff, intercept_dof = Fitting._gcv_sample_count(
		length(dataset.y_E), length(dataset.y_T), torque_weight)
	return X, y, n_eff, intercept_dof
end

"""
	gcv(f::SCEFit) -> Float64

Combined energy+torque generalized cross-validation score for the fitted model
`f`, evaluated on its training dataset and the weighted objective `f` was fit
with (same `torque_weight` and estimator).

GCV estimates the out-of-sample prediction error from a single fit via the hat
matrix `H` (`ŷ = H y`):

```
GCV = (‖r‖² / N) / (1 − tr(H)/N)²,
```

where `r` is the augmented weighted residual, `tr(H)` the effective degrees of
freedom, and `N` the number of *live* rows — energy plus torque, minus any block
zeroed by the weighting (`torque_weight = 1` drops the energy block,
`torque_weight = 0` drops the torque block). The eliminated reference energy
`j0` counts one degree of freedom only when the energy block is live
(`torque_weight < 1`). The score is in the weighted-objective unit, not eV²;
compare scores (e.g. across penalties or data sizes), not the absolute
magnitude.

Defined only for linear estimators (`OLS`, `Ridge`, `AdaptiveRidge`).

# Arguments
- `f::SCEFit`: A fitted model whose estimator is linear.

# Returns
- `Float64`: The GCV score, or `NaN` if the model is numerically saturated
  (`tr(H) ≥ N`).

# Throws
- `ArgumentError` if `f.estimator` is non-linear (`ElasticNet` / `Lasso` /
  `AdaptiveLasso`).

# Examples
```julia
f = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
gcv(f)
```
"""
function gcv(f::SCEFit)::Float64
	Fitting._require_linear_estimator(f.estimator)
	X, y, n_eff, intercept_dof = _gcv_core(f.dataset, f.torque_weight)
	score, _ = Fitting._gcv_single(
		X, y, f.residuals, f.estimator, f.jphi, n_eff, intercept_dof)
	return score
end

"""
	gcv_r2(f::SCEFit) -> Float64

GCV-based predictive R² for the fitted model `f`: the [`gcv`](@ref) score
normalized against the null-model mean square `msy = ‖y‖² / N`, namely

```
R²_gcv = 1 − GCV / msy.
```

The null model is `β = 0` on the weighted, energy-centered augmented system
(energy predicted at its mean, torque predicted as zero), so `R²_gcv` measures
the cross-validated variance explained on a fixed scale: `1` is a perfect fit,
`0` matches the null model, and a negative value means the fit predicts worse
than the null (over-parameterized / too little data). Unlike the raw `gcv`
score — which is in the weighted-objective unit and only meaningful in relative
comparison — this value can be read in isolation.

Defined only for linear estimators (`OLS`, `Ridge`, `AdaptiveRidge`).

# Arguments
- `f::SCEFit`: A fitted model whose estimator is linear.

# Returns
- `Float64`: The predictive R², or `NaN` if the model is numerically saturated
  (`tr(H) ≥ N`).

# Throws
- `ArgumentError` if `f.estimator` is non-linear (`ElasticNet` / `Lasso` /
  `AdaptiveLasso`).

# Examples
```julia
f = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
gcv_r2(f)    # ~1 good, ~0 no better than the mean / zero-torque model
```
"""
function gcv_r2(f::SCEFit)::Float64
	Fitting._require_linear_estimator(f.estimator)
	X, y, n_eff, intercept_dof = _gcv_core(f.dataset, f.torque_weight)
	score, _ = Fitting._gcv_single(
		X, y, f.residuals, f.estimator, f.jphi, n_eff, intercept_dof)
	return Fitting._gcv_r2(score, Fitting._gcv_msy(y, n_eff))
end

"""
	gcv_lambda(dataset::SCEDataset, lambdas::AbstractVector{<:Real};
	           torque_weight::Real = 1.0) -> GCVLambdaPath

Ridge GCV penalty sweep: compute the combined energy+torque GCV score for every
`lambda` and report the minimizer. A single SVD of the weighted, energy-centered
design matrix serves the whole path, so passing a fine `lambdas` grid is cheap.

# Arguments
- `dataset::SCEDataset`: The training data (design matrices built once).
- `lambdas::AbstractVector{<:Real}`: Non-negative ridge penalties to evaluate.
- `torque_weight::Real = 1.0`: Convex energy/torque weight in `[0, 1]`, as in
  `fit`.

# Returns
- `GCVLambdaPath`: Per-`lambda` GCV, predictive R² (`gcv_r2`), and effective
  dof, plus `lambda_best`.

# Throws
- `ArgumentError` if `lambdas` is empty, contains a negative value,
  `torque_weight` is outside `[0, 1]`, or every penalty gives a non-finite GCV.

# Examples
```julia
path = gcv_lambda(dataset, 10.0 .^ (-6:0.5:0))
f    = fit(SCEFit, dataset, Ridge(lambda = path.lambda_best))
```
"""
function gcv_lambda(
	dataset::SCEDataset,
	lambdas::AbstractVector{<:Real};
	torque_weight::Real = 1.0,
)::GCVLambdaPath
	isempty(lambdas) && throw(ArgumentError("gcv_lambda: lambdas must be non-empty."))
	any(<(0), lambdas) && throw(ArgumentError(
		"gcv_lambda: all lambdas must be non-negative; got minimum $(minimum(lambdas))"))
	_check_gcv_torque_weight("gcv_lambda", torque_weight)
	X, y, n_eff, intercept_dof = _gcv_core(dataset, torque_weight)
	gcvs, dofs = Fitting._gcv_lambda_path(X, y, lambdas, n_eff, intercept_dof)
	msy = Fitting._gcv_msy(y, n_eff)
	r2s = Fitting._gcv_r2.(gcvs, msy)
	best_idx = _argmin_ignore_nan(gcvs)
	return GCVLambdaPath(
		collect(Float64, lambdas), gcvs, r2s, dofs,
		Float64(lambdas[best_idx]), Float64(torque_weight))
end

"""
	gcv_learning_curve(dataset::SCEDataset, estimator::AbstractEstimator = OLS();
	                   sizes::AbstractVector{<:Integer} = <auto grid>,
	                   repeats::Integer = 5, seed::Integer = 0,
	                   torque_weight::Real = 1.0) -> GCVSizeCurve

Data-sufficiency GCV learning curve. At each training-set size, draw `repeats`
random config subsets, fit `estimator` to each, and average their combined GCV
scores; a curve that flattens with size indicates enough training data.

Subsets are drawn from the full `dataset` with a seeded RNG (reproducible). Each
draw reuses the prebuilt design matrices via row slicing — the heavy
design-matrix construction is not repeated. A draw that yields a rank-deficient
`OLS` solve or a saturated model (`tr(H) ≥ N`) contributes `NaN` and is dropped
from that size's statistics (with a warning); if every draw at a size fails, the
size reports `NaN`.

# Arguments
- `dataset::SCEDataset`: The full training data.
- `estimator::AbstractEstimator = OLS()`: Linear estimator fit on each subset.
- `sizes::AbstractVector{<:Integer}`: Training-set sizes. Defaults to six points
  spanning `max(p + 2, 10)` to `length(dataset)` (`p` = number of SALCs).
- `repeats::Integer = 5`: Random draws averaged per size.
- `seed::Integer = 0`: RNG seed.
- `torque_weight::Real = 1.0`: Convex energy/torque weight in `[0, 1]`, as in
  `fit`.

# Returns
- `GCVSizeCurve`: `sizes`, `gcv_mean`, `gcv_std`, the predictive-R² summary
  (`gcv_r2_mean`, `gcv_r2_std`), and the sweep settings.

# Throws
- `ArgumentError` if `estimator` is non-linear, `repeats < 1`, `torque_weight`
  is outside `[0, 1]`, or any size is outside `1:length(dataset)`.

# Examples
```julia
curve = gcv_learning_curve(dataset, Ridge(lambda = 1e-4); repeats = 8)
curve.sizes        # training-set sizes
curve.gcv_mean     # mean GCV at each size — look for a plateau
curve.gcv_r2_mean  # mean predictive R² — ~1 good, plateau ⇒ enough data
curve.gcv_std      # spread across draws — large spread suggests more data needed
```
"""
function gcv_learning_curve(
	dataset::SCEDataset,
	estimator::AbstractEstimator = OLS();
	sizes::AbstractVector{<:Integer} = _default_gcv_sizes(dataset),
	repeats::Integer = 5,
	seed::Integer = 0,
	torque_weight::Real = 1.0,
)::GCVSizeCurve
	Fitting._require_linear_estimator(estimator)
	repeats >= 1 ||
		throw(ArgumentError("gcv_learning_curve: repeats must be >= 1; got $repeats"))
	_check_gcv_torque_weight("gcv_learning_curve", torque_weight)
	n_total = length(dataset)
	sorted_sizes = sort(unique(Int.(sizes)))
	(!isempty(sorted_sizes) && all(s -> 1 <= s <= n_total, sorted_sizes)) ||
		throw(ArgumentError(
			"gcv_learning_curve: every size must be in 1:$n_total; got $(collect(sizes))"))
	rng = MersenneTwister(seed)
	means = Vector{Float64}(undef, length(sorted_sizes))
	stds = Vector{Float64}(undef, length(sorted_sizes))
	r2_means = Vector{Float64}(undef, length(sorted_sizes))
	r2_stds = Vector{Float64}(undef, length(sorted_sizes))
	for (i, n) in enumerate(sorted_sizes)
		scores = Vector{Float64}(undef, repeats)
		r2_scores = Vector{Float64}(undef, repeats)
		for r = 1:repeats
			idx = randperm(rng, n_total)[1:n]
			scores[r], r2_scores[r] = _gcv_subset(dataset[idx], estimator, torque_weight)
		end
		keep = .!isnan.(scores)
		n_valid = count(keep)
		if n_valid == 0
			@warn "gcv_learning_curve: all $repeats draws at size $n were numerically " *
				"saturated (effective dof >= sample count); reporting NaN."
			means[i] = NaN
			stds[i] = NaN
			r2_means[i] = NaN
			r2_stds[i] = NaN
		else
			n_valid < repeats && @warn "gcv_learning_curve: " *
				"$(repeats - n_valid) of $repeats draws at size $n were " *
				"saturated and dropped."
			valid = scores[keep]
			valid_r2 = r2_scores[keep]
			means[i] = mean(valid)
			stds[i] = n_valid > 1 ? std(valid) : 0.0
			r2_means[i] = mean(valid_r2)
			r2_stds[i] = n_valid > 1 ? std(valid_r2) : 0.0
		end
	end
	return GCVSizeCurve(
		sorted_sizes, means, stds, r2_means, r2_stds, Int(repeats), Int(seed),
		estimator, Float64(torque_weight))
end

# GCV score and predictive R² for one (already-sliced) subset: assemble the
# weighted problem, solve, and score. A rank-deficient OLS solve (small draws)
# surfaces as the Ridge-pointing `ArgumentError` from `solve_coefficients`; we
# map it to `(NaN, NaN)` so the sweep continues rather than aborting.
function _gcv_subset(
	sub::SCEDataset,
	estimator::AbstractEstimator,
	torque_weight::Real,
)::Tuple{Float64, Float64}
	X, y, n_eff, intercept_dof = _gcv_core(sub, torque_weight)
	jvals = try
		Fitting.solve_coefficients(estimator, X, y)
	catch err
		err isa ArgumentError && return (NaN, NaN)
		rethrow()
	end
	residuals = y .- X * jvals
	score, _ = Fitting._gcv_single(
		X, y, residuals, estimator, jvals, n_eff, intercept_dof)
	return score, Fitting._gcv_r2(score, Fitting._gcv_msy(y, n_eff))
end

# Default size grid for `gcv_learning_curve`: six ascending points from a safe
# lower bound up to the full dataset. The lower bound `max(p + 2, 10)` keeps the
# design overdetermined (torque adds `3·n_atoms` rows per config); for datasets
# smaller than that bound the grid is limited to the available range, with a
# warning, since the small-size GCV may be unstable.
function _default_gcv_sizes(dataset::SCEDataset)::Vector{Int}
	n_total = length(dataset)
	p = size(dataset.X_E, 2)
	lo = max(p + 2, 10)
	if lo >= n_total
		@warn "gcv_learning_curve: dataset has only $n_total configs, below the safe " *
			"lower bound max(p+2, 10) = $lo for p = $p SALCs; the size grid is " *
			"limited and GCV may be unstable at the small end."
		lo = max(2, cld(n_total, 2))
	end
	pts = round.(Int, range(lo, n_total; length = 6))
	return sort(unique(clamp.(pts, 1, n_total)))
end

# argmin ignoring NaN entries; errors only if every entry is NaN.
function _argmin_ignore_nan(v::AbstractVector{<:Real})::Int
	best_i = 0
	best_v = Inf
	for (i, x) in enumerate(v)
		if !isnan(x) && x < best_v
			best_v = x
			best_i = i
		end
	end
	best_i == 0 && throw(ArgumentError(
		"_argmin_ignore_nan: every entry is NaN; cannot select a minimizer."))
	return best_i
end


# --- save / load --------------------------------------------------------

function _require_xml_extension(path::AbstractString)
	endswith(lowercase(path), ".xml") || throw(ArgumentError(
		"save / load handle XML only; the path must end in '.xml', got: $path"))
	return nothing
end

"""
	save(obj::SCEBasis, path::AbstractString)
	save(obj::SCEModel, path::AbstractString)
	save(f::SCEFit, path::AbstractString)

Write `obj` to an XML file at `path`. An `SCEBasis` is written as
structure, symmetry parameters, and SALC basis; an `SCEModel` adds the
fitted reference energy and SCE coefficients in a `<JPhi>` block. An
`SCEFit` is serialized as the corresponding `SCEModel(f)` — the
fit-time dataset and estimator are not persisted.

The path must end in `.xml`; any other extension is an error. Use
`load` to read the file back.

# Arguments
- `obj::Union{SCEBasis, SCEModel}`: The object to serialize.
- `f::SCEFit`: Trained fit; saved as `SCEModel(f)`.
- `path::AbstractString`: Output file path; must end in `.xml`.

# Returns
- `nothing`. The XML file is written as a side effect.

# Throws
- `ArgumentError` if `path` does not end in `.xml`.

# Examples
```julia
Magesty.save(basis, "basis.xml")
Magesty.save(model, "model.xml")
```
"""
function save(obj::SCEBasis, path::AbstractString)
	_require_xml_extension(path)
	XMLIO.write_basis_xml(
		obj.structure,
		obj.symmetry,
		obj.salcbasis,
		obj.isotropy,
		path,
	)
	return nothing
end

function save(obj::SCEModel, path::AbstractString)
	_require_xml_extension(path)
	XMLIO.write_model_xml(
		obj.basis.structure,
		obj.basis.symmetry,
		obj.basis.salcbasis,
		obj.basis.isotropy,
		obj.j0,
		obj.jphi,
		path,
	)
	return nothing
end

save(f::SCEFit, path::AbstractString) = save(SCEModel(f), path)

"""
	load(::Type{SCEBasis}, path::AbstractString) -> SCEBasis
	load(::Type{SCEModel}, path::AbstractString) -> SCEModel

Read an `SCEBasis` or `SCEModel` from the XML file at `path`. `structure`
is parsed from the file, `symmetry` is recomputed from `structure` and
the stored `tolerance_sym`, and the SALC basis is reconstructed from the
stored SALC data (the expensive SALC computation is skipped).

`load(SCEModel, path)` additionally reads the `<JPhi>` block;
`load(SCEBasis, path)` accepts an `SCEModel` XML as well and simply
ignores the `<JPhi>` block.

The path must end in `.xml`; any other extension is an error.

# Arguments
- `T::Type`: `SCEBasis` or `SCEModel`.
- `path::AbstractString`: Input file path; must end in `.xml`.

# Returns
- `SCEBasis` or `SCEModel`: the deserialized object, matching the
  requested type `T`.

# Throws
- `ArgumentError` if `path` does not end in `.xml`, if a required schema
  attribute is missing, or (for `SCEModel`) if the `<JPhi>` block is
  absent.

# Examples
```julia
basis = Magesty.load(SCEBasis, "basis.xml")
model = Magesty.load(SCEModel, "model.xml")
```
"""
function load(::Type{SCEBasis}, path::AbstractString)::SCEBasis
	_require_xml_extension(path)
	structure, symmetry, salcbasis, isotropy =
		XMLIO.read_basis_components_from_xml(path)
	return SCEBasis(structure, symmetry, salcbasis, isotropy)
end

function load(::Type{SCEModel}, path::AbstractString)::SCEModel
	_require_xml_extension(path)
	structure, symmetry, salcbasis, isotropy, j0, jphi =
		XMLIO.read_model_components_from_xml(path)
	return SCEModel(SCEBasis(structure, symmetry, salcbasis, isotropy), j0, jphi)
end

# --- Compact show methods for the main public types ---------------------
# Default Julia display unfolds every field; for SCEBasis / SCEDataset /
# SCEFit / SCEModel that means dumping the SALCBasis (thousands of entries),
# the design matrices, and the residuals. These compact forms keep REPL
# output usable.

function _print_fit_summary(f::SCEFit, elapsed_time::Real)
	n_E = length(f.dataset.y_E)
	n_T = length(f.dataset.y_T)
	# Evaluate the energy and torque predictions once each; the six
	# in-sample metrics below are all derived from these two vectors,
	# avoiding three redundant `X_E * jphi` / `X_T * jphi` products.
	obs_E, pred_E = _eval_energy(f, f.dataset)
	obs_T, pred_T = _eval_torque(f, f.dataset)
	rss_E = sum(abs2, obs_E .- pred_E)
	rss_T = sum(abs2, obs_T .- pred_T)
	rmse_E = Fitting._calc_rmse(obs_E, pred_E)
	rmse_T = Fitting._calc_rmse(obs_T, pred_T)
	r2_E = Fitting._calc_r2score(obs_E, pred_E)
	r2_T = Fitting._calc_r2score(obs_T, pred_T)
	println(
		"""

		FIT
		===
		""",
	)
	println("Estimator      : ", f.estimator)
	println("torque_weight  : ", f.torque_weight)
	println("n_configs (n_E): ", n_E)
	println("n_torques (n_T): ", n_T)
	println("num_coefs      : ", length(f.jphi), " (+ j0)")
	println("nonzero coefs  : ", count(!iszero, f.jphi), " / ", length(f.jphi))
	println(@sprintf("j0             : %+.6e eV", f.j0))
	println(@sprintf("RSS  (energy)  :  %.6e eV²    (Σ residuals², no 1/n_E)", rss_E))
	println(@sprintf("RSS  (torque)  :  %.6e eV²    (Σ residuals², no 1/n_T)", rss_T))
	println(@sprintf("RMSE (energy)  :  %.6e eV     (= √(RSS_E / n_E))", rmse_E))
	println(@sprintf("RMSE (torque)  :  %.6e eV     (= √(RSS_T / n_T))", rmse_T))
	println(@sprintf("R²   (energy)  :  %.6f", r2_E))
	println(@sprintf("R²   (torque)  :  %.6f", r2_T))
	println(@sprintf(" Time Elapsed: %.6f sec.", elapsed_time))
	println("-------------------------------------------------------------------")
end

function _summarize_jphi(io::IO, jphi::AbstractVector{<:Real})
	n = length(jphi)
	preview = view(jphi, 1:min(3, n))
	print(io, "[")
	for (i, v) in enumerate(preview)
		i == 1 || print(io, ", ")
		print(io, v)
	end
	n > 3 && print(io, ", ...")
	print(io, "] (", n, " elements)")
end

function Base.show(io::IO, basis::SCEBasis)
	print(io, "SCEBasis(num_atoms=", basis.structure.supercell.num_atoms,
		", num_salcs=", length(basis.salcbasis.salc_list),
		", isotropy=", basis.isotropy, ")")
end

function Base.show(io::IO, dataset::SCEDataset)
	print(io, "SCEDataset(num_configs=", length(dataset.spinconfigs),
		", num_atoms=", dataset.basis.structure.supercell.num_atoms, ")")
end

function Base.show(io::IO, f::SCEFit)
	print(io, "SCEFit(estimator=", f.estimator,
		", num_configs=", length(f.dataset.spinconfigs),
		", jphi=")
	_summarize_jphi(io, f.jphi)
	print(io, ", j0=", f.j0, ", torque_weight=", f.torque_weight, ")")
end

function Base.show(io::IO, m::SCEModel)
	print(io, "SCEModel(num_atoms=", m.basis.structure.supercell.num_atoms,
		", jphi=")
	_summarize_jphi(io, m.jphi)
	print(io, ", j0=", m.j0, ")")
end

# Fit-quality text writers (consumed by tools/FitCheck_*.py); included
# last as they depend on the predictor types and the predict_* verbs.
include("FitCheckIO.jl")

# VASP-to-extxyz conversion API; depends on the VaspIO / ExtXYZ modules.
include("VaspConvert.jl")

# MFA spin-sampling API; depends on the IncarIO / MfaSampling modules.
include("VaspSampling.jl")

# Sunny.jl LSWT script export; depends on the SCEModel / SALCBasis types and
# the Fitting normalization. Emits text only (no Sunny.jl dependency).
include("SunnyExport.jl")


end # module Magesty
