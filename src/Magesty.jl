"""
	Magesty

The main module of the Magesty package, providing the entry point for spin cluster expansion analysis and optimization.

This module provides the following main features:
- Magnetic structure setup and management
- Symmetry analysis
- Cluster expansion
- Basis function generation
- Spin configuration optimization

# Usage
```julia
using Magesty

# Load configuration from a TOML file
sc = SpinCluster("config.toml")

```

# Main Types
- `System`: Collection of structure, symmetry, cluster, and basis set
- `SpinCluster`: Extension of System with optimization capabilities

# Submodules
- `Structures`: Crystal structure processing
- `Symmetries`: Symmetry operations processing
- `Clusters`: Cluster expansion processing
- `SALCBases`: Basis function generation
- `Optimize`: Spin configuration optimization
"""
module Magesty

using Printf
using TOML
import AtomsBase
import StatsAPI                          # for the StatsAPI.RegressionModel supertype
import StatsAPI: fit, coef, nobs, dof    # extended with SCEFit / SCEModel methods

include("common/version.jl")
using .Version

include("common/SortedCounter.jl")
include("types/AtomCells.jl")
include("SpinConfigs.jl")
using .SpinConfigs

include("utils/SphericalHarmonicsTransforms.jl")
include("utils/AngularMomentumCoupling.jl")
include("types/Basis.jl")
using .AngularMomentumCoupling
using .Basis

include("utils/ConfigParser.jl")
include("utils/atomsbase_adapter.jl")
include("utils/RotationMatrix.jl")
include("utils/MySphericalHarmonics.jl")
using .ConfigParser
using .AtomsBaseAdapter

include("Structures.jl")
include("Symmetries.jl")
include("Clusters.jl")
include("SALCBases.jl")
include("Optimize.jl")

using .Structures
using .Symmetries
using .Clusters
using .SALCBases
using .Optimize
# extended in Magesty with SCEFit / SpinConfig / SCEDataset methods
import .Optimize: SCEModel, predict_energy, predict_torque

include("utils/xml_io.jl")
using .XMLIO
include("utils/EnergyTorque.jl")
using .EnergyTorque

export System, SpinCluster, SCEBasis, SCEDataset, SCEFit, VERSION, install_tools
export SCEModel, fit_sce_model, predict_energy, predict_torque, AbstractEstimator, OLS, Ridge
export fit, coef, intercept, nobs, dof
export r2_energy, r2_torque, rss_energy, rss_torque
export residuals_energy, residuals_torque, rmse_energy, rmse_torque
export build_sce_basis, build_sce_basis_from_xml
export write_xml

# Re-export read_embset from SpinConfigs for user convenience
const read_embset = SpinConfigs.read_embset
export read_embset

"""
	System

A collection of structure, symmetry, cluster, and basis set.

# Fields
- `config::Parser`: Configuration parser
- `structure::Structure`: Crystal structure information
- `symmetry::Symmetry`: Symmetry operations
- `cluster::Cluster`: Cluster information
- `basisset::SALCBasis`: Basis set information
"""
struct System
	structure::Structure
	symmetry::Symmetry
	cluster::Cluster
	basisset::SALCBasis
end

# Shared skeleton for all input-driven constructors below.
# Returns the (structure, symmetry, cluster) triplet; callers append
# the SALCBasis — either computed via `SALCBasis(...)` or loaded from XML.
function _build_structure_skeleton(
	config::Config4System;
	verbosity::Bool = true,
)
	structure::Structure = Structure(config, verbosity = verbosity)
	symmetry::Symmetry = Symmetry(structure, config, verbosity = verbosity)
	cluster::Cluster = Cluster(structure, symmetry, config, verbosity = verbosity)
	return structure, symmetry, cluster
end

"""
	System

Create a `System` instance from either a dictionary of input parameters or a TOML configuration file.

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing input parameters.
- `toml_file::AbstractString`: Path to the TOML configuration file.
- `verbosity::Bool=true`: Whether to print detailed information during initialization.

# Returns
- `System`: A new `System` instance containing structure, symmetry, cluster, and basis set.

# Throws
- `SystemError`: If the TOML file cannot be read.
- `ErrorException`: If required parameters are missing, invalid, or the TOML parsing fails.

# Examples
```julia
# Create a System from a dictionary
input_dict = Dict("key" => "value")
system = System(input_dict)

# Create a System from a TOML file
system = System("config.toml")
```
"""
function System(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
	config::Config4System = Config4System(input_dict)
	structure, symmetry, cluster = _build_structure_skeleton(config; verbosity = verbosity)
	basisset::SALCBasis = SALCBasis(structure, symmetry, cluster, config, verbosity = verbosity)
	return System(structure, symmetry, cluster, basisset)
end

function System(toml_file::AbstractString; verbosity::Bool = true)
	try
		open(toml_file) do io
			toml = read(io, String)
			input_dict = TOML.parse(toml)
			return System(input_dict; verbosity = verbosity)
		end
	catch e
		if isa(e, SystemError)
			throw(SystemError("Failed to read file: $toml_file"))
		else
			throw(ErrorException("Failed to parse TOML file: $toml_file"))
		end
	end
end


"""
	build_sce_basis(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool=false) -> System

Build a `System` (structure + symmetry + cluster + basis set) from a parsed input
dictionary. This is the headless alternative to `System(input_dict)` and skips the
header banner by default; the SALC basis is computed from scratch (use
`build_sce_basis_from_xml` to load a precomputed basis instead).

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing the parsed
  TOML input parameters (typically the result of `TOML.parsefile`).
- `verbosity::Bool=false`: Whether to print detailed progress information.

# Returns
- `System`: A `System` instance containing structure, symmetry, cluster, and basis set.

# Throws
- `ErrorException`: If required parameters are missing or invalid.

# Examples
```julia
using TOML
input_dict = TOML.parsefile("input.toml")
system = build_sce_basis(input_dict)
```
"""
function build_sce_basis(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = false)::System
	config::Config4System = Config4System(input_dict)
	structure, symmetry, cluster = _build_structure_skeleton(config; verbosity = verbosity)
	basisset::SALCBasis = SALCBasis(structure, symmetry, cluster, config, verbosity = verbosity)
	return System(structure, symmetry, cluster, basisset)
end

"""
	build_sce_basis_from_xml(input_dict::Dict{<:AbstractString, <:Any}, xml_file::AbstractString; verbosity::Bool = true)::System

Build System from input.toml dictionary and XML file. This function constructs structure, symmetry, and cluster
from the input dictionary, but loads the basis set from the XML file to avoid expensive SALC computation.

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing input parameters (parsed from input.toml)
- `xml_file::AbstractString`: Path to XML file containing basis set information
- `verbosity::Bool=true`: Whether to print detailed information during initialization

# Returns
- `System`: A new `System` instance with basis set loaded from XML file

# Throws
- `ErrorException` if required parameters are missing or XML file format is invalid

# Examples
```julia
using TOML
input_dict = TOML.parsefile("input.toml")
system = build_sce_basis_from_xml(input_dict, "scecoeffs.xml")
```
"""
function build_sce_basis_from_xml(
	input_dict::Dict{<:AbstractString, <:Any},
	xml_file::AbstractString;
	verbosity::Bool = true,
)::System
	config::Config4System = Config4System(input_dict)
	structure, symmetry, cluster = _build_structure_skeleton(config; verbosity = verbosity)

	# Load basis set from XML file instead of computing it
	if verbosity
		println("Loading basis set from XML file: $xml_file")
	end
	basisset::SALCBasis = XMLIO.read_salcbasis_from_xml(xml_file)
	if verbosity
		println("Successfully loaded basis set from XML file")
	end

	return System(structure, symmetry, cluster, basisset)
end


"""
	SCEBasis

Material + basis: crystal structure, symmetry, cluster, and SALC basis.
The heavy part is the SALC construction; an `SCEBasis` can be persisted
and reused.

# Fields
- `structure::Structure`: Crystal structure information.
- `symmetry::Symmetry`: Symmetry operations.
- `cluster::Cluster`: Cluster information.
- `salcbasis::SALCBasis`: SALC basis functions.
"""
struct SCEBasis
	structure::Structure
	symmetry::Symmetry
	cluster::Cluster
	salcbasis::SALCBasis
end

"""
	SCEBasis(input_dict::AbstractDict; verbosity = true) -> SCEBasis
	SCEBasis(toml_path::AbstractString; verbosity = true) -> SCEBasis
	SCEBasis(system::AtomsBase.AbstractSystem; interaction, name, tolerance_sym, isotropy, verbosity) -> SCEBasis
	SCEBasis(; lattice, kd, kd_list, positions, periodicity, interaction, name, tolerance_sym, isotropy, verbosity) -> SCEBasis

Construct an `SCEBasis` from one of four input paths: a parsed TOML
dictionary, a TOML file, an `AtomsBase.AbstractSystem`, or raw Julia
keyword arguments. All paths build a `Config4System` internally and run
the same structure / symmetry / cluster / SALC construction.

The `interaction` argument is a nested `NamedTuple` keyed `body1`,
`body2`, ...; see `docs/specs/260514-sce-public-api/` for its format
and the per-path Unitful requirements.
"""
function SCEBasis(
	input_dict::AbstractDict{<:AbstractString, <:Any};
	verbosity::Bool = true,
)
	config::Config4System = Config4System(input_dict)
	structure, symmetry, cluster = _build_structure_skeleton(config; verbosity = verbosity)
	salcbasis::SALCBasis =
		SALCBasis(structure, symmetry, cluster, config, verbosity = verbosity)
	return SCEBasis(structure, symmetry, cluster, salcbasis)
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
		else
			throw(ErrorException("Failed to parse TOML file: $toml_path"))
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
	input_dict = AtomsBaseAdapter.system_to_input_dict(
		system,
		interaction;
		name = name,
		tolerance_sym = tolerance_sym,
		isotropy = isotropy,
	)
	return SCEBasis(input_dict; verbosity = verbosity)
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
	input_dict = AtomsBaseAdapter.kwargs_to_input_dict(
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
	return SCEBasis(input_dict; verbosity = verbosity)
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
- `X_E::Matrix{Float64}`: Energy design matrix, bias column at column 1,
  one row per configuration. Unweighted.
- `X_T::Matrix{Float64}`: Torque design matrix, no bias column, a
  `3 * num_atoms` block of rows per configuration. Unweighted.
- `y_E::Vector{Float64}`: Observed energies, length `n_configs`.
- `y_T::Vector{Float64}`: Observed torques, flattened, length
  `3 * num_atoms * n_configs`.
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
	SCEDataset(system::AtomsBase.AbstractSystem, spinconfigs; interaction, ...) -> SCEDataset
	SCEDataset(toml_path::AbstractString, spinconfigs::AbstractVector{SpinConfig})

Build a `SCEDataset` from a basis and training data. The base method
takes an explicit `SCEBasis` and a vector of `SpinConfig`; it builds the
unweighted energy and torque design matrices once. The `embset_path`
method reads the configurations from an EMBSET file first. The two sugar
methods build a throwaway `SCEBasis` internally (from an
`AtomsBase.AbstractSystem` or a TOML file) and embed it in the dataset;
for workflows that share one basis across several datasets, construct the
`SCEBasis` explicitly and pass it to the base method.
"""
function SCEDataset(basis::SCEBasis, spinconfigs::AbstractVector{SpinConfig})
	X_E = Optimize.build_design_matrix_energy(
		basis.salcbasis.salc_list,
		spinconfigs,
		basis.symmetry,
	)
	X_T = Optimize.build_design_matrix_torque(
		basis.salcbasis.salc_list,
		spinconfigs,
		basis.structure.supercell.num_atoms,
		basis.symmetry,
	)
	y_E::Vector{Float64} = [sc.energy for sc in spinconfigs]
	y_T::Vector{Float64} =
		isempty(spinconfigs) ? Float64[] :
		reduce(vcat, (vec(sc.torques) for sc in spinconfigs))
	return SCEDataset(basis, collect(SpinConfig, spinconfigs), X_E, X_T, y_E, y_T)
end

function SCEDataset(basis::SCEBasis, embset_path::AbstractString)
	return SCEDataset(basis, SpinConfigs.read_embset(embset_path))
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
	return SCEDataset(basis, spinconfigs)
end

function SCEDataset(
	toml_path::AbstractString,
	spinconfigs::AbstractVector{SpinConfig};
	verbosity::Bool = true,
)
	return SCEDataset(SCEBasis(toml_path; verbosity = verbosity), spinconfigs)
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

function Base.vcat(d1::SCEDataset, ds::SCEDataset...)::SCEDataset
	for d in ds
		d.basis === d1.basis || throw(ArgumentError(
			"vcat requires all SCEDataset arguments to share the same SCEBasis"))
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
coefficients, the estimator and torque weight used, the residuals of the
augmented (weighted) least-squares system, and cached in-sample metrics.

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
  use `metrics` for interpretable in-sample errors.
- `metrics::Dict{Symbol, Any}`: In-sample RMSE / R² for energy and torque
  (`:rmse_energy`, `:rmse_torque`, `:r2score_energy`, `:r2score_torque`),
  computed on the unweighted predictions.
"""
struct SCEFit <: StatsAPI.RegressionModel
	dataset::SCEDataset
	j0::Float64
	jphi::Vector{Float64}
	estimator::AbstractEstimator
	torque_weight::Float64
	residuals::Vector{Float64}
	metrics::Dict{Symbol, Any}
end

"""
	fit(::Type{SCEFit}, dataset::SCEDataset, estimator::AbstractEstimator;
	    torque_weight::Real = 0.5) -> SCEFit

Fit SCE coefficients on `dataset` with `estimator`, returning a `SCEFit`.

`torque_weight` in `[0, 1]` sets the convex combination of the per-sample
energy and torque mean squared errors that the augmented least-squares
problem minimizes: `(1 - torque_weight) * MSE_energy + torque_weight *
MSE_torque`. The design matrices stored in `dataset` are unweighted, so a
`torque_weight` sweep reuses one `SCEDataset` without rebuilding them.

This is the StatsAPI `fit` verb; `using Magesty` re-exports it.
"""
function fit(
	::Type{SCEFit},
	dataset::SCEDataset,
	estimator::AbstractEstimator;
	torque_weight::Real = 0.5,
)::SCEFit
	X, y, bias_col = Optimize.assemble_weighted_problem(
		dataset.X_E,
		dataset.X_T,
		dataset.y_E,
		dataset.y_T,
		torque_weight,
	)
	j_values = Optimize.solve_coefficients(estimator, X, y; bias_col = bias_col)
	j0, jphi = Optimize.extract_j0_jphi(j_values, dataset.X_E, dataset.y_E)
	residuals::Vector{Float64} = y .- X * j_values
	predicted_energy::Vector{Float64} = dataset.X_E[:, 2:end] * jphi .+ j0
	predicted_torque::Vector{Float64} = dataset.X_T * jphi
	metrics = Dict{Symbol, Any}(
		:rmse_energy => Optimize.calc_rmse(dataset.y_E, predicted_energy),
		:rmse_torque => Optimize.calc_rmse(dataset.y_T, predicted_torque),
		:r2score_energy => Optimize.calc_r2score(dataset.y_E, predicted_energy),
		:r2score_torque => Optimize.calc_r2score(dataset.y_T, predicted_torque),
	)
	return SCEFit(
		dataset,
		j0,
		jphi,
		estimator,
		Float64(torque_weight),
		residuals,
		metrics,
	)
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
	intercept(f::SCEFit) -> Float64
	intercept(m::SCEModel) -> Float64

The fitted reference energy `j0` (bias term). `intercept` is a
Magesty-native verb: StatsAPI has no intercept concept, and the SCE model
keeps `j0` separate from the coefficients `jphi`.
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

Extract a lightweight `SCEModel` from a fitted `SCEFit`: the fitted
`(j0, jphi)` together with the SALC basis, symmetry, and atom count
needed for prediction. The training dataset, residuals, and cached
metrics held by the `SCEFit` are dropped.
"""
function SCEModel(f::SCEFit)::SCEModel
	return SCEModel(
		f.j0,
		f.jphi,
		f.dataset.basis.salcbasis,
		f.dataset.basis.symmetry,
		f.dataset.basis.structure.supercell.num_atoms,
	)
end


# --- Prediction ---------------------------------------------------------

"""
	predict_energy(model::SCEModel, sc::SpinConfig) -> Float64
	predict_energy(f::SCEFit, spin_directions::AbstractMatrix{<:Real}) -> Float64
	predict_energy(f::SCEFit, sc::SpinConfig) -> Float64
	predict_energy(predictor, dataset::SCEDataset) -> Vector{Float64}

Predict SCE energies. The single-configuration forms return one
`Float64`; the `SCEDataset` form returns one energy per configuration.
`SCEFit` predictors delegate through `SCEModel(f)`. The `SCEDataset` form
requires the dataset to share the predictor's SALC basis.
"""
predict_energy(model::SCEModel, sc::SpinConfig)::Float64 =
	predict_energy(model, sc.spin_directions)
predict_energy(f::SCEFit, spin_directions::AbstractMatrix{<:Real})::Float64 =
	predict_energy(SCEModel(f), spin_directions)
predict_energy(f::SCEFit, sc::SpinConfig)::Float64 =
	predict_energy(SCEModel(f), sc)

function predict_energy(model::SCEModel, dataset::SCEDataset)::Vector{Float64}
	_check_basis(model, dataset)
	return dataset.X_E[:, 2:end] * model.jphi .+ model.j0
end
predict_energy(f::SCEFit, dataset::SCEDataset)::Vector{Float64} =
	predict_energy(SCEModel(f), dataset)

"""
	predict_torque(model::SCEModel, sc::SpinConfig) -> Matrix{Float64}
	predict_torque(f::SCEFit, spin_directions::AbstractMatrix{<:Real}) -> Matrix{Float64}
	predict_torque(f::SCEFit, sc::SpinConfig) -> Matrix{Float64}
	predict_torque(predictor, dataset::SCEDataset) -> Vector{Matrix{Float64}}

Predict per-atom SCE torques. The single-configuration forms return one
`3×num_atoms` matrix; the `SCEDataset` form returns one such matrix per
configuration. `SCEFit` predictors delegate through `SCEModel(f)`. The
`SCEDataset` form requires the dataset to share the predictor's SALC
basis.
"""
predict_torque(model::SCEModel, sc::SpinConfig)::Matrix{Float64} =
	predict_torque(model, sc.spin_directions)
predict_torque(f::SCEFit, spin_directions::AbstractMatrix{<:Real})::Matrix{Float64} =
	predict_torque(SCEModel(f), spin_directions)
predict_torque(f::SCEFit, sc::SpinConfig)::Matrix{Float64} =
	predict_torque(SCEModel(f), sc)

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
	model.salcbasis === dataset.basis.salcbasis || throw(ArgumentError(
		"evaluation SCEDataset was built from a different SCEBasis than the predictor"))
	return nothing
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
`predictor` is a `SCEModel` or `SCEFit`; `data` is a `SCEDataset`, a
vector of `SpinConfig`, or an EMBSET file path. The single-argument form
evaluates `f` in-sample on its own training dataset.
"""
function r2_energy(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_energy(predictor, data)
	return Optimize.calc_r2score(observed, predicted)
end
r2_energy(f::SCEFit)::Float64 = r2_energy(f, f.dataset)

"""
	r2_torque(predictor, data) -> Float64
	r2_torque(f::SCEFit) -> Float64

Coefficient of determination (R²) of the SCE torque predictions,
computed over the flattened torque components. See [`r2_energy`](@ref)
for the argument forms.
"""
function r2_torque(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_torque(predictor, data)
	return Optimize.calc_r2score(observed, predicted)
end
r2_torque(f::SCEFit)::Float64 = r2_torque(f, f.dataset)

"""
	rss_energy(predictor, data) -> Float64
	rss_energy(f::SCEFit) -> Float64

Residual sum of squares of the SCE energy predictions. See
[`r2_energy`](@ref) for the argument forms.
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
flattened torque components. See [`r2_energy`](@ref) for the argument
forms.
"""
function rss_torque(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_torque(predictor, data)
	return sum(abs2, observed .- predicted)
end
rss_torque(f::SCEFit)::Float64 = rss_torque(f, f.dataset)

"""
	residuals_energy(predictor, data) -> Vector{Float64}
	residuals_energy(f::SCEFit) -> Vector{Float64}

Per-configuration energy residuals `observed - predicted` (eV). See
[`r2_energy`](@ref) for the argument forms.
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

Flattened torque residuals `observed - predicted` (eV), length
`3 * num_atoms * n_configs`. See [`r2_energy`](@ref) for the argument
forms.
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

Root mean squared error of the SCE energy predictions (eV). See
[`r2_energy`](@ref) for the argument forms.
"""
function rmse_energy(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_energy(predictor, data)
	return Optimize.calc_rmse(observed, predicted)
end
rmse_energy(f::SCEFit)::Float64 = rmse_energy(f, f.dataset)

"""
	rmse_torque(predictor, data) -> Float64
	rmse_torque(f::SCEFit) -> Float64

Root mean squared error of the SCE torque predictions (eV), over the
flattened torque components. See [`r2_energy`](@ref) for the argument
forms.
"""
function rmse_torque(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData)::Float64
	observed, predicted = _eval_torque(predictor, data)
	return Optimize.calc_rmse(observed, predicted)
end
rmse_torque(f::SCEFit)::Float64 = rmse_torque(f, f.dataset)


"""
	SpinCluster

An extension of System with optimization capabilities.

# Fields
- `config::Parser`: Configuration parser
- `structure::Structure`: Crystal structure information
- `symmetry::Symmetry`: Symmetry operations
- `cluster::Cluster`: Cluster information
- `basisset::SALCBasis`: Basis set information
- `optimize::Optimizer`: Optimizer instance
"""
struct SpinCluster
	structure::Structure
	symmetry::Symmetry
	cluster::Cluster
	basisset::SALCBasis
	optimize::Optimizer
end

"""
	SpinCluster

Create a `SpinCluster` instance from either a dictionary of input parameters, a TOML configuration file, or an existing `System` instance. This is an extension of `System` that includes optimization capabilities.

# Arguments
- `input_dict::Dict{<:AbstractString, <:Any}`: Dictionary containing input parameters.
- `toml_file::AbstractString`: Path to the TOML configuration file.
- `system::System`: An existing `System` instance.
- `verbosity::Bool=true`: Whether to print detailed information during initialization.

# Returns
- `SpinCluster`: A new `SpinCluster` instance containing structure, symmetry, cluster, basis set, and optimizer.

# Throws
- `SystemError`: If the TOML file cannot be read.
- `ErrorException`: If required parameters are missing, invalid, or the TOML parsing fails.

# Examples
```julia
# Create a SpinCluster from a dictionary
input_dict = Dict("key" => "value")
spin_cluster = SpinCluster(input_dict)

# Create a SpinCluster from a TOML file
spin_cluster = SpinCluster("config.toml")

# Create a SpinCluster from an existing System
system = System("config.toml")
spin_cluster = SpinCluster(system)
```
"""
function SpinCluster(input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)
	config_system::Config4System = Config4System(input_dict)
	structure, symmetry, cluster = _build_structure_skeleton(config_system; verbosity = verbosity)
	basisset::SALCBasis =
		SALCBasis(structure, symmetry, cluster, config_system, verbosity = verbosity)

	config_optimize::Config4Optimize = Config4Optimize(input_dict)
	optimize::Optimizer =
		Optimizer(structure, symmetry, basisset, config_optimize, verbosity = verbosity)

	return SpinCluster(structure, symmetry, cluster, basisset, optimize)
end

function SpinCluster(toml_file::AbstractString; verbosity::Bool = true)
	try
		open(toml_file) do io
			toml = read(io, String)
			input_dict = TOML.parse(toml)
			return SpinCluster(input_dict, verbosity = verbosity)
		end
	catch e
		if isa(e, SystemError)
			throw(SystemError("Failed to read file: $toml_file"))
		else
			throw(ErrorException("Failed to parse TOML file: $toml_file"))
		end
	end
end

"""
	SpinCluster(system::System, input_dict::Dict{<:AbstractString, <:Any}; verbosity::Bool = true)

Creates a `SpinCluster` instance by extending an existing `System` object with optimization capabilities. 
This constructor uses a dictionary of input parameters to configure the optimization process.

# Arguments
- `system::System`: An existing `System` instance containing structure, symmetry, cluster, and basis set information.
- `input_dict::Dict{<:AbstractString, <:Any}`: A dictionary containing input parameters for optimization.
- `verbosity::Bool=true`: Whether to print detailed information during initialization.

# Returns
- `SpinCluster`: A new `SpinCluster` instance containing structure, symmetry, cluster, basis set, and optimizer.

# Throws
- `ErrorException`: If required parameters are missing or invalid.

# Example
```julia
input_dict = Dict("key" => "value")
system = System(input_dict)
spin_cluster = SpinCluster(system, input_dict)
```
"""
function SpinCluster(
	system::System,
	input_dict::Dict{<:AbstractString, <:Any},
	;
	verbosity::Bool = true,
)
	config::Config4Optimize = Config4Optimize(input_dict)
	optimize =
		Optimizer(system.structure, system.symmetry, system.basisset, config, verbosity = verbosity)
	return SpinCluster(
		system.structure,
		system.symmetry,
		system.cluster,
		system.basisset,
		optimize,
	)
end

"""
	SpinCluster(system::System,
	            input_dict::AbstractDict{<:AbstractString, <:Any},
	            spinconfig_list::AbstractVector{SpinConfig};
	            verbosity::Bool=true)

Create a `SpinCluster` from an existing `System`, an input dictionary, and a
caller-supplied list of `SpinConfig` objects. Use this overload when the spin
configurations come from somewhere other than the EMBSET path declared in
`input_dict` (for example, configurations generated programmatically or pulled
from a non-standard format).

# Arguments
- `system::System`: An existing `System` instance.
- `input_dict::AbstractDict{<:AbstractString, <:Any}`: Dictionary with
  optimization parameters (`[regression]` section).
- `spinconfig_list::AbstractVector{SpinConfig}`: Training spin configurations
  used in place of those that would be loaded from disk.
- `verbosity::Bool=true`: Whether to print detailed information during fitting.

# Returns
- `SpinCluster`: A fitted `SpinCluster` instance.

# Throws
- `ErrorException`: If required `[regression]` parameters are missing or invalid.

# Examples
```julia
system = System("input.toml")
configs = [SpinConfig(...), SpinConfig(...)]
sc = SpinCluster(system, input_dict, configs)
```
"""
function SpinCluster(
	system::System,
	input_dict::AbstractDict{<:AbstractString, <:Any},
	spinconfig_list::AbstractVector{SpinConfig},
	;
	verbosity::Bool = true,
)
	config::Config4Optimize = Config4Optimize(input_dict)
	optimize = Optimizer(
		system.structure,
		system.symmetry,
		system.basisset,
		config.alpha,
		config.lambda,
		config.weight,
		spinconfig_list,
		verbosity = verbosity,
	)
	return SpinCluster(
		system.structure,
		system.symmetry,
		system.cluster,
		system.basisset,
		optimize,
	)

end


"""
	calc_energy(sc::SpinCluster, spin_config::AbstractMatrix{<:Real})

Calculate the energy of a spin configuration using the spin cluster expansion.

# Arguments
- `sc::SpinCluster`: A `SpinCluster` instance containing structure, symmetry, basis set, and optimization information.
- `spin_config::AbstractMatrix{<:Real}`: A 3xN matrix representing the spin configuration, where N is the number of atoms in the supercell.

# Returns
- `Float64`: The calculated energy of the given spin configuration.

# Throws
- `ArgumentError`: If the number of columns in `spin_config` does not match the number of atoms in the supercell.

# Example
```julia
spin_config = rand(3, sc.structure.supercell.num_atoms) # Random spin configuration
energy = calc_energy(sc, spin_config)
```
"""
function calc_energy(spincluster::SpinCluster, spin_config::AbstractMatrix{<:Real})::Float64
	if spincluster.structure.supercell.num_atoms != size(spin_config, 2)
		num_atoms = spincluster.structure.supercell.num_atoms
		throw(
			ArgumentError(
				"spin_config must be 3xN matrix where N is the number of atoms in the supercell. $num_atoms",
			),
		)
	end
	return EnergyTorque.calc_energy(
		spincluster.basisset.salc_list,
		spin_config,
		spincluster.symmetry,
		spincluster.optimize,
	)
end

"""
	calc_torque(sc::SpinCluster, spin_config::AbstractMatrix{<:Real})

Calculate the torque (local magnetic field) for each atom in a spin configuration using the spin cluster expansion.

# Arguments
- `sc::SpinCluster`: A `SpinCluster` instance containing structure, symmetry, basis set, and optimization information.
- `spin_config::AbstractMatrix{<:Real}`: A 3xN matrix representing the spin configuration, where N is the number of atoms in the supercell.

# Returns
- `Matrix{Float64}`: A 3xN matrix representing the torque for each atom, where each column corresponds to the torque vector (in eV) for one atom.

# Throws
- `ArgumentError`: If the number of columns in `spin_config` does not match the number of atoms in the supercell.

# Example
```julia
spin_config = rand(3, sc.structure.supercell.num_atoms) # Random spin configuration
torque = calc_torque(sc, spin_config)
```
"""
function calc_torque(spincluster::SpinCluster, spin_config::AbstractMatrix{<:Real})::Matrix{Float64}
	if spincluster.structure.supercell.num_atoms != size(spin_config, 2)
		num_atoms = spincluster.structure.supercell.num_atoms
		throw(
			ArgumentError(
				"spin_config must be 3xN matrix where N is the number of atoms in the supercell. $num_atoms",
			),
		)
	end
	return EnergyTorque.calc_torque(
		spincluster.basisset.salc_list,
		spin_config,
		spincluster.symmetry,
		spincluster.optimize,
	)
end

"""
	write_xml(sc::SpinCluster, filename::AbstractString="jphi.xml"; write_jphi::Bool=true)

Write the structure, symmetry, basis set, and fitted SCE coefficients
held by a `SpinCluster` to an XML file in SCE format.

# Arguments
- `sc::SpinCluster`: Spin cluster (structure + symmetry + basis set + optimizer)
- `filename::AbstractString="jphi.xml"`: Output XML file name
- `write_jphi::Bool=true`: Whether to embed the J_ij parameters
  (set to `false` to write only the structure/symmetry/basis section)

# Examples
```julia
write_xml(spin_cluster)
write_xml(spin_cluster, "output.xml", write_jphi=false)
```

To export a `System` (no fitted coefficients), use the
`write_xml(::System, ...)` overload instead. To export from
standalone components (`structure`, `symmetry`, `salc_basis`,
`optimizer`), wrap them with `SpinCluster(structure, symmetry,
cluster, salc_basis, optimizer)` first.
"""
function write_xml(
	sc::SpinCluster,
	filename::AbstractString = "jphi.xml";
	write_jphi::Bool = true,
)
	XMLIO.write_xml(
		sc.structure,
		sc.symmetry,
		sc.basisset,
		sc.optimize,
		filename;
		write_jphi = write_jphi,
	)
end

"""
	write_xml(system::System, filename::AbstractString="jphi.xml")

Write System information to an XML file. This saves the structure, symmetry, and basis set
information without optimization results (JPhi).

# Arguments
- `system::System`: The System to write
- `filename::AbstractString`: Output XML file name (default: "jphi.xml")

# Examples
```julia
system = build_sce_basis(input_dict)
write_xml(system, "system.xml")
```
"""
function write_xml(
	system::System,
	filename::AbstractString = "jphi.xml",
)
	XMLIO.write_xml(system.structure, system.symmetry, system.basisset, filename)
end

"""
	write_energies(sc::SpinCluster, filename::AbstractString="energy_list.txt")

Write the observed (DFT) and predicted (SCE) energies to a text file.

# Arguments
- `sc::SpinCluster`: A `SpinCluster` instance containing optimization results.
- `filename::AbstractString="energy_list.txt"`: Output file name.

# Output Format
The file contains:
- Header line: `# data index,    DFT_Energy,    SCE_Energy`
- Data lines: index, observed energy (eV), predicted energy (eV)

# Example
```julia
write_energies(sc, "my_energies.txt")
```
"""
function write_energies(
	sc::SpinCluster,
	filename::AbstractString = "energy_list.txt",
)

	observed_energy_list = [spinconfig.energy for spinconfig in sc.optimize.spinconfig_list]
	predicted_energy_list = sc.optimize.predicted_energy_list
	# Write to file
	try
		open(filename, "w") do f
			# Write header
			println(f, "# data index,    DFT_Energy,    SCE_Energy\n# unit of energy is eV")

			# Write data
			idx_width = ndigits(length(observed_energy_list))
			for i in eachindex(sc.optimize.spinconfig_list)
				str = @sprintf(
					" %*d    % 15.10e    % 15.10e\n",
					idx_width,
					i,
					observed_energy_list[i],
					predicted_energy_list[i],
				)
				write(f, str)
			end
		end
	catch e
		@error "Failed to write lists to file" exception = (e, catch_backtrace())
		rethrow(e)
	end
end

"""
	write_torques(sc::SpinCluster, filename::AbstractString="torque_list.txt")

Write the observed (DFT) and predicted (SCE) torques for each atom to a text file.

# Arguments
- `sc::SpinCluster`: A `SpinCluster` instance containing optimization results.
- `filename::AbstractString="torque_list.txt"`: Output file name.

# Output Format
The file contains:
- Header line: `# atom index,    element,   DFT_torque_x,    DFT_torque_y,    DFT_torque_z,    SCE_torque_x,    SCE_torque_y,    SCE_torque_z`
- Data lines: atom index, element symbol, observed torque components (eV), predicted torque components (eV)
- Data is grouped by configuration index

# Example
```julia
write_torques(sc, "my_torques.txt")
```
"""
function write_torques(
	sc::SpinCluster,
	filename::AbstractString = "torque_list.txt",
)
	predicted_torque_list::Vector{Matrix{Float64}} = sc.optimize.predicted_torque_list
	observed_torque_list::Vector{Matrix{Float64}} =
		[spinconfig.torques for spinconfig in sc.optimize.spinconfig_list]

	# Write to file
	open(filename, "w") do f
		# Write header
		println(
			f,
			"# atom index,    element,   DFT_torque_x,    DFT_torque_y,    DFT_torque_z,    SCE_torque_x,    SCE_torque_y,    SCE_torque_z\n# unit of torque is eV",
		)

		# Write data
		idx_width = ndigits(length(sc.optimize.spinconfig_list))
		element_string_list =
			[sc.structure.kd_name[elm_idx] for elm_idx in sc.structure.supercell.kd_int_list]
		element_width = maximum(length.(element_string_list))

		for (ndata, (obs_torque_matrix, pred_torque_matrix)) in
			enumerate(zip(observed_torque_list, predicted_torque_list))
			println(f, "# data index: $ndata")
			for (iatom, (obs_torque, pred_torque)) in
				enumerate(zip(eachcol(obs_torque_matrix), eachcol(pred_torque_matrix)))
				# obs_torque and pred_torque are length-3 vectors (x, y, z)
				str = @sprintf(
					" %*d %*s  % 15.10e   % 15.10e   % 15.10e    % 15.10e   % 15.10e   % 15.10e\n",
					idx_width,
					iatom,
					element_width,
					element_string_list[iatom],
					obs_torque[1],
					obs_torque[2],
					obs_torque[3],
					pred_torque[1],
					pred_torque[2],
					pred_torque[3],
				)
				write(f, str)
			end
		end
	end
end

"""
	get_j0(sc::SpinCluster) -> Float64

Return the fitted reference energy (bias term `j0`) from a `SpinCluster`.

# Arguments
- `sc::SpinCluster`: Spin cluster with a fitted optimizer.

# Returns
- `Float64`: Reference energy in eV (the constant offset of the SCE expansion).

# Examples
```julia
sc = SpinCluster("input.toml")
j0 = get_j0(sc)
```
"""
function get_j0(sc::SpinCluster)::Float64
	return sc.optimize.reference_energy
end

"""
	get_jphi(sc::SpinCluster) -> Vector{Float64}

Return the fitted SCE coefficients (`jphi`) from a `SpinCluster`.

# Arguments
- `sc::SpinCluster`: Spin cluster with a fitted optimizer.

# Returns
- `Vector{Float64}`: SCE coefficients, one entry per SALC basis function.

# Examples
```julia
sc = SpinCluster("input.toml")
jphi = get_jphi(sc)
```
"""
function get_jphi(sc::SpinCluster)::Vector{Float64}
	return sc.optimize.SCE
end

"""
	get_j0_jphi(sc::SpinCluster) -> Tuple{Float64, Vector{Float64}}

Return the reference energy and SCE coefficients as a tuple.

# Arguments
- `sc::SpinCluster`: Spin cluster with a fitted optimizer.

# Returns
- `Tuple{Float64, Vector{Float64}}`: `(j0, jphi)` where `j0` is the bias term in eV
  and `jphi` is the vector of SCE coefficients.

# Examples
```julia
sc = SpinCluster("input.toml")
j0, jphi = get_j0_jphi(sc)
```
"""
function get_j0_jphi(sc::SpinCluster)::Tuple{Float64, Vector{Float64}}
	return sc.optimize.reference_energy, sc.optimize.SCE
end

"""
    install_tools(; bindir::AbstractString=joinpath(homedir(), ".julia", "bin"))

Install CLI wrapper scripts for Magesty tools into `bindir` (default: `~/.julia/bin`).
Each wrapper calls `julia /path/to/script.jl "\$@"`, so the correct package version
is always used regardless of where the package is installed.

After running this once, add `~/.julia/bin` to your PATH:
    export PATH="\$HOME/.julia/bin:\$PATH"

# Arguments
- `bindir::AbstractString=joinpath(homedir(), ".julia", "bin")`: Target directory
  for the wrapper scripts. The directory is created if it does not exist.

# Returns
- `Nothing`: Writes wrapper scripts as a side effect and prints their paths.

# Throws
- `ErrorException`: If the Magesty package directory cannot be determined.

# Available commands after installation
- `vasp2extxyz`           — convert a single VASP output to extxyz
- `vasp2extxyz_recursive` — recursively convert VASP outputs under a directory

# Examples
```julia
using Magesty
Magesty.install_tools()
# or to install into a custom location
Magesty.install_tools(bindir="/usr/local/bin")
```
"""
function install_tools(; bindir::AbstractString = joinpath(homedir(), ".julia", "bin"))
    pkg_dir = pkgdir(Magesty)
    if isnothing(pkg_dir)
        error("Cannot determine Magesty package directory")
    end

    mkpath(bindir)

    tools = [
        "vasp2extxyz"           => joinpath("tools", "vasp", "vasp2extxyz.jl"),
        "vasp2extxyz_recursive" => joinpath("tools", "vasp", "vasp2extxyz_recursive.jl"),
    ]

    global_env = "@v$(Sys.VERSION.major).$(Sys.VERSION.minor)"
    for (name, rel_path) in tools
        script = joinpath(pkg_dir, rel_path)
        wrapper = joinpath(bindir, name)
        write(wrapper, "#!/bin/sh\nexec julia --project=$global_env \"$script\" \"\$@\"\n")
        chmod(wrapper, 0o755)
        println("Installed: $wrapper")
    end

    println("\nMake sure ~/.julia/bin is in your PATH:")
    println("  export PATH=\"\$HOME/.julia/bin:\$PATH\"")
end

"""
	versioninfo(io::IO = stdout)

Print the Magesty version and the active Julia version to `io`, following
the same minimal style as `Base.versioninfo()`. Useful for recording the
runtime context in scripts and notebooks.

# Examples
```julia
julia> Magesty.versioninfo()
Magesty Version 0.1.0
Julia Version 1.11.0
```
"""
function versioninfo(io::IO = stdout)
	println(io, "Magesty Version ", Version.version_string())
	println(io, "Julia Version ", VERSION)
end

end # module Magesty
