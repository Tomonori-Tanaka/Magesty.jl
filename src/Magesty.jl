"""
	Magesty

The main module of the Magesty package, providing the entry point for
spin-cluster expansion (SCE) analysis and fitting.

# Usage
```julia
using Magesty

basis   = SCEBasis("input.toml")
dataset = SCEDataset(basis, "EMBSET.dat")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)
save(SCEModel(f), "model.xml")
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
- `Optimize`: design-matrix construction and regression internals
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
using .SpinConfigs: SpinConfig, read_embset

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

include("utils/xml_io.jl")
using .XMLIO

export SCEBasis, SCEDataset, SCEFit, SCEModel
export AbstractEstimator, OLS, Ridge
export predict_energy, predict_torque
export fit, coef, intercept, nobs, dof
export r2_energy, r2_torque, rss_energy, rss_torque
export residuals_energy, residuals_torque, rmse_energy, rmse_torque
export save, load
export read_embset
export VERSION, install_tools

# Shared skeleton for the SCEBasis input-driven constructors.
# Returns the (structure, symmetry, cluster) triplet; callers append the
# SALCBasis — either computed via `SALCBasis(...)` or loaded from XML.
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
"""
struct SCEBasis
	structure::Structure
	symmetry::Symmetry
	salcbasis::SALCBasis
	isotropy::Bool
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
`body2`, ...; see the API documentation for the accepted format and
the per-path Unitful requirements.
"""
function SCEBasis(
	input_dict::AbstractDict{<:AbstractString, <:Any};
	verbosity::Bool = true,
)
	config::Config4System = Config4System(input_dict)
	structure, symmetry, cluster = _build_structure_skeleton(config; verbosity = verbosity)
	salcbasis::SALCBasis =
		SALCBasis(structure, symmetry, cluster, config, verbosity = verbosity)
	# `cluster` is a construction step only; it is not stored in SCEBasis.
	return SCEBasis(structure, symmetry, salcbasis, config.isotropy)
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

Extract a lightweight `SCEModel` from a fitted `SCEFit`: the `SCEBasis`
together with the fitted `(j0, jphi)`. The training dataset, residuals,
and cached metrics held by the `SCEFit` are dropped.
"""
SCEModel(f::SCEFit)::SCEModel = SCEModel(f.dataset.basis, f.j0, f.jphi)


# --- Prediction ---------------------------------------------------------

"""
	predict_energy(model::SCEModel, spin_directions::AbstractMatrix{<:Real}) -> Float64
	predict_energy(model::SCEModel, sc::SpinConfig) -> Float64
	predict_energy(f::SCEFit, spin_directions::AbstractMatrix{<:Real}) -> Float64
	predict_energy(f::SCEFit, sc::SpinConfig) -> Float64
	predict_energy(predictor, dataset::SCEDataset) -> Vector{Float64}

Predict SCE energies. The single-configuration forms return one
`Float64`; the `SCEDataset` form returns one energy per configuration.
`SCEFit` predictors delegate through `SCEModel(f)`. The `SCEDataset` form
requires the dataset to share the predictor's SALC basis.
"""
predict_energy(model::SCEModel, spin_directions::AbstractMatrix{<:Real})::Float64 =
	Optimize._predict_energy(
		model.j0, model.jphi,
		model.basis.salcbasis.salc_list, model.basis.symmetry, spin_directions)
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
	predict_torque(model::SCEModel, spin_directions::AbstractMatrix{<:Real}) -> Matrix{Float64}
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
predict_torque(model::SCEModel, spin_directions::AbstractMatrix{<:Real})::Matrix{Float64} =
	Optimize._predict_torque(
		model.jphi, model.basis.salcbasis.salc_list, model.basis.symmetry, spin_directions)
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
	model.basis === dataset.basis || throw(ArgumentError(
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


# --- save / load --------------------------------------------------------

function _require_xml_extension(path::AbstractString)
	endswith(lowercase(path), ".xml") || throw(ArgumentError(
		"save / load handle XML only; the path must end in '.xml', got: $path"))
	return nothing
end

"""
	save(obj::SCEBasis, path::AbstractString)
	save(obj::SCEModel, path::AbstractString)

Write `obj` to an XML file at `path`. An `SCEBasis` is written as
structure, symmetry parameters, and SALC basis; an `SCEModel` adds the
fitted reference energy and SCE coefficients in a `<JPhi>` block.

The path must end in `.xml`; any other extension is an error. Use
`load` to read the file back.

# Arguments
- `obj::Union{SCEBasis, SCEModel}`: The object to serialize.
- `path::AbstractString`: Output file path; must end in `.xml`.

# Throws
- `ArgumentError` if `path` does not end in `.xml`.
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

# Throws
- `ArgumentError` if `path` does not end in `.xml`, if a required schema
  attribute is missing, or (for `SCEModel`) if the `<JPhi>` block is
  absent.
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
