# Evaluation.jl -- included into module Magesty (not a standalone module).
#
# Prediction verbs and the accuracy metrics over a (predictor, data) pair.
# This file is spliced into `Magesty` via `include`, so every name below is a
# binding of `Magesty`, exactly as if defined inline.
#
# Required in scope at include time (all provided by Magesty.jl above the
# include site):
#   - core types: SCEBasis, SCEModel, SCEDataset, SCEFit
#   - SpinConfig (referenced by `const SCEEvalData`, evaluated at parse time)
#   - Fitting internals: _predict_energy, _predict_torque, _calc_r2score,
#     _calc_rmse
#   - SpinConfigs.read_embset
# Include after the core type declarations and `using .Fitting`, and before
# the trailing domain includes (FitCheckIO.jl etc.) that call predict_*.

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
