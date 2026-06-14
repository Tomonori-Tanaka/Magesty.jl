# GCV.jl -- included into module Magesty (not a standalone module).
#
# Combined energy+torque generalized cross-validation diagnostics. This file
# is spliced into `Magesty` via `include`, so every name below is a binding of
# `Magesty`, exactly as if defined inline.
#
# Parse-time dependencies (must be defined before this include, all provided
# by Magesty.jl above the include site):
#   - SCEFit, SCEDataset (used in method return-type annotations)
#   - estimator types: AbstractEstimator (the GCVSizeCurve `estimator` field)
#     and OLS (the gcv_learning_curve keyword default)
# Call-time dependencies (resolved when a verb runs, not at parse; all in
# scope via Magesty.jl's `using` lines):
#   - SCEDataset integer-vector slicing
#   - Fitting internals: assemble_weighted_problem, _gcv_sample_count,
#     _gcv_single, _gcv_lambda_path, _gcv_msy, _gcv_r2,
#     _require_linear_estimator, solve_coefficients
#   - mean, std (Statistics); MersenneTwister, randperm (Random)
# Include after `using .Fitting` and the core type declarations.

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
