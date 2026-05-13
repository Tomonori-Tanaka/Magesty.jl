# Design: Optimize.jl estimator dispatch refactor

Status: draft (2026-05-13)
Spec: [`requirements.md`](requirements.md)
Source design memo: `DESIGN_NOTES.md` §"設計案: Optimize.jl の estimator
dispatch リファクタリング（2026-05-13、未着手）"

## Module layout

Everything lives inside `src/Optimize.jl`. No new files.

Before (current):

```
_fit_sce_model_internal(...)        # isa-branches on estimator
  ├─ fit_sce_model_ols(...)         # delegates to elastic_net with α=λ=0
  └─ fit_sce_model_elastic_net(...) # delegates to elastic_net_regression
       └─ elastic_net_regression    # assemble + solve + extract (3 jobs)
```

After (proposed):

```
_fit_sce_model_internal(...)
  ├─ assemble_weighted_problem(...)            # estimator-agnostic
  ├─ solve_coefficients(estimator, X, y; ...)  # dispatch boundary
  └─ extract_j0_jphi(...)                      # estimator-agnostic
```

`fit_sce_model_ols` and `fit_sce_model_elastic_net` become thin shims
that build an estimator struct and call `_fit_sce_model_internal`. They
are NOT exported (already not exported today) but stay as internal
helpers so existing `test/examples/*/test.jl` and `benchmark_optimize.jl`
do not need to change.

## Types

```julia
abstract type AbstractEstimator end

struct OLS <: AbstractEstimator end

struct Ridge <: AbstractEstimator
    lambda::Float64
end
Ridge(; lambda::Real = 0.0) = Ridge(Float64(lambda))
```

`ElasticNet` is **removed**. Today's `ElasticNet(alpha, lambda)`
struct is functionally a pure L2 ridge (its `alpha` field is
silently ignored by `elastic_net_regression`'s solve block). The
new `Ridge` reflects that reality and frees the `ElasticNet` name
for a future spec that adds a real mixed-norm estimator with an
honored `α ∈ (0, 1)`.

### Migration of `ElasticNet` call sites

| Site | Change |
|------|--------|
| `src/Optimize.jl` struct + ctor | Replace `ElasticNet` definition with `Ridge`. Remove `ElasticNet` from `export`. |
| `src/Magesty.jl` re-export | Replace `ElasticNet` with `Ridge`. |
| `src/Optimize.jl` `Optimizer(...)` default | `ElasticNet(alpha=alpha, lambda=lambda)` → `Ridge(lambda=lambda)`. Non-zero `alpha` triggers `@warn` (`"alpha is deprecated and ignored — pass lambda only"`) then is dropped. |
| `src/utils/ConfigParser.jl` `Config4Optimize` | Keep `alpha::Float64` field + `:alpha => 0.0` default so existing `input.toml` files still parse. After parse, if `alpha != 0.0`, emit `@warn` once. The constructor that builds the estimator now uses `Ridge(lambda=lambda)`. |
| `test/examples/fept_tetragonal_2x2x2/test.jl` | `ElasticNet(alpha=..., lambda=...)` → `Ridge(lambda=...)`. |
| `test/benchmark_optimize.jl` | Same. |
| `tools/check_convergence_embset.jl` | Same. |
| `SPEC.md`, `docs/src/{api,tutorial,examples}.md` | Prose + code samples updated. |

No deprecated `const ElasticNet = Ridge` alias is introduced — the
spec's intent is to free the `ElasticNet` name. A `MethodError` at
call sites is the explicit migration signal.

## Three-layer API

### Layer 1: `assemble_weighted_problem` (estimator-agnostic)

```julia
"""
    assemble_weighted_problem(Xe, Xt, ye, yt, weight)
        -> (X, y, bias_col)

Build the augmented `(X, y)` system that all estimators solve.

# Arguments
- `Xe::AbstractMatrix{<:Real}`: Energy design matrix (bias column at
  column 1).
- `Xt::AbstractMatrix{<:Real}`: Torque design matrix (no bias column).
- `ye::AbstractVector{<:Real}`: Observed energies.
- `yt::AbstractVector{<:AbstractMatrix{<:Real}}`: Observed torques as
  `3×n_atoms` matrices per spin configuration.
- `weight::Real`: Trade-off; energy rows are scaled by √(1 - weight),
  torque rows by √weight.

# Returns
- `X::Matrix{Float64}`: Augmented design matrix (rows = energies
  stacked above flattened torques; column 1 is the bias column).
- `y::Vector{Float64}`: Augmented observation vector.
- `bias_col::Int`: Column index of the bias term (= 1).
"""
```

Implementation is exactly the body of the current
`elastic_net_regression` up to (but not including) the
`lambda_vec`/solver block:

- `w_e = 1 - weight`, `w_m = weight`.
- Flatten `yt` into a column vector.
- Scale `Xe`, `Xt`, `ye`, `yt_flat` by √w.
- Reset `Xe[:, 1] .= 1.0` (bias column kept at 1, not √w).
- Prepend a zero bias column to scaled torque matrix.
- `vcat` energy and torque blocks.
- Return `(X, y, 1)`.

### Layer 2: `solve_coefficients` (dispatch boundary)

```julia
"""
    solve_coefficients(estimator::AbstractEstimator, X, y; bias_col)
        -> Vector{Float64}

Solve for the augmented coefficient vector. The element at index
`bias_col` is the bias (j0 component); estimators that apply
regularization MUST exclude that column from the penalty.

Each `AbstractEstimator` subtype defines exactly one method of this
function. Unknown types raise `MethodError`.
"""
function solve_coefficients end

solve_coefficients(::OLS, X, y; bias_col::Int = 1) = X \ y

function solve_coefficients(e::Ridge, X, y; bias_col::Int = 1)
    if e.lambda ≈ 0.0
        return X \ y
    end
    lambda_vec = fill(e.lambda, size(X, 2))
    lambda_vec[bias_col] = 0.0
    return ridge(X, y, lambda_vec; bias = false)
end
```

Behavior of `solve_coefficients(::Ridge, ...)` is identical to the
current `elastic_net_regression`'s solve block when `alpha = 0`
(same `λ ≈ 0` fallback, same `lambda_vec[1] = 0.0`, same
`MultivariateStats.ridge` call). Since the current code silently
ignored `alpha` anyway, this is a numerical no-op for all existing
test data.

`solve_coefficients(::OLS, ...)` is `X \ y`. Today's OLS path goes
through `fit_sce_model_ols → fit_sce_model_elastic_net(α=0, λ=0) →
elastic_net_regression`, which also ends at `X \ y` because the
`λ ≈ 0` short-circuit fires. So numerical output is identical.

### Layer 3: `extract_j0_jphi` (estimator-agnostic)

```julia
"""
    extract_j0_jphi(j_values, Xe, ye) -> (j0::Float64, jphi::Vector{Float64})

Split the augmented coefficient vector into the bias (j0) and the
SCE coefficients (jphi). j0 is re-estimated from the energy residual
to undo the √weight scaling on the bias column.

# Notes
- `Xe` is the **unscaled** energy design matrix; `ye` is the
  **unscaled** observed energy vector.
- `j_values[2:end]` are taken as `jphi`; `j_values[1]` is discarded
  in favor of the residual-based estimate
  `mean(ye .- Xe[:, 2:end] * jphi)`, matching the current
  `elastic_net_regression` behavior.
"""
```

This is the existing `j0 = mean(observed_energy_list .- design_matrix_energy[:, 2:end] * jphi)`
moved out of the solver. Bias estimation logic is unchanged.

## Orchestrator

```julia
function _fit_sce_model_internal(
    design_matrix_energy::AbstractMatrix{<:Real},
    design_matrix_torque::AbstractMatrix{<:Real},
    observed_energy_list::AbstractVector{<:Real},
    observed_torque_list::AbstractVector{<:AbstractMatrix{<:Real}},
    estimator::AbstractEstimator,
    weight::Real,
)
    X, y, bias_col = assemble_weighted_problem(
        design_matrix_energy,
        design_matrix_torque,
        observed_energy_list,
        observed_torque_list,
        weight,
    )
    j_values = solve_coefficients(estimator, X, y; bias_col = bias_col)
    return extract_j0_jphi(j_values, design_matrix_energy, observed_energy_list)
end
```

Three lines. No `isa`. No `else throw(ArgumentError(...))`.

## Backward-compat shims

```julia
function fit_sce_model_ols(Xe, Xt, ye, yt, weight)
    return _fit_sce_model_internal(Xe, Xt, ye, yt, OLS(), weight)
end

function fit_sce_model_elastic_net(Xe, Xt, ye, yt, alpha, lambda, weight)
    if alpha != 0.0
        @warn "alpha is deprecated and ignored; use Ridge(lambda)"
    end
    return _fit_sce_model_internal(
        Xe, Xt, ye, yt,
        Ridge(lambda = lambda),
        weight,
    )
end
```

`fit_sce_model_elastic_net` is renamed to `fit_sce_model_ridge`
(internal, not exported). The old name is kept as a thin deprecated
alias for one release so the TOML/Config4Optimize path that calls
into it does not break in a single jump. Both keep their current
docstrings (lightly edited) so JET/Aqua do not regress.

`elastic_net_regression` itself: removed. `grep` (recorded in this
spec's research phase) confirms no caller outside `src/Optimize.jl`,
and it is not exported.

## `Optimizer` outer constructor

Signature stays:

```julia
function Optimizer(
    structure, symmetry, basisset,
    alpha::Real, lambda::Real, weight::Real,
    spinconfig_list::AbstractVector{SpinConfig};
    verbosity::Bool = true,
    estimator::AbstractEstimator = (alpha != 0.0
        ? (@warn("alpha is deprecated and ignored; use Ridge(lambda)");
           Ridge(lambda = lambda))
        : Ridge(lambda = lambda)),
)
```

Today this signature is ambiguous: if the caller passes `estimator`,
`alpha` / `lambda` are silently ignored. We leave that behavior
unchanged (still out of scope for this spec) and add a one-line
doc note:

> If `estimator` is passed explicitly, the positional `alpha` and
> `lambda` arguments are ignored. The positional `alpha` is
> deprecated and is only consulted to emit a warning when non-zero;
> drive estimator choice via the `estimator` keyword.

A cleaner constructor (`Optimizer(structure, symmetry, basisset,
estimator, weight, spinconfig_list)`) is left for the user-facing
API spec.

## `Config4Optimize` (TOML)

`src/utils/ConfigParser.jl` keeps the `alpha::Float64` field and
the `:alpha => 0.0` default so existing `input.toml` files parse
unchanged. After parsing, if `alpha != 0.0`, emit a single
`@warn` per `Config4Optimize` instance:

> `[regression].alpha` is deprecated and currently has no effect.
> Use `Ridge`-style L2 regularization via `lambda` only. Remove
> this key from your TOML; future releases may reject it.

The constructor that builds the estimator from
`Config4Optimize` switches from `ElasticNet(alpha, lambda)` to
`Ridge(lambda)`. The TOML key itself is **not** removed in this
spec (would break every existing input file in `test/examples/`);
removal is a follow-up after a deprecation window.

## Test plan

1. **Golden regression test** (new), in `test/component_test/test_Optimize.jl`
   (create the file if missing — currently empty per grep):
   - Build a small synthetic `Xe`, `Xt`, `ye`, `yt`.
   - Call `_fit_sce_model_internal(...)` with `OLS()` and check
     `(j0, jphi)` matches a precomputed reference (captured from
     `main` at branch start, hard-coded in the test).
   - Same with `Ridge(lambda=0.1)`.
   - Assert `bias_col = 1` and `assemble_weighted_problem` reset
     of `Xe[:, 1] .= 1.0`.
2. **Smoke**: `test/examples/dimer`, `chain`, `2d_fcc_2x2x2`,
   `fept_tetragonal_2x2x2` — updated to `Ridge(lambda=...)` where
   they construct an estimator directly; otherwise unchanged.
   Run through `make test-integration`.
3. **TOML deprecation path**: confirm that a `[regression].alpha = 0.0`
   key (the existing test inputs) parses without warning, and that
   a `[regression].alpha = 0.5` (synthetic) emits exactly one
   `@warn` and still produces identical numerics to today's run
   (since `alpha` was always ignored at the solver level).
4. **Static analysis**: `make test-jet`, `make test-aqua` stay green.
5. **Benchmark**: `julia --project test/benchmark_optimize.jl --with-fit --samples 20` before vs after, recorded in
   `.claude/bench_log.md`.

## Linked sections (must update when this lands)

- `src/Magesty.jl`: re-export `Ridge` instead of `ElasticNet`.
- `DESIGN_NOTES.md`: append a "完了 → PR #XXX" line to the
  estimator-dispatch section after merge; note the `ElasticNet → Ridge`
  rename in the same line.
- The user-facing API design memo references this refactor as a
  prerequisite — its example code uses `ElasticNet(...)` symbolically;
  update that example to `Ridge(...)` to stay consistent with reality.

## Open questions (resolved)

- ~~**Q1**: `_fit_sce_model_internal` positional vs kwarg estimator?~~
  → Positional (diff minimality).
- ~~**Q2**: Keep `elastic_net_regression` as deprecated wrapper or
  delete?~~ → Delete. Grep confirmed no external callers.
- ~~**Q3**: `bias_col` trait?~~ → Not needed (postpone).
- ~~**Q4**: `ElasticNet` rename direction?~~ → Option A: rename to
  `Ridge`, no alias. Discussed and approved 2026-05-13.
