# Design: energy-centered design matrix

Status: draft (2026-05-18)

## Summary

Replace the "bias column inside `X_E`" representation with energy-only
centering. The bias term `j0` is eliminated analytically before the
solver sees the problem, so `assemble_weighted_problem` returns an
augmented `(X, y)` whose energy block is mean-subtracted and whose
torque block is unchanged. `solve_coefficients` returns `jphi`
directly (no `j_values[1]` slot to discard); `extract_j0_jphi`
becomes a one-line post-processing step that computes
`j0 = mean(y_E - X_E * jphi)` from the unscaled inputs (same formula
as today, with the `[:, 2:end]` indexing dropped because there is no
longer a bias column).

This is the canonical setup that GLMNet expects in its
`intercept=false` mode, so the same code path is shared between
`OLS`, `Ridge`, and the forthcoming `ElasticNet` without any
estimator-specific branching for bias handling.

## Mathematical derivation

The mixed energy / torque objective on the per-sample-normalized
problem is

```
L(j0, jphi) = (1 - w) / n_E * sum_i [y_E[i] - j0 - X_E[i, :] * jphi]^2
            +       w / n_T * sum_k [y_T[k]      - X_T[k, :] * jphi]^2
```

`j0` enters only the first sum. Setting the partial derivative with
respect to `j0` to zero gives

```
j0_star(jphi) = (1 / n_E) * sum_i (y_E[i] - X_E[i, :] * jphi)
              = mean(y_E - X_E * jphi).
```

Substituting `j0 = j0_star(jphi)` back into `L`:

```
L(jphi) = (1 - w) / n_E
            * sum_i [(y_E[i] - mean(y_E))
                     - (X_E[i, :] - mean(X_E, dims = 1)) * jphi]^2
        + w / n_T * sum_k [y_T[k] - X_T[k, :] * jphi]^2
       = (1 - w) / n_E * ||y_tilde_E - X_tilde_E * jphi||^2
        +       w / n_T * ||y_T       - X_T       * jphi||^2.
```

The torque block is unchanged. The energy block becomes ordinary
least squares on the centered inputs. After the standard
`sqrt(weight / n)` row-whitening, the augmented system

```
X = [sqrt((1 - w) / n_E) * X_tilde_E ;
     sqrt(      w / n_T) * X_T      ]
y = [sqrt((1 - w) / n_E) * y_tilde_E ;
     sqrt(      w / n_T) * y_T      ]
```

is exactly what every estimator solves for `jphi`. Recovering `j0`
afterward uses the same closed form on the unscaled inputs:

```julia
j0 = mean(observed_energy_list .- design_matrix_energy * jphi)
```

This is the formula `extract_j0_jphi` already uses today; the
refactor only removes the `[:, 2:end]` slicing.

For `Ridge`, the L2 penalty is `lambda * ||jphi||^2`. The penalty
applies uniformly to every SCE coefficient (it never reached `j0`
under the existing `bias = false` flag plus zeroed `lambda_vec`
entry), so dropping the per-column lambda mask is a pure
simplification: every column of the new `X` is an SCE coefficient
that should be penalized.

For `OLS` (no regularization), the augmented `X \ y` returns the
same `jphi` as today by construction.

For the forthcoming `ElasticNet`, GLMNet with `intercept = false`
solves

```
(1 / 2N) * ||y - X * jphi||^2
    + lambda * [(1 - alpha) / 2 * ||jphi||_2^2 + alpha * ||jphi||_1]
```

on exactly the centered augmented system, with no need for a
`_center_energy_block!` helper or `n_energy` keyword — those appear
only because the current 260517 spec was written against the
bias-in-`X` representation.

## Module layout

| Target | Change |
|---|---|
| `src/Fitting.jl` `build_design_matrix_energy` | Allocate `(num_spinconfigs, num_salcs)`. Remove the `[:, 1] .= 1.0` init. |
| `src/Fitting.jl` `assemble_weighted_problem` | Center the energy block before row scaling. Drop the `[:, 1] .= 1.0` reset and the torque-side `hcat(zeros(...), ...)`. Return `(X, y)`. |
| `src/Fitting.jl` `solve_coefficients` | Drop `bias_col::Int = 1` from the contract and from `OLS` / `Ridge` methods. |
| `src/Fitting.jl` `extract_j0_jphi` | `jphi = collect(j_values)`; `j0 = mean(y_E .- X_E * jphi)`. |
| `src/Magesty.jl` `SCEDataset` field docstring | Drop the "bias column at column 1" sentence for `X_E`. |
| `src/Magesty.jl` body of `fit(::Type{SCEFit}, ...)` (around line 581) | Destructure `(X, y)`. Drop the `bias_col = bias_col` keyword on the `solve_coefficients` call. |
| `src/Magesty.jl` `predict_energy(::SCEModel, ::SCEDataset)` | `dataset.X_E[:, 2:end]` -> `dataset.X_E`. |
| `src/Magesty.jl` `fit` docstring | Update the `j0 = mean(y_E - X_E[:, 2:end] * jphi)` reference to drop the slice. |
| `test/component/test_SCEDataset.jl` | Drop `all(d.X_E[:, 1] .== 1.0)`. Adjust column-count assertions. |
| `test/component/test_Fitting_dispatch.jl` | Drop `bias_col` destructuring and solver kwarg. Adjust the energy / torque block slicing in the row assertions. |
| `test/component/test_SCEFit.jl` | `dataset.X_E[:, 2:end] * f.jphi` -> `dataset.X_E * f.jphi`. |
| `test/integration/febcc_2x2x2_pm/test.jl` | `size(...) == num_salcs + 1` -> `== num_salcs`. `X_E[:, 2:end] * jphi_true` -> `X_E * jphi_true`. |
| `test/integration/fege_2x2x2/test.jl` | Same two changes as above. |
| `test/component/test_energy_centered_design_matrix.jl` (new) | Regression test pinning `(j0, jphi)` to pre-refactor reference literals. |
| `docs/src/api.md` | Update the `SCEDataset` description for the new `X_E` shape (if the page references column count or bias). |
| `SPEC.md` | Same. |
| `CHANGELOG.md` | `[Unreleased]` -> `### Internal` entry. |

## API

```julia
# build_design_matrix_energy returns (num_spinconfigs, num_salcs).
function build_design_matrix_energy(
    salc_list::AbstractVector{Vector{CoupledBases.CoupledBasis_with_coefficient}},
    spinconfig_list::AbstractVector{SpinConfig},
    symmetry::Symmetry,
)::Matrix{Float64}

# assemble_weighted_problem returns (X, y) only.
function assemble_weighted_problem(
    design_matrix_energy::AbstractMatrix{<:Real},
    design_matrix_torque::AbstractMatrix{<:Real},
    observed_energy_list::AbstractVector{<:Real},
    observed_torque,                                  # matrices-per-config or flattened
    weight::Real,
)::Tuple{Matrix{Float64}, Vector{Float64}}

# solve_coefficients contract drops bias_col.
solve_coefficients(estimator::AbstractEstimator,
                   X::AbstractMatrix{<:Real},
                   y::AbstractVector{<:Real})::Vector{Float64}
# OLS:   X \ y  (unchanged)
# Ridge: lambda_vec = fill(e.lambda, size(X, 2));
#        ridge(X, y, lambda_vec; bias = false)

# extract_j0_jphi: signature unchanged, body simplifies.
function extract_j0_jphi(
    j_values::AbstractVector{<:Real},
    design_matrix_energy::AbstractMatrix{<:Real},
    observed_energy_list::AbstractVector{<:Real},
)
    jphi = collect(j_values)
    j0 = mean(observed_energy_list .- design_matrix_energy * jphi)
    return j0, jphi
end
```

## Types and conventions

- `SCEDataset.X_E` shape changes from `(n_E, num_salcs + 1)` to
  `(n_E, num_salcs)`. The "bias column at column 1" line in the
  field docstring is removed.
- No change to `_cluster_scaling`, `Z_lm`, SALC ordering, or any
  unit / sign convention.
- The numerical penalty `Ridge(lambda = lambda)` applies is
  unchanged: the `bias = false` flag plus the absence of a bias
  column means `Ridge` continues to penalize every SCE coefficient
  uniformly, exactly as before.
- The per-column-mask `lambda_vec[bias_col] = 0.0` line in
  `Ridge` disappears: every column of the new `X` is an SCE
  coefficient that should be penalized.
- `extract_j0_jphi` keeps its three-argument signature, but the
  expected length of `j_values` changes from `num_salcs + 1` to
  `num_salcs` (no bias slot to discard). The only caller is the
  body of `fit(::Type{SCEFit}, ...)` in `Magesty.jl`, updated in
  the same commit. The function docstring must call out the new
  expected length so misuse is caught on first read.

## Impact on linked sites

- [x] Spherical-harmonics convention (`TesseralHarmonics`):
      unchanged.
- [x] SCE coefficient XML (`save` / `load`): unchanged. `X_E` is
      not serialized; `jphi` / `j0` round-trip exactly as today.
- [x] `Fitting` <-> `SALCBasis`: column ordering of `X_E` / `X_T`
      is preserved; only the leading bias column of `X_E` is
      removed.
- [x] `.claude/agents/` references: none of the agent prompts
      mention the bias column or the `bias_col` keyword. No sweep
      required.
- [x] `SPEC.md` / `docs/src/api.md`: update the `SCEDataset`
      description for the new `X_E` shape.

## Test strategy

1. **Pre-refactor reference values (M1)**. Before any source edit,
   `fit(SCEFit, dataset, OLS())` at `weight ∈ {0.0, 0.5, 1.0}`
   and `fit(SCEFit, dataset, Ridge(lambda = 0.1))` at
   `weight = 0.5`. Two datasets are used (`febcc_2x2x2_pm` and
   `fege_2x2x2`). `(j0, jphi)` are captured at full `Float64`
   precision via `repr(::Float64)` (round-trips exactly) and
   placed as literal arrays in a draft regression test file. These
   literals are the *reference output of the pre-refactor code* —
   not a bit-for-bit target, since the new code uses a different
   floating-point operation sequence.
2. **Refactor regression test (M3)**. The file from M1 becomes
   `test/component/test_energy_centered_design_matrix.jl`. For
   each `(estimator, weight)` pair, run `fit(...)` after the
   refactor and assert `(j0, jphi)` matches the reference. The
   tolerance is measured during M3: run the test once, record the
   maximum observed `rtol` across all pairs, then set the
   regression tolerance to that value padded by 10x (rounded up
   to the next power of 10). Expected order of magnitude is
   `rtol ≈ 1e-10` based on the condition number of typical SCE
   design matrices; any observed disagreement substantially
   larger than that suggests a transcription error in the
   refactor and should be investigated, not absorbed. The chosen
   tolerance is recorded in the test docstring and in the
   `CHANGELOG.md` `[Unreleased]` entry.
3. **Existing test updates**. Every test currently checking the
   `X_E` shape or `[:, 2:end]` slicing is updated to the new
   layout. No new test logic is added; only the assertions about
   layout change.
4. **`assemble_weighted_problem` block test**. The block test in
   `test_Fitting_dispatch.jl` is adjusted to verify that the
   centered energy block reproduces
   `sqrt((1 - w) / n_E) * (X_E .- mean(X_E, dims = 1))` within
   `rtol = 1e-10`, and that the torque block is
   `sqrt(w / n_T) * X_T` within the same tolerance. Bit-for-bit
   equality is *not* asserted: the test and the production code
   may compute the centering with slightly different summation
   orders, and FP-rounding-noise-level differences are tolerated.

## Risks and open items

- **Risk: integration-test target values.** Both
  `febcc_2x2x2_pm` and `fege_2x2x2` carry hand-computed or
  fitted coefficient targets. Only the shape / slicing
  assertions should need updating; the numeric targets should
  not move. The M1 reference capture validates this implicitly
  before the source is edited.
- **Risk: synthetic target generation.** Both integration tests
  build `y_E_synth = X_E[:, 2:end] * jphi_true .+ j0_true`. The
  refactor changes this to `y_E_synth = X_E * jphi_true .+ j0_true`.
  This is a one-line change per test, but if missed, the test
  will silently fail with a shape mismatch.
- **Open: one extra allocation.** The new
  `assemble_weighted_problem` allocates one `(n_E, num_salcs)`
  centered matrix. Negligible at fit time (which already
  allocates the augmented `X`), but recorded here so the
  decision is on file. If this ever shows up in a profile, the
  centering can be folded into the row-scaling loop.
- **Out of scope: follow-up edit to the ElasticNet spec.** The
  existing
  [`260517-elasticnet-estimator`](../260517-elasticnet-estimator/)
  spec assumes the bias-in-`X` representation in several places.
  Those edits belong on the `refactor/lasso-estimator` branch
  (where the 260517 spec lives) — once this refactor lands on
  `main`, the lasso branch rebases onto `main` and the spec is
  simplified there. This refactor's tasklist does *not* carry
  those edits; tracking them on this branch would require
  cherry-picking the 260517 spec commit, which adds coupling
  between two otherwise-independent units.
