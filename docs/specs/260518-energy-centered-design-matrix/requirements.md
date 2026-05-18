# Requirements: energy-centered design matrix

Status: draft (2026-05-18)

## Goal

Remove the bias column from `SCEDataset.X_E` and
`build_design_matrix_energy`, and replace the post-scaling
`[:, 1] .= 1.0` reset in `assemble_weighted_problem` with explicit
energy-only centering of the energy block. The bias term `j0` is no
longer a column of the augmented design matrix; it is recovered
analytically by `j0 = mean(y_E - X_E * jphi)` from the unscaled
energy data after the solver returns `jphi` directly.

Net effect: `X_E` and `X_T` have the same column count
(`num_salcs`); `assemble_weighted_problem` returns `(X, y)` with no
`bias_col`; `solve_coefficients` no longer threads `bias_col`;
`extract_j0_jphi` collapses to a one-line `j0` computation; the
coefficient vector returned by the solver is exactly `jphi` (no slot
to discard).

This is an internal refactor only. The two formulations are
mathematically equivalent — `j0` is eliminated analytically (block
elimination using `∂L/∂j0 = 0`) and the residual problem in `jphi`
alone is what the solver sees. The floating-point operation
sequence is different, however: different matrix shape, different
LAPACK / QR path, different summation order inside `mean(...)`.
`(j0, jphi)` from `OLS()` and `Ridge(lambda = …)` are therefore
expected to agree within floating-point rounding noise from those
order changes, not bit-for-bit. The actual disagreement is
measured in M3 and the regression-test tolerance is set from that
measurement (see `tasklist.md`).

## Background

Two pieces of `Fitting.jl` are currently load-bearing but hard to
read on first contact:

1. `assemble_weighted_problem` scales every column of the energy
   block by `√((1 - weight) / n_E)` and then immediately rewrites
   the bias column back to `1.0`. The "scale then partially
   un-scale" step exists only to prevent `j_values[1]` from being
   absorbed into the weighting factor. A reader who does not know
   the convention assumes the reset is a bug.
2. `extract_j0_jphi` then discards `j_values[1]` and re-fits `j0`
   from the unscaled energy data via
   `mean(y_E - X_E[:, 2:end] * jphi)`. The augmented intercept the
   solver computes is unused.

Both pieces become unnecessary if `j0` is eliminated analytically
before the solve. The mixed energy / torque objective

```
L(j0, jphi) = (1 - w) / n_E * ||y_E - 1 * j0 - X_E * jphi||^2
            +       w / n_T * ||y_T -        X_T * jphi||^2
```

has the closed-form minimizer `j0(jphi) = mean(y_E - X_E * jphi)`
(torque rows are insensitive to `j0`). Substituting back gives the
centered objective

```
L(jphi) = (1 - w) / n_E * ||y_tilde_E - X_tilde_E * jphi||^2
        +       w / n_T * ||y_T       - X_T       * jphi||^2
```

with `y_tilde_E = y_E - mean(y_E)` and
`X_tilde_E = X_E - mean(X_E, dims=1)`. This is mathematically
equivalent to the current setup for all `0 <= weight <= 1` and for
every estimator that does not penalize `j0` — currently `OLS` and
`Ridge`. The derivation is reproduced in `design.md`.

Forward link: the upcoming ElasticNet spec
([`260517-elasticnet-estimator`](../260517-elasticnet-estimator/))
originally proposed energy-only centering as a special path for the
L1 case only, with bias-handling helpers and an `n_energy` keyword
threaded through `solve_coefficients`. With this refactor in place,
energy-only centering is the default representation; ElasticNet
plumbing simplifies to a straight `glmnet(X, y; intercept=false, ...)`
call. Editing the 260517 spec to reflect that simplification is
a follow-up on the `refactor/lasso-estimator` branch (where the
spec lives), not part of this branch — see the "Follow-up" note
in `tasklist.md`.

## Scope

Includes:

- `build_design_matrix_energy`: returns `Matrix{Float64}` of shape
  `(num_spinconfigs, num_salcs)`. The leading `1.0` bias column is
  removed; the function no longer initializes
  `design_matrix[:, 1] .= 1.0`.
- `SCEDataset.X_E::Matrix{Float64}` shape is `(n_E, num_salcs)`,
  matching `X_T`'s width. The field's docstring is updated.
- `assemble_weighted_problem`: rewritten to center the energy block
  (`y_E .- mean(y_E)` and `X_E .- mean(X_E, dims = 1)`) before row
  scaling. Drops the `[:, 1] .= 1.0` reset and the
  `hcat(zeros(...), X_T)` prepend on the torque block. Returns
  `(X::Matrix{Float64}, y::Vector{Float64})` — no `bias_col`.
- `solve_coefficients`: drops the `bias_col::Int = 1` keyword from
  the contract and from the `OLS` / `Ridge` methods. `Ridge` no
  longer needs the `lambda_vec[bias_col] = 0.0` line — every column
  of the new `X` is an SCE coefficient that should be penalized.
- `extract_j0_jphi`: simplifies to
  `j0 = mean(observed_energy_list .- design_matrix_energy * jphi)`
  (no `[:, 2:end]` indexing). Signature unchanged.
- `predict_energy(::SCEModel, ::SCEDataset)`:
  `dataset.X_E[:, 2:end] * model.jphi .+ model.j0` becomes
  `dataset.X_E * model.jphi .+ model.j0`.
- Docstring updates in `Fitting.jl` and `Magesty.jl` for the new
  `X_E` shape, the new `assemble_weighted_problem` return tuple,
  and the new `solve_coefficients` contract.
- Test updates: `test/component/test_SCEDataset.jl`,
  `test/component/test_SCEFit.jl`,
  `test/component/test_Fitting_dispatch.jl`,
  `test/integration/febcc_2x2x2_pm/test.jl`,
  `test/integration/fege_2x2x2/test.jl` — drop the `+ 1` on
  expected column counts, drop the `[:, 2:end]` slicing, drop the
  `bias_col` destructuring.
- A new regression test
  (`test/component/test_energy_centered_design_matrix.jl`) that
  pins `(j0, jphi)` to pre-refactor reference literals for
  representative `(estimator, weight)` pairs.

Excludes:

- New estimators (the ElasticNet / Lasso work is covered by spec
  [260517-elasticnet-estimator](../260517-elasticnet-estimator/)).
- Changes to XML on-disk format. `X_E` is not serialized; only
  `jphi` / `j0` are. Confirmed during scoping.
- Public API signature changes: `fit`, `predict_energy`,
  `predict_torque`, `coef`, `intercept`, `Magesty.save`,
  `Magesty.load` are signature-unchanged.
- Hot-path optimization. `assemble_weighted_problem` is called once
  per fit and is not a hot path; the one extra
  `(n_E, num_salcs)` centering allocation is negligible. No
  `make bench-fitting` entry needed.

## Invariants

- For every existing fixture and synthetic example in the test
  suite, `(j0, jphi)` produced by `OLS()` and `Ridge(lambda=…)`
  after the refactor must agree with the pre-refactor values
  within floating-point rounding noise. The two formulations are
  mathematically equivalent (analytic elimination of `j0` via
  `∂L/∂j0 = 0`), but their floating-point operation sequences
  differ, so a non-zero residual is expected. The regression test
  pins the tolerance to a value measured in M3 (typically
  `rtol ≈ 1e-10`, padded by 10x over the worst observed
  disagreement). Any disagreement substantially larger than the
  measured baseline must be traced to a transcription error in
  the refactor, not absorbed by relaxing the tolerance.
- All other physics conventions are untouched: `3 x n_atoms` spin
  direction layout, tesseral spherical-harmonics normalization,
  SALC and Clebsch-Gordan conventions, `_cluster_scaling(N) =
  (4 * pi)^(N / 2)`.
- XML round-trip for `SCEBasis` / `SCEModel` continues byte-for-byte;
  this refactor does not touch anything that is serialized.
- `Fitting` <-> `SALCBasis` link (the key-group order of `X_E` /
  `X_T` columns) is untouched. Column count drops by 1; the
  remaining columns retain their order.
- Public API surface (`fit(SCEFit, ...)`, `predict_energy`,
  `predict_torque`, `coef`, `intercept`, the estimator types) is
  signature-unchanged. Users do not see the refactor.

## Completion criteria

- [ ] `make test-all` green.
- [ ] `make test-jet`, `make test-aqua` clean (no new warnings).
- [ ] `test/component/test_energy_centered_design_matrix.jl`
      compares `(j0, jphi)` from `fit(SCEFit, ...)` against
      pre-refactor reference literals for the four pairs:
      `(OLS, weight = 0.0)`, `(OLS, weight = 0.5)`,
      `(OLS, weight = 1.0)`, `(Ridge(lambda = 0.1), weight = 0.5)`.
      Tolerance is measured in M3 and recorded in the test
      docstring and in `CHANGELOG.md` (expected
      `rtol ≈ 1e-10`).
- [ ] `extract_j0_jphi` returns `(j0, jphi)` consistent with the
      new contract (no `j_values[1]` discard, no `[:, 2:end]`
      slice) for every input used in
      `test/component/test_Fitting_dispatch.jl`.
- [ ] `docs/src/api.md`, `SPEC.md`, and the relevant docstrings in
      `Fitting.jl` / `Magesty.jl` are updated to reflect the new
      `X_E` shape and the new `assemble_weighted_problem` return
      tuple.
- [ ] `CHANGELOG.md` `[Unreleased]` records the internal refactor
      under "Internal" with an explicit "no behavior change for OLS
      / Ridge" line.
- [ ] `Status:` line in this spec's `tasklist.md` and the matching
      row in [`docs/specs/README.md`](../README.md) updated to
      `complete (YYYY-MM-DD)` together.

## References

- Spec affected by this refactor (will be simplified afterward):
  [`260517-elasticnet-estimator`](../260517-elasticnet-estimator/)
- Earlier estimator-dispatch spec (sets the
  `solve_coefficients(::AbstractEstimator, …)` contract):
  [`260513-estimator-dispatch`](../260513-estimator-dispatch/)
- CLAUDE.md "Linked sites" section — `Fitting` <-> `SALCBasis`
  ordering invariant.
