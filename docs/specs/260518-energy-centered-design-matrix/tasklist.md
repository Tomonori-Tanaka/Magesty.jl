# Tasklist: energy-centered design matrix

Status: draft (2026-05-18)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Capture pre-refactor reference values

- [ ] Pick two integration datasets (`febcc_2x2x2_pm`,
      `fege_2x2x2`). For each, fit `OLS()` at
      `weight ∈ {0.0, 0.5, 1.0}` and `Ridge(lambda = 0.1)` at
      `weight = 0.5`. Print `(j0, jphi)` at full `Float64`
      precision via `repr(::Float64)` (which round-trips
      exactly).
- [ ] Drop the captured values into a new test file
      `test/component/test_energy_centered_design_matrix.jl`,
      with the actual `@test` calls commented out / wrapped in
      `@test_skip` for now. The file lives in-tree so the
      literals survive into M3.
- [ ] No source code changes in M1. This milestone exists to
      lock in the reference output of the current code before
      anything in `Fitting.jl` or `Magesty.jl` moves.

### M2 — Refactor (one commit)

The `X_E` shape change, the `assemble_weighted_problem`
rewrite, the `solve_coefficients` contract change, and the
`extract_j0_jphi` simplification are tightly coupled —
committing between them leaves the test suite red. Do them
together in one commit.

- [ ] `build_design_matrix_energy`: shape
      `(num_spinconfigs, num_salcs)`; remove the
      `[:, 1] .= 1.0` initialization.
- [ ] `SCEDataset` field docstring: remove the "bias column at
      column 1" line.
- [ ] `assemble_weighted_problem`: compute
      `mean_y = mean(observed_energy_list)` and
      `mean_X = mean(design_matrix_energy, dims = 1)` once;
      build the centered blocks; row-scale each block; vcat.
      Drop the `[:, 1] .= 1.0` reset and the
      `hcat(zeros(...), normalized_design_matrix_torque)`
      prepend. Return `(X, y)`.
- [ ] `solve_coefficients`: drop `bias_col::Int = 1` from the
      contract and from `OLS` / `Ridge` methods. `Ridge` no
      longer needs `lambda_vec[bias_col] = 0.0`.
- [ ] `extract_j0_jphi`: `jphi = collect(j_values)`;
      `j0 = mean(observed_energy_list .- design_matrix_energy * jphi)`.
      Update the docstring.
- [ ] Body of `fit(::Type{SCEFit}, ...)` (around
      `src/Magesty.jl:581`): destructure `(X, y)`; drop the
      `bias_col` kwarg on the `solve_coefficients` call.
- [ ] `predict_energy(::SCEModel, ::SCEDataset)`:
      `dataset.X_E[:, 2:end]` -> `dataset.X_E`.
- [ ] `fit` docstring in `Magesty.jl`: update the closed-form
      reference-energy recovery line to drop the slice.
- [ ] Update test layout assertions:
      - `test/component/test_SCEDataset.jl`: drop
        `all(d.X_E[:, 1] .== 1.0)`; column-count expectations
        align to `num_salcs`.
      - `test/component/test_Fitting_dispatch.jl`: drop
        `bias_col` destructuring and solver kwarg; adjust the
        energy / torque block slicing in row assertions.
      - `test/component/test_SCEFit.jl`:
        `X_E[:, 2:end] * f.jphi` -> `X_E * f.jphi`.
      - `test/integration/{febcc_2x2x2_pm,fege_2x2x2}/test.jl`:
        `size(...) == num_salcs + 1` -> `== num_salcs`;
        `X_E[:, 2:end] * jphi_true` -> `X_E * jphi_true`.

### M3 — Enable the regression test and measure tolerance

- [ ] Un-skip the regression test from M1.
- [ ] Run once with a placeholder loose tolerance (e.g.,
      `rtol = 1e-6`) and record the maximum observed `rtol`
      across all four `(estimator, weight)` pairs and both
      datasets. Expected order of magnitude is `rtol ≈ 1e-10`.
- [ ] Set the test's `rtol` to that maximum, padded by 10x and
      rounded up to the next power of 10. Record the chosen
      value and the worst observed value in the test docstring
      so the next reviewer can verify the choice.
- [ ] If the observed disagreement is substantially larger than
      `~1e-10` for well-conditioned cases, stop and investigate
      — likely a transcription error in M2, not a real
      floating-point limitation.
- [ ] `make test-all` green.
- [ ] `make test-jet`, `make test-aqua` clean.

### M4 — Docs, comment, and CHANGELOG

- [ ] `docs/src/api.md`: update any `SCEDataset` `X_E`
      reference that mentions a bias column or `num_salcs + 1`
      width.
- [ ] `SPEC.md`: same sweep.
- [ ] **Write the post-refactor explanatory comment** above the
      `solve_coefficients` / `extract_j0_jphi` call pair in the
      body of `fit(::Type{SCEFit}, ...)` in `src/Magesty.jl`.
      The comment should answer "why two steps?" from the
      perspective of the *new* code:
      - `j0` is eliminated analytically before the solve (the
        energy block is mean-centered inside
        `assemble_weighted_problem`).
      - The solver therefore returns only `jphi` for the
        centered, weighted augmented system.
      - `extract_j0_jphi` recovers `j0` from the *unscaled,
        un-centered* energy data via
        `mean(y_E - X_E * jphi)`, keeping `j0` in the input
        energy unit and independent of `torque_weight` /
        estimator choice.
      Keep it short (4–6 lines).
- [ ] `CHANGELOG.md` `[Unreleased]` `### Internal`: record the
      refactor with an explicit "no observable behavior change
      for OLS / Ridge beyond floating-point rounding noise"
      line and the measured `rtol` value from M3.

### M5 — Wrap-up

- [ ] Update `Status:` in this file and in
      [`docs/specs/README.md`](../README.md) to
      `complete (YYYY-MM-DD)` together.
- [ ] Append implementation commit SHAs below.

## Follow-up (not in this branch)

The ElasticNet spec
([`260517-elasticnet-estimator`](../260517-elasticnet-estimator/))
was drafted against the bias-in-`X` representation and assumes a
special energy-only-centering path for the L1 case. Once this
refactor merges to `main`, the `refactor/lasso-estimator` branch
should rebase onto `main` and the 260517 spec should be
simplified there:

- Drop the `_center_energy_block!` helper, the `n_energy`
  keyword on `solve_coefficients`, and the `bias_col` plumbing
  from the ElasticNet `solve_coefficients` path.
- `solve_coefficients(::ElasticNet, X, y)` becomes a straight
  `glmnet(X, y; alpha, lambda, standardize, intercept = false)`
  call.
- Adjust the `docs/specs/README.md` one-line summary for 260517
  to drop the now-redundant centering qualifier.

These edits are tracked on the lasso branch (not on this
branch's tasklist) because the 260517 spec files do not exist
here.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] Regression test
      (`test/component/test_energy_centered_design_matrix.jl`)
      passes at the `rtol` measured in M3, for all four
      `(estimator, weight)` pairs captured in M1.
- [ ] Post-refactor explanatory comment is in place above the
      two-step coefficient recovery in `src/Magesty.jl`.
- [ ] ~~If public API changed: `SPEC.md` and `docs/src/api.md`
      updated.~~ No public API change; `SPEC.md` /
      `docs/src/api.md` are still touched for the internal
      `X_E` shape description.
- [ ] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ `assemble_weighted_problem` is
      called once per fit (not hot); no bench entry needed.
- [ ] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ None changed.
- [ ] `CHANGELOG.md` `[Unreleased]` updated, with the chosen
      `rtol` recorded.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.

## Implementation commits

<!-- Append `<short SHA>  <subject>` lines here as commits land. -->
