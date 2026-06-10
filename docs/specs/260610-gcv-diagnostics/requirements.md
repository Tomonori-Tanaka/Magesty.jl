# Requirements: GCV diagnostics (hat-matrix generalized cross-validation)

Status: draft (2026-06-10)

## Goal

Add a generalized cross-validation (GCV) diagnostic for SCE fitting, usable two
ways: (1) selecting the ridge penalty `lambda` by minimizing GCV, and (2)
checking data sufficiency by sweeping the number of training spin configurations
and watching whether the GCV score plateaus.

## Background

There is no cross-validation or model-selection diagnostic in the package today.
Users currently have only in-sample metrics (`rmse_*`, `r2_*`). GCV estimates
out-of-sample prediction error from a single fit via the hat (influence) matrix
`H` (`ŷ = H y`):

```
GCV = (‖r‖² / N) / (1 − tr(H)/N)²
```

The fit already minimizes a single weighted least-squares objective stacking
energy and torque rows (`assemble_weighted_problem`), so the combined
energy+torque GCV is the natural, exactly-defined quantity for linear
estimators. The user wants both penalty selection and a data-sufficiency check,
with simple table output that feeds the existing FitCheck Python plotting flow.

## Scope

Includes:

- New exported core API in `Magesty`: `gcv(f::SCEFit)`,
  `gcv_lambda(dataset, lambdas; torque_weight)`,
  `gcv_learning_curve(dataset, estimator; sizes, repeats, seed, torque_weight)`.
- GCV for the linear estimators `OLS` / `Ridge` (exact), and `AdaptiveRidge`
  (conditional: dof from the converged reweighting treated as fixed; the weights
  themselves are recovered exactly from `coef(f)`).
- Result types `GCVLambdaPath` / `GCVSizeCurve`.
- Table writers `write_gcv_lambda` / `write_gcv_learning_curve` in `FitCheckIO`.
- Python plot scripts under `tools/`.
- One example script; docs (narrative page + `api.md` `@docs`).

Excludes:

- Non-linear estimators — `ElasticNet` (and its `Lasso` alias, which constructs
  an `ElasticNet`) and `AdaptiveLasso`: no exact hat matrix; these raise
  `ArgumentError`.
- k-fold / leave-one-out CV, automatic convergence detection, automatic
  `lambda` grid generation beyond a simple default.
- Any change to the fitting objective, solvers, or numerical conventions.

## Invariants

- No change to existing numerical results: `fit`, `solve_coefficients`,
  `assemble_weighted_problem`, the design matrices, and `j0` recovery are
  untouched. GCV only reads them.
- Existing XML round-trips byte-for-byte (no I/O format change).
- Spin direction stays unit-vector with `3 × n_atoms` layout.
- GCV target is always the combined, weighted energy+torque objective (the same
  thing the fit minimizes), consistent with `torque_weight`. The sample count
  `N` and the dof count only *live* rows: a block zeroed by the weighting
  (`torque_weight = 1` drops energy, `0` drops torque) does not inflate `N`.
- `j0` (eliminated by energy mean-centering) contributes one effective degree of
  freedom only when the energy block is live (`torque_weight < 1`).

## Completion criteria

- [ ] `make test-all` passes, including new `test/component/test_gcv.jl`.
- [ ] `make test-jet` / `make test-aqua` clean (no new warnings).
- [ ] GCV reduces to the closed-form OLS value on a synthetic design (analytic
      test, not tailored to output).
- [ ] SVD-based ridge dof and λ-path RSS match dense recomputation to ~1e-10.
- [ ] `gcv_lambda` at a given λ equals `gcv(fit(..., Ridge(lambda=λ)))`.
- [ ] Non-linear estimators raise `ArgumentError`.
- [ ] New API in `docs/src/api.md` and a narrative diagnostics page; one example.

## References

- Related specs: `260518-energy-centered-design-matrix` (the `j0`-by-centering
  design GCV relies on), `260520-fitcheck-io-writers` (table-writer convention),
  `260526-solver-unification-and-memory-log` (Cholesky-on-`X'X` solvers).
- Golub, Heath, Wahba (1979), "Generalized Cross-Validation as a Method for
  Choosing a Good Ridge Parameter", Technometrics 21(2).
