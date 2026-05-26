# Requirements: Solver unification (Cholesky) and design-matrix memory logging

Status: draft (2026-05-26)

## Goal

Unify the OLS and Ridge linear solvers around Cholesky on the normal
equations `X'X` (matching the existing `AdaptiveRidge` path), and surface
the design-matrix memory footprint to stdout so users can judge whether a
given problem fits on their local machine before they wait through the
build.

## Background

`solve_coefficients(::OLS, X, y)` currently uses non-pivoted QR with a
pivoted-QR fallback (commit 0c10182). QR is `κ(X)`-stable but ~2× slower
than Cholesky on the normal equations and needs to hold the full Householder
factorization in memory. Ridge already moves through `X'X` (line 1374:
`Symmetric(XtX + λI) \ Xty`, which dispatches to Bunch-Kaufman), and
`AdaptiveRidge` does the same per iteration.

The design philosophy supporting this unification: if `X` is so
ill-conditioned that squaring `κ(X)` matters, the OLS fit is itself
meaningless — the right answer in that regime is `Ridge(λ)`, not a more
defensive OLS solver. Routing both paths through Cholesky makes the
`PosDefException` an explicit user signal ("use Ridge") rather than a silent
loss of accuracy.

Separately, users have no way to know up front that a particular
`(num_spinconfigs, num_salcs)` pair will require, say, 40 GB of design
matrix. A two-line `stdout` log answers the "will this run on my laptop?"
question without changing any numerics.

## Scope

Includes:

- `solve_coefficients(::OLS, X, y)`: switch to `cholesky(Symmetric(X'X)) \ (X'y)`
  with `PosDefException` rethrown as a clearer error pointing to `Ridge`.
- `solve_coefficients(::Ridge, X, y)`: switch the regular path to
  `cholesky(Symmetric(X'X + λI)) \ (X'y)` (replacing the implicit
  Bunch-Kaufman dispatch via `Symmetric \`); route the `lambda ≈ 0` fast
  path to the OLS Cholesky path so the same `PosDefException` behavior
  applies (a user who asked for unregularized fit should see the same
  failure mode).
- `solve_coefficients(::AdaptiveRidge, X, y)`: convert the two `Symmetric \`
  calls to explicit `cholesky` for consistency. Route the `lambda ≈ 0` fast
  path to the OLS Cholesky path for the same reason.
- Remove the `MultivariateStats` dependency from `Project.toml` and the
  `using MultivariateStats` line from `Fitting.jl`. `MultivariateStats.ridge`
  was the only call site; no other module imports it.
- `build_design_matrix_energy` / `build_design_matrix_torque`: emit a
  pre-construction memory estimate (energy block, torque block, Gram matrix
  `num_salcs²`) and a post-construction `Base.summarysize` line, under the
  existing `verbosity = true` gate.
- OLS / Ridge docstring updates reflecting the new solver.
- Regression test (Cholesky-OLS agrees with the prior QR-OLS on a
  well-conditioned design) and a `PosDefException` rethrow test on a
  rank-deficient design.

Excludes:

- Changes to `ElasticNet` / `Lasso` / `AdaptiveLasso` / `PrecomputedPilot`
  (these dispatch through GLMNet or return precomputed values).
- Iterative solvers (LSQR / CG) — deferred.
- Truncated-SVD regularization — deferred.
- Sparse-matrix support — deferred.
- Logging memory for caches outside `Fitting` (e.g. SH cache, SALC
  storage).

## Invariants

- Numerical results agree with the prior QR solver to ~`1e-10` on
  well-conditioned designs (the same algebraic problem, different
  factorization).
- Bias term `j0` is still recovered separately by `extract_j0_jphi`; it
  does not enter the solve.
- Public API of `solve_coefficients`, `OLS`, `Ridge`, `AdaptiveRidge`
  unchanged (same constructors, same return type).
- `MultivariateStats.ridge` is no longer called from `Fitting.jl`. The
  dependency is dropped from `Project.toml`. (A one-off comparison against
  `MultivariateStats.ridge` is welcome in a test-time `using` block during
  the migration, but not retained.)
- The pre-construction memory log is computed from `(num_spinconfigs,
  num_salcs, n_atoms)` and the `Float64` element size; it does not require
  allocating the matrix first.

## Completion criteria

- [ ] `make test-all` passes (including new tests).
- [ ] `make test-jet` / `make test-aqua` show no new warnings.
- [ ] Pre-construction memory log shows the correct byte count on a small
      integration case (`test/integration/`).
- [ ] `PosDefException` on rank-deficient `X` is rethrown with a message
      that names `Ridge` and explains the likely cause.
- [ ] Cholesky-Ridge agrees with `MultivariateStats.ridge` (loaded once at
      test time only) to `1e-10` on a representative case.
- [ ] `.claude/bench_log.md` updated with before/after `solve_coefficients`
      timings on a representative case.
- [ ] OLS and Ridge docstrings updated to describe Cholesky.

## References

- Related commits: `0c10182` (the prior QR switch this spec partially
  reverses for OLS).
- Related code: `src/Fitting.jl:1217-1402` (the four `solve_coefficients`
  methods touched).
