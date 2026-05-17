# Tasklist: Optimize.jl estimator dispatch refactor

Status: **complete (2026-05-13)** — branch `refactor/estimator-dispatch`,
commits `5efbf1b..03168dc`.
Spec: [`requirements.md`](requirements.md) · [`design.md`](design.md)

Milestones (coarse-grained; daily TODOs use `TaskCreate`). Items below
are kept as-is for historical reference; the actual per-step progress
was tracked in the session-scoped task list.

## M0 — Baseline capture (before any code change)  ✅

- [x] Run `make test-all` on `main` HEAD; record pass/fail.
- [x] Run `julia --project test/benchmark_optimize.jl --input test/examples/fept_tetragonal_2x2x2/input.toml --with-fit --samples 20` on `main` HEAD; append to
      `.claude/bench_log.md` with the tag `pre-estimator-dispatch`.
      (Initially captured at `--samples 5` in commit `6ab646e`;
      re-measured at `--samples 20` in M6 for fair comparison.)
- [x] On a tiny synthetic fixture, capture golden `(j_values, j0, jphi)`
      from the current `elastic_net_regression` for both OLS-equivalent
      (α=0, λ=0) and ElasticNet (α=0, λ=0.1) paths. Hard-code these
      into the new regression test.

## M1 — Extract `assemble_weighted_problem` + `extract_j0_jphi`  ✅

- [x] Move the "scale + concat + bias-column reset" block out of
      `elastic_net_regression` into a new `assemble_weighted_problem`
      function. Same body, same outputs.
- [x] Move the `j0 = mean(...)` block into `extract_j0_jphi`.
- [x] Keep `elastic_net_regression` calling these two helpers
      internally for one commit so the diff is reviewable. (commit
      `92923fc`.)
- [x] `make test-all` green.

## M2 — Introduce `solve_coefficients` dispatch (still `ElasticNet`)  ✅

- [x] Define `solve_coefficients(::OLS, ...)` and
      `solve_coefficients(::ElasticNet, ...)` — keep the current name
      for one milestone so the diff is bounded.
- [x] Switch `elastic_net_regression` body to: `assemble` →
      `solve_coefficients(ElasticNet(α, λ), ...)` → `extract`.
- [x] `make test-all` green; golden regression test passes.

## M3 — Replace `isa` orchestrator  ✅

- [x] Rewrite `_fit_sce_model_internal` as the 3-line orchestrator
      (`design.md` "Orchestrator" section).
- [x] Convert `fit_sce_model_ols` and `fit_sce_model_elastic_net` to
      thin shims that build the estimator and call
      `_fit_sce_model_internal`.
- [x] Delete `elastic_net_regression`. **Spec correction**: the initial
      grep claim "no external callers" was wrong; two test files
      (`test/examples/{fege_2x2x2,febcc_2x2x2_pm}/test.jl`) used it
      directly. Both were updated to `_fit_sce_model_internal(...,
      OLS(), weight)` since they passed `alpha=lambda=0`.
- [x] `make test-all` green.

## M4 — Rename `ElasticNet` → `Ridge` (breaking)  ✅

- [x] In `src/Optimize.jl`:
      - [x] Replace the `ElasticNet` struct + ctor with
            `struct Ridge <: AbstractEstimator; lambda::Float64 end`
            and the `Ridge(; lambda=0.0)` kwarg ctor.
      - [x] Update `solve_coefficients(::ElasticNet, ...)` signature
            to `solve_coefficients(::Ridge, ...)`.
      - [x] Rename `fit_sce_model_elastic_net` → `fit_sce_model_ridge`.
            **Spec deviation**: the `alpha` positional arg was dropped
            entirely (not kept with a warn) because grep showed zero
            in-repo callers after M3 — keeping it would have been
            dead-arg pollution. New signature:
            `fit_sce_model_ridge(Xe, Xt, ye, yt, lambda, weight)`.
      - [x] Update `Optimizer` outer constructor default to
            `Ridge(lambda=lambda)`. Emit a `@warn maxlog=1` if
            `alpha != 0.0` via the new `_default_estimator(alpha,
            lambda)` helper.
      - [x] Remove `ElasticNet` from `export`. Add `Ridge`.
- [x] In `src/Magesty.jl`: update the re-export line
      (`ElasticNet` → `Ridge`).
- [x] In `src/utils/ConfigParser.jl`:
      - [x] Keep `Config4Optimize.alpha`, `:alpha => 0.0` default,
            and the validator. The comment on
            `DEFAULT_VALUES_OPTIMIZE[:alpha]` was also updated to
            note the deprecation.
      - [x] In the constructor, after parsing, emit a one-shot
            `@warn maxlog=1` when `alpha != 0.0`.
      - [x] Config4Optimize does NOT itself build an estimator; the
            estimator construction happens in `Optimizer`'s outer
            ctor and was switched there (see above).
- [x] Update in-repo call sites:
      - [x] `test/examples/fept_tetragonal_2x2x2/test.jl`
      - [x] `test/benchmark_optimize.jl`
      - [x] `tools/check_convergence_embset.jl`
- [x] Update docs:
      - [x] `SPEC.md`
      - [x] `docs/src/api.md`
      - [x] `docs/src/tutorial.md`
      - [x] `docs/src/examples.md`
      - [x] `.claude/agents/test-runner.md` (bonus catch via grep)
- [x] `make test-all` green; new TOML-depwarn unit test
      (`alpha = 0.5` emits exactly one warning) passes (35 → 36
      ConfigParser tests).

## M5 — Docstrings  ✅

- [x] Add docstrings to `Ridge`, `assemble_weighted_problem`,
      `solve_coefficients`, `extract_j0_jphi` (English, Julia
      standard format `# Arguments` / `# Returns`). Per-method
      docstring also added to `solve_coefficients(::Ridge, ...)` to
      document the `λ ≈ 0` short-circuit; `::OLS` left undocumented
      (trivial body).
- [x] Update `Optimizer` outer constructor doc note re: `estimator`
      vs positional `alpha`/`lambda` and the deprecation path.
      (Existing `Optimizer` struct had no docstring; added a full
      `# Fields` / `# Arguments` / `# Keyword arguments` block.)

## M6 — Validation + benchmark  ✅

- [x] `make test-unit` / `make test-integration` / `make test-jet` /
      `make test-aqua` all green.
- [x] Bench command appended to `.claude/bench_log.md` with tag
      `post-estimator-dispatch`. Within ±1% of baseline at
      `--samples 20` (240.023 ms vs 242.456 ms median). The
      `--samples 5` trial ran into inter-run thermal noise — see
      bench_log caveat.
- [x] Numerical equivalence on integration examples confirmed via
      `make test-integration` (155/155 unchanged). Golden values on
      the synthetic dispatch fixture are bit-identical.

## M7 — Wrap-up  ✅

- [x] Update `DESIGN_NOTES.md` estimator-dispatch section: added a
      "complete" status header noting the rename and bench result.
- [x] Update the user-facing API design memo in `DESIGN_NOTES.md`
      to swap `ElasticNet(lambda=...)` for `Ridge(lambda=...)` in the
      code samples. The "other estimators" placeholder list still
      mentions `ElasticNet(alpha=0.5, lambda=0.1)` deliberately as
      a reference to the *future* genuine mixed-norm estimator the
      name is now reserved for.
- [x] Push `refactor/estimator-dispatch` to `origin`.
- [x] Open PR `refactor(optimize): 3-layer estimator dispatch +
      ElasticNet→Ridge rename` against `main` →
      [PR #1](https://github.com/Tomonori-Tanaka/Magesty.jl/pull/1).
      Ubuntu CI green (`Julia 1 - ubuntu-latest`, `Documentation`).
- [x] PR merged via merge commit
      [`4f0eefd`](https://github.com/Tomonori-Tanaka/Magesty.jl/commit/4f0eefd).
      Local branch deleted.

### Follow-up (not part of this spec)

After this merge, new estimators (`Lasso` / true mixed-norm
`ElasticNet` / `Bayesian` / NNLS / `RegularizationPath`, ...) can be
added each as a separate spec under
`docs/specs/[YYMMDD]-<estimator-name>/`. Each is a `<: AbstractEstimator`
struct plus one `solve_coefficients` method.

## Out-of-scope reminders

- Do NOT touch `SALCBases.jl`, `Symmetry.jl`, `Structure.jl`, SALC
  ordering, or `MySphericalHarmonics.jl` in this branch.
- Do NOT add genuinely new estimator types (`Lasso`, true
  `ElasticNet`, `Bayesian`, ...) in this branch — they each get
  their own spec + PR.
- Do NOT remove `[regression].alpha` from the TOML schema in this
  branch (kept as a deprecated-but-accepted key; full removal is a
  follow-up).
- Do NOT introduce StatsAPI / GLM-style user-facing `fit` API —
  that is the user-facing API spec, a separate effort.
