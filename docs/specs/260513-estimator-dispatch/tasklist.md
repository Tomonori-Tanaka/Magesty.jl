# Tasklist: Optimize.jl estimator dispatch refactor

Status: **complete (2026-05-13)** — branch `refactor/estimator-dispatch`,
commits `5efbf1b..03168dc`.
Spec: [`requirements.md`](requirements.md) · [`design.md`](design.md)

Milestones (coarse-grained; daily TODOs use `TaskCreate`). Items below
are kept as-is for historical reference; the actual per-step progress
was tracked in the session-scoped task list.

## M0 — Baseline capture (before any code change)

- [ ] Run `make test-all` on `main` HEAD; record pass/fail.
- [ ] Run `julia --project test/benchmark_optimize.jl --input test/examples/fept_tetragonal_2x2x2/input.toml --with-fit --samples 20` on `main` HEAD; append to
      `.claude/bench_log.md` with the tag `pre-estimator-dispatch`.
- [ ] On a tiny synthetic fixture, capture golden `(j_values, j0, jphi)`
      from the current `elastic_net_regression` for both OLS-equivalent
      (α=0, λ=0) and ElasticNet (α=0, λ=0.1) paths. Hard-code these
      into the new regression test.

## M1 — Extract `assemble_weighted_problem` + `extract_j0_jphi`

- [ ] Move the "scale + concat + bias-column reset" block out of
      `elastic_net_regression` into a new `assemble_weighted_problem`
      function. Same body, same outputs.
- [ ] Move the `j0 = mean(...)` block into `extract_j0_jphi`.
- [ ] Keep `elastic_net_regression` calling these two helpers
      internally for one commit so the diff is reviewable.
- [ ] `make test-all` green.

## M2 — Introduce `solve_coefficients` dispatch (still `ElasticNet`)

- [ ] Define `solve_coefficients(::OLS, ...)` and
      `solve_coefficients(::ElasticNet, ...)` — keep the current name
      for one milestone so the diff is bounded.
- [ ] Switch `elastic_net_regression` body to: `assemble` →
      `solve_coefficients(ElasticNet(α, λ), ...)` → `extract`.
- [ ] `make test-all` green; golden regression test passes.

## M3 — Replace `isa` orchestrator

- [ ] Rewrite `_fit_sce_model_internal` as the 3-line orchestrator
      (`design.md` "Orchestrator" section).
- [ ] Convert `fit_sce_model_ols` and `fit_sce_model_elastic_net` to
      thin shims that build the estimator and call
      `_fit_sce_model_internal`.
- [ ] Delete `elastic_net_regression` (grep already confirmed no
      external callers).
- [ ] `make test-all` green.

## M4 — Rename `ElasticNet` → `Ridge` (breaking)

- [ ] In `src/Optimize.jl`:
      - [ ] Replace the `ElasticNet` struct + ctor with
            `struct Ridge <: AbstractEstimator; lambda::Float64 end`
            and the `Ridge(; lambda=0.0)` kwarg ctor.
      - [ ] Update `solve_coefficients(::ElasticNet, ...)` signature
            to `solve_coefficients(::Ridge, ...)`.
      - [ ] Rename `fit_sce_model_elastic_net` → `fit_sce_model_ridge`.
            The shim now builds `Ridge(lambda)` and emits a
            `@warn` if a non-zero `alpha` is still passed by an
            internal caller.
      - [ ] Update `Optimizer` outer constructor default to
            `Ridge(lambda=lambda)`. Emit a `@warn` once if
            `alpha != 0.0`.
      - [ ] Remove `ElasticNet` from `export`. Add `Ridge`.
- [ ] In `src/Magesty.jl`: update the re-export line
      (`ElasticNet` → `Ridge`).
- [ ] In `src/utils/ConfigParser.jl`:
      - [ ] Keep `Config4Optimize.alpha`, `:alpha => 0.0` default,
            and the validator.
      - [ ] In the constructor, after parsing, emit a one-shot
            `@warn` when `alpha != 0.0`.
      - [ ] Change the estimator construction at the end of the
            constructor from `ElasticNet(alpha, lambda)` to
            `Ridge(lambda)` (verify by reading `ConfigParser.jl`
            around L430).
- [ ] Update in-repo call sites:
      - [ ] `test/examples/fept_tetragonal_2x2x2/test.jl`
      - [ ] `test/benchmark_optimize.jl`
      - [ ] `tools/check_convergence_embset.jl`
- [ ] Update docs:
      - [ ] `SPEC.md` (lines 14, 79, 105 per grep).
      - [ ] `docs/src/api.md` (L43).
      - [ ] `docs/src/tutorial.md` (L118).
      - [ ] `docs/src/examples.md` (L42, L183).
- [ ] `make test-all` green; new TOML-depwarn unit test
      (`alpha = 0.5` emits exactly one warning) passes.

## M5 — Docstrings

- [ ] Add docstrings to `Ridge`, `assemble_weighted_problem`,
      `solve_coefficients`, `extract_j0_jphi` (English, Julia
      standard format `# Arguments` / `# Returns`).
- [ ] Update `Optimizer` outer constructor doc note re: `estimator`
      vs positional `alpha`/`lambda` and the deprecation path.

## M6 — Validation + benchmark

- [ ] `make test-unit` / `make test-integration` / `make test-jet` /
      `make test-aqua` all green.
- [ ] `julia --project test/benchmark_optimize.jl --input test/examples/fept_tetragonal_2x2x2/input.toml --with-fit --samples 20`; append result to `.claude/bench_log.md`
      with tag `post-estimator-dispatch`. Confirm no regression
      (within run-to-run noise).
- [ ] Compare `(j0, jphi)` from `test/examples/dimer` before vs after.
      Expected: bit-identical or within `atol=1e-12 rtol=1e-10`.

## M7 — Wrap-up

- [ ] Update `DESIGN_NOTES.md` estimator-dispatch section: add a
      "完了 → PR #..." line and note `ElasticNet → Ridge` rename.
- [ ] Update the user-facing API design memo in `DESIGN_NOTES.md`
      to swap `ElasticNet(...)` for `Ridge(...)` in code samples.
- [ ] Open PR `refactor: 3-layer estimator dispatch + ElasticNet→Ridge`
      against `main`. Body: link to this spec; describe before/after;
      summarize numerical-equivalence evidence; call out the breaking
      `ElasticNet` removal explicitly.
- [ ] After merge, follow-up specs may add `Lasso` / true mixed-norm
      `ElasticNet` / `Bayesian` — each as a separate spec under
      `docs/specs/[YYMMDD]-<estimator-name>/`.

## Out-of-scope reminders

- Do NOT touch `BasisSets.jl`, `Symmetry.jl`, `Structure.jl`, SALC
  ordering, or `MySphericalHarmonics.jl` in this branch.
- Do NOT add genuinely new estimator types (`Lasso`, true
  `ElasticNet`, `Bayesian`, ...) in this branch — they each get
  their own spec + PR.
- Do NOT remove `[regression].alpha` from the TOML schema in this
  branch (kept as a deprecated-but-accepted key; full removal is a
  follow-up).
- Do NOT introduce StatsAPI / GLM-style user-facing `fit` API —
  that is the user-facing API spec, a separate effort.
