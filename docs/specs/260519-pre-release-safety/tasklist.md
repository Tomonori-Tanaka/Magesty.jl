# Tasklist: pre-release safety net

Status: draft (2026-05-19)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Spin unit-vector validation

- [ ] Add per-column norm check + `atol_unit_norm` keyword to the
      `SpinConfig` inner constructor (`src/SpinConfigs.jl`).
- [ ] Add a component test verifying rejection of a non-unit column and
      acceptance under a loose tolerance.
- [ ] `make test-unit` clean.

### M2 — Basis fingerprint

- [ ] Implement `salc_fingerprint(::SALCBasis)::UInt64` in
      `src/SALCBases.jl` using `Base.hash` over the wire tuple
      pinned in `design.md` "Types and conventions". The tuple per
      `coupled::CoupledBasis_with_coefficient` is
      `(coupled.ls, coupled.Lf, coupled.Lseq, coupled.atoms,
      coupled.multiplicity, length(coupled.coefficient))`; iteration
      order is outer key-group index, then inner SALC index;
      seed with `hash(:Magesty_SALC_fingerprint_v1)` so the recipe
      can be versioned later without colliding with old saved
      values.
- [ ] Add `salc_fingerprint::UInt64` as the 5th field of `SCEBasis`
      (`src/Magesty.jl:125-130`). Add a 4-arg outer convenience
      constructor that computes the fingerprint via
      `salc_fingerprint(salcbasis)` and forwards to the default
      5-arg constructor; existing call sites
      (`src/Magesty.jl:161, 1234, 1241`) keep their current shape.
- [ ] Restructure `_check_basis` (`src/Magesty.jl:894-904`; the
      current single-expression `model.basis === dataset.basis \|\|
      throw(...)` becomes two short-circuits: identity, then
      fingerprint, then throw) so that
      `model.basis.salc_fingerprint == dataset.basis.salc_fingerprint`
      is checked after the existing `===` fast-path. Signature
      unchanged. The `SCEFit` delegate at line 905 is left alone.
- [ ] Add a component test (new file
      `test/component/test_fit_basis_check.jl`) that builds a second
      basis with a forced (mismatching) fingerprint via the default
      5-arg constructor and verifies `predict_energy` raises.

### M3 — `torque_weight` persistence — Deferred

Dropped from this spec (decision recorded 2026-05-20). Rationale:
adding `torque_weight` to `SCEModel` and the model XML schema starts a
slippery slope of estimator-specific audit fields (Ridge `lambda`,
AdaptiveLasso `gamma` / `lambda_grid` / `pilot_strategy`, ...) that
would inevitably want the same treatment, and the model XML would grow
estimator-specific nodes without bound. For v0.1.0, users who need the
weighting can keep the `SCEFit` object alongside the model (it already
carries `torque_weight`) or record the value in the TOML input. A
future spec may introduce a dedicated provenance side-channel if this
becomes a recurring need. See `requirements.md` "Excludes".

### M4 — Integration test gaps

- [ ] **Build an `SCEDataset` in `square_lattice/test.jl`** —
      currently the file has no `SCEDataset` and no `SpinConfig`
      list; it only calls `SCEModel(basis, 0.0, jphi_list)` directly.
      Synthesize a small list of `SpinConfig`s for the FM + AFM
      configurations (and a handful of intermediate ones to make the
      `fit` well-posed), with per-config energies computed from
      `E = Σ_{<ij>} Jij * Si·Sj` for the existing `J = -1.0 eV`
      Heisenberg model and per-atom `local_magfield_i = -m_i *
      Σ_j J_ij * S_j`. Follow the `dimer` / `chain` /
      `febcc_2x2x2_pm` patterns under `test/integration/` for the
      `SCEDataset(basis, configs)` call shape. The new fixture is
      contained inside the new `@testset` (no new files; no EMBSET
      file required for this short list).
- [ ] Add a new `@testset` to
      `test/integration/square_lattice/test.jl` (next to the
      existing "Hypothetical SCE model (isotropic)" block) that
      consumes the dataset above, calls
      `fit(SCEFit, dataset, OLS(); torque_weight = 0.5,
      verbosity = false)`, and verifies
      `predict_energy(SCEModel(fit), dataset)` agrees with
      `predict_energy(fit, dataset)` to machine precision. Do not
      touch the existing direct-`SCEModel`-construction testset.
- [ ] In the same testset, exercise the batched `predict_*` overloads
      (`AbstractVector{SpinConfig}` and
      `AbstractVector{<:AbstractMatrix}`) against the dataset path.
- [ ] In the same testset, verify `Magesty.save(SCEModel(fit), ...)`
      followed by `Magesty.load(SCEModel, ...)` produces a model whose
      `basis.salc_fingerprint` (recomputed on load) equals the
      original, and whose `predict_energy(reloaded_model, dataset)`
      still passes `_check_basis` and reproduces the in-memory
      prediction to machine precision (covers the load-then-predict
      path that motivated the fingerprint).

### M5 — Docs, changelog, design-note alignment

- [ ] Update the `SPEC.md` `SCEBasis` type diagram to show the new
      `salc_fingerprint::UInt64` field (`SCEModel` is unchanged by
      this spec).
- [ ] Update `docs/src/api.md` with the new field and the strict
      spin-norm invariant.
- [ ] `CHANGELOG.md` `[Unreleased]` entry: stricter `SpinConfig`
      validation and in-memory fingerprint-based basis check. The
      on-disk XML schema is unchanged.
- [ ] Remove the `SCEModel(fit)` round-trip item from the
      pre-release cleanup design note and add a one-line
      cross-reference into this spec, so that note keeps single
      ownership of broken examples / Documenter / docstring
      tidy-up only.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes (component + integration).
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
      These are not part of `make test-all`; run them as separate
      targets (also covered by `make ci-local`).
- [ ] If results changed: regression or validation test added.
      (Results should NOT change for valid inputs; this spec is
      strictly safety-net work.)
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ — no hot path touched.
- [ ] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ — no such changes.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
