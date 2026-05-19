# Requirements: pre-release safety net

Status: draft (2026-05-19)

## Goal

Close three pre-v0.1.0 safety gaps surfaced by a package-wide review:
silent basis mismatch on reload, missing spin unit-vector validation,
and missing integration coverage for `SCEModel`, batch `predict_*`, and
cross-basis mismatch detection.

## Background

A package-wide extensibility review on 2026-05-19 surfaced a small set of
items that should land before v0.1.0:

- `_check_basis` (`src/Magesty.jl:894-904`, with the `SCEFit` delegate at line 905) compares `SCEBasis` by `===`
  identity, so a basis reloaded from disk that happens to have the same
  SALC count as the one used to fit (but a different SALC ordering) will
  pass the check and produce silently wrong predictions. The XML
  round-trip is supposed to preserve ordering, but the only safeguard
  today is the linked-site discipline documented in `CLAUDE.md`.
- `SpinConfig` only validates the `3 × n_atoms` shape; the documented
  unit-vector invariant is enforced nowhere. A non-normalized direction
  silently breaks `Zₗₘ` semantics.
- The square_lattice integration test (`test/integration/square_lattice/test.jl`)
  currently exercises `SCEModel` only via direct construction
  (`SCEModel(basis, j0, jphi_list)`); the `SCEModel(fit::SCEFit)` path
  that production code uses is not covered by any integration test.
  The `AbstractVector{SpinConfig}` / `AbstractVector{<:AbstractMatrix}`
  `predict_*` overloads are covered at the component level
  (`test/component/test_SCEFit.jl`) but not at integration level. The
  cross-basis mismatch error in `_check_basis`
  (`src/Magesty.jl:894-904`) has no test at any level.

Other pre-release findings tracked in the design-note index under the
"pre-release cleanup" note — broken example paths, Documenter
`checkdocs` tightening, `predict_*` docstring section normalization,
`CITATION.cff`, `Project.toml` version / compat tidy-up — remain in
that note and ship in a separate PR. The one item in that note that
overlaps with this spec (an `SCEModel(fit)` round-trip in
`square_lattice/test.jl`) is **moved into this spec** and removed from
the cleanup note as part of M5, so the two documents do not
double-own it.

Larger-scope items (I/O backend abstraction, observable trait,
SOC / complex harmonics support, CG cache eviction, `InputSpecs`
extension hatch) are out of scope and will get their own specs /
notes if and when they are triggered.

## Scope

Includes:

- Add a `salc_fingerprint::UInt64` (or equivalent) field to `SCEBasis`
  computed deterministically from the key-group ordering of
  `salcbasis.salc_list`. The fingerprint is in-memory metadata: it is
  recomputed from the reconstructed `salcbasis` on every XML load and
  is not part of the on-disk schema.
- Compare the fingerprint (alongside the existing `===` fast-path)
  inside `_check_basis`.
- Validate unit-norm spin directions in `SpinConfig`'s inner constructor.
- Add integration tests for `SCEModel` round-trip prediction, the
  batched `predict_*` overloads, and the cross-basis mismatch error.
- Update `CHANGELOG.md`, `SPEC.md`, and `docs/src/api.md` for the new
  fields / behavior.

Excludes:

- Persisting `torque_weight` on `SCEModel` or in model XML. Audit /
  reproducibility metadata is fit-time hyperparameter information; if
  v0.1.0 starts encoding it on `SCEModel`, every future estimator
  (Ridge `lambda`, AdaptiveLasso `gamma` / `lambda_grid` /
  `pilot_strategy`, ...) will inevitably want the same treatment, and
  the model XML schema would grow estimator-specific fields without
  bound. Users who need the weighting can keep the `SCEFit` object
  alongside the model (it already carries `torque_weight`) or record
  the value in the TOML input. A future spec may introduce a dedicated
  provenance side-channel if this becomes a recurring need.
- Any new fields on the model / basis XML schema. The existing payload
  is unchanged by this spec; `salc_fingerprint` lives only in memory
  (recomputed on load).
- `examples/01_basic_flow.jl` / `examples/03_save_load.jl` path fixes
  (trivial; tracked in the pre-release cleanup design-note).
- Documenter `checkdocs = :exports`, `predict_*` docstring section
  normalization, `CITATION.cff`, `Project.toml` version / compat
  tidy-up (all tracked in the pre-release cleanup design-note).
- I/O backend abstraction, observable trait, CG cache eviction,
  `InputSpecs` extension hatch, SOC / complex-harmonics support.
- Any change to SALC ordering itself.

## Invariants

- Existing XML round-trips remain byte-identical: this spec adds no
  new on-disk fields. The `salc_fingerprint` lives in memory only and
  is recomputed on every load.
- SALC ordering inside `salcbasis.salc_list` is unchanged.
- Numerical results of `predict_energy` / `predict_torque` are unchanged
  for any basis + model pair that was valid before this change.
- Spin direction layout stays `3 × n_atoms`; the only new constraint is
  per-column `‖·‖ ≈ 1`.
- Public API signatures are not removed; new keyword arguments and new
  read-only fields only.

## Completion criteria

- [ ] `make test-all` passes (component + integration).
- [ ] New component test: `SpinConfig` rejects non-unit columns.
- [ ] New component test: `_check_basis` raises on a basis with the
      same SALC count but a different fingerprint.
- [ ] New `@testset` added to `test/integration/square_lattice/test.jl`
      that goes through `SCEFit` and verifies
      `predict_energy(SCEModel(fit), dataset)` agrees with
      `predict_energy(fit, dataset)` to machine precision.
- [ ] At integration level, the same `@testset` also exercises the
      batched `predict_*` overloads (`AbstractVector{SpinConfig}` and
      `AbstractVector{<:AbstractMatrix}`) against the per-config /
      dataset paths. Component-level coverage already exists in
      `test/component/test_SCEFit.jl:106-141`; the integration test
      verifies the same overloads on a realistic basis.
- [ ] XML round-trip test verifies that the in-memory
      `salc_fingerprint` of a reloaded basis equals the original (i.e.
      the recompute on load reproduces the same value), and that the
      reloaded `SCEModel` still passes `_check_basis` against the
      original dataset.
- [ ] `CHANGELOG.md` `[Unreleased]` updated. The XML schema is
      unchanged; note the fingerprint-based safety check and the
      stricter `SpinConfig` validation only.
- [ ] `SPEC.md` and `docs/src/api.md` reflect the new fields.

## References

- Related: the pre-release cleanup design note (broken examples,
  documenter strictness, predict_* docstring normalization,
  `CITATION.cff`, `Project.toml` cleanup) — separate PR; this spec's
  M5 removes the `SCEModel(fit)` item from that note to keep
  ownership single.
- Linked sites (see CLAUDE.md): SCE coefficient I/O, Fitting ↔ SALCBasis
- Source touchpoints:
  - `src/Magesty.jl:125-130` (`SCEBasis` struct — no inner constructor)
  - `src/Magesty.jl:161, 1234, 1241` (existing 4-arg `SCEBasis(...)`
    construction sites)
  - `src/Magesty.jl:894-904` (`_check_basis(::SCEModel, ::SCEDataset)`
    method body; the `SCEFit` delegate at line 905 stays untouched)
  - `src/SpinConfigs.jl:101-157` (`SpinConfig` inner constructor —
    dimension and non-negative-moment checks only)
  - `src/SALCBases.jl` (`SALCBasis.salc_list` key-group ordering,
    source of fingerprint)
