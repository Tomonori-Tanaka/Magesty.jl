# Requirements: src architecture and hot-path cleanup

Status: draft (2026-06-01)

## Goal

Remove two architectural smells and a set of hot-path allocations surfaced by
a multi-axis review of `src/`, without changing any numerical result or
on-disk format. Specifically: encapsulate the angular-momentum coupling cache
(drop the cross-module reach into `CoupledBases._angular_momentum_cache`),
de-duplicate the folded-tensor contraction skeleton shared by the energy and
torque design-matrix kernels, and cut avoidable allocations in the SALC and
SH-cache hot paths.

## Background

A review panel (maintainability + performance axes) over the whole `src/`
tree flagged these items. Lower-risk findings (docstrings, magic numbers,
internal-helper naming, the Sunny SALC-classification duplication) were
already handled directly on `main`. The remaining items are mid-sized:
they cross module boundaries (cache) or touch documented hot paths
(Fitting / SALCBases / TesseralHarmonics), so they go through this spec on a
`refactor/src-cleanup` branch.

## Scope

Includes:

- **A1 — cache encapsulation.** Replace `SALCBases.jl`'s direct access to
  `CoupledBases._angular_momentum_cache` (and the `:none` magic symbol /
  key tuple structure) with a public-internal `CoupledBases` accessor.
- **A2 — Fitting kernel de-duplication.** Factor the shared
  "folded-tensor index expansion + leave-one-out product" skeleton used by
  `design_matrix_energy_element` and the torque accumulators
  (`_accumulate_grad_torque_*`), *only if* it can be done without a hot-path
  performance regression.
- **C — hot-path allocations.** `build_sh_cache_torque` double Legendre
  recursion; `projection_matrix_coupled_basis` temporary-matrix allocation;
  `is_translationally_equivalent_coupled_basis` / `find_translation_atoms`
  per-call vectors; `[false, true]` → `(false, true)` literal.

Excludes:

- Public-API changes (e.g. removing `*_unsafe` from `export`, the `isotropy`
  positional→keyword move, the `kd_name` XML `<SpeciesOrder>` schema change).
  These are the review's "E" items and are proposed individually, each with
  explicit user confirmation, outside this spec.

## Invariants

- **Numerical results are bit-for-bit unchanged** for energy and torque
  design matrices, predictions, and the Sunny export.
- Existing XML `save` / `load` round-trips byte-for-byte.
- Real-tesseral `Zₗₘ` normalization and signs are untouched.
- Spin directions stay unit vectors with `3 × n_atoms` layout.
- The SALC key-group order inside `SALCBasis` (the outer index that
  `Fitting` and the XML I/O depend on) is preserved.

## Completion criteria

- [ ] `make test-all` passes.
- [ ] `make test-jet` / `make test-aqua` clean (no new findings).
- [ ] A1: no remaining reference to `CoupledBases._angular_momentum_cache`
      outside `CoupledBases`.
- [ ] C (and A2 if landed): before/after `@btime` medians recorded in
      `.claude/bench_log.md`; no measurable regression, ideally a measurable
      win on the targeted benchmarks.
- [ ] A2, if deferred instead of landed: the decision and reason are recorded
      in `tasklist.md` M3.
- [ ] Tier 2 four-axis review panel run and findings resolved.

## References

- Related issues / PRs: —
- Related specs / design notes: review panel output (this session).
