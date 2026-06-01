# Tasklist: src architecture and hot-path cleanup

Status: merged (2026-06-01, PR #33) — M1 + M2 landed; M3 (A2) deferred by decision

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — A1: cache encapsulation

- [x] Add a peek accessor `CoupledBases.cached_coupling_results(ls, isotropy)`
      (returns built results or `nothing`; does not build — preserves the
      original `haskey(...) || continue` semantics exactly).
- [x] Repoint the `SALCBases` lookup to the accessor; drop the direct
      `_angular_momentum_cache` / `:none` reference.
- [x] `make test-all` green; committed (`6a7e359`).

### M2 — C: hot-path allocations

- [x] C5 `(false, true)` literal.
- [x] C3 / C4 `is_translationally_equivalent_coupled_basis` allocation removal
      (the reviewer's `find_translation_atoms` label was a misattribution; the
      flagged lines live in `is_translationally_equivalent_coupled_basis`).
- [x] C2 `projection_matrix_coupled_basis` in-place matrix write.
- [x] C1 `build_sh_cache_torque` single Legendre recursion via the new
      `TesseralHarmonics.Zₗₘ_grad_unsafe` combined value+gradient entry point
      (~2.4× faster per-atom SH fill; bit-identity verified directly).
- [x] `bench-salcbasis` + a direct `Zₗₘ_grad_unsafe` micro-benchmark recorded
      in `.claude/bench_log.md`; `make test-all` green; committed (`90f8755`).

### M3 — A2: Fitting kernel de-duplication (deferred)

- [x] Decision: **deferred, no commit.** Reading the three kernels showed the
      shared "skeleton" the review imagined does not cleanly exist:
  - `design_matrix_energy_element` uses a deliberately different optimized
    structure (last-site split with a reused partial product, scalar
    accumulation, `Z`-only). It cannot share the torque loop.
  - The two torque accumulators (`_accumulate_grad_torque_cluster!` /
    `_accumulate_grad_torque_scaled!`) are genuine near-twins (differ only in
    output buffer and the `mult` vs `mult*coeff` factor). They *could* be
    unified bit-identically (delegate `cluster!` to `scaled!` with `coeff=1.0`
    and a 2-D `view`, since `mult*1.0 === mult`).
  - But `_accumulate_grad_torque_cluster!` is on the package's hottest path
    (the threaded torque design-matrix build). Routing it through a `view`
    adds indirection there for a maintainability-only gain (~30 duplicated
    lines), which is not worth the regression risk. This is the
    `[contention: performance]` the performance reviewer flagged.
  - Revisit only if a future change already touches these kernels or a
    benchmark shows the `view` path is neutral.

## Exit checklist

- [x] `make test-all` passes.
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] ~~If results changed: regression or validation test added.~~ (results
      do not change; guarded by the existing suite + direct bit-identity
      checks for C1 and the SALC C-items)
- [x] ~~If public API changed: `SPEC.md` and `docs/src/api.md` updated.~~
      (internal-only; `Zₗₘ_grad_unsafe` is exported alongside the other
      hot-path `*_unsafe` variants but is not part of the documented public
      surface)
- [x] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved. Numerical: clean (0 findings, all
      bit-preserving). Two consensus majors (maintainability + API: the new
      `Zₗₘ_grad_unsafe` was missing from the module-docstring buffer list, and
      the export line exceeded 92 chars) plus docstring minors were applied.
      Performance minors (complex-power `^` → division; cold-path `collect`)
      left as out-of-scope follow-ups.
- [x] ~~If module names or Makefile targets changed: `.claude/agents/`
      swept.~~ (none)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hashes appended below.

## Implementation commits

- `f354976` docs(specs): add 260601-src-refactor spec
- `6a7e359` refactor(salc): encapsulate angular-momentum coupling cache access (M1 / A1)
- `90f8755` perf(salc,sh): cut hot-path allocations in SALC build and torque SH cache (M2 / C1–C5)
- M3 / A2: deferred (no commit; see M3 above)
