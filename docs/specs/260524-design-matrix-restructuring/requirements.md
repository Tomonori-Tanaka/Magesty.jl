# Requirements: Design-matrix algorithmic restructuring

Status: draft (2026-05-24)

## Goal

Reduce the wall-time of `build_design_matrix_energy` and
`build_design_matrix_torque` in `src/Fitting.jl` — currently the
dominant cost of `fit` — through four staged algorithmic
restructurings: (A) folding the SALC `coefficient` vector into the
coupled-basis tensor at construction time, (C) pre-enumerating orbit
clusters at SALC-build time, (B) caching tesseral-harmonic values per
spinconfig, and (D) recasting the torque computation as a reverse-mode
Jacobian over clusters. Each stage is gated by a benchmark on the
FeGe B20 2x2x2 fixture before proceeding to the next.

## Background

Previous specs `260516-coupled-basis-typeparam`,
`260516-optimize-workspace`, and `260523-design-matrix-3body-perf`
exhausted micro-optimizations within the current loop structure
(rank type-parameterization, scratch pooling, dedup data-structure
swap, SubArray removal). The remaining wins are structural:
reordering loops and lifting work out of the inner kernel. See the
companion design note
[`docs/design-notes/design-matrix-algorithmic-restructuring.md`](../../design-notes/design-matrix-algorithmic-restructuring.md)
for the full analysis, reviewer findings, and the per-stage
tolerance table.

Baseline (FeGe B20 2x2x2, 100 spinconfigs, 4 threads, commit
`4809a5a`, from `.claude/bench_log.md`):

| function | time | allocations |
|---|---:|---:|
| `build_design_matrix_energy` | 1.23 s | 7.02e6 |
| `build_design_matrix_torque` | 10.11 s | 3.36e8 |

The 10:1 torque:energy ratio is structural; D removes the dominant
`num_atoms / N` excess factor while C removes the historical ~20–25%
dedup overhead.

## Scope

Includes:

- (A) `Mf` folding into `CoupledBasis_with_coefficient` at SALC
  build time, dropping the `mf_idx` loop from the hot kernel.
- (C) Pre-enumeration of orbit clusters at SALC-build time, dropping
  `searched_pairs` / `_atoms_hash_key` / per-translation sort
  buffers from the hot kernel.
- (B) Per-spinconfig Zₗₘ / ∂ᵢZₗₘ table cache built in
  `build_design_matrix_*`.
- (D) Cluster-major loop reorder for `build_design_matrix_torque`
  with reverse-mode Jacobian accumulation per `(sc_idx, atom,
  salc_idx)` and a single `cross` per cell.
- Updates to `src/XMLIO.jl` save/load to persist the new on-disk
  shape (and schema-version bump if needed; see §A in design.md).
- Per-stage benchmark entries in `.claude/bench_log.md` (FeGe 2x2x2
  fixture, 5-trial median, before/after).

Excludes:

- (E) SIMD-friendly contraction (`SMatrix` / Tullio). Deferred to a
  follow-up spec once kernel shape stabilizes after A + B.
- (F) Concrete-type `salc_list` (eliminate dispatch boxing).
  Deferred; D drops the boxing frequency by `O(num_atoms / N)`, so
  F's priority is re-evaluated post-D.
- GPU offload, batched contraction across spinconfigs, sparsity-
  aware contraction.

## Invariants

- **Spherical-harmonics convention** (CLAUDE.md "Linked sites"):
  `Zₗₘ` / `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` normalizations and signs
  unchanged. The `test/component/test_sphericart_agreement.jl`
  agreement holds throughout.
- **SCE design-feature formula** in `docs/src/technical_notes.md`:
  the value of every design-matrix entry is mathematically
  identical to the current implementation, up to the per-stage
  floating-point tolerance documented in design.md.
- **Spin convention**: directions stay unit-vector with `3 × n_atoms`
  layout. The unit-vector constraint and the `cross(spin, ∇ₑu)`
  tangent-projection are preserved by D's reordering (cross
  linearity argument).
- **Pure-translation symmetry in hot path**: the design-matrix hot
  path's symmetry application is pure translation (verified at
  `src/Symmetries.jl` `symnum_translation` and at the
  `@inbounds for itrans in symmetry.symnum_translation` loops in
  `src/Fitting.jl` `design_matrix_energy_element` and
  `calc_∇ₑu!`). B's cache is only valid under this invariant.
- **Cluster-atom distinctness in D**: every translated cluster has
  pairwise-distinct atom indices. Enforced today by
  `Clusters.jl`'s `allunique(atom_list_all)` filter and by
  symmetry-operation bijectivity. D's gradient derivation
  depends on this; the spec must assert it in a debug build.
- **XML round-trip**: existing `SCEBasis` save/load files must
  either continue to round-trip byte-for-byte (strategy a), or
  the schema-version bump must make mismatched files reject loudly
  (strategy b — silent misinterpretation must be impossible).
  See design.md §A for the chosen strategy.
- **SALC key-group order** (CLAUDE.md "Linked sites: Fitting and
  SALCBasis"): the outer index of
  `Vector{Vector{CoupledBasis_with_coefficient}}` is preserved; the
  physical meaning of `Jφ` does not change.
- **Energy units**: no change to the eV / `j0` split or to the
  `(4π)^(N/2)` `_cluster_scaling` factor.

## Completion criteria

Per-stage gates (must hold before the next stage is started):

- [ ] **After A**: `make test-all` passes. Pre/post design matrices
      match at `rtol = 1e-14, atol = 1e-15` on FeGe + FePt
      fixtures (bit-identical if `Mf` accumulation order is
      preserved). Bench entry in `.claude/bench_log.md`.
- [ ] **After C**: `make test-all` passes. Pre/post match at
      `rtol = 1e-13, atol = 1e-14`. Torque speedup ≥ 1.15× over
      baseline (the dedup-removal floor; spec `260523` M2 measured
      ~9% from `Set` → `Vector` alone). Bench entry.
- [ ] **After B**: `make test-all` passes. Pre/post match at
      `rtol = 1e-14, atol = 1e-15` (no contraction-order change;
      caching is exact). Energy speedup ≥ 2× over post-C baseline.
      Bench entry.
- [ ] **After D**: `make test-all` passes. Pre/post match at
      `rtol = 1e-12, atol = 1e-13` (energy block) and
      `rtol = 1e-10, atol = 1e-12` (torque block). Torque speedup
      ≥ 5× over post-B baseline. Bench entry.

Project-level gates:

- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] End-to-end fit results on FeGe + FePt integration tests
      unchanged within the per-stage tolerance.
- [ ] Tier 2 review panel (numerical / maintainability /
      performance / API) run after D lands, findings resolved.
- [ ] `SPEC.md` updated if `CoupledBasis_with_coefficient` field
      layout changes (strategy b in §A).
- [ ] `CHANGELOG.md` `[Unreleased]` updated.

## References

- Design note:
  [`docs/design-notes/design-matrix-algorithmic-restructuring.md`](../../design-notes/design-matrix-algorithmic-restructuring.md)
- Prior specs: `260516-coupled-basis-typeparam`,
  `260516-optimize-workspace`, `260523-design-matrix-3body-perf`.
- Math: `docs/src/technical_notes.md`.
- Baseline benchmarks: `.claude/bench_log.md`.
