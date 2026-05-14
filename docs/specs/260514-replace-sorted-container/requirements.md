# Requirements: replace custom SortedContainer with plain + sort-once pattern

Status: draft (2026-05-14)
Owner: T. Tanaka
Branch: `refactor/replace-sorted-container`

## Goal

Remove the custom `src/common/SortedContainer.jl` module
(`SortedVector` / `SortedUniqueVector` / `SortedCountingUniqueVector`)
and migrate its call sites to a simpler "plain container during build
+ sort once at finalize" pattern, backed by `Base.Vector` and
`Base.Dict`. No tree-based data structure from `DataStructures.jl`
(`SortedDict`, `SortedMultiDict`) is adopted.

This direction was chosen after micro-benchmarking
(`test/benchmark_sorted_container.jl`) showed the plain pattern is
1.4×–47× faster than the current eagerly-sorted custom containers at
N = 1k–10k, while tree-based alternatives offer at best a 1.4–1.6×
gain at N = 10k and are slower at smaller N. The earlier
`docs/design-notes/replace-sorted-container.md` proposal assumed
`DataStructures.SortedMultiSet` exists; it does not.

## Scope

In scope:

- `src/common/SortedContainer.jl` — delete the entire module.
- `src/Clusters.jl`
  - `SortedVector{AtomCell}` in `interaction_cutoff_dict` and
    `interaction_clusters` → plain `Vector{AtomCell}` with a single
    `sort!` after build.
  - `SortedCountingUniqueVector{Vector{Int}}` in
    `irreducible_cluster_dict` → new `SortedCountedDict{T}` helper
    (or callers use bare `Dict{T,Int}` + cached sorted keys; final
    naming decided in design.md).
  - `Cluster.irreducible_cluster_dict` field type changes accordingly.
- `src/SALCBases.jl`
  - `SortedCountingUniqueVector{Basis.CoupledBasis}` in
    `SALCBasis.coupled_basislist` and the many intermediate
    `OrderedDict{Int, SortedCountingUniqueVector{...}}` carriers →
    same new helper type.
  - `SALCBasis.coupled_basislist` field type changes accordingly.
- `src/utils/xml_io.jl` — `build_sce_basis_from_xml` reconstructs the
  basis container; uses the new helper.
- `test/component_test/test_SortedContainer.jl` — delete.
- `test/component_test/test_Basis.jl` — drop the
  `using .SortedContainer: SortedCountingUniqueVector` reference if
  it's still needed and pick a substitute path.
- `test/benchmark_sorted_container.jl` — retained as the micro-bench
  record; updated to reflect the chosen design.

Out of scope:

- `SALCBasis.coupled_basislist` is technically a public field type
  (the struct is exported), but no documented external code reads
  `.coupled_basislist` directly outside Magesty itself. We accept
  this incidental type change; renaming the struct field or the
  basis container is part of the SCE public API spec
  (`docs/design-notes/sce-public-api.md`).
- Renaming or restructuring any module beyond what the migration
  forces.
- Optimization passes beyond replacing the container type (e.g.,
  reducing the number of intermediate dicts in `SALCBases.jl` is
  out of scope here).

## Invariants

1. **Numerical results unchanged.** `j_phi`, `j0`, design matrix,
   torque output, RMSE — all byte-identical with baseline.
2. **XML output byte-identical** on every integration example
   (`dimer`, `chain`, `fept`, `fege_2x2x2`, `square_lattice`,
   `febcc_pm`, `2d_fcc`). Iteration order of the basis container
   must match baseline.
3. **Sort order convention**: the new helper iterates keys in the
   same `isless` ordering the current `SortedCountingUniqueVector`
   uses (i.e. `Base.isless` on the key type). The
   `SortedCountingUniqueVector{Basis.CoupledBasis}` ordering relies
   on `Basis.CoupledBasis`'s `isless` definition (sorted by
   `(atoms, Lseq, Lf, ...)` — confirmed in `src/types/Basis.jl`).
4. **Multiplicity semantics unchanged.** `getcount(scuv, k)`
   returning `0` for absent keys is preserved on the new helper
   under whatever name we pick.
5. **`Cluster` and `SALCBasis` public constructors keep the same
   call signatures.** Only internal field types change.
6. **Physics conventions unchanged.** No SH / SALC / CG / sign /
   normalization touched.

## Non-goals

- No new `DataStructures.jl` types are pulled in for this work
  (Dict/Vector are sufficient). The `DataStructures` dependency is
  still needed elsewhere (`OrderedDict` in `SALCBases.jl`,
  `OrderedDict` introduced in `Clusters.jl` from the prior
  CountingContainer refactor).
- No deprecation shim. `SortedVector` / `SortedUniqueVector` /
  `SortedCountingUniqueVector` are not exported from `Magesty`, so
  removing them does not break documented public surface.
- No new exports.

## Completion criteria

- [ ] `src/common/SortedContainer.jl` deleted; `Magesty.jl` and
      `test/runtests.jl` no longer include it.
- [ ] All caller files migrated; `grep -r "SortedVector\\|SortedUniqueVector\\|SortedCountingUniqueVector"`
      in `src/` and `test/` returns no hits.
- [ ] `make test-unit`, `make test-integration`, `make test-jet`,
      `make test-aqua` all green.
- [ ] Integration test XML diff vs baseline (current `main`):
      byte-identical on the seven example directories listed above.
- [ ] `test/benchmark_sorted_container.jl` updated to bench the
      chosen design (old vs new) and the speedup recorded in
      `.claude/bench_log.md` and `docs/design-notes/replace-sorted-container.md`.
- [ ] `DESIGN_NOTES.md` index and `docs/design-notes/replace-sorted-container.md`
      marked complete; refactor-sweep R8 Plan C dropped (resolved
      by this refactor).

## Risk

| Risk | Mitigation |
|------|------------|
| Iteration order drift breaks XML byte-identical | The new helper sorts keys lazily via `Base.sort` on the same key type as before. Add an explicit round-trip XML diff test in CI for `fege_2x2x2` and `fept` before merge. |
| `SALCBasis.coupled_basislist` field type change breaks an out-of-repo caller | The field is not documented as a stable accessor (`SPEC.md` describes `SALCBasis` at the constructor level). The diff is documented in the commit message and DESIGN_NOTES. SCE public API spec will codify accessor methods separately. |
| Performance regression on small inputs (overhead of finalize step) | Micro-bench results already show wins down to N=100. Re-confirm on `dimer` / `chain` integration tests during PR. |
| Loss of `findfirst` / `findall` / `delete!` / `clear!` semantics | These methods are used only inside `SortedContainer.jl` self and inside the tests we are deleting. Verified by grep before merge. |
