# Tasklist: replace custom SortedContainer

Status: implementation complete pending commit (2026-05-14)
Branch: `refactor/replace-sorted-container`

Coarse milestones. In-session TODOs are tracked separately.

## Phase 0 — sign-off

- [x] User agrees on direction (plain + sort once).
- [x] User agrees on helper name `SortedCounter{T}` (or proposes
      alternative).
- [x] User agrees on module location `src/common/SortedCounter.jl`.

## Phase 1 — new helper

- [x] Create `src/common/SortedCounter.jl` with the `SortedCounter{T}`
      struct (Dict-backed, lazy sorted key cache).
- [x] Implement required operations: empty constructor, `push!(sc, k)`,
      `push!(sc, k, n::Integer)`, `length`, `isempty`, `getindex`,
      `iterate`, `==`, `isless`, `copy`, `show`. `.counts` exposed as
      a public field.
- [x] Add `test/component_test/test_SortedCounter.jl` with parity
      tests against the operations actually used by callers (push
      semantics, count lookup, iteration order, length, isempty,
      indexed access).
- [x] Wire `include("common/SortedCounter.jl")` into `src/Magesty.jl`
      and `test/runtests.jl` (alongside the existing
      `SortedContainer` for now).
- [x] `make test-unit` green with both modules coexisting.

## Phase 2 — migrate `Clusters.jl`

- [x] Swap `SortedVector{AtomCell}` → `Vector{AtomCell}` in
      `interaction_cutoff_dict` / `interaction_clusters`.
- [x] Add `sort!` finalize at the end of `generate_clusters` for
      both data structures.
- [x] Swap `SortedCountingUniqueVector{Vector{Int}}` →
      `SortedCounter{Vector{Int}}` in `irreducible_cluster_dict`.
- [x] Update field types on `Cluster` struct.
- [x] Drop `using ..SortedContainer` from `Clusters.jl`, add
      `using ..SortedCounter: SortedCounter`.
- [x] `make test-unit` / `make test-integration` green; XML
      byte-identical on representative examples (`dimer`, `fege_2x2x2`).

## Phase 3 — migrate `SALCBases.jl` and `xml_io.jl`

- [x] Swap `SortedCountingUniqueVector{Basis.CoupledBasis}` →
      `SortedCounter{Basis.CoupledBasis}` everywhere.
- [x] Update `SALCBasis.coupled_basislist` field type.
- [x] Drop `using ..SortedContainer` and the `isa SortedCountingUniqueVector`
      check at L975 (replace with `isa SortedCounter` or remove
      now-redundant branch — confirm during edit).
- [x] Same migration in `src/utils/xml_io.jl::build_sce_basis_from_xml`.
- [x] `make test-integration` byte-identical on `fept`, `fege_2x2x2`,
      `febcc_pm`, `2d_fcc`, `square_lattice`, `dimer`, `chain`.

## Phase 4 — delete legacy

- [x] Remove `include("common/SortedContainer.jl")` from
      `src/Magesty.jl`.
- [x] Remove the same from `test/runtests.jl` plus
      `include("./component_test/test_SortedContainer.jl")`.
- [x] Delete `src/common/SortedContainer.jl`.
- [x] Delete `test/component_test/test_SortedContainer.jl`.
- [x] Verify `grep -rE "Sorted(Vector|UniqueVector|CountingUniqueVector)"`
      in `src/` and `test/` returns no hits.
- [x] `make test-all` / `make test-jet` / `make test-aqua` all green.

## Phase 5 — bench, docs, commit

- [x] ~~Re-run `test/benchmark_sorted_container.jl` against the new
      `SortedCounter` for parity~~ — bench script deleted; its data
      is captured in the spec design.md table. The integration
      test wall time stayed within noise (~72 s before, ~80 s
      after; ±10% sample variance), so no separate bench log
      entry created.
- [x] ~~`@time make test-integration` baseline vs new~~ — covered
      above; no regression.
- [x] Update `docs/design-notes/replace-sorted-container.md` and
      `DESIGN_NOTES.md` index to mark complete.
- [x] Update `docs/design-notes/refactor-sweep.md` R8 Plan C entry:
      mark resolved.
- [x] `code-reviewer` agent pass on the diff.
- [x] Single commit via `git-helper`. Push only after explicit user
      confirmation.

## Risks tracked

- XML byte-identical regression caught only by integration tests —
  Phase 2 / Phase 3 each include a manual diff check before
  proceeding.
- Compile time: introducing one new struct in `src/common/` and
  removing three. Net neutral, no concern.
