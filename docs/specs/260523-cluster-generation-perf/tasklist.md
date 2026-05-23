# Tasklist: Cluster-generation performance

Status: complete (2026-05-23)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Canonical-form helper

- [x] Add `_translation_canonical_form(cluster, symmetry)` to
      `src/Clusters.jl`.
- [x] Add a testset to `test/component/test_Clusters.jl` covering:
      sort-invariance, equivalent-orbit consistency, distinct-rep
      distinctness.
- [x] `make test-unit` green.

### M2 — Dict-keyed irreducible_clusters

- [x] Reimplement `irreducible_clusters` to use a per-body
      `Set{Vector{Int}}` keyed by canonical form; preserve
      "first-seen-wins" by skipping clusters whose canonical form is
      already in the map.
- [x] Leave `is_translationally_equiv_cluster` in place — existing tests
      reference it; inline comment notes its role as the test-only
      predicate form.
- [x] `make test-all` green; existing component tests unchanged.

### M3 — Remove `set_mindist_pairs` double computation

- [x] Change `generate_clusters` signature to take `min_distance_pairs`
      as the last positional argument; remove the internal recompute.
- [x] Compute `min_distance_pairs` once in the `Cluster` constructor,
      pass it to `generate_clusters`, and store it in the struct field
      as before.
- [x] Update `bench/benchmark_cluster.jl` to construct
      `min_distance_pairs` once and pass it into the timed
      `generate_clusters` call. The standalone
      `set_mindist_pairs` stage timing stays.
- [x] Ship `bench/fixtures/fege_2x2x2_3body_all_open/` (all body-3
      cutoffs = -1) as a third default benchmark fixture — sub-second
      after the optimization, the original 4000 s pain point.
- [x] `make test-all` and `make bench-cluster` green.

### M4 — Benchmark and record

- [x] Run `make bench-cluster` on all three shipped 3-body fixtures.
- [x] Append a before/after stage-table comparison to
      `.claude/bench_log.md`, citing this spec slug. Baseline measured
      at `bb1af45`, post-optimization at `da83dda`.
- [x] Add an `Internal`-section entry to `CHANGELOG.md` `[Unreleased]`
      summarizing the speedup.

### M5 — Tier 2 review and merge

- [x] Run the four-axis review panel
      (numerical / maintainability / performance / API) on the diff.
      Result: numerical 0 findings; remaining axes 7 majors + selected
      minors, all addressed in `71d57ff`.
- [x] Resolve blockers / major findings.
- [x] Update `Status:` in this file and the table row in
      `docs/specs/README.md` together when complete.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes.
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] If results changed: regression or validation test added.
      (results are bit-identical; the new
      `_translation_canonical_form` testset covers the helper directly,
      and the existing component / integration / save-load round-trip
      tests guard the bit-identity claim.)
- [ ] ~~If public API changed: `SPEC.md` and `docs/src/api.md` updated.~~
      (no public-API change)
- [x] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [ ] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ (no agent-facing change)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [x] Implementation commit hashes appended below.

## Implementation commits

- `3f33232` — docs(specs): draft cluster-generation-perf spec
- `724a292` — perf(clusters): add _translation_canonical_form helper (M1)
- `6cf10b5` — perf(clusters): replace irreducible_clusters O(N^2) scan
  with canonical-form lookup (M2)
- `8dcb80d` — Merge branch 'bench/cluster-generation-benchmark' into
  refactor/cluster-generation-perf
- `da83dda` — perf(clusters): thread min_distance_pairs through
  generate_clusters; ship all-open fixture (M3)
- `0f46cac` — docs(changelog): record cluster-generation perf wins (M4)
- `71d57ff` — refactor(clusters): apply review-panel findings on
  cluster-generation-perf (M5)
