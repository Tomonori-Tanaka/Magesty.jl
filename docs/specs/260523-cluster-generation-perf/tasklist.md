# Tasklist: Cluster-generation performance

Status: draft (2026-05-23)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Canonical-form helper

- [ ] Add `_translation_canonical_form(cluster, symmetry)` to
      `src/Clusters.jl`.
- [ ] Add a testset to `test/component/test_Clusters.jl` covering:
      sort-invariance, equivalent-orbit consistency, distinct-rep
      distinctness.
- [ ] `make test-unit` green.

### M2 — Dict-keyed irreducible_clusters

- [ ] Reimplement `irreducible_clusters` to use a per-body
      `Dict{Vector{Int}, Nothing}` (or `Set`) keyed by canonical form;
      preserve "first-seen-wins" by skipping clusters whose canonical
      form is already in the map.
- [ ] Leave `is_translationally_equiv_cluster` in place — existing tests
      reference it; add an inline comment noting it is now used only by
      tests.
- [ ] `make test-all` green; existing component tests unchanged.

### M3 — Remove `set_mindist_pairs` double computation

- [ ] Change `generate_clusters` signature to take `min_distance_pairs`
      as the last positional argument; remove the internal recompute.
- [ ] Compute `min_distance_pairs` once in the `Cluster` constructor,
      pass it to `generate_clusters`, and store it in the struct field
      as before.
- [ ] Update `bench/benchmark_cluster.jl` to construct
      `min_distance_pairs` once and pass it into the timed
      `generate_clusters` call. The standalone
      `set_mindist_pairs` stage timing stays.
- [ ] `make test-all` and `make bench-cluster` green.

### M4 — Benchmark and record

- [ ] Run `make bench-cluster` on both shipped 3-body fixtures.
- [ ] Append a before/after stage-table comparison to
      `.claude/bench_log.md`, citing this spec slug.
- [ ] Add a `Performance` entry to `CHANGELOG.md` `[Unreleased]` linking
      to the bench-log entry.

### M5 — Tier 2 review and merge

- [ ] Run the four-axis review panel
      (numerical / maintainability / performance / API) on the diff.
- [ ] Resolve blockers / major findings.
- [ ] Update `Status:` in this file and the table row in
      `docs/specs/README.md` together when complete.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] If results changed: regression or validation test added.
- [ ] ~~If public API changed: `SPEC.md` and `docs/src/api.md` updated.~~
      (no public-API change)
- [ ] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [ ] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ (no agent-facing change)
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.

## Implementation commits

<!-- Appended after each milestone commit, e.g.:
- `abcdef1` — perf(clusters): switch irreducible_clusters to canonical-form Dict lookup
-->
