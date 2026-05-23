# Requirements: Cluster-generation performance

Status: draft (2026-05-23)

## Goal

Cut `Cluster` construction time on three-body systems from minutes/hours
down to seconds by removing the two leading bottlenecks identified by the
companion benchmark spec (`260522-cluster-generation-benchmark`):

1. `irreducible_clusters` runs an O(N_clusters^2) linear scan against all
   previously accepted representatives, each comparison itself iterating
   over all pure-translation symmetry operations.
2. `set_mindist_pairs` is computed twice — once inside `generate_clusters`
   and once again in the `Cluster` constructor — for the same structure.

## Background

The benchmark established that on
`bench/fixtures/fege_2x2x2_3body_fefe_open/input.toml` (body-3 Fe-Fe pair
opened to all neighbors, all other body-3 cutoffs finite),
`irreducible_clusters` accounts for ~85% of construction time and ~2 GiB
of allocation: body-3 raw clusters jump 240 -> 1788 and the quadratic scan
dominates. The 4000 s run reported by the user opens *all* body-3 pairs;
the linear-scan / pure-translation product grows as
O(N_clusters^2 · N_translations · cluster_size).

`set_mindist_pairs` is a smaller win in absolute terms but is structurally
wasted work — the constructor already computes it for the struct field and
`generate_clusters` recomputes the same matrix internally. Both are
read-only consumers; passing one matrix through is straightforward.

## Scope

Includes:

- `src/Clusters.jl`: replace the O(N^2) linear scan in
  `irreducible_clusters` with a canonical-form lookup against a `Dict`
  keyed by the lex-minimum translation image of each cluster.
- `src/Clusters.jl`: compute `min_distance_pairs` once in the `Cluster`
  constructor and thread it into `generate_clusters` so it is not
  recomputed.
- New unit tests in `test/component/test_Clusters.jl` covering the
  canonical-form helper.
- Updated `bench/benchmark_cluster.jl` call site for the new
  `generate_clusters` signature, plus a before/after entry in
  `.claude/bench_log.md`.
- New benchmark fixture
  `bench/fixtures/fege_2x2x2_3body_all_open/input.toml` (all body-3
  cutoffs at -1). With the optimization in place this case runs in
  sub-second time and is cheap enough to ship as a regression-sized
  benchmark, validating the win on the original 4000 s pain point.

Excludes:

- `cluster_orbits` algorithmic changes (it is BFS-bounded by orbit size,
  not the N_clusters^2 scan; leave alone for now). Becomes the dominant
  stage after this spec lands; addressed in a separate follow-up if
  needed.
- `is_within_cutoff` allocation cleanup (`combinations` + `collect` inside
  the inner loop). Track as a follow-up if benchmarks show it surfaces.
- Threading any stage. Single-threaded wins should already drop runtime
  below the user's pain threshold; threading is a separate spec if needed.

## Invariants

Numerical / structural invariants that MUST hold across the change.

- The number of irreducible representatives per body and per primitive
  atom (`length(cluster.irreducible_cluster_dict[body])`) is unchanged on
  every existing fixture (dimer, chain, fept_tetragonal_2x2x2, fege_2x2x2,
  the two new 3-body fixtures).
- The orbit partition (`cluster_orbits_dict`) is unchanged: same number of
  orbits per body, same multiset of orbit sizes, same set of cluster
  members (orderings of `Vector{Vector{Int}}` may differ but the set of
  underlying atom lists must match).
- The representative atom-list chosen for each equivalence class is the
  first cluster encountered in the existing iteration order over
  `(body, prim_atom_sc, cluster_dict key)`. This preserves the current
  "first-seen-wins" semantics so that downstream XML I/O, SALC ordering,
  and `Jφ` indexing remain bit-identical.
- `min_distance_pairs` (the matrix stored on `Cluster`) is bit-identical
  before and after.
- Public API surface is unchanged: `Cluster(structure, symmetry, nbody,
  cutoff_radii)` and `Cluster(structure, symmetry, interaction)` keep the
  same signatures and the same struct fields.

## Completion criteria

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean.
- [ ] New unit tests for the canonical-form helper pass, including the
      "equivalent clusters share canonical form" and "distinct
      irreducible representatives have distinct canonical forms"
      properties.
- [ ] On the `fefe_open` benchmark fixture, `irreducible_clusters` stage
      time drops by at least 10x (target: from ~0.44 s to ~40 ms or less)
      and the dominant-stage line no longer points at
      `irreducible_clusters`.
- [ ] `make bench-cluster` overall wall-clock on `fefe_open` drops by at
      least 5x compared with the baseline recorded on
      `bench/cluster-generation-benchmark`.
- [ ] Before/after `@btime` (or `bench-cluster` stage table) recorded in
      `.claude/bench_log.md`.
- [ ] Integration test results on existing fixtures byte-identical
      (`make test-integration` passes with no diffs in regression-pinned
      outputs).

## References

- Companion spec:
  [`260522-cluster-generation-benchmark`](../260522-cluster-generation-benchmark/)
  (benchmark fixtures, `bench/benchmark_cluster.jl`, `bench-cluster`).
- Affected code: `src/Clusters.jl` (`irreducible_clusters`,
  `is_translationally_equiv_cluster`, `generate_clusters`,
  `set_mindist_pairs`, `Cluster` constructor).
- Existing tests: `test/component/test_Clusters.jl`.
