# Requirements: Cluster-generation benchmark

Status: complete (2026-05-22)

## Goal

Build a `bench/` benchmark that measures `Cluster` construction stage by
stage, so the bottleneck in three-body (`nbody = 3`) cluster generation can
be identified with numbers rather than guesses.

## Background

Extending the `fege_2x2x2` integration system to three-body interactions
(`nbody = 3`, every pair `cutoff = -1`, i.e. all neighbors) made `Cluster`
construction take roughly 4000 seconds. The constructor runs
`generate_clusters` -> `irreducible_clusters` -> `cluster_orbits` ->
`set_mindist_pairs` in sequence (`src/Clusters.jl`), and it is currently
unknown which stage dominates.

This spec covers the **measurement infrastructure only**. Algorithm
optimization is deferred to a follow-up spec, decided after the numbers are
in. Per the requester, work starts with a **lightweight fixture** (finite
cutoff, runs in seconds) so the benchmark script can be developed and
iterated quickly; the full `cutoff = -1` fixture is added later.

## Scope

Includes:

- New three-body benchmark fixtures under `bench/fixtures/` (`input.toml`
  only; structure is defined inline): a lightweight all-finite-cutoff
  fixture, and a `fefe_open` fixture with only the body-3 Fe-Fe cutoff
  opened to `-1`.
- New benchmark script `bench/benchmark_cluster.jl`: per-stage wall-time
  measurement plus a `Profile` flat/tree report; `--input` is repeatable.
- New `bench-cluster` Makefile target.

Excludes:

- The full `cutoff = -1` four-thousand-second fixture (added later, once the
  script is stable).
- Any change to `src/Clusters.jl` or other source modules — including the
  known double computation of `set_mindist_pairs` (once inside
  `generate_clusters`, once in the `Cluster` constructor). Optimization is a
  separate spec.
- Threading or algorithmic changes to cluster generation.

## Invariants

- No source module changes: `src/Clusters.jl` and all of `src/` are
  untouched. Numerical results are unaffected by definition.
- Physics conventions, public API, and XML format are not touched.
- The benchmark only *reads* via existing constructors and the (non-exported
  but reachable) `Clusters` stage functions; it introduces no new public API.

## Completion criteria

- [x] `make bench-cluster` runs on the lightweight fixture and completes in
      a few seconds.
- [x] Output includes a per-stage wall-time table (`generate_clusters` /
      `irreducible_clusters` / `cluster_orbits` / `set_mindist_pairs`),
      cluster-count metrics, a `Profile` flat/tree report, and a one-line
      "dominant stage" conclusion.
- [x] The script also runs on the existing two-body
      `test/integration/fege_2x2x2/input.toml` via `--input` without error.
- [x] On the lightweight fixture, a compiled full `Cluster` construction
      (0.0614 s) is within 2x of the per-stage total (0.0484 s), confirming
      the stage functions reproduce the constructor path.

## References

- Related issues / PRs: none.
- Related specs / design notes: follow-up optimization spec to be created
  after this benchmark produces results.
