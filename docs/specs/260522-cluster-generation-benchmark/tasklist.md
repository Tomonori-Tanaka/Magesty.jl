# Tasklist: Cluster-generation benchmark

Status: complete (2026-05-22)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Three-body fixtures

- [x] Create the `bench/fixtures/` directory and add
      `bench/fixtures/fege_2x2x2_3body_light/input.toml`: `fege_2x2x2`
      structure inline + three-body interaction with a finite cutoff tuned
      for a few-second run. Header comment records the cutoff choice.
- [x] Add `bench/fixtures/fege_2x2x2_3body_fefe_open/input.toml`: identical
      to the light fixture except `[interaction.body3]` `cutoff."Fe-Fe" = -1`,
      isolating the cost of opening a single species pair.

### M2 — Benchmark script

- [x] Add `bench/benchmark_cluster.jl`: argument parsing, input/skeleton
      construction, per-stage `@timed` timing with cluster-count metrics,
      `Profile` flat/tree report, dominant-stage conclusion line.
- [x] Add the `bench-cluster` Makefile target.

### M3 — Verify and report

- [x] Run `make bench-cluster`; confirms seconds-scale completion and the
      full output (stage table / counts / profile / conclusion).
- [x] Run with `--input test/integration/fege_2x2x2/input.toml`; runs with
      no error.
- [x] Cross-check: compiled full construction (0.0614 s) is within 2x of
      the per-stage total (0.0484 s), confirming the stage functions
      reproduce the constructor path.
- [x] Done: the script prints a dominant-stage line per fixture. On the
      light fixture `set_mindist_pairs` dominates (~42 %). On the
      `fefe_open` fixture, opening just the body-3 Fe-Fe pair grows body-3
      raw clusters 240 -> 1788 and makes `irreducible_clusters` dominant
      (~85 %, ~2 GiB allocated) — its O(N^2) linear scan over clusters is
      the prime suspect for the ~4000 s all-neighbor case. This is the
      input to the follow-up optimization spec.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] ~~`make test-all` passes.~~ (no source change)
- [ ] ~~`make test-aqua` / `make test-jet` clean.~~ (no source change)
- [ ] ~~If results changed: regression or validation test added.~~ (read-only)
- [ ] ~~If public API changed: `SPEC.md` and `docs/src/api.md` updated.~~
- [ ] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ (no source change; new tooling)
- [ ] ~~Tier 2 review panel run.~~ (benchmark tooling, no source change)
- [x] `.claude/agents/profiler.md` swept: `bench-cluster` added to both
      benchmark tables.
- [ ] ~~`CHANGELOG.md` `[Unreleased]` updated.~~ (bench tooling is not
      tracked in the changelog, consistent with other `bench-*` targets)
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [x] Implementation commit hash appended below.

## Implementation commits

- `0538918` — perf(bench): add cluster-generation benchmark with
  stage-level timing (two three-body fixtures, `bench/benchmark_cluster.jl`,
  `bench-cluster` Makefile target).
