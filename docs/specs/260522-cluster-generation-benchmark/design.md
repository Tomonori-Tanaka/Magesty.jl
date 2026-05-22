# Design: Cluster-generation benchmark

Status: complete (2026-05-22)

## Summary

`Cluster` construction (`src/Clusters.jl:117-167`) is a fixed four-stage
pipeline: `generate_clusters` -> `irreducible_clusters` -> `cluster_orbits`
-> `set_mindist_pairs`. To attribute the ~4000 s three-body cost to a stage,
the benchmark calls these stage functions directly (they live in module
`Clusters`; only `Cluster` is exported, but `Magesty.Clusters.<fn>` reaches
them) and times each with `@elapsed`, then runs a single `Profile`-sampled
`Cluster(...)` call for a flat/tree breakdown.

`@benchmark` with multiple samples is not used: a 4000 s call cannot be
sampled repeatedly. The benchmark does one warmup call (to pay JIT cost)
followed by one measured pass per stage.

Work starts with a **lightweight fixture** — the `fege_2x2x2` structure with
a three-body interaction but a *finite* cutoff, so the combinatorics stay
small and the script runs in seconds. The full `cutoff = -1` fixture is a
later addition; the script takes `--input`, so no code change is needed to
point it at the heavy case.

## Module layout

| Target | Change |
|---|---|
| `bench/fixtures/fege_2x2x2_3body_light/input.toml` | New (creates the `bench/fixtures/` directory). `fege_2x2x2` structure inline + three-body interaction with finite cutoff. |
| `bench/fixtures/fege_2x2x2_3body_fefe_open/input.toml` | New. Same as the light fixture except `[interaction.body3]` `cutoff."Fe-Fe" = -1`; isolates the cost of opening one species pair. |
| `bench/benchmark_cluster.jl` | New. Per-stage timing + `Profile` report. `--input` is repeatable, so several fixtures run in one invocation. |
| `Makefile` | New `bench-cluster` target (benchmarks both fixtures by default). |

No `src/` files change. No new `bench/Project.toml` dependencies are
required: `TOML` and `Profile` are stdlib, and the script uses `@elapsed`
rather than `@benchmark`.

## API

No public or internal API is added or changed. The benchmark consumes
existing entry points:

```julia
# Input parsing and skeleton construction (existing)
input  = TOML.parsefile(path)
system_spec, interaction, options = Magesty.parse_toml_inputs(input)
structure = Structure(system_spec; verbosity = false)
symmetry  = Symmetry(structure, options; verbosity = false)

# Stage functions, reached through the Clusters module (non-exported)
cluster_dict = Magesty.Clusters.generate_clusters(
    structure, symmetry, interaction.bodyn_cutoff, interaction.nbody)
irr = Magesty.Clusters.irreducible_clusters(cluster_dict, symmetry)
orb = Magesty.Clusters.cluster_orbits(irr, symmetry)
mdp = Magesty.Clusters.set_mindist_pairs(
    structure.supercell.num_atoms, structure.x_image_cart,
    structure.exist_image; tol = symmetry.tol)
```

`interaction.bodyn_cutoff` is an `OffsetArray{Float64, 3}` with body-axis
`2:nbody` (and species axes `1:nkd`). The benchmark must derive it from
`parse_toml_inputs` and pass it straight to `generate_clusters` — do not
reshape, convert, or rebuild it as a plain `Array`, or the offset axes will
mismatch.

Benchmark CLI (hand-rolled `parse_args`, matching `benchmark_salcbasis_hotspots.jl`):

- `--input <path>` — default: the lightweight fixture.
- `--profile-delay <sec>` — `Profile.init` delay; default small, raise for
  long runs to avoid buffer overflow.
- `--no-profile` — per-stage timing only, skip the `Profile` pass.

## Types and conventions

No impact on physics conventions, units, or numerical conventions. The
benchmark is read-only with respect to source behavior.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): no.
- [ ] SCE coefficient XML (`save` / `load`): no.
- [ ] `Fitting` <-> `SALCBasis`: no.
- [ ] `.claude/agents/` references: new Makefile target `bench-cluster`;
      sweep agent docs that enumerate bench targets if any.
- [ ] `SPEC.md` / `docs/src/api.md` updates: none (no API change).

## Test strategy

No unit/integration tests are added — this is benchmark tooling, not source
behavior. Verification is manual:

1. `make bench-cluster` completes in seconds on the lightweight fixture and
   prints the four-stage table, cluster counts, `Profile` flat/tree, and the
   dominant-stage conclusion line.
2. `--input test/integration/fege_2x2x2/input.toml` (existing two-body)
   runs without error.
3. Sum of per-stage `@elapsed` is consistent with the constructor's printed
   `Time Elapsed`.

## Risks and open items

- The lightweight fixture's finite cutoff value is chosen empirically so the
  run takes a few seconds; the exact value is tuned during implementation
  and recorded in the fixture file's header comment.
- `set_mindist_pairs` is computed twice in the real construction path (once
  in `generate_clusters`, once in the constructor). The benchmark reports
  both timings as-is; it does not deduplicate them. The finding is noted for
  the follow-up optimization spec.
- `Profile` buffer can overflow on very long runs; `--profile-delay`
  mitigates this, and the heavy fixture is out of scope here anyway.
