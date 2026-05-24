# Design: Design-matrix algorithmic restructuring

Status: draft (2026-05-24)

## Summary

Four staged changes to `src/Fitting.jl` and adjacent SALC / I/O code,
applied in order A → C → B → D, each gated by a FeGe 2x2x2 benchmark.
The full algorithmic analysis lives in the design note
[`docs/design-notes/design-matrix-algorithmic-restructuring.md`](../../design-notes/design-matrix-algorithmic-restructuring.md);
this document records the decisions specific to this spec.

The progression is intentional: A is bit-identical and validates the
benchmark / regression-test harness; C removes hot-path dedup
overhead (~20–25% of torque per historical measurement) and
simplifies the kernel for B and D; B caches Zₗₘ values per
spinconfig, removing orders of magnitude of redundant SH calls; D
rewrites the torque kernel as cluster-major reverse-mode AD,
eliminating the `num_atoms / N` excess factor.

Stages B and D can in principle be developed in parallel after C
lands (D does not strictly require B); we keep the linear order for
benchmark attribution.

## Module layout

| Target | Change |
|---|---|
| `src/CoupledBases.jl` | A: add `folded_tensor` field (strategy a); B: no change; C: add `clusters::Vector{Vector{Int}}` field per `cbc`; D: no change |
| `src/SALCBases.jl` | C: compute `clusters` during basis build (translation enumeration + canonical sort + dedup, run once at construction) |
| `src/Fitting.jl` | A: drop `mf_idx` loop in `design_matrix_energy_element` / `calc_∇ₑu!`; C: replace `for itrans ...` with `for atoms in cbc.clusters`; B: add SH cache build at start of each spinconfig iteration in `build_design_matrix_*`, pass cache through; D: rewrite `build_design_matrix_torque` + replace `calc_∇ₑu!` with cluster-major Jacobian kernel |
| `src/XMLIO.jl` | A: serialize `folded_tensor` alongside `coeff_tensor` / `coefficient` (strategy a); C: do **not** serialize the `clusters` list — recompute lazily on load from `cbc.atoms` + `symmetry.map_sym` |
| `test/component/test_fitting.jl` | Add per-stage pre/post design-matrix equivalence tests with the staged tolerances |
| `test/regression_fixtures/` | New directory holding plain-text reference design matrices (`writedlm` format, no extra dependency) captured at the start of each stage |
| `bench/bench_b1_3body_fege.jl` | Reuse as the bench fixture for every stage; record before/after in `.claude/bench_log.md` |

## API

No public-API signature changes are planned. All changes are
internal to the SALC build / fitting hot path.

The `CoupledBasis_with_coefficient` field set changes (additive
under strategy a):

```julia
struct CoupledBasis_with_coefficient{R}
    ls::Vector{Int}
    Lf::Int
    Lseq::Vector{Int}
    atoms::Vector{Int}
    coeff_tensor::Array{Float64, R}        # kept (strategy a)
    coefficient::Vector{Float64}           # kept (strategy a)
    multiplicity::Int
    folded_tensor::Array{Float64, R - 1}   # NEW (A)
    clusters::Vector{Vector{Int}}          # NEW (C); each inner Vector has length R - 1
end
```

The `clusters` field deliberately uses `Vector{Vector{Int}}` rather
than a parametric `Vector{NTuple{N, Int}}` (which would force a
second type parameter `{R, N}` on the struct). Spec
`260516-coupled-basis-atoms-svector` previously showed that adding
an `N` type parameter to a field stored in `Vector{CoupledBasis...}`
containers regresses performance due to UnionAll erasure on
iteration. `Vector{Vector{Int}}` keeps the element type concrete and
matches the existing `cbc.atoms::Vector{Int}` precedent; the
per-cluster length-`R-1` invariant is enforced in the inner
constructor rather than in the type.

Strategy choice for A: **(a) additive — keep `coeff_tensor` and
`coefficient`, add `folded_tensor`.** Doubles per-`cbc` tensor
memory; preserves XML round-trip identity at the field level and
keeps existing introspection call sites
(`test/integration/fept_tetragonal_2x2x2/test.jl` field-equality
round-trip; `test/integration/square_lattice/test.jl`
`norm(salc[1].coefficient)`). Strategy (b) (replace) is rejected for
this spec to minimize the test / I/O surface; reconsider if
per-`cbc` memory becomes an issue.

Hot-path function signatures gain a cache parameter under B:

```julia
function design_matrix_energy_element(
    cbc, spin_directions, symmetry, ws,
    sh_cache,                                 # NEW (B)
)::Float64
```

`sh_cache` shape and lifetime are documented in the function
docstring per CLAUDE.md exported-API rules (the cache type itself is
not exported).

## Types and conventions

### Per-stage floating-point tolerances

Each stage relaxes the floating-point reduction order differently;
regression tests use staged tolerances:

| After | rtol | atol | reason |
|---|---:|---:|---|
| A | 1e-14 | 1e-15 | Bit-identical if `Mf` accumulation order preserved (FMA-stable) |
| A + C | 1e-13 | 1e-14 | C may iterate clusters in build-time canonical order vs runtime translation order |
| A + C + B | 1e-14 | 1e-15 | B is exact value caching; no order change |
| A + C + B + D (energy block) | 1e-12 | 1e-13 | D doesn't touch energy, but full-pipeline harness picks up combined drift |
| A + C + B + D (torque block) | 1e-10 | 1e-12 | Per-(atom, salc) accumulation across clusters before `cross` |

The `atol` floor is required so design-matrix entries that are zero
by symmetry are not falsely flagged.

### Invariants introduced

- **Cluster-atom distinctness (D)**: every entry in
  `cbc.clusters[k]` is pairwise distinct. Today's enforcement chain
  (`Clusters.jl` `allunique(atom_list_all)` filter,
  `symnum_translation` bijection) makes this hold by construction.
  D's gradient derivation depends on it. Asserted in a debug build
  as `@assert allunique(translated_atoms)`.
- **Pure-translation hot-path symmetry (B)**: see requirements.md.
  Any future change admitting non-trivial rotations in the
  design-matrix hot-path symmetry set invalidates the SH cache.

## Impact on linked sites

- [x] **Spherical-harmonics convention** (`TesseralHarmonics`): no
      change. B caches the existing `Zₗₘ_unsafe` /
      `∂ᵢZlm_unsafe` values; values themselves are unchanged.
- [x] **SCE coefficient XML** (`save` / `load`): A adds
      `folded_tensor` to the on-disk record; C adds the `clusters`
      list (or documents recompute-on-load). Round-trip identity is
      preserved (strategy a). Save format is additive only —
      readers built before this spec must continue to load files
      written before this spec; readers built after must accept old
      files and recompute the new fields lazily on load. The
      version field is bumped to flag the new fields, but
      backward compatibility is preserved.
- [x] **`Fitting` ↔ `SALCBasis`**: SALCBasis gains responsibility
      for computing `folded_tensor` (A) and `clusters` (C). Hot
      path reads only the precomputed fields; the outer-index
      ordering of `salc_list` is unchanged.
- [x] **`.claude/agents/` references**: N/A — none of A–D rename
      modules or Makefile targets; agent files unchanged.
- [x] **`SPEC.md` / `docs/src/api.md` updates**: N/A —
      `CoupledBasis_with_coefficient` is internal (not part of the
      public API listed in `docs/src/api.md`), and no exported
      signatures change.

## Test strategy

### Per-stage equivalence harness

A new test in `test/component/test_fitting.jl`:

```julia
@testset "design-matrix equivalence — stage X" begin
    # Build baseline matrices from a captured pre-stage state
    # (loaded from a JLD2 / on-disk reference file).
    # Build new matrices from the current implementation.
    # Compare with the staged (rtol, atol) from design.md table.
end
```

Reference matrices are captured once at the start of each stage
(via `tools/personal/` capture scripts) and stored in
`test/regression_fixtures/` as plain text via `writedlm` (no extra
dependency required; `DelimitedFiles` is in stdlib). Fixtures are
checked into the repo (small — FeGe 2x2x2 design matrices are
~hundreds of KB) and regenerated only when the prior stage's gate
is met. Each fixture file carries a header comment with the
generating commit hash so regression diffs are traceable.

End-to-end integration tests (`test/integration/dimer/`,
`test/integration/fept_tetragonal_2x2x2/`,
`test/integration/fege_2x2x2/`) run unchanged as the behavioral
regression gate; their built-in numerical tolerances are loose
enough to accept the staged drift.

### Benchmark gates

Each stage records a before/after entry in `.claude/bench_log.md`
using the FeGe B20 2x2x2 fixture (`bench/bench_b1_3body_fege.jl`),
5-trial median via `@timed`, `JULIA_NUM_THREADS=4`. Per the
completion criteria in requirements.md:

| Stage | Required speedup gate |
|---|---|
| A | (no minimum; clean-up step) |
| C | torque ≥ 1.15× over baseline |
| B | energy ≥ 2× over post-C |
| D | torque ≥ 5× over post-B |

A stage that fails its gate must be analyzed (or reverted) before
the next stage starts.

### Tier 2 review

After D lands, run the four-axis Tier 2 review panel (numerical /
maintainability / performance / API) per CLAUDE.md and resolve
findings before merge.

## Risks and open items

- **Strategy (a) memory overhead in A**: per-`cbc` tensor memory
  roughly doubles. For typical SCE basis sizes this is negligible
  (megabytes), but if the basis grows to gigabyte scale we may need
  to revisit strategy (b). Measure peak `SALCBasis` memory after A
  lands.
- **Energy threading axis in B**: design note §2.B / §4 question 2.
  Two layouts:
  - Keep `@threads for i = 1:num_salcs` and rebuild the cache per
    spinconfig inside the inner loop (cheap rebuild, no shared
    state).
  - Switch to `@threads for j = 1:num_spinconfigs` with one cache
    per spinconfig built up front.
  Benchmark both before committing to one. Default plan: keep the
  current axis, rebuild per spinconfig (lower-risk change).
- **C's serialization choice** (resolved): the `clusters` list is
  **not** serialized to XML; it is recomputed lazily on load from
  `cbc.atoms` + `symmetry.map_sym`. Keeps the file format
  unchanged for C. Revisit only if load time becomes a measurable
  bottleneck on large basis sets.
- **D's accumulation buffer**: per-thread, per-spinconfig buffer of
  shape `(3, num_atoms, num_salcs)`. For FeGe-class sizes this is
  ~150 KB per thread; large fixtures (`num_atoms ~ 256`,
  `num_salcs ~ 500`) push toward 3 MB per thread, still well within
  L2/L3. Flag if a future fixture pushes beyond this.
- **Sparsity of folded `T̃` (A)**: how sparse is the folded tensor?
  If many entries are zero by symmetry, contraction can skip them.
  Profile after A lands; do not add complexity speculatively.
- **Dispatch boxing (F deferred)**: post-D, F's contribution drops
  by `O(num_atoms / N)`. Decide whether to schedule F as a
  follow-up spec based on post-D allocation counts.
