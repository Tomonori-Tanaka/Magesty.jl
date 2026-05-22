# Design: Cluster-generation performance

Status: draft (2026-05-23)

## Summary

Replace `irreducible_clusters`'s O(N_clusters^2) linear scan with a
`Dict`-keyed lookup against the **lex-minimum translation image** of each
cluster (its canonical form under the pure-translation subgroup). Two
clusters are translationally equivalent iff they share a canonical form,
so the existing `is_translationally_equiv_cluster` predicate is replaced
by one `==` on the canonical form — O(1) lookup instead of O(N) scan.
Computing the canonical form for a single cluster is
O(N_translations · cluster_size), the same factor each `equiv?` call
already pays, so the total cost drops from
O(N_clusters^2 · N_translations · cluster_size) to
O(N_clusters · N_translations · cluster_size).

The representative atom-list stored in the `SortedCounter` is still the
first cluster encountered in the existing iteration order over
`(body, prim_atom_sc, cluster_dict key)`, exactly matching today's
"first-seen-wins" behavior — the canonical form is only used as a lookup
key, not as the stored representative. This keeps downstream XML I/O,
SALC ordering, and `Jφ` indexing bit-identical.

In the same change, the `Cluster` constructor computes
`min_distance_pairs` once and threads it into `generate_clusters` as an
explicit positional argument. The matrix is read-only in both call sites,
so passing it through removes a duplicated O(N_atoms^2 · 27) computation
without changing any output.

Alternatives considered:

- **Hashing the canonical form via `hash` + `IdDict`**: rejected. The
  canonical form is already a small `Vector{Int}` and a `Dict{Vector{Int},
  Vector{Int}}` lookup is fine; `IdDict` would not buy anything.
- **Storing the canonical form as the representative**: rejected. It
  would change the atom-list that lands in `irreducible_cluster_dict`,
  changing XML output and basis ordering. The brief explicitly forbids
  this in the Invariants section of `requirements.md`.
- **Bucket per `(body, prim_atom_sc)` rather than per `body`**: rejected.
  Today's `is_translationally_equiv_cluster` already checks across all
  primitive atoms (translations can map one prim atom to another), so the
  canonical-form `Dict` is correspondingly keyed at the `body` level
  only.

## Module layout

| Target | Change |
|---|---|
| `src/Clusters.jl` | Add internal `_translation_canonical_form(cluster, symmetry)::Vector{Int}`. Reimplement `irreducible_clusters` using a `Dict{Vector{Int}, Nothing}` keyed by canonical form. Keep `is_translationally_equiv_cluster` for the existing component tests but mark it as test-only via comment. |
| `src/Clusters.jl` | Change `generate_clusters` to take `min_distance_pairs` as the last positional argument; remove the internal recompute. Update the `Cluster` constructor to compute `min_distance_pairs` once before calling `generate_clusters`. |
| `bench/benchmark_cluster.jl` | Update the direct `generate_clusters` call to pass `min_distance_pairs`. The per-stage measurement of `set_mindist_pairs` stays — it still measures the matrix construction in isolation. |
| `test/component/test_Clusters.jl` | Add a new testset for `_translation_canonical_form`: (i) idempotent on a representative, (ii) equivalent clusters share canonical form, (iii) distinct irreducible representatives have distinct canonical forms. |
| `.claude/bench_log.md` | Append before/after entry from `make bench-cluster` on both shipped 3-body fixtures. |

## API

Internal helper (not exported):

```julia
"""
Return the lex-minimum atom-list reachable from `cluster` under any pure
translation in `symmetry.symnum_translation` (or its inverse). Two
clusters in the same translation orbit share this canonical form, so it
acts as a `Dict` key for grouping translation-equivalent clusters.

The result is always a sorted `Vector{Int}` of the same length as
`cluster`.
"""
function _translation_canonical_form(
    cluster::AbstractVector{<:Integer},
    symmetry::Symmetry,
)::Vector{Int}
```

Reimplemented internal:

```julia
function irreducible_clusters(
    cluster_dict::Dict{Int, Dict{Int, OrderedDict{Vector{Int}, Int}}},
    symmetry::Symmetry,
)::Dict{Int, SortedCounter{Vector{Int}}}
    # Per-body Dict{canonical_key => representative_atom_list} replaces the
    # quadratic linear scan against accepted representatives. Iteration
    # order over (prim_atom_sc, cluster_dict key) is preserved, so the
    # first-seen cluster in each translation orbit is the stored
    # representative — bit-identical to the previous behavior.
end
```

Changed positional signature:

```julia
function generate_clusters(
    structure::Structure,
    symmetry::Symmetry,
    cutoff_radii::AbstractArray{<:Real, 3},
    nbody::Integer,
    min_distance_pairs::AbstractMatrix{<:AbstractVector{<:DistInfo}},
)::Dict{Int, Dict{Int, OrderedDict{Vector{Int}, Int}}}
```

The `Cluster` constructor:

```julia
function Cluster(structure, symmetry, nbody, cutoff_radii; verbosity=true)
    min_distance_pairs = set_mindist_pairs(
        structure.supercell.num_atoms,
        structure.x_image_cart,
        structure.exist_image;
        tol = symmetry.tol,
    )
    cluster_dict = generate_clusters(
        structure, symmetry, cutoff_radii, nbody, min_distance_pairs,
    )
    irreducible_cluster_dict = irreducible_clusters(cluster_dict, symmetry)
    cluster_orbits_dict = cluster_orbits(irreducible_cluster_dict, symmetry)
    # ...same downstream as today...
end
```

Public-facing signatures (`Cluster(structure, symmetry, nbody,
cutoff_radii)` and `Cluster(structure, symmetry, interaction)`) and the
`Cluster` struct fields are unchanged.

## Types and conventions

The canonical form is defined as

```
canonical(c) = min { sort([map_sym[a, T]      for a in c]) : T in trans }
             ∪ min { sort([map_sym_inv[a, T]  for a in c]) : T in trans }
```

under lexicographic ordering on `Vector{Int}`. The inclusion of both
`map_sym` and `map_sym_inv` exactly mirrors today's
`is_translationally_equiv_cluster`. The pure-translation operations form
a group, so the inverse iteration is in principle redundant; we still
include it to keep the equivalence relation pointwise identical to the
existing test-pinned behavior.

No physics conventions, units, signs, or normalizations are touched.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): no change.
- [ ] SCE coefficient XML (`save` / `load`): no change. The first-seen
      cluster is preserved as the representative, so the XML key order is
      bit-identical. Verified by `make test-integration` (existing
      round-trip tests).
- [ ] `Fitting` <-> `SALCBasis`: no change. The outer index of
      `Vector{Vector{CoupledBasis_with_coefficient}}` is driven by the
      representative atom-list, which is unchanged.
- [ ] `.claude/agents/` references: no change. Bench tables already list
      `bench-cluster`. No new Makefile targets.
- [ ] `SPEC.md` / `docs/src/api.md` updates: no change (no public-API
      surface change).
- [ ] `CHANGELOG.md` `[Unreleased]`: add a `Performance` entry pointing
      at the speedup figure recorded in `.claude/bench_log.md`.

## Test strategy

- **Property tests** (new) on `_translation_canonical_form`:
  - `canonical(c) == canonical(sort(c))` (sort-invariance).
  - For every cluster reachable from `c` by translation, the canonical
    form matches `canonical(c)`.
  - Distinct irreducible representatives have distinct canonical forms.
- **Cross-check against current code** (new, scoped): build `Cluster`
  using the existing fixtures (dimer, chain, fept_tetragonal_2x2x2),
  compare `irreducible_cluster_dict` keys and counts against a baseline
  captured before the refactor. Run once in the spec branch, then keep
  the count assertions as permanent tests; the dict-equality assertion
  can stay if cheap.
- **Existing tests** in `test/component/test_Clusters.jl` (length / orbit
  assertions on dimer / chain / fept_tetragonal_2x2x2) must continue to
  pass unchanged.
- **Benchmark**: rerun `make bench-cluster` on both 3-body fixtures.
  Record before/after stage tables in `.claude/bench_log.md`.

## Risks and open items

- **Canonical-form definition mismatch.** If `map_sym` and `map_sym_inv`
  reference the same group element under a different convention than
  assumed, the canonical form could disagree with the existing equivalence
  relation. The property test "every translation of `c` has the same
  canonical form as `c`" mitigates this — it asserts the canonical form is
  in fact a representative of the same orbit.
- **Representative drift.** If the iteration order over the outer loops
  changes for any reason in the future (e.g., a Dict key ordering change),
  the chosen representative could drift even though the equivalence
  relation is unchanged. The "first-seen-wins" invariant is documented in
  `requirements.md`; if it ever becomes load-bearing for downstream
  reproducibility, we should make the choice deterministic by sorting on
  e.g. the lex-min atom-list per orbit. Out of scope for this spec.
- **`is_within_cutoff` allocation.** The benchmark still has `collect(
  combinations(...))` allocations inside the inner loop. Not addressed
  here; flagged as a follow-up in the Excludes section of
  `requirements.md`.
