# Design: replace custom SortedContainer with plain + sort-once pattern

Status: draft (2026-05-14)

## Benchmark evidence (chose this direction)

Source: `test/benchmark_sorted_container.jl`.

| Workload | N | Custom (current) | DataStructures tree | Plain + sort! once |
|---|---:|---:|---:|---:|
| **[A]** `push!`+iterate `Int` (duplicates) | 100 | 1.08 µs | 2.46 µs (`SortedMultiDict`) | **0.38 µs (2.9× ↑)** |
| | 1 000 | 25.8 µs | 27.8 µs | **2.9 µs (8.9× ↑)** |
| | 10 000 | 1 180 µs | 820 µs (1.4× ↑) | **25 µs (47× ↑)** |
| **[B]** counted-unique, max_val=8 | 100 | 4.1 µs | 4.7 µs (`SortedDict`) | **3.0 µs (1.4× ↑)** |
| | 1 000 | 25 µs | 37 µs (0.7×) | **17 µs (1.5× ↑)** |
| | 10 000 | 219 µs | 433 µs (0.5×) | **156 µs (1.4× ↑)** |
| **[B']** counted-unique, max_val=64 | 1 000 | 114 µs | 80 µs (1.4× ↑) | **51 µs (2.3× ↑)** |
| | 10 000 | 3 920 µs | 2 460 µs (1.6× ↑) | **716 µs (5.5× ↑)** |
| **[C]** `in` lookup | 10 000 | 172 µs | 166 µs (1.0×) | n/a (different op) |

Key reading:

- `DataStructures.SortedMultiSet` **does not exist**; only
  `SortedSet` / `SortedDict` / `SortedMultiDict` are available. The
  earlier design note `docs/design-notes/replace-sorted-container.md`
  was authored on the false premise that it did.
- Tree-based containers help only at N ≥ 10 000 with mostly-unique
  keys, and even then by ≤ 1.6×. Below that they are slower than
  the current `Vector`-backed custom types.
- The clear winner is **plain `Vector` / `Dict` during construction
  and a single `sort!` at finalize/iteration**. Wins are 1.4× at
  small N and 47× at large N for the SortedVector case; 1.4×–5.5×
  for the SortedCountingUniqueVector case.

## Caller operation audit

`grep` over `src/` and `test/` (excluding `SortedContainer.jl` and
`test_SortedContainer.jl`):

### Operations actually used by external callers

`SortedVector{AtomCell}` (only in `Clusters.jl`):

- `SortedVector{AtomCell}()` empty constructor.
- `push!(sv, ac::AtomCell)` — sorted insertion (duplicates allowed).
- `push!(sv, SortedVector([ac]))` — pushing a vector wrapper.
- Iteration in sorted order: `for cluster::SortedVector{AtomCell} in interaction_clusters[body][prim_atom_sc]`.

`SortedCountingUniqueVector{T}` (Clusters.jl, SALCBases.jl,
xml_io.jl, where T is `Vector{Int}` or `Basis.CoupledBasis`):

- Empty constructor `SortedCountingUniqueVector{T}()`.
- `push!(scv, k)` — increment count or insert with count = 1.
- `push!(scv, k, n::Integer)` — insert/increment by `n`
  (SALCBases.jl L513 only).
- Iteration in sorted unique-key order.
- `length(scv)`, `isempty(scv)`, `scv[1]` (first element by index).
- **Direct field access `scv.counts[k]`** (very common):
  - `haskey(scv.counts, k)` (Clusters.jl L638).
  - `scv.counts[k] += n` (Clusters.jl L639).
  - `scv.counts[k] = n` (Clusters.jl L642).
  - `scv.counts[k]` lookup, assumed present (SALCBases.jl L112,
    L204, L362, L499, L1118).
  - `get(scv.counts, k, 1)` (SALCBases.jl L975).

### Operations defined on the custom types but never used externally

- `findfirst`, `findall`, `delete!`, `deleteat!`, `deleteall!`,
  `clear!`, `copy`, `==`, `isless`, `in` (membership), `append!`,
  `addcount!`, `getcount` (replaced everywhere by `.counts[k]`
  access).
- `SortedUniqueVector{T}` has **zero external callers** — it exists
  only as the backing storage of `SortedCountingUniqueVector.data`.

This narrows the migration surface considerably: we do not need to
reproduce every method, only the operations actually used.

## Replacement plan per type

### `SortedVector{AtomCell}` → `Vector{AtomCell}` + `sort!` once

`Clusters.jl` builds two structures of this type, both populated
during `generate_clusters` and then iterated downstream.

- Replace the field type with `Vector{AtomCell}`.
- Replace `SortedVector{AtomCell}()` with `AtomCell[]`.
- Replace `push!(sv, ac)` (sorted insertion) with plain `push!`,
  then `sort!(vec)` at the end of `generate_clusters` before
  returning the `interaction_clusters` / `interaction_cutoff_dict`
  structures.
- The `SortedVector([single_element])` construction (Clusters.jl
  L387) becomes `[single_element]` (a one-element Vector). Likewise
  `SortedVector(atom_combination)` becomes `sort(atom_combination)`
  for the comparable-to-baseline ordering.
- Downstream iteration code is unchanged (`for x in vec`).

Sort happens once per `(body, prim_atom_sc)` bucket at finalize.
Bucket sizes are bounded by neighbor-list cardinality, well within
the regime where the bench shows 9–47× speedup.

### `SortedCountingUniqueVector{T}` → `SortedCounter{T}` (new helper)

We introduce one private helper struct in
`src/common/SortedCounter.jl`. The struct preserves the external API
the callers already use (`.counts`, `length`, `isempty`, `getindex`,
iteration, `push!(scv, k)`, `push!(scv, k, n)`) but uses
`Dict{T, Int}` as the canonical store with a lazily-sorted view
cache.

Sketch:

```julia
mutable struct SortedCounter{T}
    counts::Dict{T, Int}        # canonical: insertion + count lookup
    _sorted_keys::Vector{T}     # lazy cache of sorted unique keys
    _dirty::Bool                # cache invalidation flag
end

SortedCounter{T}() where T = SortedCounter{T}(Dict{T,Int}(), T[], false)

function Base.push!(sc::SortedCounter{T}, k::T) where T
    if haskey(sc.counts, k)
        sc.counts[k] += 1
    else
        sc.counts[k] = 1
        sc._dirty = true
    end
    return sc
end

function Base.push!(sc::SortedCounter{T}, k::T, n::Integer) where T
    if haskey(sc.counts, k)
        sc.counts[k] += n
    else
        sc.counts[k] = n
        sc._dirty = true
    end
    return sc
end

function _sorted_keys(sc::SortedCounter{T})::Vector{T} where T
    if sc._dirty
        resize!(sc._sorted_keys, length(sc.counts))
        i = 0
        for k in keys(sc.counts)
            sc._sorted_keys[i += 1] = k
        end
        sort!(sc._sorted_keys)
        sc._dirty = false
    end
    return sc._sorted_keys
end

Base.length(sc::SortedCounter) = length(sc.counts)
Base.isempty(sc::SortedCounter) = isempty(sc.counts)
Base.getindex(sc::SortedCounter, i::Int) = _sorted_keys(sc)[i]
Base.iterate(sc::SortedCounter) = iterate(_sorted_keys(sc))
Base.iterate(sc::SortedCounter, state) = iterate(_sorted_keys(sc), state)
```

Design notes:

- **Not `<: AbstractVector{T}`** — current `SortedCountingUniqueVector`
  inherits from `AbstractVector{T}` (transitively via
  `AbstractSortedVector`) but no caller relies on the full
  `AbstractVector` interface (`size`, broadcast, slicing, etc.).
  Keeping it as a plain struct is simpler and avoids method
  ambiguities with `Base.push!(::AbstractVector, ...)`.
- **`counts` stays as a public-ish field**: callers already read
  `scv.counts[k]` directly (~11 sites in `SALCBases.jl` /
  `Clusters.jl`). Replicating that access is the lowest-churn
  migration. We document the field contract: `counts` is the
  source of truth, `_sorted_keys` is a derived cache.
- **Cache invalidation is monotonic**: setting `_dirty = true` on
  new insert is sufficient because we never delete from the
  Dict in caller code (audit confirmed `delete!`/`deleteat!`/`clear!`
  unused externally). If future code adds deletion paths, they must
  also flip `_dirty`.
- **`SortedCounter` is internal**: not exported from `Magesty`.
  Lives in `src/common/SortedCounter.jl` as a sibling module to
  the (deleted) `SortedContainer`. Callers `using ..SortedCounter`
  the same way they did `using ..SortedContainer` before.
- **Cost model**: `push!` O(1) amortized. First `iterate` /
  `getindex` after any push: O(N log N) for the cached sort. Until
  the next push, all iteration is O(N) walk over the cached Vector,
  and `.counts[k]` lookup stays O(1). This matches the build/read
  bimodal pattern observed in callers.

### `SortedUniqueVector{T}` → deleted

Zero external callers. The container exists solely as a backing
storage for `SortedCountingUniqueVector.data`, which itself is
being replaced. Drop the type alongside `SortedContainer.jl`.

## Module layout after migration

| File | Change |
|---|---|
| `src/common/SortedContainer.jl` | **deleted** |
| `src/common/SortedCounter.jl` | **new** (private helper module) |
| `src/Magesty.jl` | swap `include("common/SortedContainer.jl")` → `include("common/SortedCounter.jl")` |
| `src/Clusters.jl` | `using ..SortedContainer` removed; new `using ..SortedCounter: SortedCounter`. `SortedVector` types replaced with plain `Vector{...}` and `sort!` finalize. `irreducible_cluster_dict` field type → `Dict{Int, SortedCounter{Vector{Int}}}` |
| `src/SALCBases.jl` | `using ..SortedContainer` removed; new `using ..SortedCounter`. All `SortedCountingUniqueVector{Basis.CoupledBasis}` type annotations and field types → `SortedCounter{Basis.CoupledBasis}` |
| `src/utils/xml_io.jl` | Same field-type swap |
| `test/runtests.jl` | swap include and testset references |
| `test/component_test/test_SortedContainer.jl` | **deleted** |
| `test/component_test/test_SortedCounter.jl` | **new** (tests for the helper) |
| `test/component_test/test_Basis.jl` | drop unused `using .SortedContainer:` import or migrate as needed |

## Linked changes

- `SALCBasis.coupled_basislist` field type changes from
  `SortedCountingUniqueVector{Basis.CoupledBasis}` to
  `SortedCounter{Basis.CoupledBasis}`. The `SALCBasis` constructor
  signature is unchanged; callers that read the field directly will
  see a type change. Within Magesty no caller relies on
  `AbstractVector`-only methods of the old type, so this is safe.
- `Cluster.irreducible_cluster_dict` field type changes from
  `Dict{Int, SortedCountingUniqueVector{Vector{Int}}}` to
  `Dict{Int, SortedCounter{Vector{Int}}}`. Same reasoning.
- The R8 Plan C bug-fix backlog
  (`docs/design-notes/refactor-sweep.md`) is fully resolved by this
  refactor: the bugs (`clear!` / `deleteat!` undefined, `delete!`
  O(N), `copy` redundant sort) are in deleted code.

## Physics & numerics

No change. The migration touches container abstractions only; the
basis ordering convention (`Basis.CoupledBasis` `isless`), the
cluster ordering convention (`Vector{Int}` lexicographic), and the
multiplicity counting are all preserved. XML output is required to
be byte-identical (verified per integration example).

## Open design questions for sign-off

1. **Helper name**: `SortedCounter{T}` proposed. Alternatives:
   `OrderedCounter`, `SortedCountedKeys`, `BasisSetCounter` (too
   domain-specific). Decision needed.
2. **Module location**: `src/common/SortedCounter.jl` proposed
   (parallels old layout). Alternative: inline as a private struct
   inside `Clusters.jl` or `SALCBases.jl`, but it's used by both
   so a shared module is cleaner.
3. **Whether to also bench end-to-end** on `fege_2x2x2` before
   merge, or rely on micro-bench + integration test correctness.
   Recommended: take a one-shot `@time make test-integration`
   before/after and append to `.claude/bench_log.md`.
