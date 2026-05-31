# Design: src architecture and hot-path cleanup

Status: draft (2026-06-01)

## Summary

Three independent workstreams, ordered by risk. A1 (cache encapsulation) and
C (allocations) are mechanical and behavior-preserving. A2 (Fitting kernel
de-duplication) is the risky one: the energy and torque kernels share a
contraction skeleton, but they are documented hot paths and the shared helper
must stay `@inline`-able. A2 lands only if a `@btime` comparison shows no
regression; otherwise it is deferred (recorded as such) rather than forced.

## Module layout

| Target | Change |
|---|---|
| `src/CoupledBases.jl` | Add a public-internal accessor that returns the angular-momentum coupling results for `(ls, isotropy)`, building and caching internally. Keep the cache `const` but stop exposing its key structure / `:none` symbol. |
| `src/SALCBases.jl` | Replace the `CoupledBases._angular_momentum_cache[(ls_vec, :none, isotropy)]` lookup (≈ lines 387–406) with a call to the new accessor. C items: `projection_matrix_coupled_basis` in-place matrix write; `is_translationally_equivalent_coupled_basis` / `find_translation_atoms` allocation removal; `(false, true)` literal. |
| `src/Fitting.jl` | C1: `build_sh_cache_torque` reuses the `plm_raw` from one `_legendre_pair_unsafe!` call to fill both `Z` and `∂Z`, instead of rebuilding the Legendre recursion once for `Zₗₘ_unsafe` and again for `∂ᵢZlm_unsafe`. A2 (conditional): extract the shared folded-tensor contraction skeleton. |
| `src/TesseralHarmonics.jl` | Possibly add a low-level combined value+gradient entry point to support C1 (only if needed; prefer reusing `_legendre_pair_unsafe!`). |

## API

Internal API only (nothing exported, nothing in `docs/src/api.md`). Proposed
`CoupledBases` accessor:

```julia
# Returns the per-Lf coupled-tensor bases and their coupling paths for the
# given per-site angular momenta, building and caching on first use.
function angular_momentum_coupling_results(
    ls::AbstractVector{<:Integer},
    isotropy::Bool,
)::Tuple{Dict{Int, Vector{Array{Float64}}}, Dict{Int, Vector{Vector{Int}}}}
```

`SALCBases` consumes it as `bases_by_L, paths_by_L = CoupledBases.angular_momentum_coupling_results(ls_vec, isotropy)`,
with the same `normalize = :none` semantics the call site relies on today
(the accessor fixes `normalize = :none` internally, matching the existing key).

The `:none` fix is safe: the only place that reads the cache by key is
`SALCBases.jl:388`, and it always uses `:none`. Other `normalize` values reach
the cache exclusively through `tesseral_coupled_bases_from_tesseral_bases`
(which keeps its `normalize` keyword and continues to populate the same cache),
so no caller depends on retrieving a non-`:none` entry through this accessor.

## Types and conventions

No change to physics conventions, units, or normalization. The cache key
currently encodes `(ls_vec, normalize, isotropy)`; the accessor keeps that
internally. `tesseral_coupled_bases_from_tesseral_bases` continues to populate
the same cache, so a mixed call sequence stays consistent.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): only if C1 adds a
      combined entry point — it must produce bit-identical `Z` / `∂Z` to the
      current `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe`. Guarded by the SpheriCart
      agreement test.
- [ ] SCE coefficient XML (`save` / `load`): not touched; round-trip test must
      still pass.
- [ ] `Fitting` ↔ `SALCBasis`: SALC key-group order preserved (A1 does not
      reorder; A2 does not change which feature maps to which `Jφ`).
- N/A `.claude/agents/` references: none touched.
- N/A `SPEC.md` / `docs/src/api.md` updates: internal-only, no public API change.

## Test strategy

- Rely on the existing suite (`make test-all`) to guard bit-for-bit results:
  energy/torque design matrices, predictions, save/load round trip, SpheriCart
  agreement, Sunny round trip.
- Benchmarks: `make bench-salcbasis` (and `bench-sphericart` for C1) before and
  after; record medians in `.claude/bench_log.md`.
- A2 gate: micro-benchmark `design_matrix_energy_element` and the torque
  accumulator before/after; if the shared helper regresses either, revert A2.

## Risks and open items

- **A2 is the main risk.** A non-inlined helper can defeat the current
  hot-path inlining. Mitigation: mark the helper `@inline`, compare `@btime`
  medians and `@code_warntype`/`@allocated`, and defer A2 if not clearly neutral
  or better.
- C1 must not change floating-point summation order in a way that perturbs
  `Z` / `∂Z`; if a combined recursion changes any digit, keep the two-call form
  and instead only avoid the redundant Legendre rebuild via the existing
  `_legendre_pair_unsafe!`.
- A1 changes process-global cache access patterns but not its contents;
  thread-safety characteristics are unchanged (build-time, single-threaded
  SALC construction).
