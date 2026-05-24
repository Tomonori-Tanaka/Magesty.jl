# Design-matrix algorithmic restructuring

**Status**: not started (2026-05-24)

Investigation note. `build_design_matrix_energy` and especially
`build_design_matrix_torque` in `src/Fitting.jl` remain the dominant
cost of `fit`. Previous specs (`260516-coupled-basis-typeparam`,
`260516-optimize-workspace`, `260523-design-matrix-3body-perf`) have
exhausted micro-optimizations within the current loop structure
(type-parameterize on rank, pool scratch, drop `Set` / `SubArray`,
buffer Legendre). The remaining wins require **structural** changes —
reordering the loops and lifting work out of the inner kernel.

This note enumerates those structural ideas, with rough impact
estimates, so we can pick one (or stage them as a spec).

Current baseline (FeGe B20 2x2x2, 100 spinconfigs, 4 threads, commit
`4809a5a`, from `.claude/bench_log.md`):

| function | time | allocations |
|---|---:|---:|
| `build_design_matrix_energy` | 1.23 s | 7.02e6 |
| `build_design_matrix_torque` | 10.11 s | 3.36e8 |

The 10:1 torque:energy ratio is the most striking number and turns out
to be structural (see §1.3).

## 1. Structural problems with the current implementation

### 1.1 Zₗₘ is recomputed an order of magnitude more often than necessary

In `design_matrix_energy_element`, `Zₗₘ_unsafe(l, m, dir, legendre_buf)`
is called inside

```
for itrans in symnum_translation
    for site_idx = 1:(R - 1)
        for m_idx = 1:(2*l + 1)
            sh_values[...] = Zₗₘ_unsafe(l, m, dir, legendre_buf)
```

and the surrounding `build_design_matrix_energy` wraps it in
`(num_salcs, num_spinconfigs, cbc-in-group)`. Total calls scale as

```
num_spinconfigs × num_salcs × Σ_cbc × n_translations × N_sites × (2l+1).
```

But the set of values that can ever appear is

```
num_spinconfigs × num_atoms × (l_max + 1)²
```

— for the FeGe benchmark that is `100 × 64 × 9 ≈ 5.8e4` values total,
versus orders-of-magnitude more `Zₗₘ_unsafe` calls today. Caching
`Zₗₘ_unsafe` allocations (the workspace work, completed under spec
`260516-optimize-workspace`) reduced cost-per-call; **the call count
itself has not yet been touched**.

### 1.2 The `Mf` loop holds compile-time constants

Inside the inner contraction we do

```julia
for mf_idx = 1:Mf_size
    for other_tuple in ...
        for m_idx_last = ...
            mf_contribution += cbc.coeff_tensor[idx_buf..., mf_idx] * (...)
        end
    end
    tensor_result += cbc.coefficient[mf_idx] * mf_contribution
end
```

`cbc.coefficient` is fixed at SALC construction time. The sum is
linear, so

```
T̃[m₁,…,m_N] := Σ_{Mf} coeff_tensor[m₁,…,m_N, Mf] · coefficient[Mf]
```

can be precomputed once per `cbc`. The inner kernel then contracts a
single rank-`N` tensor with the SH product and drops the entire `Mf`
axis. For typical `Mf_size ∈ [3, 7]`, this is a direct contraction-cost
reduction of the same factor. Numerically: linear recombination, only
the floating-point reduction order differs.

### 1.3 Torque does `num_atoms` × redundant work, dedup overhead included

`build_design_matrix_torque` is structured as

```
for sc_idx = 1:num_spinconfigs
    for iatom = 1:num_atoms                       # 64
        for (salc_idx, key_group) in salc_list    # ~hundreds
            for cbc in key_group
                calc_∇ₑu!(...)                    # walks all translations
                    │
                    ├─ translate + sort + hash + searched_pairs membership
                    │   (runs on every translation, regardless of match)
                    └─ if atom_site_idx == 0; continue; end
```

`calc_∇ₑu!` loops over every translation. The dedup machinery
(`copyto!` + `sort!` + `_atoms_hash_key` + linear-scan `in`) runs on
**every** translation — the `atom_site_idx == 0` early-return fires
only *after* the dedup block (`src/Fitting.jl:854–865`). For a cluster
of `N = 2..4` sites and `num_atoms = 64`, only `N` of those
`num_atoms` iterations carry real SH/contraction work; the other ~60
still pay the full dedup cost before returning.

Spec `260523-design-matrix-3body-perf` measured that switching dedup
from `Set{UInt}` to `Vector{UInt}` alone shaved ~9% of torque wall
time. The underlying dedup work is not "modest single-digit-percent"
overhead — historical evidence puts it at roughly **20–25% of current
torque wall time** for this fixture. Two structural changes
contribute to fixing it independently:

- §2.C eliminates dedup from the hot path entirely (orbit
  enumeration moves to SALC-build time).
- §2.D restructures the loop to `(salc, cbc, cluster)`-major, so
  `(iatom, itrans)` is no longer a cross product and the
  `atom_site_idx == 0` iterations vanish.

The "10:1 torque:energy" ratio is therefore the sum of two
independent excess factors (dedup overhead, plus `num_atoms / N`
iteration multiplication), and §2.C / §2.D should be treated as
**co-equal wins** for the torque path, not as a precondition and a
follow-up.

### 1.4 Translation dedup is re-done per `(spinconfig, cbc)` call

`for itrans in symnum_translation` + `_atoms_hash_key` + the
`searched_pairs` linear-scan is, in effect, enumerating the distinct
clusters in the orbit. That set depends only on
`(cbc.atoms, symmetry.map_sym)` — **not** on spin directions and not
on `spinconfig`. Today we redo this enumeration on every call.

## 2. Restructuring proposals

Listed in order of (estimated impact / implementation cost). Several
compose; rough numbers assume FeGe-scale benchmarks.

### A. Fold `Mf` into the coefficient tensor at construction time

`CoupledBasis_with_coefficient` carries both `coeff_tensor` (rank `R`)
and `coefficient` (length `Mf_size`). Precompute a folded tensor

```
T̃[m₁,…,m_N] := Σ_{Mf} coeff_tensor[m₁,…,m_N, Mf] · coefficient[Mf]
```

of shape `(2l₁+1, …, 2l_N+1)` (rank `R-1`) once during basis build,
and read `T̃` instead of the `(coeff_tensor, coefficient)` pair in the
hot kernel.

**Field-layout decision** (must be made up front; both strategies are
viable):

- **(a) Add a third field `folded_tensor`, keep `coeff_tensor` and
  `coefficient`.** Memory roughly doubles per `cbc`. Preserves XML
  round-trip identity at the field level, all existing introspection
  call sites still work. Recommended default — safest for
  numerical-reproducibility audits.
- **(b) Replace `coeff_tensor` + `coefficient` with `folded_tensor`
  only.** Minimal memory. Requires schema-version bump in
  `src/XMLIO.jl` (new readers must reject old files; old readers must
  reject new files); requires updates to integration tests that
  introspect `cbc.coefficient` (`test/integration/fept_tetragonal_2x2x2/test.jl`
  field-equality round-trip; `test/integration/square_lattice/test.jl`
  `norm(salc[1].coefficient)` analysis metric); requires reviewing
  any callers in `tools/` that walk individual `cbc` fields.

The spec must pick one strategy and enumerate the touched test/I/O
surface.

- **Scope**: `CoupledBases.jl` constructor; `design_matrix_energy_element`
  and `calc_∇ₑu!` drop the `mf_idx` loop; `src/XMLIO.jl` save/load
  (depends on strategy); integration tests above (only if strategy b).
- **Expected impact**: contraction cost ×(1/`Mf_size`) on the inner
  loop. Per-`cbc` `Mf_size` is `2 Lf + 1`; for FeGe body-2 with
  `Lf ∈ {0, 1, 2}` it is 1, 3, or 5 (the "3–7×" range in earlier
  drafts assumed mixed body / higher `Lf`). End-to-end wall-time
  reduction depends on what fraction of `design_matrix_*` time is in
  the contraction vs SH evaluation; A alone is unlikely to dominate
  the overall speedup, but it's a near-free prerequisite for B and D.
- **Numerical risk**: if `T̃` is accumulated in the same `Mf` order
  the current code uses, A alone is **bit-identical** to the existing
  output (FMA-stability assumed). Roundoff tolerance enters only when
  the SH product / cluster sum order also changes (B and D).
- **Cost**: ~1 day. Cleanest first step.

### B. Precompute Zₗₘ / ∂ᵢZₗₘ tables per spinconfig

Hoist SH evaluation out of `design_matrix_energy_element` /
`calc_∇ₑu!` and into `build_design_matrix_*`.

Sketch:

```
for j in spinconfigs:
    # Once per spinconfig:
    #   Z[atom, lm_index]           — Float64    (energy + torque)
    #   ∂Z[atom, lm_index]          — SVector{3} (torque only)
    # sized num_atoms × Σ_l (2l+1) where l ranges over the
    # union of cbc.ls values present in salc_list.

    for salc_idx, cbc, cluster:
        # SH values become pure table lookups; inner kernel is just
        # the (now Mf-folded) tensor contraction.
```

**Physical invariant the cache relies on** — the symmetry application
in the design-matrix hot path is **pure translation**, not rotation.
Verified at `src/Symmetries.jl` (`symnum_translation` is built from
identity-rotation operations) and at `src/Fitting.jl:660-664`
(energy) / `:886-890` (grad) which read `spin_directions[:, translated_atom]`
without rotating the spin vector. So `Zₗₘ(spin[a])` depends only on
`(atom_index, l, m)` for a fixed spinconfig, and caching is exact.
Any future extension that admits non-trivial rotations in the
hot-path symmetry set (e.g. magnetic point groups with spin-orbit
coupling) would invalidate this cache and must be flagged.

**Cache layout, per function** (the two paths thread differently and
need different placements):

- `build_design_matrix_torque` threads over `sc_idx`. The cache is
  per-thread, built once at the top of each `sc_idx` iteration before
  the `iatom`/`salc_idx` loops. Memory: one cache per active thread.
- `build_design_matrix_energy` threads over `i = 1:num_salcs`, with a
  serial `for j = 1:num_spinconfigs` inside. The natural placement is
  *inside* the inner `j` loop, recomputed per spinconfig — meaning
  the cache is rebuilt `num_threads × num_spinconfigs` times instead
  of `num_spinconfigs` times. To get the single-build amortization
  the sketch above implies, swap energy's threading axis to be over
  spinconfigs (or build the per-spinconfig caches once up front into
  a shared read-only `Vector{Cache}` of length `num_spinconfigs`).
  The spec must pick.

**Cache size & footprint**: `num_atoms × Σ_l (2l+1)` Float64 entries
(plus 3× that for `∂Z` stored as `SVector{3, Float64}`). Worked
examples:

- FeGe 2×2×2, body-2, only `l = 1` per site:
  `64 × 3 = 192` entries · 8 B = 1.5 KB (energy);
  + 192 × 24 B = 4.5 KB (grad). Both fit easily in L1.
- Hypothetical 3-body with `l_max = 2`:
  `64 × (1+3+5) = 576` entries · 8 B = 4.5 KB; grad 14 KB. Still L1.

False sharing is not a concern at these sizes.

- **Scope**: per-`(l_max, atom)` SH layout decision (flat `Vector`
  with offsets vs `Matrix`); cache lifetime & threading-axis choice
  for the energy path; update `design_matrix_energy_element` /
  `calc_∇ₑu!` to accept and read from the cache instead of calling
  `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` directly.
- **Expected impact**: `Zₗₘ_unsafe` call count drops by orders of
  magnitude (cluster/translation/cbc iteration multiplicity is
  removed). Most of energy's 1.23 s today is SH evaluation +
  contraction; this targets the SH half directly.
- **Numerical risk**: none (we are caching identical values).
- **Cost**: ~2 days.
- **Synergy**: lets D (below) drop right in. (Note: D is *not*
  strictly dependent on B — D can call `Zₗₘ_unsafe` directly and
  still get the dominant loop-count reduction. B provides incremental
  gain on top.)

### C. Pre-enumerate clusters-in-orbit at SALC construction time

Promote the orbit enumeration ("walk all translations, sort, hash,
dedup") to a basis-construction step. Attach the result to each
`cbc` (or to its orbit container) as
`clusters::Vector{NTuple{N, Int}}`.

Then `design_matrix_energy_element` becomes

```julia
for atoms in cbc.clusters
    # SH lookup or eval, then contraction
end
```

- **Scope**: extend `SALCBasis` build to compute and persist the
  cluster list per `cbc` (or share across `cbc`s with the same
  `cbc.atoms`). Drop `searched_pairs` / `_atoms_hash_key` /
  `translated_atoms` sort buffers / `MVector` of translated atoms from
  the hot path entirely. XML I/O must serialize the cluster list (or
  recompute on load — see §4).
- **Expected impact**: removes dedup overhead that historical
  measurements (spec `260523-design-matrix-3body-perf` M2) attribute
  ~20–25% of torque wall time to. For the energy path, dedup is a
  smaller fraction but still measurable. C is a **co-equal torque
  win with D**, not a precondition. Either alone removes a different
  excess factor (see §1.3); together they target the full structural
  gap.
- **Numerical risk**: none (deterministic enumeration of the same
  set, with construction-time sort).
- **Cost**: ~2–3 days. Touches the SALC-build / XML I/O surface area.

### D. Compute torque via the energy Jacobian (one pass per cluster)

Reverse the current `(sc_idx, iatom, salc_idx, cbc, itrans)` loop
order to `(sc_idx, salc_idx, cbc, cluster, site_in_cluster → iatom)`.

For one cluster touching atoms `a₁, …, a_N`, we already evaluate the
SH product. With the folded tensor `T̃` and the SH vectors `u^{(k)}`
for each site `k`, the per-cluster, per-`cbc` energy contribution is

```
Φ_cbc(cluster) = multiplicity · Σ_{m₁,…,m_N} T̃[m₁,…,m_N] · ∏_k u^{(k)}_{m_k}
```

(the `_cluster_scaling((4π)^(N/2))` factor is applied once after
summing all clusters, exactly as today.) The gradient with respect to
spin direction at site `k`'s atom is

```
∂Φ_cbc / ∂u^{(k)}_{m_k}  = multiplicity · Σ_{m≠k indices}
                              T̃[…] · ∏_{j≠k} u^{(j)}_{m_j}
∂Φ_cbc / ∂spin[a_k]      = Σ_m (∂Φ_cbc/∂u^{(k)}_m) · ∂ᵢZ_{l_k,m}(spin[a_k])
```

i.e. ordinary reverse-mode AD applied to a tensor product. The `N`
partial-product slices needed for the `N` sites can be assembled at
the same time as the forward `Φ`.

**Invariant the derivation depends on**: within a single translated
cluster, the atom indices `a₁, …, a_N` are pairwise distinct, so each
`∂Φ_cbc / ∂spin[a_k]` is a single-site partial (not a sum of
site-partials over multiple matches). Today's enforcement chain:
`Clusters.jl` filters input clusters via `allunique(atom_list_all)`,
and `symnum_translation` operations are bijections of the supercell,
so distinctness propagates per translation. This invariant must be
called out in `requirements.md` (any future on-site / self-cluster
extension would need a separate code path that sums site-partials),
and asserted in the prototype as `@assert allunique(translated_atoms)`
inside a debug build. The current `calc_∇ₑu!` masks a violation
silently by keeping only the last matching `atom_site_idx`.

**`cross` placement (linearity argument)**: today the cross product
`cross(spin[iatom], grad_u)` happens once per
`(sc_idx, iatom, salc_idx)` after summing all `cbc` contributions
into `group_grad`. Under D, we accumulate `∂Φ/∂spin[a_k]` into a
per-`(sc_idx, atom, salc_idx)` slot across clusters, and then take
`cross(spin[a_k], ·)` once at the end. This is exact because
`cross(s, ·)` is linear in its second argument and `spin[a_k]` does
not change inside the spinconfig.

**Accumulation buffer & threading**: per-thread per-spinconfig
buffer of shape `(3, num_atoms, num_salcs)` Float64. For FeGe
(`num_atoms = 64`, `num_salcs ≈ 100`) one buffer is
`3 × 64 × 100 × 8 B ≈ 150 KB`; per-thread at 4 threads ≈ 600 KB
total — fits comfortably in L2. Threading axis stays as today
(parallel over `sc_idx`); the inner loop walks `(salc_idx, cbc,
cluster)` and writes into the local buffer; at the end of the
spinconfig the buffer is reduced into the global design matrix with a
single `cross` per `(atom, salc_idx)` cell. No cross-thread
synchronization required.

(Alternative: parallel over `(sc_idx, salc_idx)` blocks for finer
granularity on machines with many cores. Adds bookkeeping; revisit
only if `num_spinconfigs < num_threads`.)

- **Scope**: substantial rewrite of `calc_∇ₑu!` and the surrounding
  `build_design_matrix_torque` loop.
- **Expected impact**: torque drops from `O(num_atoms × n_clusters × N)`
  to `O(n_clusters × N)`, removing the `num_atoms / N` factor
  (≈ 16–32× for N=2..4 on a 64-atom cell). Realistic end-to-end gain
  ~5–15× after accumulation-buffer traffic, the deferred per-spinconfig
  reduction, and the still-present dispatch boxing (§2.F).
- **Numerical risk**: low (mathematically identical to the current
  formula, only the loop order changes). Float-roundoff order changes;
  regression tolerance needed (see §4).
- **Cost**: ~3–5 days. Biggest single win. Best done after A + C; B
  is incremental on top, not a strict prerequisite.

### E. SIMD-friendly contraction (`SMatrix` / Tullio)

After A and B, the inner kernel is a small tensor inner product with
known-at-compile-time dimensions (for fixed `R`). For 2-body, `l ≤ 4`,
this is `u' * T̃ * v` with `T̃::SMatrix{9, 9, Float64}`. For 3-body,
Tullio.jl or a hand-written `SArray{Tuple{d, d, d}}` contraction.

- **Scope**: limited to the contraction kernel.
- **Expected impact**: 1.5–3× on the contraction step alone (highly
  size- and platform-dependent).
- **Cost**: ~1–2 days, deferred until A + B settle the kernel shape.

### F. Concrete-type `salc_list` (eliminate dispatch boxing)

The deferred item from spec `260523`: `salc_list` carries the
rank-erased `Vector{Vector{CoupledBasis_with_coefficient}}`, so
`for cbc in key_group` incurs boxing (~110M of 336M torque
allocations). Independent from A–E.

After D the `(iatom, itrans)` cross product is gone, so the
`for cbc` loop fires at `O(num_atoms / N) ≈ 16–32×` lower frequency.
The remaining boxing allocations therefore drop to roughly `3–7M`,
which is no longer a leading allocation source. **F's priority falls
significantly once D lands**; consider whether it is still worth
pursuing on its own after benchmarking post-D.

## 3. Recommended sequencing

1. **A** (Mf folding) — small diff, clean, contraction-step
   speedup ×(1/Mf_size). Warm-up; doesn't touch the loop structure.
2. **C** (orbit pre-enumeration) — moves dedup out of the hot path.
   For torque this alone removes ~20–25% of wall time (per spec
   `260523` historical data). Co-equal torque win with D.
3. **B** (per-spinconfig Z table) — biggest energy-side win; near-
   free numerically.
4. **D** (torque as Jacobian) — biggest torque-side win. Removes the
   `num_atoms / N` excess factor. B is incremental on top, not a
   prerequisite — D can be developed in parallel once C lands.
5. **E** (SIMD contraction) — last-mile polish.
6. **F** (concrete-type `salc_list`) — reassess priority after D
   benchmarks land; absolute allocation contribution drops
   proportionally with D's loop-count reduction, so F may no longer
   be worth its structural cost.

Order-of-magnitude target after A + B + C + D (FeGe benchmark, very
rough):

| function | now | after | factor |
|---|---:|---:|---:|
| `build_design_matrix_energy` | 1.23 s | ~0.1–0.2 s | 5–10× |
| `build_design_matrix_torque` | 10.11 s | ~0.3–0.5 s | 20–30× |

These are speculative. Each step needs a benchmark gate before
proceeding to the next (see `.claude/bench_log.md`).

## 4. Open questions before promoting to a spec

1. **XML I/O surface (A and C)**: changes what `SCEBasis` persists.
   Two options per change:
   - **A**: choose field-layout strategy (a) keep `coeff_tensor` +
     `coefficient` and add `folded_tensor`, or (b) replace both with
     `folded_tensor` only. (b) requires an `XMLIO.jl` schema-version
     bump so old/new readers reject mismatched files (silent
     misinterpretation must be impossible — e.g. a folded tensor
     stored with `Mf_size = 1` would silently deserialize as a
     singleton `coefficient = [1.0]` and double-multiply on load).
   - **C**: serialize the per-`cbc` cluster list, or recompute on
     load. Recomputing is simpler but shifts SALC-build cost to
     deserialization. Decide based on typical `load` use cases.
2. **B cache layout & lifetime**: see §2.B "Cache layout, per
   function". The energy threading axis is the open decision: keep
   `@threads for i = 1:num_salcs` and rebuild per spinconfig inside
   the inner loop (cheap rebuild, no shared state), or switch to
   `@threads for j = 1:num_spinconfigs` (single build per spinconfig,
   but changes the threading granularity and may interact poorly with
   small `num_spinconfigs` fixtures). Benchmark both before commit.
3. **Folded-tensor sparsity (A)**: how sparse is `T̃` after folding?
   If many entries are zero by symmetry, we can skip them in the
   contraction. Profile before adding complexity.
4. **Validation harness — staged tolerances**: end-to-end fit on
   FePt + FeGe with pre/post comparison of the design matrices.
   Tolerance must be staged per proposal because each step relaxes
   the floating-point reduction order differently:

   | After step(s) | rtol | atol | notes |
   |---|---:|---:|---|
   | A only (same Mf order at build) | 1e-14 | 1e-15 | Bit-identical if FMA-stable |
   | A + B (SH cached, contraction order unchanged) | 1e-14 | 1e-15 | Still no order change |
   | A + B + C (orbit enumeration order) | 1e-13 | 1e-14 | C may iterate clusters in build-time order vs runtime-translation order |
   | A + B + C + D (loop reorder, energy) | 1e-12 | 1e-13 | Per-cluster contraction order changes |
   | A + B + C + D (torque) | 1e-10 | 1e-12 | Per-(atom, salc) accumulation across clusters before `cross` |

   `atol` floor is required so design-matrix entries that are zero by
   symmetry aren't falsely flagged. Existing dimer / FePt integration
   tests run as regression gates throughout.

## 5. References

- `src/Fitting.jl` — current implementation of
  `build_design_matrix_energy`, `design_matrix_energy_element`,
  `build_design_matrix_torque`, `calc_∇ₑu!`.
- `src/SALCBases.jl`, `src/CoupledBases.jl` —
  `CoupledBasis_with_coefficient{R}` and SALC construction.
- `src/TesseralHarmonics.jl` — `Zₗₘ_unsafe`, `∂ᵢZlm_unsafe`, buffered
  overloads.
- `docs/src/technical_notes.md` — design-feature formula.
- `.claude/bench_log.md` — historical benchmarks, current baseline.
- Completed prior work: `docs/specs/260516-coupled-basis-typeparam/`,
  `docs/specs/260516-optimize-workspace/`,
  `docs/specs/260523-design-matrix-3body-perf/`.
