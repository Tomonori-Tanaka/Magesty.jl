# SpheriCart.jl adoption assessment

**Status**: complete (2026-05-11)

**Conclusion**: do not adopt SpheriCart in Magesty.jl (opposite of the
SpinClusterMC decision).

## Background

A SpinClusterMC-side investigation (same day) confirmed that Magesty's
`Zₗₘ_unsafe` is bit-exact equivalent to SpheriCart's `Y_l^m_real`.
We checked whether Magesty itself would benefit from adoption.

## Functional check

- `make test-sphericart` -> `SpheriCart vs TesseralHarmonics agreement`
  passes 700/700 (`test/component/test_sphericart_agreement.jl`,
  lmax=4, all `(l, m)`, 7 unit vectors).
- Numerical identity is complete (same as the SpinClusterMC
  investigation).

## Measured performance (`bench-sphericart` family)

`test/benchmark_optimize_sphericart.jl` implements both v1 (drop-in:
per-atom `compute`) and v2 (per-spinconfig precomputed cache). Median
comparison:

### fege_2x2x2 (64 atoms, lmax = 1, 146 SALCs, 100 spinconfigs)

| Function | baseline (Magesty) | v1 (SpheriCart drop-in) | v2 (SpheriCart precomputed) |
|---|---|---|---|
| `design_matrix_energy_element` | **54.0 μs** | 133 μs (2.46x slower) | 113 μs (2.09x slower) |
| `calc_∇ₑu` | **8.4 μs** | 23.5 μs (2.80x slower) | 21.5 μs (2.55x slower) |
| `build_design_matrix_torque` (1 config) | **208 ms** | — | 503 ms (2.42x slower) |

### fept_tetragonal_2x2x2 (16 atoms, lmax = 2, 31 SALCs, 30 spinconfigs)

| Function | baseline | v1 | v2 |
|---|---|---|---|
| `design_matrix_energy_element` | **27.7 μs** | 2.40x slower | 2.07x slower |
| `calc_∇ₑu` | **7.7 μs** | 3.11x slower | 2.71x slower |
| `build_design_matrix_torque` | **6.56 ms** | — | 2.39x slower |

Memory in v2 also regresses: torque matrix at **785 MiB** vs the
baseline **293 MiB** (2.7x).

## Why the conclusion is opposite to SpinClusterMC

Access pattern difference:

- **SpinClusterMC**: per-atom cache during a sweep wants "all `(l, m)`
  up to `L_max` for one atom". SpheriCart returns an SVector batch that
  matches exactly -> **10–20x speedup**.
- **Magesty `design_matrix_energy_element` / `calc_∇ₑu`**: the loop is
  `for (site_idx, atom) in translated_atoms; l = cbc.ls[site_idx]; for m_idx in 1:2l+1`,
  i.e., **a fixed `l` per site**, only `2l+1` `m` values needed.
  SpheriCart still computes every `(l, m)` up to `L_max`, wasting the
  rest. Even at `L_max = 2`, of `(L+1)² = 9` values only ~1–3 are
  needed, so the overhead of computing the unused `(l', m')` outweighs
  any SIMD parallelism gain.

v2 (precomputed cache) does not help because:
- `Y_cache, dY_cache` per spinconfig is cheap, but each `itrans × site`
  iteration builds
  `sh_values = [[Y_cache[...] for mi in 1:2l+1] for si in 1:N]`,
  allocating an inner `Vector` on the heap every time. That allocation
  cost beats the value-computation cost of the built-in `Zₗₘ_unsafe`.

## Conditions where SpheriCart could win in Magesty (future reference)

In principle, SpheriCart could win when:
1. All SALC `cbc.ls` values cluster near the maximum (i.e., we
   actually use every `(l, m)` per call).
2. The inner `Vector{Float64}` allocations are removed entirely
   (`MMatrix` / linear indexing).
3. Per-spinconfig caches stay as `Matrix{Float64}` and are read
   row-by-row.

Today, Magesty SALCs mix `ls` across clusters (anywhere from 1 to
lmax), and high-lmax SALCs are a small fraction of the total. The
restructure cost is unlikely to be repaid.

## Future re-evaluation triggers

Re-investigate when these conditions all hold:
- 50%+ of SALCs have `L_max >= 4`.
- The hot path changes from `design_matrix_*` to a Monte-Carlo-style
  "per-atom, all `(l, m)`" access pattern.
- We can structurally drop the inner `sh_values` allocation (e.g.,
  switch to `MMatrix{N, 2L+1}`).

Until then, Magesty keeps **`TesseralHarmonics`**. SpheriCart stays
adopted only on the SpinClusterMC side; Magesty does not need the
dependency.

## Artifacts retained

- `test/component/test_sphericart_agreement.jl` stays as a drift
  detector (catches divergence whenever SpheriCart upgrades).
- `bench/benchmark_sphericart.jl` stays as a baseline for re-evaluation
  if the access pattern changes later.
- `Project.toml [compat] SpheriCart = "0.2.3"` is loaded only in tests
  via `[extras]`, so it can stay.
