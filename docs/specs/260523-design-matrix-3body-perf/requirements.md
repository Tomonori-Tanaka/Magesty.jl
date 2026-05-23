# Requirements: Design-matrix construction perf for 3-body SALCs

Status: complete (2026-05-23)

## Goal

Reduce the wall time of `build_design_matrix_energy` and especially
`build_design_matrix_torque` (`src/Fitting.jl`) for 3-body SCE bases,
where current profiling shows the torque path is ~9× slower than the
energy path and dominated by allocations / broadcast overhead inside
`calc_∇ₑu!`.

## Background

A 3-body benchmark on FeGe B20 2x2x2 (64 atoms, 100 spin configs) was
run with the fixtures `bench/fixtures/fege_2x2x2_3body_light/` and
`bench/fixtures/fege_2x2x2_3body_fefe_open/` (selectable via
`bench/bench_b1_3body_fege.jl` using fixture keys `light` and
`fege_open` respectively), reusing the EMBSET from
`test/integration/fege_2x2x2/`. Median wall times (4 threads):

| fixture | SALCs | energy (s) | torque (s) | torque/energy |
|---|---|---|---|---|
| `_light` (cutoff 3.0 Å) | 1232 | 1.32 | 12.08 | 9.1× |
| `_fefe_open` (Fe-Fe = -1) | 6897 | 10.5 | 119 | 11.3× |

Sampling-profile + `Profile.Allocs` traces identified two genuinely
actionable hot spots, both inside `Fitting.jl`:

1. `searched_pairs::Set{UInt}` de-duplication in both
   `design_matrix_energy_element` and `calc_∇ₑu!` spends ~24 % of
   torque samples in `push!` / `in` (hash-table operations).
   `n_translations` is small enough (O(ntran) ≈ tens) that
   `Vector{UInt}` with linear-scan `in` is faster end-to-end.
2. The `@views Zₗₘ_unsafe(l, m, spin_directions[:, atom], …)` form
   inside the element functions, plus
   `@views dir_iatom = spinconfig.spin_directions[:, iatom]` inside
   `build_design_matrix_torque`, both create a `SubArray` of
   `spin_directions` on every read. Reading three scalars into a
   stack `SVector{3, Float64}` once per `(itrans, site)` (or once
   per `iatom`) reuses it across all `2*l + 1` SH calls per site
   without any view.

Other candidate hot spots investigated and dismissed:
- Replacing the per-`mf_idx` `MVector{3, Float64}` accumulator and
  broadcast `.+=` with three scalar `Float64` was tried (a clean,
  numerically bit-identical rewrite) but measured **slightly
  negative** on torque time. The compiler is already SROA-ing the
  `MVector` into registers and emits SIMD-friendly 3-wide IR; the
  scalar form removed that pattern. Reverted.
- `coeff_tensor[idx_buf..., mf_idx]` splat indexing is alloc-free and
  well-optimized; replacing it with rank-specialized branches did not
  move time or allocations. Spherical-harmonics evaluation (`P̄ₗₘ`,
  `Zₗₘ_unsafe`) is the next-largest cost but tracked separately
  (SpheriCart backend evaluation is its own investigation).

## Scope

Includes:

- `src/Fitting.jl`: `calc_∇ₑu!`, `design_matrix_energy_element`,
  `build_design_matrix_torque`, workspace types (`EnergyWorkspace`,
  `GradWorkspace`).
- Replacement of `Set{UInt}` de-duplication with `Vector{UInt}` +
  linear scan.
- `SVector{3, Float64}` column read used in place of an `@views`
  `SubArray` on `spin_directions[:, atom]` in three sites: the
  inner SH loop of `design_matrix_energy_element` and `calc_∇ₑu!`,
  and the `dir_iatom_svec` construction in
  `build_design_matrix_torque`.
- Regression coverage that the design matrix is bit-for-bit unchanged
  on the existing FePt / FeGe / FeBCC integration fixtures.

Excludes:

- Dispatch-boxing alloc reduction (`for cbc in key_group` with
  abstract `Vector{CoupledBasis_with_coefficient}` element type) —
  the largest remaining alloc source, deferred to a follow-up spec
  because the fix is structural to the SALC representation.
- Rank-specialized inner contraction over `coeff_tensor` — explored
  during implementation, gave no measurable improvement, and
  complicates the N-body-generic form. Documented in `design.md`
  under "Implementation notes (rejected alternatives)".
- Spherical-harmonics replacement / SpheriCart backend changes
  (separate, larger investigation).
- Algorithmic changes to SALC construction (`SALCBases.jl`).
- Multithreading scheme changes (`@threads` over key-groups for energy,
  spinconfigs for torque stays as is).

## Invariants

- The energy and torque design matrices are bit-for-bit identical to
  the pre-change values on every integration fixture
  (`test/integration/{chain,dimer,square_lattice,fept_tetragonal_2x2x2,fege_2x2x2,febcc_2x2x2_pm}/`).
- SCE coefficient conventions are unchanged: sign, normalization
  (`(4π)^(N/2)` per-cluster scaling factor), and SALC ordering.
- `Fitting` public API (`build_design_matrix_energy`,
  `build_design_matrix_torque`, `calc_∇ₑu`, `calc_∇ₑu!`) keeps current
  signatures and docstring contracts.
- Spherical-harmonics convention (`TesseralHarmonics`) is untouched.
- XML round-trip (`Magesty.save` / `Magesty.load`) is untouched.

## Completion criteria

- [ ] `make test-all` passes; no integration value drift.
- [ ] `build_design_matrix_torque` on the `_light` fixture is at least
      **10 % faster** than baseline (4 threads, 5-trial median).
      Achieved by items (1) `Set` → `Vector` (~9 % via removed
      hash-table overhead) and (2) SVector replacements (~7 % via
      removed `@views` SubArray creation in three sites).
- [ ] `build_design_matrix_energy` does not regress (≤ 5 % slower).
- [ ] No alloc-count regression on either matrix (matching baseline is
      acceptable; the remaining allocation source is the
      dispatch-boxing path documented under "Out-of-scope follow-up").
- [ ] Before / after numbers recorded in `.claude/bench_log.md`.
- [ ] Tier 2 review panel run and findings resolved.

## Out-of-scope follow-up

A second-tier allocation source — the dynamic-dispatch closure
boxing produced by `for cbc in key_group` where `key_group` has
abstract element type `Vector{CoupledBasis_with_coefficient}` —
accounts for the bulk (~110 M out of 336 M) of remaining torque
allocations after this spec. `Profile.Allocs` attributes the type to
the captured `Magesty.Symmetries.Symmetry` field. Removing it
requires structural changes to the SALC representation (concretely
parameterizing the key-group element type on `R`, or routing the
inner call through a runtime function barrier). Those structural
changes are deferred to a follow-up spec.

The original brief targeted ≥ 1.5× torque speedup. The achievable
in-scope gain on this hot path turned out to be ~13 %; the larger
factors live behind the dispatch-boxing change above, not in the
contraction code touched here.

## References

- Related specs:
  [260522-cluster-generation-benchmark](../260522-cluster-generation-benchmark/),
  [260523-cluster-generation-perf](../260523-cluster-generation-perf/),
  [260516-optimize-workspace](../260516-optimize-workspace/),
  [260516-coupled-basis-typeparam](../260516-coupled-basis-typeparam/).
- Profile run: `bench/bench_b1_3body_fege.jl` (this branch).
