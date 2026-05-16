# Tasklist: CoupledBasis.atoms as SVector

**Status: ABANDONED (2026-05-16).** Spec retained as a record so the
same trap is not re-attempted blindly. See "Outcome" below.

## Step 1 — baseline benchmark (done)

- [x] Baseline captured on fept_tetragonal_2x2x2 (JULIA_NUM_THREADS=4):
  - `SALCBasis` full build: 0.043 s / 108 MB / 2.73 M allocs
  - `build_design_matrix_energy`: 0.0015 s / 478 KB / 9784 allocs
  - `build_design_matrix_torque`: 0.0032 s / 6.95 MB / 112827 allocs

## Step 2 — field type change (done, later reverted)

- [x] `CoupledBasis{R, N}.atoms` → `SVector{N, Int}` (added `N` type
      parameter; `R = N + 1` invariant).
- [x] `CoupledBasis_with_coefficient{R, N}.atoms` → `SVector{N, Int}`.
- [x] `AngularMomentumCouplingResult` unchanged (no `atoms` field).
- [x] Inner ctors construct `SVector{N, Int}(ntuple(i -> Int(atoms[i]), N))`.
- [x] `reorder_atoms` forwards via inner ctor.

## Step 3 — caller adjustments (done, later reverted)

- [x] `SALCBases.jl` hot loop tried `map(a -> symmetry.map_sym[a, n], cb1.atoms)`;
      this regressed (closure capture of `n` / `map_sym`). Reverted to the
      original hoisted `Vector{Int}` buffer.
- [x] `Optimize.jl` `where {R}` dispatch continued to work transparently
      with the new `{R, N}` struct (N becomes a UnionAll-bound free
      parameter at the signature level, specialized at call sites).
- [x] `xml_io.jl`, `Basis.jl:isless / show` transparent under `SVector`.

## Step 4 — tests (done)

- [x] `make test-unit` (20607/20607) green.
- [x] `make test-integration` (813/813) green; XML byte-identical.
- [x] `make test-aqua` (10/10) green.
- [x] `make test-jet` green.

## Step 5 — benchmark (done; result: regression)

After (atoms→SVector, no caller rewrite that helped):

| Target | Baseline | After | Delta |
|---|---|---|---|
| `SALCBasis` full build | 0.043 s / 108 MB / 2.73 M allocs | 0.048 s / 114 MB / 2.87 M allocs | **+12% / +6 MB / +140K allocs** |
| `build_design_matrix_energy` | 0.0015 s / 478 KB / 9784 allocs | 0.0015 s / 479 KB / 9815 allocs | ≈no change |
| `build_design_matrix_torque` | 0.0032 s / 6.95 MB / 112827 allocs | 0.0033 s / 6.95 MB / 112858 allocs | ≈no change |

## Step 6 — close out (abandoned)

- [ ] ~~code-reviewer + commit + PR~~ — refactor abandoned; src reverted
      to `main` state, topic branch deleted.

## Outcome and rationale

Hot design-matrix kernels (`design_matrix_energy_element`, `calc_∇ₑu!`)
are insulated by function barriers (`where {R}` specialization) and
saw no measurable change.

`projection_matrix_coupled_basis` iterates the coupled basis list as a
UnionAll-typed container (`SortedCounter{CoupledBasis}` /
`Vector{CoupledBasis}`). With `atoms::SVector{N, Int}`, the per-element
length `N` is type-erased at the iteration site, so the SVector
indexing path cannot be specialized. `Vector{Int}` indexing remains
fast in the same UnionAll context because its element type is fixed
and its length is dynamic by design. Net effect: SALC build regressed
~12% in time, ~6 MB / +140K in allocations, with no offsetting win on
the design-matrix path.

Re-attempt prerequisites (if revisited later):

- Container specialization: store `CoupledBasis` instances in a
  container whose element type is concrete (`Vector{CoupledBasis{R, N}}`
  per (R, N) group, or grouped iteration via function barrier). Only
  then does inline `SVector{N, Int}` indexing become specializable in
  the loop body.
- Eliminate the `map(closure, atoms)` hot-loop rewrite, or rewrite it
  with `let`-bound captures to avoid boxing.

Until those prerequisites are met, this refactor is net negative and
should not be re-attempted.
