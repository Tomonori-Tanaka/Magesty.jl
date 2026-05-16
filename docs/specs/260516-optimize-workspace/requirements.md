# Requirements: Optimize hot-path scratch workspace

Created 2026-05-16. Topic branch: `refactor/optimize-workspace`.

## Purpose

Eliminate per-call heap allocations in `design_matrix_energy_element`
and `calc_∇ₑu!` by hoisting the remaining scratch state into a reusable
per-thread workspace passed by the caller (`build_design_matrix_energy`,
`build_design_matrix_torque`, `_predict_energy`).

This is the natural follow-up to B1+B3: those refactors made the
indexing path type-stable and converted small scratch buffers to
`MVector`, but per-call heap allocations from the remaining scratch
state (`Set{UInt}`, `Vector{Vector{Float64}}` for spherical harmonics)
still dominate the remaining 262K / 1.17M allocations.

## Scope

In scope:

- Define a per-thread workspace struct (or pair of structs) that owns
  the resizable scratch:
  - `searched_pairs::Set{UInt}` (cleared via `empty!` per call)
  - `sh_values::Vector{Vector{Float64}}` (outer resized to N=R-1, each
    inner resized to `2lᵢ+1`)
  - `atom_grad_values::Vector{Vector{Float64}}` (grad only)
- Thread the workspace through `design_matrix_energy_element` and
  `calc_∇ₑu!`.
- Caller (`build_design_matrix_energy`, `build_design_matrix_torque`,
  `_predict_energy`) creates one workspace per thread inside its
  `@threads` loop.
- For the non-bang `calc_∇ₑu` convenience wrapper, allocate a workspace
  internally (it's not a hot path).

Out of scope:

- 旧 #10 (`atoms::SVector`).
- B2 architectural rewrite of `sh_values` to a flat storage layout.
  The current `Vector{Vector{Float64}}` is kept but **reused**.
- API surface change to public `SCEDataset` / `SCEFit` constructors.

## Invariants

1. **Numerical equivalence**: every existing test passes bit-stable.
2. **Thread safety**: each `@threads` worker holds its own workspace.
   No mutable state shared between threads.
3. **Public API unchanged**: `SCEDataset`, `SCEFit`, `predict_energy`,
   `predict_torque` signatures unchanged. The workspace is purely an
   internal optimization concern.
4. **No regression on cold path**: top-level callers without
   `@threads` (e.g. one-shot `predict_energy` for inference) continue
   to work; the internal helper that allocates a workspace on first
   use is sufficient.

## Completion criteria

- `make test-all` green; `make test-aqua`, `make test-jet` green.
- Per-call allocations for `design_matrix_energy_element` drop from
  ~74-110/call to ≤ 5/call (target: 0-1 once workspace is warm).
- `build_design_matrix_energy` total alloc count drops by 5×+;
  `build_design_matrix_torque` similarly.
- Wall-time improvement recorded in `.claude/bench_log.md`. Even if
  wall-time delta is modest (current 2.1ms energy is small), the
  per-call alloc reduction translates to GC pressure relief in
  multi-config fits.

## Risks

- **Workspace lifecycle**: Set / Vector{Vector{Float64}} need correct
  reset semantics. `empty!(::Set)` keeps capacity (good).
  `Vector{Vector{Float64}}` resize: if outer resized smaller, inner
  vectors are reused; if larger, new inner vectors must be allocated
  lazily. Helper `_ensure_sh_buffers!(ws, R, ls)` handles this.
- **Threading correctness**: must verify each `@threads` iteration
  uses its thread-local workspace, not a shared one. Pattern:
  `ws = workspaces[threadid()]` or `ws = EnergyWorkspace()` inside
  the threaded loop (the latter is simpler and only allocates
  `nthreads()` workspaces).
- **`calc_∇ₑu` non-bang wrapper**: must allocate a workspace
  internally; cannot share with caller. This is fine since it's not
  hot-path.
