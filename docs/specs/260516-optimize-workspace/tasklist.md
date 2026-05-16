# Tasklist: Optimize hot-path scratch workspace

## Step 1 — baseline (done)

- [x] Per-call alloc profile captured via
      `test/develop_tmp/profile_remaining_allocs.jl`:
      `design_matrix_energy_element` ≈ 74-110 allocs/call;
      `calc_∇ₑu!` ≈ 10-50 allocs/call.

## Step 2 — workspace structs (done)

- [x] `EnergyWorkspace` / `GradWorkspace` mutable structs added in
      `src/Optimize.jl` (not exported).
- [x] `_ensure_sh_buffer!(buf, ls)` helper added for `sh_values`.
- [x] Discovery during Step 2: `∂ᵢZlm_unsafe` returns
      `SVector{3, Float64}`. The pre-refactor
      `atom_grad_values::Vector{Vector{Float64}}` was forcing a
      `convert` on every assignment. `GradWorkspace.atom_grad_values`
      is now `Vector{SVector{3, Float64}}`, which only needs
      `resize!` (no inner alloc) — that subsumes the "grad analog
      helper" originally planned.

## Step 3 — thread workspace through kernels (done, combined with Step 4)

- [x] `design_matrix_energy_element`: takes `ws::EnergyWorkspace`.
- [x] `calc_∇ₑu!`: takes `ws::GradWorkspace`.
- [x] `calc_∇ₑu` non-bang wrapper: allocates `GradWorkspace`
      internally.
- [x] Discovery: the dominant remaining allocation was *not* the Set
      or Vector{Vector} containers — it was the unbuffered
      `Zₗₘ_unsafe(l, m, uvec)` / `∂ᵢZlm_unsafe(l, m, uvec)` paths,
      which each reallocate a Legendre cache per `(l, m)` call. Added
      `legendre_buf::Vector{Float64}` to both workspaces and switched
      the kernels to the 4-argument buffered overloads. This single
      change accounts for most of the remaining reduction.

## Step 4 — caller integration (done with Step 3)

- [x] `build_design_matrix_energy`: `EnergyWorkspace` allocated inside
      `@threads` loop (one per column iteration).
- [x] `build_design_matrix_torque`: `GradWorkspace` inside `@threads`.
- [x] `_predict_energy`: one `EnergyWorkspace` at the top.
- [x] `_predict_torque`: one `GradWorkspace` at the top.

## Step 5 — tests (done)

- [x] `make test-unit` (20607/20607).
- [x] `make test-integration` (813/813); XML byte-stability passes.
- [x] `make test-aqua` (10/10).
- [x] `make test-jet`: no issues.

## Step 6 — benchmark (done)

- [x] Per-call allocs (`profile_remaining_allocs.jl`):
      energy 74-110 → **2**, grad 10-50 → **6**.
- [x] `bench_b1_design_matrix.jl` on fept fixture:
      energy 2.1 ms → 1.5 ms (×1.4), torque 6.4 ms → 4.1 ms (×1.6);
      total allocs 262K → 9.8K (energy, ×27), 1.17M → 113K (torque, ×10).
      Recorded under "C — workspace" in `.claude/bench_log.md`.

## Step 7 — close out

- [x] `code-reviewer` agent pass: 1 Critical (sequential
      `legendre_buf` safety) addressed via clarifying comment; 1
      Minor (`@inbounds` scope) tightened.
- [x] Single commit via `git-helper`: `aca0c9a perf(optimize): pool
      hot-path scratch in workspace structs`.
- [x] PR #8 squash-merged to `main` as `e1a4326` on 2026-05-16.
      Topic branch deleted.
