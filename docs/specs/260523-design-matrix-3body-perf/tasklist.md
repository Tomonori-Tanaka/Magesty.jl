# Tasklist: Design-matrix construction perf for 3-body SALCs

Status: complete (2026-05-23)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Replace `searched_pairs::Set{UInt}` with `Vector{UInt}`

- [ ] Change the field type on both `EnergyWorkspace` and
      `GradWorkspace`. Initializers updated.
- [ ] In-function usage unchanged (`empty!` / `in` / `push!` all
      apply to `Vector{UInt}` with the same semantics, with `in`
      switching from O(1) hash to O(N) linear scan; N is small).
- [ ] Bench on `_light` confirms ~5–10 % torque time reduction from
      M1 alone.

### M2 — `SVector{3, Float64}` column reads instead of `@views` SubArrays

- [ ] In `design_matrix_energy_element`: read
      `(spin_directions[1, atom], spin_directions[2, atom],
      spin_directions[3, atom])` into an `SVector{3, Float64}` once
      per `(itrans, site)`; drop the `@views` SubArray on
      `spin_directions[:, atom]`.
- [ ] In `calc_∇ₑu!`: same replacement for the inner SH loop
      (`Zₗₘ_unsafe` and `∂ᵢZlm_unsafe` both see the `SVector`).
- [ ] In `build_design_matrix_torque`: replace
      `@views dir_iatom = spinconfig.spin_directions[:, iatom];
      dir_iatom_svec = SVector{3, Float64}(dir_iatom)` with three
      direct scalar reads into the `SVector` constructor.
- [ ] Bench confirms additional ~5–7 % torque-time reduction on
      `_light` and that `make test-all` is still green.

### M3 — Bench bookkeeping and review

- [ ] Append before / after to `.claude/bench_log.md` for both
      `build_design_matrix_energy` and `build_design_matrix_torque`
      on the `_light` fixture, 4-thread medians (5 trials).
- [ ] Run Tier 2 review panel (numerical / maintainability /
      performance / API).
- [ ] Resolve every blocker / major finding; record contention
      escalations to the user as required.

## Rejected during implementation

- **Scalar `Float64` gradient accumulators** in place of
  `mf_grad_contribution` / `grad_result` `MVector{3, Float64}`
  broadcast `.+=`. Bit-identical, but slightly slower; the
  compiler is already SROA-ing the `MVector` and emitting
  3-wide IR. See `design.md` "Implementation notes".
- **Rank-specialized `if R == 2 / 3 / 4` contractions** of
  `coeff_tensor[idx_buf..., mf_idx]`. No measurable gain.
- **Generic `CartesianIndices(dims_t)` contraction with
  workspace-backed `Vector{Int}` scratch.** -67 % allocs but
  +8 % torque time; the remaining allocs were dispatch-boxing
  outside this spec's scope.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- ~~If results changed: regression or validation test added.~~
  (No numerical change; the same integration fixtures cover the
  rewritten paths.)
- ~~If public API changed: `SPEC.md` and `docs/src/api.md` updated.~~
- [ ] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~
- ~~`CHANGELOG.md` `[Unreleased]` updated.~~ (Internal perf only; no
  user-visible behavior change.)
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
