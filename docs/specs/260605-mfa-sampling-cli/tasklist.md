# Tasklist: MFA spin-sampling CLI (`magesty vasp mfa`)

Status: complete (2026-06-05)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Layer 1: code-agnostic MFA sampler

- [x] `src/MfaSampling.jl`: `thermal_averaged_m`, `tau_from_magnetization`,
      hand-rolled `sample_vmf_direction`, `mfa_sample`, rotation helpers,
      `parse_atom_index_spec`, `mfa_sweep` driver.
- [x] Add `Roots` to `Project.toml` `[deps]` / `[compat]`.
- [x] `include` in `src/Magesty.jl`.

### M2 — Layer 2: INCAR I/O + orchestration

- [x] `src/IncarIO.jl`: port `parse_incar` / `write_incar`
      (`DataStructures.OrderedDict`).
- [x] `src/VaspSampling.jl`: `sample_mfa_incar`; `export`.

### M3 — Layer 3: CLI

- [x] `mfa` leaf in `@cast module Vasp` (`cli/src/MagestyCLI.jl`).

### M4 — Tests

- [x] `test_MfaSampling.jl` (analytic vMF validation), `test_IncarIO.jl`,
      CLI regression in `cli/test/runtests.jl`, INCAR fixture.
- Exit: `make test-all` and `make test-cli` pass with the new tests wired in.

### M5 — Docs + cleanup + review

- [x] `CHANGELOG.md`, `docs/src/api.md`, CLI narrative docs, `SPEC.md`.
- [x] Remove `tools/sampling_mfa.jl` (keep `tools/vasp/vasptools.jl` —
      used by `tools/personal/`).
- [x] Tier-2 review panel.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes.
- [x] `make test-cli` passes.
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] If results changed: regression or validation test added (vMF analytic
      validation; behavior otherwise preserved).
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [x] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md` (sampler hoist optimization recorded).
- [x] Tier 2 review panel run and findings resolved.
- [x] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ (no agent module enumerations
      needed updating; `test-runner.md` reference is illustrative)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [x] Implementation commit hash appended below.

## Implementation

- `38fb9d8` — feat(cli): promote MFA sampler to exported API and
  `magesty vasp mfa` command (branch `feat/mfa-sampling-cli`).
