# Tasklist: `write_energies` / `write_torques` for the current API

Status: draft (2026-05-20)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Implement the writers

- [x] Add `src/FitCheckIO.jl` with `write_energies` / `write_torques`, the
      two signature families, and private resolution helpers.
- [x] Wire `include("FitCheckIO.jl")` and `export write_energies,
      write_torques` into `src/Magesty.jl`.

### M2 — Test

- [x] Add `test/component/test_fitcheck_io.jl`; register it in
      `test/runtests.jl`.
- [x] `make test-unit` / `make test-aqua` / `make test-jet` clean.
- [x] End-to-end: `FitCheck_energy.py` / `FitCheck_torque.py` parse the
      generated files.

### M3 — Document

- [x] Update `docs/src/api.md`, `docs/src/tools.md`, `SPEC.md`,
      `CHANGELOG.md` `[Unreleased]`.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-unit` passes (`make test-all` to be run with integration).
- [x] `make test-aqua` / `make test-jet` clean (no new warnings).
- [x] If results changed: regression or validation test added.
      (No numerical results change; format regression test added.)
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] ~~If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.~~ (No hot path touched.)
- [ ] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ (No such change.)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
