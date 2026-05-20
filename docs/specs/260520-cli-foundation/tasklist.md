# Tasklist: CLI foundation (Comonicon) with VASP-to-extxyz pilot

Status: draft (2026-05-20)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — Comonicon integration and CLI skeleton (done)

- [x] Add `Comonicon` to `Project.toml` `[deps]`; `[compat] Comonicon = "1"`.
- [x] Add `[tool.comonicon]` (`name = "magesty"`).
- [x] Create `src/CLI.jl` with a minimal `Comonicon.@main`; `include` it
      last in `src/Magesty.jl`.
- [x] Install path: `deps/build.jl` calls `Magesty.comonicon_install()`;
      documented in `installation.md`.
- [x] `magesty`, `magesty --help`, `magesty version` work after install.
- [x] Load-time delta measured (~0.20 s upper bound); Comonicon kept as a
      plain `[deps]` entry (see design.md "Risks and open items").

### M2 — Move conversion modules into `src/` (done)

- [x] `tools/ExtXYZ.jl` -> `src/ExtXYZ.jl` (module `ExtXYZWriter` ->
      `ExtXYZ`); `tools/vasp/VaspParser.jl` -> `src/VaspIO.jl` (module
      `VaspParser` -> `VaspIO`); both `include`d in `src/Magesty.jl`
      after `XMLIO.jl`.
- [x] Relocate `test_ExtXYZ.jl` / `test_VaspParser.jl` (-> `test_VaspIO.jl`)
      and `fixtures/` to `test/component/`; retargeted at `Magesty.ExtXYZ`
      / `Magesty.VaspIO`; added to `test/runtests.jl`.
- [x] Update `tools/vasp/vasp2extxyz.jl` and `vasp2extxyz_recursive.jl` to
      `using Magesty.ExtXYZ` / `Magesty.VaspIO` (both kept as scripts;
      `vasp2extxyz.jl` is removed in M3, recursive deferred).
- [x] Remove `tools/ExtXYZ.jl`, `tools/vasp/VaspParser.jl`, `tools/test/`;
      drop the `make test-tools` target and its `ci-local` reference;
      drop the CI `Tools` job; update the `test-tools` / `tools/test/`
      references in `CLAUDE.md`, `CONTRIBUTING.md`, and `SPEC.md`.

### M3 — `vasp_to_extxyz` API and `vasp extxyz` subcommand

- [ ] Implement `vasp_to_extxyz` in `src/CLI.jl` (logic copied verbatim
      from `tools/vasp/vasp2extxyz.jl`); export it.
- [ ] Add the `@cast module Vasp` command group with an `extxyz` leaf
      subcommand wrapping it (`magesty vasp extxyz`).
- [ ] Regression test: `vasp_to_extxyz` byte-identical output vs. a
      stored reference.
- [ ] Smoke test: `magesty vasp extxyz` converts a fixture and exits 0
      with output matching the reference.
- [ ] Remove `tools/vasp/vasp2extxyz.jl`.

### M4 — Documentation and cleanup

- [ ] Polish the `installation.md` "Command-Line Interface" section
      (drafted in M1) once the `vasp extxyz` subcommand exists.
- [ ] Update `tools.md` (drop the removed `vasp2extxyz.jl` entry),
      `api.md`, `SPEC.md`, `CHANGELOG.md`.
- [ ] Sweep `.claude/agents/` for module / test-layout references
      (`test-runner.md` test-target table, `spec-reviewer.md` `tools/test`
      mention) and the `tools/test` reference in the pre-release-cleanup
      design note.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] If results changed: regression or validation test added.
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`. (N/A — no hot path touched.)
- [ ] If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
