# Tasklist: CLI foundation (Comonicon) with VASP-to-extxyz pilot

Status: complete (2026-05-21)

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

### M3 — `vasp_to_extxyz` API and `vasp extxyz` subcommand (done)

- [x] Implement `vasp_to_extxyz` in `src/CLI.jl` (logic copied verbatim
      from `tools/vasp/vasp2extxyz.jl`); exported from `Magesty`.
- [x] Add the `Comonicon.@cast module Vasp` command group with an
      `extxyz` leaf subcommand wrapping it (`magesty vasp extxyz`).
- [x] Regenerate the golden `FeRh.extxyz` / `IrMn3.extxyz` fixtures from
      the current `vasp2extxyz.jl` — the committed copies were stale
      (`MAGMOM_smoothed` vs. the current `magmom_smoothed`).
- [x] Regression test (`test/component/test_vasp_to_extxyz.jl`):
      `vasp_to_extxyz`, its `output` keyword, and the `magesty vasp
      extxyz` subcommand all produce byte-identical output vs. the golden
      fixtures; subcommand exits 0.
- [x] Remove `tools/vasp/vasp2extxyz.jl`.

### M4 — Documentation and cleanup (done)

- [x] Polished the `installation.md` "Command-Line Interface" section
      with a `magesty vasp extxyz` example.
- [x] Updated `tools.md` (the `vasp2extxyz.jl` entry became
      `magesty vasp extxyz`), `api.md` (new "VASP conversion" section),
      `SPEC.md` (`ExtXYZ` / `VaspIO` / `CLI` modules, `vasp_to_extxyz`,
      `Comonicon`), and `CHANGELOG.md`.
- [x] Swept `.claude/agents/` (`test-runner.md` test-target table,
      `spec-reviewer.md` `tools/test` mention) and the
      pre-release-cleanup design note.
- [x] `test/jet.jl` marked `@test_broken`: JET's `report_package` cannot
      macro-expand Comonicon's `@main` under its loader. Tracked in
      `DESIGN_NOTES.md` (extract the CLI into a separate package); the
      guard reverts to a real assertion once `report_package` succeeds.

## Exit checklist

- [x] `make test-all` passes.
- [x] `make test-aqua` clean; `make test-jet` exits clean — JET is
      `@test_broken` (Comonicon/JET incompatibility, see M4 and
      `DESIGN_NOTES.md`).
- [x] Regression test added (`test/component/test_vasp_to_extxyz.jl`);
      behavior-preserving — no numerical results changed.
- [x] Public API changed (`vasp_to_extxyz` exported): `SPEC.md` and
      `docs/src/api.md` updated.
- ~~If a hot path was touched: bench_log entry.~~ N/A — no hot path.
- [x] Module / Makefile-target changes: `.claude/agents/` swept.
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line here and the `docs/specs/README.md` row synced.
- [x] Implementation commit hashes appended below.

## Commits

- `f6d34dc` — spec
- `94d8238` — M1 (Comonicon integration and CLI skeleton)
- `874e9ac` — M2 (VASP/extxyz converters promoted into `src/`)
- `36b67a0` — M3 (`vasp_to_extxyz` API and `magesty vasp extxyz`)
- M4 (documentation and cleanup) — this commit
