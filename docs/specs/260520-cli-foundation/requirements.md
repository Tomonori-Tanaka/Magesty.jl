# Requirements: CLI foundation (Comonicon) with VASP-to-extxyz pilot

Status: draft (2026-05-20)

## Goal

Establish a proper, installable command-line interface for Magesty.jl
built on Comonicon.jl, and validate it end-to-end by promoting one tool —
the VASP-to-extxyz converter — into a package API function plus a CLI
subcommand. Remaining `tools/` scripts are migrated in later specs that
follow the pattern established here.

## Background

`docs/src/installation.md`'s "CLI Tools" section currently tells users to
run `tools/*.jl` scripts by filesystem path. That is fragile: a
`Pkg.add`-installed package lives under a version-hashed directory that
moves on every update, and distributing loosely-runnable scripts is not
standard Julia practice — functionality belongs in the package API and/or
a real installable command.

A previous `Magesty.install_tools()` function attempted this by writing
shell shims that baked an absolute `pkgdir` path; it broke on package
updates and was removed. A first cleanup step has already landed
(commit `bbf0220`): five `tools/` scripts that loaded the package via
`include("../src/Magesty.jl")` now use `using Magesty`, so those tools
depend on the real precompiled package rather than an anonymous source
copy. `vasp2extxyz_recursive.jl` was not among the five — it `include`s
the local converter modules directly and is updated in M2 of this spec.

This spec is the next step: a Comonicon-based CLI whose installed launcher
resolves the package by name through an environment, so it follows package
updates instead of breaking.

## Scope

Includes:

- Add `Comonicon` to `Project.toml` `[deps]` and `[compat]`; add the
  `[tool.comonicon]` configuration (command name `magesty`).
- A CLI entrypoint in `src/` (Comonicon `@cast` / `@main`) and a documented
  install path (plain install only; no sysimage/app build).
- Move the pure-library conversion modules into the package:
  `tools/ExtXYZ.jl` -> `src/ExtXYZ.jl`, `tools/vasp/VaspParser.jl` ->
  `src/VaspIO.jl`.
- A high-level API function `vasp_to_extxyz` and a `magesty vasp extxyz`
  CLI subcommand — the pilot. Subcommands are grouped per DFT code
  (`vasp` is a Comonicon `@cast module` group).
- Relocate `tools/test/test_ExtXYZ.jl` / `test_VaspParser.jl` into
  `test/component/`, retargeted at the new `src/` modules.
- Update `vasp2extxyz_recursive.jl` to consume the package modules via
  `using Magesty` (kept as a script; full promotion deferred).
- Documentation: rewrite the installation "CLI Tools" section; update
  `tools.md`, `api.md`, `SPEC.md`, `CHANGELOG.md`.

Excludes:

- Promotion of any other tool (`vasp2extxyz_recursive`, `pos2toml`,
  `oszicar2magmom`, plotting / MFA / sampling tools) — separate specs.
- Comonicon `[sysimg]` / `[application]` build modes (deferred; plain
  install only).

## Invariants

- The extxyz output of `vasp_to_extxyz` / the `magesty vasp extxyz`
  subcommand is byte-identical to the current `tools/vasp/vasp2extxyz.jl`
  for the same inputs.
- No physics-convention, sign, unit, or numerical change. The conversion
  logic is moved verbatim.
- SCE XML `save` / `load`, spherical-harmonics conventions, and the
  `SALCBasis` <-> `Fitting` ordering are untouched.
- `using Magesty` continues to work (it now also loads Comonicon).

## Completion criteria

- [ ] `make test-all`, `make test-aqua`, `make test-jet` pass cleanly.
- [ ] The `magesty` command installs via the documented path; `magesty`,
      `magesty --help`, and `magesty vasp extxyz --help` work.
- [ ] Regression test: `vasp_to_extxyz` output is byte-identical to a
      stored reference from the existing `tools/test/fixtures/`.
- [ ] End-to-end: `magesty vasp extxyz` converts a fixture input, exits 0,
      and produces output matching the reference.
- [ ] `ExtXYZ` / `VaspIO` component tests pass under `test/component/`.
- [ ] `installation.md` "CLI Tools", `tools.md`, `api.md`, `SPEC.md`
      updated; `CHANGELOG.md` `[Unreleased]` updated.

## References

- Related commit: `bbf0220` (`refactor(tools): replace source include with
  package import`).
- Removed predecessor: `Magesty.install_tools()` (added `c67581a`, removed
  `53476a5`).
- Comonicon.jl: https://comonicon.org
