# Requirements: extract the CLI into a `MagestyCLI` subdirectory package

Status: complete (2026-05-21)

## Goal

Move the Comonicon-based `magesty` command-line interface out of the
`Magesty` core package into a sibling `MagestyCLI` package located in the
`cli/` subdirectory of this repository. The core package then carries no
`Comonicon` dependency, which restores `JET.report_package(Magesty)` and
lets `test/jet.jl` return to a plain assertion.

## Background

Spec `260520-cli-foundation` added `src/CLI.jl`, which uses Comonicon's
`@cast` / `@main` macros. `JET.report_package` macro-expands every
top-level form in the package through a JuliaInterpreter-based loader, and
Comonicon's command-tree code generation throws under that loader
(`Comonicon.AST.Entry has no field name`). No JET configuration avoids
this — the crash is during macro expansion, which JET performs
unconditionally. `test/jet.jl` is therefore currently `@test_broken`,
leaving the numerical core without static analysis.

The fix is structural: remove the Comonicon code from the `Magesty`
source tree that `report_package` walks. The chosen layout — a separate
package in a `cli/` subdirectory rather than a separate repository — is
recorded in the design note on CLI package extraction (indexed in
`DESIGN_NOTES.md`).

## Scope

Includes:

- New package `MagestyCLI` under `cli/` (`cli/Project.toml`,
  `cli/src/MagestyCLI.jl`, `cli/deps/build.jl`, `cli/test/`).
- Move the Comonicon command tree (`version`, the `Vasp` command group,
  the `extxyz` subcommand, `@main`) from `src/CLI.jl` into `MagestyCLI`.
- Keep the `vasp_to_extxyz` API and the `VaspIO` / `ExtXYZ` modules in the
  core package; relocate them out of `src/CLI.jl` into a Comonicon-free
  core file.
- Remove `Comonicon` from the core `[deps]` / `[compat]` and drop the
  `[tool.comonicon]` section and `deps/build.jl` from the core.
- Revert `test/jet.jl` to a plain assertion.
- Update CI and the Makefile to build and test `MagestyCLI`.
- Update installation / CLI documentation for the new install path.

Excludes:

- Promoting the remaining `tools/` scripts to package API / CLI
  subcommands (future spec).
- Registering `Magesty` or `MagestyCLI` in any registry.
- Adding new CLI subcommands or changing existing CLI behavior.

## Invariants

- `vasp_to_extxyz` stays exported from `Magesty` with an unchanged
  signature and byte-identical output.
- `magesty vasp extxyz` produces byte-identical output to the current
  implementation.
- No physics / numerical convention changes — the CLI and the converter
  are pure I/O.
- The core public API (everything in `Magesty`'s `export` list) is
  unchanged.

## Completion criteria

- [x] `make test-all` passes.
- [x] `make test-jet` passes with a plain `@test` (no `@test_broken`);
      `report_package(Magesty)` returns and reports 0 issues.
- [x] `MagestyCLI` instantiates against the dev-ed core, and the installed
      `magesty vasp extxyz` produces byte-identical output to a recorded
      golden file.
- [x] CI builds and tests `MagestyCLI` (the `cli` job runs `make test-cli`).
- [x] `docs/src/installation.md` reflects the new install path; `SPEC.md`
      and `CHANGELOG.md` updated. (`docs/src/tools.md` needed no change.)

## References

- Related specs / design notes: spec `260520-cli-foundation` (predecessor);
  the CLI-package-extraction design note (indexed in `DESIGN_NOTES.md`).
