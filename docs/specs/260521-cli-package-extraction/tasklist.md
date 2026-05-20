# Tasklist: extract the CLI into a `MagestyCLI` subdirectory package

Status: complete (2026-05-21)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

Branch: `refactor/cli-package-extraction`.

## Milestones

### M1 — De-Comonicon the core

- [x] Add `src/VaspConvert.jl` with `vasp_to_extxyz` and its two helpers
      (moved verbatim from `src/CLI.jl`).
- [x] Update `src/Magesty.jl`: `include("VaspConvert.jl")` in place of
      `include("CLI.jl")`; delete `src/CLI.jl`.
- [x] Remove `Comonicon` from `Project.toml` `[deps]` / `[compat]` and the
      `[tool.comonicon]` section; delete `deps/build.jl`.
- [x] Revert `test/jet.jl` to a plain `report_package` assertion.
- [x] Remove the `Magesty.command_main` CLI block from
      `test/component/test_vasp_to_extxyz.jl` (keep the API and
      `output`-keyword tests).
- [x] Verify `make test-all` and `make test-jet` pass.

### M2 — Create the `MagestyCLI` package

- [x] Create `cli/Project.toml`, `cli/Comonicon.toml`,
      `cli/src/MagestyCLI.jl`, `cli/deps/build.jl` (command tree moved
      from the old `src/CLI.jl`).
- [x] Add `cli/test/runtests.jl`: a `MagestyCLI.command_main` byte-identical
      `extxyz` check, reusing the core fixtures in
      `test/component/fixtures/` (resolved relative to the repo root).
- [x] Verify the installed `magesty vasp extxyz` produces byte-identical
      output to the golden fixtures.

### M3 — CI and Makefile integration

- [x] Add a Makefile target that instantiates and tests `cli/`
      (`test-cli`; also wired into `ci-local`).
- [x] Add a `MagestyCLI` job/step to `.github/workflows/CI.yml`.
- [x] Sweep `.claude/agents/` if a new Makefile target was added
      (`test-runner.md` target table updated).

### M4 — Documentation and close-out

- [x] Update `docs/src/installation.md` install path (`MagestyCLI` /
      `Pkg.build("MagestyCLI")`). `docs/src/tools.md` needed no change —
      it only links to the install instructions.
- [x] Update `SPEC.md` (module table, directory layout, external-library
      table) and `CHANGELOG.md`.
- [x] Update the CLI-package-extraction design note `Status:` line and
      the `DESIGN_NOTES.md` index row.
- [x] Update this `Status:` line and the `docs/specs/README.md` table.
- [x] Re-resolve `docs/Manifest.toml` (the core dropped `Comonicon`).

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [x] `make test-all` passes (22797 passed).
- [x] `make test-aqua` / `make test-jet` clean (`test/jet.jl` is a plain
      assertion again, 0 issues); `make test-cli` passes.
- [x] ~~If results changed: regression or validation test added.~~ (no
      numerical change)
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
      (core API unchanged; `SPEC.md` updated for the layout change;
      `api.md`'s `vasp_to_extxyz` `@docs` block is unaffected)
- [x] ~~If a hot path was touched: `.claude/bench_log.md`.~~ (no hot path)
- [x] If module names or Makefile targets changed: `.claude/agents/`
      swept and updated (`test-runner.md` target table).
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [x] Implementation commit hashes appended below.

## Commits

- `8314043` docs(specs): add the CLI-package-extraction spec
- `bbfd9dd` refactor(cli): remove the Comonicon CLI from the core package (M1)
- `82b5892` feat(cli): add the MagestyCLI package (M2)
- M3 (CI / Makefile integration) and M4 (docs / close-out): see the
  branch history for `refactor/cli-package-extraction`.
