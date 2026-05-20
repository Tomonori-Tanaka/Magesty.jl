# Tasklist: extract the CLI into a `MagestyCLI` subdirectory package

Status: draft (2026-05-21)

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

- [ ] Add a Makefile target that instantiates and tests `cli/`.
- [ ] Add a `MagestyCLI` job/step to `.github/workflows/CI.yml`.
- [ ] Sweep `.claude/agents/` if a new Makefile target was added.

### M4 — Documentation and close-out

- [ ] Update `docs/src/installation.md` / `docs/src/tools.md` install path.
- [ ] Update `SPEC.md` (module table, directory layout) and `CHANGELOG.md`.
- [ ] Update the CLI-package-extraction design note `Status:` line.
- [ ] Update this `Status:` line and the `docs/specs/README.md` table.

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (`test/jet.jl` is a plain
      assertion again, 0 issues).
- [ ] ~~If results changed: regression or validation test added.~~ (no
      numerical change)
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
      (core API unchanged; `SPEC.md` still updated for the layout change)
- [ ] ~~If a hot path was touched: `.claude/bench_log.md`.~~ (no hot path)
- [ ] If module names or Makefile targets changed: `.claude/agents/`
      swept and updated.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
