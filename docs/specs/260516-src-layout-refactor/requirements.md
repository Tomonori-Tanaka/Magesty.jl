# Requirements: `src/` layout refactor for publication readiness

## Background

Ahead of publishing `v0.1.0-DEV` (intended for GitHub release and the
General registry), tidy up `src/` file names, module names, and
directory layout so a first-time reader can follow them without
friction. Implementation behavior is unchanged — names and placement
only.

## Goals

1. Make file names match the `module X` declared inside them (the
   standard Julia convention).
2. Rename names that feel awkward for a public package, such as the
   `My` prefix family.
3. Resolve the current `common/` / `types/` / `utils/` split: the
   boundaries are vague and file counts are skewed.
4. Make the `include` block at the top of `src/Magesty.jl` readable as
   a dependency-ordered sequence.

## Non-goals

- Signature changes to the public API (`export`-ed symbols). Renames
  affect **internal module names only**.
- Physics-convention changes (signs, units, normalization).
- Performance changes.
- Test-structure changes other than adding / fixing tests.
- Restructuring `tools/` (its internal refactor is a separate effort).

## Invariants

R1. **Public API stays identical**: `SCEBasis`, `SCEDataset`,
    `SCEFit`, `SCEModel`, `OLS`, `Ridge`, `AbstractEstimator`, `fit`,
    `coef`, `intercept`, `nobs`, `dof`, `predict_energy`,
    `predict_torque`, `r2_*`, `rmse_*`, `rss_*`, `residuals_*`,
    `save`, `load`, `read_embset`, `VERSION`, `install_tools` — all
    callable with the same signatures.

R2. **Tests stay green as before**: `make test-unit` /
    `make test-integration` / `make test-aqua` / `make test-jet` /
    `make test-tools` produce the same outcomes as pre-refactor.

R3. **XML schema and EMBSET format unchanged**: existing `.xml` /
    `EMBSET.dat` files continue to load.

R4. **Docs build passes**: `make build` runs cleanly. `@docs` blocks
    follow the renamed modules / functions
    (e.g., `MySphericalHarmonics.Zₗₘ` -> the new name).

R5. **`SPEC.md` reflects the post-refactor layout.**

## Scope

### In-scope

- Renames and directory restructure under `src/` (file names,
  directory layout, internal module names).
- Adjust `include` paths and `using .X` in `src/Magesty.jl`.
- Sync the directory diagram and table in `SPEC.md`.
- Sync `@docs` module names in `docs/src/api.md`.
- Sync references to names like `MySphericalHarmonics` inside
  `CLAUDE.md`.
- Update `using` statements in `test/` when a module name changes.
- **`tools/` follow-on edits only**: minimal `include` / `using` /
  path-reference sync forced by the rename or move; no structural
  change.

### Out-of-scope

- Restructuring or rewriting `tools/` scripts (separate effort).
- `docs/specs/` (past spec bodies are kept as history).
- `examples/` (uses only the public API, so unaffected).

## Completion criteria

C1. Every `.jl` file under `src/` agrees with its `module`
    declaration.
C2. Every name flagged for reconsideration (e.g., the `My` prefix) is
    resolved.
C3. `make test-all` / `make test-aqua` / `make test-jet` /
    `make test-tools` / `make build` all green.
C4. Public-API tests under `test/component/` keep working without any
    edits to call sites that start from `using Magesty`.
C5. `SPEC.md` and `docs/src/api.md` reflect the post-refactor layout.
C6. No references to old paths or old module names remain under
    `tools/` (confirm via `grep`).
C7. The PR description includes both "public API unchanged" and the
    rename-mapping table.
