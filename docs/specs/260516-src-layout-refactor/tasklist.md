# Tasklist: `src/` layout refactor

Tracked in-session with `TaskCreate`. This file records commit-sized
milestones.

## Milestone 1: file-name alignment + directory flattening

- [ ] `src/common/version.jl` -> `src/Version.jl` (`git mv`)
- [ ] `src/common/SortedCounter.jl` -> `src/SortedCounters.jl`
- [ ] `src/types/AtomCells.jl` -> `src/AtomCells.jl`
- [ ] `src/types/Basis.jl` -> `src/Basis.jl` (no module rename yet)
- [ ] `src/utils/SphericalHarmonicsTransforms.jl` ->
      `src/SphericalHarmonicsTransforms.jl`
- [ ] `src/utils/AngularMomentumCoupling.jl` ->
      `src/AngularMomentumCoupling.jl`
- [ ] `src/utils/ConfigParser.jl` -> `src/ConfigParser.jl`
- [ ] `src/utils/atomsbase_adapter.jl` -> `src/AtomsBaseAdapter.jl`
- [ ] `src/utils/RotationMatrix.jl` -> `src/RotationMatrix.jl`
- [ ] `src/utils/MySphericalHarmonics.jl` ->
      `src/MySphericalHarmonics.jl` (no module rename yet)
- [ ] `src/utils/xml_io.jl` -> `src/XMLIO.jl`
- [ ] Remove the now-empty `common/`, `types/`, `utils/`.
- [ ] Update every `include(...)` in `src/Magesty.jl` to a top-level
      path.
- [ ] Confirm `make test-all` / `make test-tools` / `make test-aqua` /
      `make build` are green.
- [ ] Commit: `refactor(src): flatten directory structure and unify file naming`.

## Milestone 2: module renames

- [ ] `MySphericalHarmonics.jl` -> `TesseralHarmonics.jl` plus
      `module` rename.
- [ ] `Basis.jl` -> `CoupledBases.jl` plus `module` rename.
- [ ] `Optimize.jl` -> `Fitting.jl` plus `module` rename.
- [ ] Replace `using ..Basis` -> `using ..CoupledBases` across `src/`.
- [ ] Replace dotted references like `Basis.CoupledBasis` with
      `CoupledBases.CoupledBasis`.
- [ ] Replace `MySphericalHarmonics` references with
      `TesseralHarmonics` across `src/`.
- [ ] Replace `Optimize` module references with `Fitting` across
      `src/`.
- [ ] Update `include` / `using .X` lines in `src/Magesty.jl` to the
      new names.
- [ ] Update `using` statements in tests and rename
      `test_MySphericalHarmonics.jl` -> `test_TesseralHarmonics.jl`.
- [ ] Sync `include` / `using` in `tools/` (grep-verify).
- [ ] Update module names in the `@docs` blocks of `docs/src/api.md`.
- [ ] Confirm `make test-all` / `make test-tools` / `make test-aqua` /
      `make test-jet` / `make build` are green.
- [ ] Commit: `refactor(src): rename internal modules for clarity`.

## Milestone 3: documentation sync

- [ ] Rewrite the directory diagram and the "Main modules" /
      "Utilities" tables in `SPEC.md` with the new layout and names.
- [ ] Replace `MySphericalHarmonics` with `TesseralHarmonics` in the
      "Spherical-harmonics convention" entry of `CLAUDE.md`.
- [ ] Minimally sync old-name references in `docs/design-notes/`.
- [ ] Confirm `git grep` returns zero for:
      - `MySphericalHarmonics`
      - `module Basis\b`
      - `module Optimize\b`
      - `common/version`, `common/SortedCounter`
      - `types/AtomCells`, `types/Basis`
      - `utils/atomsbase_adapter`, `utils/xml_io`,
        `utils/MySphericalHarmonics`,
        `utils/SphericalHarmonicsTransforms`,
        `utils/AngularMomentumCoupling`, `utils/ConfigParser`,
        `utils/RotationMatrix`
- [ ] Commit: `docs: sync SPEC.md and CLAUDE.md with src layout refactor`.

## Milestone 4: PR

- [ ] Include "public API unchanged" and the rename-mapping table in
      the PR description.
- [ ] Brief reviewers on the publication-readiness context.

## Verification commands (run after each milestone)

```bash
make test-all
make test-tools
make test-aqua
make test-jet
make build
```
