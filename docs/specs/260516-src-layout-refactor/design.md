# Design: `src/` layout refactor

## 1. Rename mapping (final)

### 1.1 File-name alignment (match the declared module)

| Old path | New path | Module |
|---|---|---|
| `src/common/version.jl` | `src/Version.jl` | `Version` |
| `src/common/SortedCounter.jl` | `src/SortedCounters.jl` | `SortedCounters` |
| `src/utils/atomsbase_adapter.jl` | `src/AtomsBaseAdapter.jl` | `AtomsBaseAdapter` |
| `src/utils/xml_io.jl` | `src/XMLIO.jl` | `XMLIO` |

### 1.2 Module renames

| Old name | New name | Old file | New file |
|---|---|---|---|
| `MySphericalHarmonics` | **`TesseralHarmonics`** | `src/utils/MySphericalHarmonics.jl` | `src/TesseralHarmonics.jl` |
| `Basis` | **`CoupledBases`** | `src/types/Basis.jl` | `src/CoupledBases.jl` |
| `Optimize` | **`Fitting`** | `src/Optimize.jl` | `src/Fitting.jl` |

`TesseralHarmonics` makes the real (tesseral) form explicit with the
proper academic term. `CoupledBases` matches the naming style of
`SALCBases` and reflects the content (couplings built from CG
coefficients). `Fitting` aligns with the StatsAPI `fit` verb.

### 1.3 Unchanged

Main modules: `Structures`, `Symmetries`, `Clusters`, `SALCBases`,
`SpinConfigs`, `AngularMomentumCoupling`, `RotationMatrix`,
`ConfigParser`, `SphericalHarmonicsTransforms`, `AtomCells`.

## 2. Directory layout (final: flat under `src/`)

```
src/
  Magesty.jl                       # main module
  AngularMomentumCoupling.jl
  AtomCells.jl                     # was types/
  AtomsBaseAdapter.jl              # was utils/atomsbase_adapter.jl
  Clusters.jl
  ConfigParser.jl                  # was utils/
  CoupledBases.jl                  # was types/Basis.jl
  Fitting.jl                       # was Optimize.jl
  RotationMatrix.jl                # was utils/
  SALCBases.jl
  SortedCounters.jl                # was common/SortedCounter.jl
  SphericalHarmonicsTransforms.jl  # was utils/
  SpinConfigs.jl
  Structures.jl
  Symmetries.jl
  TesseralHarmonics.jl             # was utils/MySphericalHarmonics.jl
  Version.jl                       # was common/version.jl
  XMLIO.jl                         # was utils/xml_io.jl
```

17 files, one level. `common/`, `types/`, `utils/` are removed.

Rationale: in a survey of comparable public packages, StatsAPI /
Spglib / MultivariateStats / AtomsBase all run with a flat or very
shallow split.

## 3. Impact (reference sites already mapped)

### 3.1 `src/Magesty.jl`

- `include("common/version.jl")` -> `include("Version.jl")`
- `include("common/SortedCounter.jl")` -> `include("SortedCounters.jl")`
- `include("types/AtomCells.jl")` -> `include("AtomCells.jl")`
- `include("types/Basis.jl")` -> `include("CoupledBases.jl")`
- `include("utils/SphericalHarmonicsTransforms.jl")` ->
  `include("SphericalHarmonicsTransforms.jl")`
- `include("utils/AngularMomentumCoupling.jl")` ->
  `include("AngularMomentumCoupling.jl")`
- `include("utils/ConfigParser.jl")` -> `include("ConfigParser.jl")`
- `include("utils/atomsbase_adapter.jl")` ->
  `include("AtomsBaseAdapter.jl")`
- `include("utils/RotationMatrix.jl")` -> `include("RotationMatrix.jl")`
- `include("utils/MySphericalHarmonics.jl")` ->
  `include("TesseralHarmonics.jl")`
- `include("utils/xml_io.jl")` -> `include("XMLIO.jl")`
- `include("Optimize.jl")` -> `include("Fitting.jl")`
- Update the corresponding `using .X` group:
  `using .Basis` -> `using .CoupledBases`, etc.

### 3.2 Cross-references inside `src/`

- `Clusters.jl:28`: `using ..SortedCounters: SortedCounter` — works
  unchanged (the module name was already correct).
- `SALCBases.jl:15`: `using ..SortedCounters: SortedCounter` — same.
- `SALCBases.jl` (many spots): `Basis.CoupledBasis` ->
  `CoupledBases.CoupledBasis`.
- `SALCBases.jl`, `Fitting.jl` (was `Optimize.jl`), and others:
  `MySphericalHarmonics` references -> `TesseralHarmonics`.
- `XMLIO.jl`: `using ..Version` — unchanged.

### 3.3 `test/`

Likely-referencing files (grep-confirmed):

- `test/runtests.jl`
- `test/component_test/test_MySphericalHarmonics.jl` -> rename to
  `test_TesseralHarmonics.jl`.
- `test/component_test/test_SortedCounter.jl` — module name is already
  `SortedCounters`; check the contents.
- `test/component_test/test_sphericart_agreement.jl`
- `test/examples/{chain,dimer}/test.jl`
- `bench/benchmark_sphericart.jl`,
  `bench/benchmark_spherical_harmonics.jl`

### 3.4 `tools/`

Apply only the minimum `include` / `using` sync forced by renames or
moves under `src/`. No structural change or internal refactor (out of
scope). When starting, confirm with:

```
grep -rln "MySphericalHarmonics\|atomsbase_adapter\|xml_io\|Optimize\|types/Basis" tools/
```

### 3.5 `docs/` / `SPEC.md` / `CLAUDE.md`

- `@docs` blocks in `docs/src/api.md`: `MySphericalHarmonics.Zₗₘ` etc.
  -> `TesseralHarmonics.Zₗₘ` etc.
- The directory diagram and main-module table in `SPEC.md`.
- The "Spherical-harmonics convention" entry in the physics-conventions
  section of `CLAUDE.md` (existing name references).
- References in `docs/design-notes/` to the old names (historical
  documents — keep edits minimal).

## 4. Commit granularity (implementation order)

After each commit, confirm that `make test-all` + `make test-tools` +
`make build` are green. Each commit is independently revertible.

### Commit 1 — File-name alignment + directory flattening

- 14 file moves (`git mv`):
  - `common/version.jl` -> `Version.jl`
  - `common/SortedCounter.jl` -> `SortedCounters.jl`
  - `types/AtomCells.jl` -> `AtomCells.jl`
  - `types/Basis.jl` -> `Basis.jl` (no module rename yet)
  - `utils/SphericalHarmonicsTransforms.jl` ->
    `SphericalHarmonicsTransforms.jl`
  - `utils/AngularMomentumCoupling.jl` -> `AngularMomentumCoupling.jl`
  - `utils/ConfigParser.jl` -> `ConfigParser.jl`
  - `utils/atomsbase_adapter.jl` -> `AtomsBaseAdapter.jl`
  - `utils/RotationMatrix.jl` -> `RotationMatrix.jl`
  - `utils/MySphericalHarmonics.jl` -> `MySphericalHarmonics.jl`
    (no module rename yet)
  - `utils/xml_io.jl` -> `XMLIO.jl`
  - `Optimize.jl` (top-level, untouched)
- Update `include` paths in `src/Magesty.jl` (do not touch module
  names yet).
- Remove the now-empty `common/`, `types/`, `utils/` directories.
- Module names are unchanged in this commit, so the diff is minimal
  and easier to follow.

### Commit 2 — Module renames

- `MySphericalHarmonics.jl` -> `TesseralHarmonics.jl`;
  `module MySphericalHarmonics` -> `module TesseralHarmonics`.
- `Basis.jl` -> `CoupledBases.jl`; `module Basis` ->
  `module CoupledBases`.
- `Optimize.jl` -> `Fitting.jl`; `module Optimize` -> `module Fitting`.
- Rewrite all `using ..Basis` -> `using ..CoupledBases`, etc. across
  `src/`.
- Rewrite dotted references: `Basis.CoupledBasis` ->
  `CoupledBases.CoupledBasis`.
- Rewrite references under `test/` and rename
  `test_MySphericalHarmonics.jl` -> `test_TesseralHarmonics.jl`.
- Rewrite references under `tools/` (minimal).
- Update `@docs` blocks in `docs/src/api.md`.
- Update docstring references inside the main `Magesty.jl` module.
- Confirm `make test-all` and `make build` are green.

### Commit 3 — Sync `SPEC.md` / `CLAUDE.md` / `docs/design-notes/`

- Replace the directory diagram in `SPEC.md` with the new layout.
- In the "Spherical-harmonics convention" entry of the physics-
  conventions section in `CLAUDE.md`, replace `MySphericalHarmonics`
  with `TesseralHarmonics`.
- Minimally sync name references in `docs/design-notes/` (these are
  historical notes — only fix broken references).
- Confirm zero matches with:
  ```
  git grep -E "MySphericalHarmonics|module Basis\b|module Optimize\b|atomsbase_adapter|xml_io.jl|common/version|types/Basis|utils/MySphericalHarmonics"
  ```

## 5. Risks and mitigations

- **Missed references**: after each commit, `git grep` for the old
  names and old paths should return zero. `make test-aqua` also
  doubles as a module-structure sanity check.
- **`@docs` resolution failure**: Documenter resolves the module path
  in `@docs` blocks. Update `docs/src/api.md` in Commit 2 and confirm
  `make build` runs without warnings.
- **`SortedCounters` plural-vs-singular drift**: keep `module
  SortedCounters` (plural) and the inner type `SortedCounter`
  (singular), following the Julia convention used by `Strings` /
  `Tests` / `LinearAlgebra`.
- **`git mv` history**: `git mv` preserves history via
  `git log --follow`. `git diff -M` detects renames.
- **JET / Aqua false positives**: after module renames, Aqua may
  flag unused exports etc. Run `make test-aqua` after each commit.

## 6. Verification checklist

After each commit, run in order:

1. `make test-unit`
2. `make test-integration`
3. `make test-tools`
4. `make test-aqua`
5. `make test-jet`
6. `make build` (docs)
7. `git grep` for old names and old paths returns zero:
   - `git grep -E "MySphericalHarmonics"` (after Commit 2)
   - `git grep -E "module Basis\b"` (after Commit 2)
   - `git grep -E "module Optimize\b"` (after Commit 2)
   - `git grep -E "common/version|common/SortedCounter|types/AtomCells|types/Basis|utils/"`
     (after Commit 1, re-check after Commit 2)
