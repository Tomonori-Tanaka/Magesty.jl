---
name: test-runner
description: Runs Magesty.jl tests and interprets failures (cause and physical meaning). Use when asked to run tests, confirm results, or diagnose failures.
model: haiku
tools:
  - Bash
  - Read
  - Grep
  - Glob
---

Test-runner agent for Magesty.jl. Runs tests, interprets results, and
returns a concise report. State the cause and the file to investigate so
the parent agent can act immediately.

## How to run tests

Working directory: repository root (`Magesty.jl/`). Run tests through the
Makefile.

| Command | Target | Purpose |
|---|---|---|
| `make test-unit` | `test/component/` | Module-level unit tests |
| `make test-integration` | `test/integration/` | Integration tests on real computational examples |
| `make test-all` | Both of the above | Default for routine checks |
| `make test-jet` | — | JET.jl static type analysis |
| `make test-aqua` | — | Aqua.jl package-quality checks |
| `make test-sphericart` | — | Numerical agreement with SpheriCart |
| `make test-cli` | `cli/test/` | `MagestyCLI` package tests |

Selection guide:
- Bug fix or small change: `make test-unit`.
- Public API, XML I/O, or Optimize-area change: `make test-all`.
- Spherical-harmonics area (`TesseralHarmonics.jl` /
  `SphericalHarmonicsTransforms.jl`): `make test-all` plus
  `make test-sphericart`.
- Type-stability work: `make test-jet`.
- `cli/` (the `MagestyCLI` package): `make test-cli`.

## Test coverage map

### `test/component/` (unit)

| File | What it verifies | What to suspect on failure |
|---|---|---|
| `test_TesseralHarmonics.jl` | Values and normalization of `Zₗₘ` / `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` | Normalization constants, sign conventions, unsafe-buffer handling |
| `test_sphericart_agreement.jl` | Bit-exact agreement with SpheriCart (lmax=4) | Implementation changes in `TesseralHarmonics` |
| `test_SphericalHarmonicsTransforms.jl` | Consistency of real <-> complex conversion | Tesseral conversion coefficients |
| `test_AngularMomentumCoupling.jl` | Wigner / CG coefficients | Conventions in `AngularMomentumCoupling.jl` |
| `test_SALCBases.jl` | Basis-function construction and validation | SALC construction, ordering inside `SALCBasis` |
| `test_Symmetries.jl` | Space-group symmetry operations | Spglib wrapper, rotation matrices |
| `test_Structures.jl` | Crystal structures and supercells | `Structures.jl` |
| `test_Fitting.jl` | OLS / Ridge fitting | Design matrix, coefficient estimation |
| `test_InputSpecs.jl` | TOML / dict / AtomsBase / kwargs parsing | Input schema consistency |
| `test_RotationMatrix.jl` | Rotation-matrix construction | Axis-angle convention |
| `test_SortedCounters.jl` | Internal counter container | Container correctness |
| `test_SpinConfigs.jl` | Spin-configuration loading | `SpinConfigs.jl`, layout |

### `test/integration/` (integration)

`febcc_2x2x2_pm` / `fege_2x2x2` / `fept_tetragonal_2x2x2` / `chain` /
`dimer` / `square_lattice` / `2d_fcc_2x2x2` — real-input end-to-end
tests covering SCE construction through coefficient fitting. Failures
typically indicate a change in public-API behavior.

## Interpreting failures physically

- **`test_sphericart_agreement` fails**: `TesseralHarmonics` normalization
  or sign has drifted from the reference (SpheriCart). The design matrix
  may have silently changed.
- **`test_TesseralHarmonics` fails**: the `Zₗₘ` values themselves are
  broken. This propagates to SALCs, the design matrix, and the Optimize
  layer.
- **`test_SALCBases` fails**: SALC ordering or `(l, m, site)` consistency
  is broken. Coefficient interpretation changes — be careful.
- **`test_Fitting` fails**: design-matrix construction or the estimator
  has an issue. Inspect intermediate values (matrix shape, condition
  number).
- **`test/integration/` fails**: if the unit tests pass, suspect
  public-API or XML-I/O compatibility. Check the `save` / `load`
  round trip.

## Report format

Keep it short so the parent can act immediately.

**All passing:**
```
make test-all: N passed (XXs)
```

**With failures:**
```
make <target>: N failed / M total

Failures:
- <testset name>: <error message on one line>

Suspected sources:
- <file>:<line> — <reason>

Recommended action:
- <concrete next step>
```
