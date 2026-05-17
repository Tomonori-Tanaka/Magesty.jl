# Julia compat policy

**Status**: operating rule (2026-05-17)

Operating notes for the `julia = "1.12"` lower bound in `Project.toml`
and the CI matrix `'1.12'` + `'1'` (latest stable).

## Current state

- `Project.toml`: `julia = "1.12"` (bumped 1.10 -> 1.12 in PR #12).
- `.github/workflows/CI.yml`: matrix over `'1.12'` (minimum) and `'1'`
  (latest stable), across ubuntu and macos.
- `Makefile`'s `ci-local`: uses the `juliaup` `release` channel
  (1.12 family at the time of writing).

## When Julia 1.13 stable ships

Decision axes:

1. **Standard library and dependency churn**: does AtomsBase / Spglib /
   StaticArrays / WignerD / ... need any breaking adjustment for 1.13?
2. **User base**: distribution of Julia versions used (HPC often holds
   on to the 1.12 LTS-style line for a long time).
3. **CI cost**: stretch the matrix axis to three (`1.12` / `1.13` / `1`)
   or keep two (`1.12` and `1`).

**Default policy**: once 1.13 stable lands,

- `1.13` is automatically picked up by the `'1'` matrix entry.
- Keep the minimum (`'1.12'`) for at least 6 months. After that, open a
  PR bumping the lower bound to 1.13 and let CI confirm no breakage.
- When changing the `julia` compat in `Project.toml`, run
  `make ci-local` and add a `BREAKING CHANGE` line in CHANGELOG.

## 1.14 and beyond

- Aim for "current stable minus 2 minor versions" as the minimum, with
  case-by-case extensions for HPC sites needing longer support.

## References

- PR #12: raise minimum Julia to 1.12
- `Makefile` `ci-local` target — local reproduction of the CI matrix
