# Pre-release cleanup (v0.1.0)

**Status**: not started (2026-05-17)

Repository-wide audit results and fixes for the first official release
v0.1.0. The plan was to land in one PR: broken examples paths (P0),
`Project.toml` version / authors / compat tidy-up, add `CITATION.cff`,
tighten Documenter checks, normalize `predict_*` docstrings, and remove
the unimplemented `square_lattice` test.

## Background

- `Project.toml` was at `version = "0.1.0-DEV"` and 60+ commits past
  the most recent tag `v0.0.2`.
- A three-agent sweep audited the package's API, documentation, tests,
  and CI for release-blocking defects and hygiene issues.
- User decision: bundle P0 fixes, metadata cleanup, doc / test
  hygiene, version bump (0.1.0), and minor-level compat relaxation in
  one PR.

## Issues found and fixes

### P0. Broken examples paths (must fix)

`FIXTURE` pointed at `test/examples/`, but the actual location is
`test/integration/`.

- `examples/01_basic_flow.jl:10` —
  `"..", "test", "examples", "fept_tetragonal_2x2x2"` ->
  `"..", "test", "integration", "fept_tetragonal_2x2x2"`.
- `examples/03_save_load.jl:12` — same fix.

`test/integration/fept_tetragonal_2x2x2/{input.toml,EMBSET}` are
present. `examples/02_cif_input.jl` has no issue.

### P1. `Project.toml` cleanup

- **version** (line 3): `"0.1.0-DEV"` -> `"0.1.0"`.
- **authors** (line 4): keep
  `"T.Tanaka <tomonori.tanaka.academic@gmail.com>"`.
- **LICENSE line 375**: `Copyright (c) 2024 TomonoriTanaka` -> normalize
  to `T.Tanaka`.
- **Relax compat patch pins to minor**:
  - `AtomsBase = "0.5.2"` -> `"0.5"`
  - `Combinat = "0.1.4"` -> `"0.1"`
  - `EzXML = "1.2.1"` -> `"1.2"`
  - `LegendrePolynomials = "0.4.5"` -> `"0.4"`
  - `MultivariateStats = "0.10.3"` -> `"0.10"`
  - `OffsetArrays = "1.17.0"` -> `"1"`
  - `StaticArrays = "1.9.10"` -> `"1.9"`
  - `Unitful = "1.28.0"` -> `"1"`
  - `WignerD = "0.1.4"` -> `"0.1"`
  - Keep `JET = "0.11"` (CI test pin).

`test/component/test_Version.jl` only checks the string match between
`Project.toml`'s version and `Magesty.VERSION`, so the bump does not
ripple into other tests.

### P2. Add `CITATION.cff` (new)

Make GitHub's "Cite this repository" work. CFF 1.2.0. Source the paper
info from the README (arXiv:2512.04458, ORCID).

```yaml
cff-version: 1.2.0
message: "If you use Magesty.jl, please cite the following work."
title: "Magesty.jl"
authors:
  - family-names: Tanaka
    given-names: Tomonori
    orcid: "https://orcid.org/0000-0001-7306-6770"
version: 0.1.0
date-released: 2026-05-17
license: MPL-2.0
repository-code: "https://github.com/Tomonori-Tanaka/Magesty.jl"
preferred-citation:
  type: article
  authors:
    - family-names: Tanaka
      given-names: Tomonori
  title: "General spin models from noncollinear spin DFT via spin-cluster expansion"
  year: 2025
  url: "https://arxiv.org/abs/2512.04458"
```

Confirm paper title / authors / arXiv ID at implementation time.

### P3. Tighten Documenter checks

- `docs/make.jl:38` — `checkdocs = :none` -> `checkdocs = :exports`.
- After the tighten, CI verifies that every exported symbol has an
  `@docs` block in `docs/src/api.md` etc. The prior survey suggests
  full coverage, but re-confirm the build passes.

### P4. Normalize `predict_energy` / `predict_torque` docstrings

- `src/Magesty.jl:529–540` (predict_energy)
- `src/Magesty.jl:559–571` (predict_torque)
- Currently description-only. Add `# Arguments` / `# Returns` per the
  CLAUDE.md rule, matching the style of `r2_energy` / `rmse_torque` in
  the same file.

### P5. Moved to `docs/specs/260519-pre-release-safety/`

The `square_lattice` `SCEModel(fit)` round-trip task — together with
the batched `predict_*` overload coverage and the basis-identity
check it motivated — was folded into M4 of the
`260519-pre-release-safety` spec and shipped from there. This note
keeps single ownership of broken examples, Documenter strictness,
docstring normalization, `CITATION.cff`, and `Project.toml` tidy-up.

## Verification

1. `julia --project=examples examples/01_basic_flow.jl` and
   `examples/03_save_load.jl` both run to completion.
2. `make test-all` — unit + integration pass (square_lattice / Version
   / Aqua included).
3. `make test-jet` / `make test-aqua` — zero warnings.
4. `julia --project=docs docs/make.jl` — builds with `:exports`, no
   undocumented exports.
5. After push, manually confirm that GitHub's "Cite this repository"
   button recognizes `CITATION.cff`.

## Not flagged (deliberately out of scope)

Pre-audit confirmed clean:

- Japanese comments / docstrings in `.jl` source.
- References to `CLAUDE.md` / `docs/specs/` / `.claude/`.
- "spin cluster expansion" hyphenation (correct at src/Magesty.jl:5).
- US English violations (only `normalisation=:L2` for SpheriCart, kept
  intentionally).
- TODO / FIXME / XXX / HACK comments (only the square_lattice one,
  handled by P5).
- Stray debug `@show` / `println`.
- Stale `@test_broken` / `@test_skip`.
- Commented-out blocks.
- `.DS_Store` in the index.
- `examples/02_cif_input.jl` API consistency.
- Math conventions in `docs/src/technical_notes.md` vs the source.
- Hot-path application of type annotations / `@views` / `@inbounds`.

## Workflow

- Agree on this note, then decide whether to lift to a spec or go
  directly to a PR.
- Split commits per concern (P0 / P1 / P2 / P3 + P4 / P5 / optional
  P6).
- `git commit` / `git push` only after the user gives an explicit
  instruction, via the `git-helper` agent.
- The `v0.1.0` release tag is created by the user after this PR is
  merged.
