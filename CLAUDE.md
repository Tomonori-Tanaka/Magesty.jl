# CLAUDE.md

## Project goal

Julia package that builds general effective spin models from noncollinear
spin DFT calculations using the spin-cluster expansion (SCE). We prioritize
**numerical correctness, reproducibility, and physical consistency** over
stylistic refactoring.

## Core rules

- Never silently change numerical conventions (signs, units, normalization).
- Before editing an algorithm, confirm the relevant equations and the
  current sign / unit conventions.
- Any change that may alter numerical results must come with:
  1. A short explanation of why the result changes.
  2. A regression or validation test.
  3. Updates to `docs/` / `examples/` if user-facing.
- Always confirm with the user before `git commit` or `git push`. Do not
  create local commits without an explicit instruction.

## Implementation rules

- Avoid hidden global state.
- Exported APIs must have explicit type annotations and docstrings.
- Use type annotations liberally (`::Float64`, `::Vector{Float64}`,
  `::SVector{3, Float64}`, etc.).
- Public-API docstrings follow the standard Julia format
  (`# Arguments` / `# Returns` / `# Examples`).
- Record performance changes (before / after) in `.claude/bench_log.md`
  (or `DESIGN_NOTES.md` for higher-level summaries).

For detailed coding style (naming, loop conventions, argument order, etc.),
see `STYLE_GUIDE.md`. **Always consult it when editing code.**

## Language and terminology

- **Conversation with the user**: Japanese is fine.
- **Everything committed to the repository**: English only. This covers
  `.jl` source, comments, docstrings, all Markdown (CLAUDE.md, SPEC.md,
  STYLE_GUIDE.md, README, `docs/**`, `.claude/agents/*.md`, etc.), shell
  scripts, TOML, commit messages, PR titles and descriptions, and issue
  templates.
- **Commit messages follow [Conventional Commits](https://www.conventionalcommits.org/).**
  Format: `<type>(<scope>): <subject>` (scope optional).
  - Types in use: `feat` / `fix` / `docs` / `test` / `refactor` / `perf` /
    `chore` / `style`.
  - Subject is imperative, lowercase, no trailing period
    (e.g., `add SCE basis loader`).
  - Breaking changes: include `BREAKING CHANGE: ...` in the body.
- **US English** throughout (`normalize` / `behavior` / `color` / `center`).
- Exception: **preserve external API spellings literally**
  (e.g., SpheriCart's `normalisation=:L2` keyword argument).
- **Japanese in any committed file is auto-blocked by the PostToolUse hook
  (`.claude/hooks/no-japanese.sh`).** The hook covers the entire repository,
  with two exemptions: `tools/personal/` (personal scripts) and
  `.claude/bench_log.md` (historical record).
- **Do not reference Claude-internal scaffolding from source code.** In `.jl`
  comments and docstrings, never name `CLAUDE.md` / `DESIGN_NOTES.md` /
  `docs/design-notes/` / `.claude/` / `docs/specs/`. These are
  Claude-collaboration scaffolding; published source must read independently.
  Summarize the relevant context inline when needed.
  - References allowed in source: `docs/src/` (Documenter), `SPEC.md`,
    `STYLE_GUIDE.md`, and
    `https://Tomonori-Tanaka.github.io/Magesty.jl/` (technical notes).
  - Exception: commit-message `Refs:` lines may cite
    `docs/design-notes/<slug>.md` style paths (commits are workflow
    artifacts).

## Tests

Always run tests via the Makefile after edits.

| Command | Target | Purpose |
|---|---|---|
| `make test-unit` | `test/component/` | Module-level unit tests |
| `make test-integration` | `test/integration/` | End-to-end integration tests |
| `make test-all` | both of the above | Default for routine checks |
| `make test-tools` | `tools/test/` | Tests for `tools/` scripts |
| `make test-jet` | — | JET.jl static type analysis |
| `make test-aqua` | — | Aqua.jl package-quality checks |
| `make test-coverage` | `src/` | Per-file coverage report (ascending). One-time `make coverage-setup` first |

See the Makefile for benchmark targets (`make bench-sphericart`,
`make bench-salcbasis`, etc.).

## Physics conventions

Easy to break silently. Confirm before touching the algorithm.

- **Spin directions are unit vectors.** `spin_directions[:, atom]` is the
  *direction* of the classical spin (`‖·‖ = 1`); the magnetic-moment
  magnitude is a separate concept.
- **Spin-matrix layout**: `3 × n_atoms` (rows = x, y, z; columns = atoms).
  Transposing breaks the entire pipeline.
- **Real (tesseral) spherical harmonics**. Internal representation is real
  `Zₗₘ`, not complex `Yₗₘ`. The number of `(l, m)` pairs is `(l_max+1)²`.
  - Hot paths use `Zₗₘ_unsafe` (buffered, no bounds check). Semantically
    identical, but faster.
- **SALC and Clebsch-Gordan conventions are defined in this repository.**
  Consult the
  [Magesty.jl technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/)
  before changing anything.
- **Energy units**: SCE coefficients `Jφ` carry the energy unit of the DFT
  input (typically eV). `j0` (reference energy) is stored separately from
  the SCE coefficients.

## Performance guidelines

Hot paths are the basis-evaluation loop in `Fitting.jl` and the SALC
construction in `SALCBases.jl`.

**StaticArrays usage**:
- Use `SVector{3, Float64}` (immutable) or `MVector{3, Float64}` (scratch
  buffer) for 3-component vectors.
- Do not allocate new `Vector`s inside loops. Pre-allocate `MVector`s and
  reuse.
- For column slices like `spin_directions[:, atom]`, use `@views` to avoid
  copies and convert to `SVector` for stack-resident processing.

```julia
# Good
dir_svec = SVector{3, Float64}(spin_directions[:, atom])
buf = MVector{3, Float64}(0.0, 0.0, 0.0)

# Bad (heap allocation)
dir_vec = spin_directions[:, atom]
buf = zeros(3)
```

**Threading**:
- The spin-configuration loop (`num_spinconfigs`) is the primary `@threads`
  target.
- The symmetry-operation loop (`isym in 1:n_operations`) is also a
  candidate for `@threads`.

**Bounds-check suppression**:
- Use `@inbounds` on inner loops whose indices are provably correct.
- For spherical harmonics in the hot path, prefer the unsafe variant
  (`Zₗₘ_unsafe`).

**Expensive operations**:
- SALC construction in `SALCBasis` is heavy. To reuse, persist with
  `Magesty.save(basis, path)` and reload with
  `Magesty.load(SCEBasis, path)` (`save` / `load` are intentionally not
  exported, to avoid clashing with the generic names from JLD2, FileIO,
  CSV.jl, etc.).

**Bench bookkeeping**:
- When touching a hot path (`Fitting` / `SALCBases` / `TesseralHarmonics` /
  `Optimize` workspaces), append a before/after entry to
  `.claude/bench_log.md`, even if the change does not affect numerical
  results. The point is systematic regression detection. A short `@btime`
  median is enough. Keep this practice after the release.

## Linked sites (change one, check all)

### Spherical-harmonics convention
The normalization and signs of `Zₗₘ` / `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe` in
`TesseralHarmonics.jl` must stay consistent with
`SphericalHarmonicsTransforms.jl`, the SALC construction in `SALCBases.jl`,
and the SpheriCart agreement test
(`test/component/test_sphericart_agreement.jl`). Changing one without the
others silently breaks the design matrix.

### SCE coefficient I/O
The round trip of `Magesty.save(basis_or_model, path)` /
`Magesty.load(SCEBasis, path)` / `Magesty.load(SCEModel, path)`
(`src/XMLIO.jl`) must agree to the last digit. When changing the on-disk
format, basis ordering, or the field layout of `SCEBasis` / `SCEModel`,
update both sides together. Do not reintroduce the old `write_xml` /
`build_sce_basis_from_xml` names — the new `save` / `load` API is
canonical, and is intentionally non-exported (call via `Magesty.save` /
`Magesty.load`).

### Fitting and SALCBasis
The design-matrix construction in `Fitting.jl` depends on the key-group
order inside `SALCBasis` (the outer index of
`Vector{Vector{CoupledBasis_with_coefficient}}`). Reordering SALCs changes
the physical meaning of `Jφ`. The XML I/O serializes coefficients in the
same order, so all three must stay synchronized.

## Managing development units

Mid-sized or larger work goes into **spec folders**. No cross-sprint
progress trackers.

- **Active development units**:
  [`docs/specs/[YYMMDD]-[slug]/`](docs/specs/), each with
  `requirements.md` / `design.md` / `tasklist.md`. Index at
  [`docs/specs/README.md`](docs/specs/README.md); template at
  [`docs/specs/_template/`](docs/specs/_template/).
- **Cross-cutting design notes, investigations, on-hold ideas**:
  [`DESIGN_NOTES.md`](DESIGN_NOTES.md) (index) with full content under
  [`docs/design-notes/`](docs/design-notes/). Operating rules at
  `docs/design-notes/README.md`.
- **Day-to-day TODOs**: `TaskCreate` (in-session only). Lightweight TODOs
  that should persist belong in `docs/design-notes/` or a spec.
- **Historical benchmark records**: `.claude/bench_log.md` and `git log`.

### Spec-folder workflow

**Always create a spec folder and agree on it before starting mid-sized or
larger work.**

Entry criteria (any of these triggers a spec):
- Multi-day effort.
- Multiple design choices (affecting API, types, or conventions).
- Mid-sized or larger change to existing behavior (new feature,
  medium-scale refactor).
- Future readers will ask "why was this done this way?".

Skip the spec for:
- Bug fixes (covered by a regression test).
- Documentation or comment fixes.
- Small refactor within a single file.
- Minor behavior tweaks already covered by existing tests.

Procedure (Claude executes):
1. Create `docs/specs/[YYMMDD]-[slug]/` (`YYMMDD` = `date +%y%m%d`; `slug`
   is English kebab-case).
2. Copy the three files from
   [`docs/specs/_template/`](docs/specs/_template/) and fill them out with
   the user:
   - `requirements.md` — goal, invariants, scope, completion criteria.
   - `design.md` — modules, API, types, conventions, impact on linked
     sites.
   - `tasklist.md` — milestones (coarse). Day-to-day tasks use
     `TaskCreate`. Includes an exit checklist.
3. Reach agreement on the spec before starting implementation.
4. Keep the folder after completion (do not delete it; it is the historical
   record). Update the `Status:` line in `tasklist.md` and the table in
   [`docs/specs/README.md`](docs/specs/README.md) together.

## Working principles for Claude

### Free to proceed without asking

- Bug fixes (minimal change plus test).
- Adding / fixing tests.
- Documentation typos.
- Adding notes to the `DESIGN_NOTES.md` index or `docs/design-notes/`.

### Sub-agent usage

- After implementing, use the `test-runner` agent to run and diagnose tests
  (`.claude/agents/test-runner.md`).
- Before committing, use the `code-reviewer` agent on the diff
  (`.claude/agents/code-reviewer.md`).
- For performance investigation, use the `profiler` agent
  (`.claude/agents/profiler.md`).
- Once the user gives an explicit commit / push instruction, hand off to
  the `git-helper` agent (`.claude/agents/git-helper.md`). It drafts the
  Conventional Commit message, runs an ASCII-only check, applies the
  commit via `Write` + `git commit -F file` (to avoid heredoc accidents),
  detects `BREAKING CHANGE`s, and adds `Refs:` lines. The main agent must
  not run `git commit -m` directly.

### Propose before implementing

- Algorithm changes (when numerical results may change).
- Refactors that cross module boundaries.
- Performance improvements (present benchmark numbers first).

### Always confirm — do not implement first

- Physics-convention changes (signs, units, normalization, SALC, CG).
- New external dependencies.
- Public-API signature changes (anything in `export`).
- Changes to spherical-harmonics normalization (`Zₗₘ` family).
- Changes to the XML I/O format.
- `git commit` / `git push` of any kind (including local commits) without
  explicit instruction.

## References

Consult as needed before working.

- `STYLE_GUIDE.md` — detailed coding-style rules.
  **Always consult when editing code.**
- `SPEC.md` — architecture, directory layout, primary types, public API.
- `docs/specs/` — active and completed specs (mid-sized or larger units).
- `DESIGN_NOTES.md` — index of design notes, investigations, and on-hold
  ideas. Bodies under `docs/design-notes/`.
- [Magesty.jl technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/)
  — mathematical conventions for SALC, CG coefficients, etc.
