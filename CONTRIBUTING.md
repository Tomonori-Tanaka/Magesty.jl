# Contributing to Magesty.jl

Thanks for your interest in contributing. Magesty.jl is a Julia package for
constructing effective spin models from noncollinear spin DFT calculations
via the spin-cluster expansion (SCE). Numerical correctness, reproducibility,
and physical consistency take precedence over stylistic refactoring.

## Reporting issues

- **Bugs**: open a GitHub issue using the *Bug report* template. Please include
  a minimal reproducer (small input file, expected vs. observed numbers).
- **Feature requests**: use the *Feature request* template. Describe the
  scientific use case so we can judge scope.
- **Security issues**: see [SECURITY.md](SECURITY.md). Do not file public
  issues for vulnerabilities.

## Development workflow

1. Fork the repository and create a topic branch from `main`. Branch naming:
   - `fix/<slug>` for bug fixes
   - `feat/<slug>` for new features
   - `refactor/<slug>` / `chore/<slug>` / `docs/<slug>` as appropriate
2. Make changes. For non-trivial work (multiple design decisions, days of
   effort, or behavioral changes) we use spec folders under `docs/specs/`.
   See [docs/specs/](docs/specs/) for examples; a template is at
   [docs/specs/_template/](docs/specs/_template/).
3. Add or update tests. Numerical changes must come with a regression or
   validation test.
4. Run the full local check before opening a PR:
   ```bash
   make test-all      # unit + integration
   make test-aqua     # package hygiene
   make test-jet      # static type analysis
   ```
   `make ci-local` runs the same matrix CI runs on Julia `release`.
5. Update documentation as needed:
   - User-facing changes → `docs/src/`
   - New public API → `SPEC.md` and `docs/src/api.md`
   - Examples (`examples/`) that exercise the changed code path

## Commit messages

We follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <subject>
```

- Types: `feat` / `fix` / `docs` / `test` / `refactor` / `perf` / `chore` / `style`
- Subject: imperative mood, lowercase, no trailing period
  (e.g., `add SCE basis loader`)
- Breaking changes: include `BREAKING CHANGE: ...` in the commit body

## Style

- US English in source, comments, docstrings, commit messages, and PRs.
  Japanese is fine for conversation and `docs/design-notes/`.
- See [STYLE_GUIDE.md](STYLE_GUIDE.md) for naming, loop conventions, argument
  order, and type annotation rules.
- Hot-path performance guidance (StaticArrays, threading, bounds-check
  suppression) lives in [CLAUDE.md](CLAUDE.md) under "Performance guidelines".

## Physics conventions

These are easy to break silently:

- Spin directions are unit vectors; the layout is `3 × n_atoms`.
- Real (tesseral) spherical harmonics `Zₗₘ`, not complex `Yₗₘ`.
- SALC / Clebsch-Gordan conventions are defined in
  [the technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/);
  consult them before changing any normalization or sign.
- Energy unit follows the DFT input (typically eV). `j0` (reference energy)
  is stored separately from the SCE coefficients `Jφ`.

## Pull requests

- One logical change per PR. Small follow-ups are preferred over giant PRs.
- Fill in the PR template; check the boxes that apply.
- CI must pass: tests, Aqua, JET, docs build.

## License

By contributing, you agree that your contributions will be licensed under the
[Mozilla Public License 2.0](LICENSE).
