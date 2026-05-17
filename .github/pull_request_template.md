<!-- Thanks for contributing to Magesty.jl. -->

## Summary

<!-- 1-3 sentences: what does this PR do and why. -->

## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would change existing behavior)
- [ ] Docs / tests / chore only

## Numerical correctness

If this PR changes any numerical result:

- [ ] I have explained *why* the result changes.
- [ ] I have added a regression or validation test.
- [ ] I have updated `docs/` / `examples/` if user-facing behavior changes.

If this PR touches physics conventions (signs, units, normalization, SALC,
Clebsch-Gordan, spherical harmonic conventions):

- [ ] I have consulted the
      [technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/).
- [ ] I have updated all linked sites listed under "Linked sites" in
      [CLAUDE.md](../CLAUDE.md).

## Checks

- [ ] `make test-all` passes locally.
- [ ] `make test-aqua` / `make test-jet` pass (or pre-existing warnings only).
- [ ] Commit messages follow
      [Conventional Commits](https://www.conventionalcommits.org/).
- [ ] Public API changes are reflected in `SPEC.md` and `docs/src/api.md`.

## Related issues

<!-- Closes #123, refs #456. -->
