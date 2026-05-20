# Requirements: `write_energies` / `write_torques` for the current API

Status: draft (2026-05-20)

## Goal

Reintroduce `write_energies` and `write_torques` — plain-text writers that
dump observed (DFT) vs. predicted (SCE) energies and torques — adapted to the
current four-type API (`SCEBasis` / `SCEDataset` / `SCEFit` / `SCEModel`), so
the existing `tools/FitCheck_energy.py` and `tools/FitCheck_torque.py`
visualization scripts can be driven without manual file assembly.

## Background

The legacy `write_energies` / `write_torques` helpers were removed together
with the old `System` / `SpinCluster` API. They consumed a monolithic
`SpinCluster` and read pre-computed predictions out of its optimizer.

The current API exposes `predict_energy`, `predict_torque`, `read_embset`, and
`Structure`, which together cover everything the old writers needed, but no
convenience writer exists. Users currently have to assemble the
`FitCheck_*.py` input files by hand. This spec restores the two writers as a
thin I/O layer over the existing prediction verbs.

## Scope

Includes:

- New file `src/FitCheckIO.jl` defining `write_energies` and `write_torques`.
- Two signature families per writer:
  - `write_*(f::SCEFit, filename)` — self-contained, evaluates the fit's own
    training dataset.
  - `write_*(predictor, data, filename)` — `predictor` is `SCEModel` or
    `SCEFit`; `data` is `SCEDataset` / `AbstractVector{SpinConfig}` /
    `embset_path::AbstractString` (validation / test sets).
- Exporting both functions from `Magesty`.
- Output file format identical to the legacy format (so `FitCheck_*.py` runs
  unchanged).
- Regression test and `docs/` updates.

Excludes:

- Any change to `tools/FitCheck_energy.py` / `tools/FitCheck_torque.py`.
- A combined `write_fitcheck` wrapper (the two Python scripts are separate).
- Changes to `predict_energy` / `predict_torque` themselves.

## Invariants

- No numerical-convention change: predictions flow through the existing
  `predict_energy` / `predict_torque` verbs unchanged. This is pure I/O.
- Legacy output format reproduced verbatim:
  - Energy file: header + 3 columns `index  DFT_Energy  SCE_Energy`, eV.
  - Torque file: header + per-config `# data index: N` blocks + 8 columns
    `atom_index  element  DFT_xyz  SCE_xyz`, eV.
- Spin direction stays unit-vector with `3 × n_atoms` layout.
- Existing XML round-trips and the public API of the four core types are
  untouched (only additive: two new exported names).

## Completion criteria

- [ ] `write_energies` / `write_torques` implemented and exported.
- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] New regression test `test/component/test_fitcheck_io.jl` registered and
      passing; covers format, signature families, and agreement between them.
- [ ] `FitCheck_energy.py` / `FitCheck_torque.py` parse the generated files
      end-to-end.
- [ ] New API reflected in `docs/src/api.md` (and `tools.md` where the
      `FitCheck_*.py` scripts are documented), `SPEC.md`, and `CHANGELOG.md`
      `[Unreleased]`.

## References

- Related issues / PRs: —
- Related specs / design notes:
  [260514-sce-public-api](../260514-sce-public-api/) (the four-type API these
  writers target).
