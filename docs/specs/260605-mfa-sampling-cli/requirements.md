# Requirements: MFA spin-sampling CLI (`magesty vasp mfa`)

Status: merged (2026-06-06, PR #34)

## Goal

Promote the standalone `tools/sampling_mfa.jl` script to a first-class
`magesty vasp mfa` CLI command, backed by a reusable, code-agnostic
Mean-Field-Approximation (MFA) spin sampler in core `Magesty`. The design
must extend naturally to non-VASP DFT codes later (`magesty qe mfa`, etc.).

## Background

`tools/sampling_mfa.jl` samples thermally conditioned spin configurations
from a VASP INCAR using the MFA with a von Mises-Fisher (vMF) direction
distribution, sweeping over `tau` (T/Tc) or `m` (magnetization) and writing
many `sample-NN.INCAR` files. It calls no `Magesty` API, ships its own INCAR
text parser (`tools/vasp/vasptools.jl`, distinct from the core `VaspIO.jl`
that parses `vasprun.xml`), and depends on `Distributions` + `Roots`.

The CLI foundation (specs `260520-cli-foundation`,
`260521-cli-package-extraction`) established the pattern: conversion/processing
logic lives in core `Magesty`; the `MagestyCLI` package holds thin Comonicon
wrappers. This spec applies that pattern to MFA sampling. This is the first
"non-converter" CLI tool (sweep producing many output files rather than one
`String`), so the API shape is a new design point.

## Scope

Includes:

- New code-agnostic MFA sampler in core (`src/MfaSampling.jl`): vMF direction
  sampling, the MFA self-consistency solver, the `tau`/`m` sweep driver.
- INCAR text I/O in core (`parse_incar` / `write_incar`), ported from
  `tools/vasp/vasptools.jl`.
- VASP orchestration + exported public API `sample_mfa_incar`.
- CLI leaf `magesty vasp mfa` in `MagestyCLI`.
- Hand-rolled exact p=3 vMF sampler to avoid a `Distributions` dependency;
  add only `Roots` as a new core dependency.
- Tests, docs, CHANGELOG; remove the promoted script.

Excludes:

- Non-VASP DFT-code adapters (future specs); the Layer1 sampler is built to
  be reused by them but no second adapter is written here.
- `tools/vasp/vasptools.jl` removal (still used by `tools/personal/`).
- Other `tools/` scripts (separate specs).

## Invariants

- **Spin layout unchanged**: directions are unit vectors; the spin matrix is
  `3 × n_atoms` (rows x/y/z, columns atoms). Magnitudes preserved per atom.
- **MFA / vMF numerical conventions preserved** relative to the script:
  - MFA equation `m = coth(3m/τ) − τ/3m`; clamps `MIN_TEMP = 1e-5`,
    `MAX_TEMP = 0.99999`.
  - `thermal_averaged_m`: `τ < MIN_TEMP ⇒ returns 1.0` (fully ordered),
    `τ > MAX_TEMP ⇒ returns 0.0` (fully disordered). (This corrects a latent
    sign mistake in the original script's boundary returns, confirmed against
    the Langevin limit `m = L(3m/τ)`; the original branches were never reached
    by the sampler, so sampling output is unchanged.)
  - `mfa_sample` (per-config sampler), concentration `κ = 3m/τ`:
    `τ > MAX_TEMP ⇒ κ = 1e-6` (near-uniform draw), `τ < MIN_TEMP ⇒ returns the
    input spin matrix unchanged`.
  - Inverse map `tau_from_magnetization`: `m → 0 ⇒ τ = 1`, `m → 1 ⇒ τ = 0`;
    magnetizations outside the bracketed range clamp to the `0.0`/`1.0` limit.
  - Zero-norm spins are skipped (left at zero).
  - `--fix` / `--uniform-atoms` / `--randomize` interaction identical to the
    script (fixed atoms get the same rotation when `--randomize`; uniform
    atoms sampled isotropically on the sphere; `M_CONSTR` written equal to
    `MAGMOM`; all other INCAR keys preserved).
- **Existing public API and XML round-trips untouched.**
- The hand-rolled vMF sampler is statistically equivalent to the original
  (validated analytically, not pinned to current output).

## Behavior change vs. the script (intentional)

- Sweep is specified by **count**, not width: `--num-points` (integer) gives
  `range(start, stop, length = num_points)`, replacing `--step` / `-w`.
- `--end` / `-e` is renamed to `--stop` (Julia reserves `end`).
- New `--outdir` (default `.`) chooses the output directory (script wrote to
  the current directory).

## Completion criteria

- [x] `make test-all` passes (new `test_MfaSampling.jl`, `test_IncarIO.jl`).
- [x] `make test-cli` passes (new `magesty vasp mfa` regression test).
- [x] `make test-jet` / `make test-aqua` clean (core stays Distributions-free).
- [x] vMF sampler validated analytically against `A₃(κ) = coth κ − 1/κ`.
- [x] `sample_mfa_incar` exported with a full Julia docstring; reflected in
      `docs/src/api.md` and the CLI narrative docs.
- [x] `tools/sampling_mfa.jl` removed; `CHANGELOG.md` updated.
- [x] Tier-2 review panel run and findings resolved.

## References

- Related specs: `260520-cli-foundation`, `260521-cli-package-extraction`,
  `260529-sce-sunny-export` (CLI wrapper precedent).
- Source: `tools/sampling_mfa.jl`, `tools/vasp/vasptools.jl`.
