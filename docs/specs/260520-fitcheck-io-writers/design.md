# Design: `write_energies` / `write_torques` for the current API

Status: draft (2026-05-20)

## Summary

`write_energies` and `write_torques` are thin text-I/O wrappers over the
existing prediction verbs. They resolve a `(model, spinconfigs)` pair from
their arguments, call `predict_energy` / `predict_torque`, read the observed
values off the `SpinConfig`s, and write a fixed-format text file matching the
legacy layout consumed by `tools/FitCheck_energy.py` /
`tools/FitCheck_torque.py`.

Each writer offers two signature families: a self-contained
`write_*(f::SCEFit, filename)` that evaluates the fit's own training dataset,
and a general `write_*(predictor, data, filename)` for evaluating any
predictor against any dataset (training, validation, or test). The
`predictor` accepts `SCEModel` or `SCEFit`; `data` accepts the same forms as
the existing `(predictor, data)` evaluation verbs (`SCEDataset` /
`AbstractVector{SpinConfig}` / EMBSET path).

The writers live in a new included file rather than a submodule, so they see
`SCEFit` / `SCEModel` / `predict_*` / `read_embset` / `Structure` directly in
the `Magesty` module scope.

## Module layout

| Target | Change |
|---|---|
| `src/FitCheckIO.jl` | New file. Defines `write_energies`, `write_torques`, and private helpers. Plain top-level definitions, no `module`. |
| `src/Magesty.jl` | Add `include("FitCheckIO.jl")` after `include("XMLIO.jl")`; add `export write_energies, write_torques`. |
| `test/component/test_fitcheck_io.jl` | New regression test. |
| `test/runtests.jl` | Register the new test in the `component` testset. |
| `docs/src/api.md`, `docs/src/tools.md`, `SPEC.md`, `CHANGELOG.md` | Document the new exported API. |

Reused unchanged: `predict_energy` / `predict_torque` batch overloads
(`src/Magesty.jl`), `read_embset` (`src/SpinConfigs.jl`), `Structure` fields
`kd_name` and `supercell.kd_int_list` (`src/Structures.jl`). `Printf` is
already imported in `Magesty.jl`, so `@sprintf` is available.

## API

```julia
# Self-contained: evaluate against the fit's own training dataset.
write_energies(f::SCEFit, filename::AbstractString = "energy_list.txt")::Nothing
write_torques(f::SCEFit,  filename::AbstractString = "torque_list.txt")::Nothing

# Explicit predictor + evaluation data. `filename` has no default here:
# `data` may itself be a string (EMBSET path), so a defaulted `filename`
# would let a 2-arg `write_*(model, "EMBSET")` call silently treat the
# EMBSET path as the output path. Requiring an explicit `filename` makes
# the 2-arg form unambiguously `(SCEFit, output_path)`.
write_energies(predictor::Union{SCEModel, SCEFit},
               data::Union{SCEDataset, AbstractVector{SpinConfig}, AbstractString},
               filename::AbstractString)::Nothing
write_torques(predictor::Union{SCEModel, SCEFit},
              data::Union{SCEDataset, AbstractVector{SpinConfig}, AbstractString},
              filename::AbstractString)::Nothing
```

- `predictor`: `SCEModel` or `SCEFit`.
- `data`: `SCEDataset` / `AbstractVector{SpinConfig}` / `embset_path::AbstractString`.
- Both functions are exported and carry standard Julia docstrings
  (`# Arguments` / `# Returns` / `# Examples`) with explicit type
  annotations.

Internal flow (shared private helpers):

1. Resolve `(model::SCEModel, spinconfigs::Vector{SpinConfig})`:
   - `f::SCEFit` alone → `model = SCEModel(f)`, `spinconfigs =
     f.dataset.spinconfigs`; predictions reuse the dataset's stored `X_E` /
     `X_T` (`predict_*(model, f.dataset)`).
   - `(predictor, data)` → `model` is `predictor` or `SCEModel(predictor)`;
     `data` is resolved: `SCEDataset` → `data.spinconfigs`;
     `AbstractVector{SpinConfig}` → as-is; `AbstractString` →
     `read_embset(data)` (read once, reuse for both observed and predicted).
2. Observed energies `[sc.energy for sc]`; predicted via `predict_energy`.
3. Observed torques `[sc.torques for sc]` (each `3 × n_atoms`); predicted via
   `predict_torque`.
4. Element symbols from the predictor's structure:
   `structure.kd_name[structure.supercell.kd_int_list[i]]`.
5. Write via `@sprintf` in the legacy layout, inside `open(...) do ... end`
   with `try` / `catch` `@error` + `rethrow`.

## Output format (legacy, reproduced verbatim)

Energy file (default `"energy_list.txt"`):

```
# data index,    DFT_Energy,    SCE_Energy
# unit of energy is eV
  1    % 15.10e    % 15.10e
  ...
```

Index width = `ndigits(n_configs)`. Energies in eV. `FitCheck_energy.py`
treats column 1 as the index and uses columns 2/3.

Torque file (default `"torque_list.txt"`):

```
# atom index,    element,   DFT_torque_x, ..., SCE_torque_z
# unit of torque is eV
# data index: 1
  1 Fe  % 15.10e   % 15.10e   % 15.10e    % 15.10e   % 15.10e   % 15.10e
  ...
# data index: 2
  ...
```

8 columns; element symbol per atom; torque components in eV.

## Types and conventions

No physics-convention, unit, or numerical-convention impact. Energies and
torques are emitted in the units of the training data (eV), matching what the
`predict_*` verbs return and what `FitCheck_*.py` expects (the Python scripts
convert eV → meV for display).

New invariant: the legacy text format is now a supported output contract;
`test/component/test_fitcheck_io.jl` pins it.

## Impact on linked sites

- ~~Spherical-harmonics convention (`TesseralHarmonics`)~~: N/A.
- ~~SCE coefficient XML (`save` / `load`)~~: N/A — separate text format.
- ~~`Fitting` <-> `SALCBasis`~~: N/A — predictions go through existing verbs.
- ~~`.claude/agents/` references~~: N/A — no module rename or Makefile target
  change.
- [x] `SPEC.md` / `docs/src/api.md` updates: add `write_energies` /
      `write_torques` to the public API listing; mention in `docs/src/tools.md`
      next to the `FitCheck_*.py` documentation.

## Test strategy

New `test/component/test_fitcheck_io.jl`, registered in the `component`
testset of `test/runtests.jl`:

- Build a small `SCEBasis` + `SCEDataset` from an existing test fixture
  (the dimer TOML + FM/AFM configs, as in `test_SCEFit.jl`), `fit` an
  `SCEFit`.
- Write to `mktempdir()` paths and assert:
  - Header lines present; energy file has `n_configs` data rows of 3 numeric
    columns; torque file has one `# data index:` block per config, each with
    `n_atoms` rows of 8 columns (the x/y/z components are columns within a
    row, not separate rows).
  - Energy columns 2/3 equal `sc.energy` and `predict_energy(f, dataset)` to
    the printed precision; torque columns equal `sc.torques` and
    `predict_torque`.
  - The `(predictor, data)` form with `data` = `SCEDataset`,
    `Vector{SpinConfig}`, and an EMBSET path all produce output identical to
    the `SCEFit`-only form on the training set.

End-to-end: run `FitCheck_energy.py` / `FitCheck_torque.py` on generated files
and confirm they parse.

No hot path touched → no `.claude/bench_log.md` entry.

## Risks and open items

- No numerical-result risk: pure I/O over existing verbs.
- Open: branch policy — feature work, not a refactor; confirm whether it lands
  on `main` or a `feat/fitcheck-io-writers` topic branch before committing.
