# Design: MFA spin-sampling CLI (`magesty vasp mfa`)

Status: complete (2026-06-05)

## Summary

Three layers, separated so future DFT codes reuse the physics:

```
Layer 1  code-agnostic MFA sampler (core)   spin_matrix + sweep params -> Vector{Matrix}
Layer 2  DFT-code I/O adapter (core)        VASP INCAR <-> spin matrix (future: QE, ...)
Layer 3  CLI wrapper (MagestyCLI)           magesty vasp mfa  (future: magesty qe mfa)
```

Layer 1 holds the von Mises-Fisher (vMF) direction sampler and the MFA
self-consistency solver and knows nothing about file formats. Layer 2 reads
the initial spin matrix from an INCAR and writes sampled matrices back as
INCARs. Layer 3 is a thin Comonicon `@cast` leaf calling the exported
orchestration function. A future `magesty qe mfa` adds only a Layer 2 adapter
plus a Layer 3 leaf and reuses Layer 1 verbatim.

**vMF sampling without `Distributions`.** The script used
`Distributions.VonMisesFisher`, whose only role is to sample a direction on
S² given a mean direction and concentration `κ`. For p = 3 this has a closed
form (Ulrich 1984 / Wood 1994): draw the cosine `w` to the mean direction via
the inverse CDF `w = 1 + (1/κ)·log(u + (1−u)·e^{−2κ})`, `u ~ U(0,1)`, pick a
uniform azimuth in the tangent plane, and rotate the tangent-normal frame onto
the mean direction. This is exact (no rejection) and ~10 lines, so we drop the
heavy `Distributions` dependency and add only the lightweight `Roots` (used by
the 1-D Brent solve of the MFA equation).

The formula degenerates as `κ → 0` (`w → 1` for all `u`), so the caller never
reaches it in that regime: `mfa_sample` follows the script's κ branches —
`τ > MAX_TEMP ⇒ κ = 1e-6` and `τ < MIN_TEMP ⇒ return input unchanged` — and
`sample_vmf_direction` additionally draws an isotropic direction when
`κ < κ_min` (a small positive threshold) instead of evaluating the inverse
CDF. For large `κ` the `log` argument is computed in a numerically safe form
(`e^{−2κ}` underflows cleanly to 0, giving `w = 1 + (1/κ)·log(u)`).

## Module layout

| Target | Change |
|---|---|
| `src/MfaSampling.jl` | New. Layer 1: code-agnostic MFA physics + sweep driver. |
| `src/IncarIO.jl` | New. Layer 2: `parse_incar` / `write_incar` ported from `tools/vasp/vasptools.jl`. |
| `src/VaspSampling.jl` | New. Layer 2 orchestration: `sample_mfa_incar`. |
| `src/Magesty.jl` | `include` the three files; `export sample_mfa_incar`. |
| `Project.toml` | Add `Roots` to `[deps]` and `[compat]`. |
| `cli/src/MagestyCLI.jl` | Add `mfa` leaf inside `@cast module Vasp`. |
| `tools/sampling_mfa.jl` | Remove (logic promoted). |

`MfaSampling.jl` uses `Random`, `LinearAlgebra`, `StaticArrays` (all existing
deps) plus `Roots`. `IncarIO.jl` uses `OrderedDict` and `Printf` (existing
deps) — no new dependency. The script ports faithfully: core already depends on
`DataStructures`, which **re-exports `OrderedCollections.OrderedDict`**
(`DataStructures.OrderedDict === OrderedCollections.OrderedDict`, verified), so
the INCAR dict type is byte-for-byte the same as `vasptools.jl`'s and no
`OrderedCollections`/`Distributions` entry is added.

## API

### Layer 1 (code-agnostic, `Magesty.`-qualified; export `mfa_sweep` TBD)

```julia
# MFA self-consistency: solve m = coth(3m/τ) − τ/3m
thermal_averaged_m(τ::Real)::Float64
tau_from_magnetization(m::Real)::Float64

# exact p=3 vMF direction sample about `mean_dir` (unit) with concentration κ
sample_vmf_direction(mean_dir::SVector{3,Float64}, κ::Real)::SVector{3,Float64}

# one sampled config (3×N); `variable` ∈ ("tau","m"), `value` the sweep value
mfa_sample(spin_matrix::AbstractMatrix{<:Real},
           variable::AbstractString, value::Real)::Matrix{Float64}

# full sweep -> configs in output order
mfa_sweep(spin_matrix::AbstractMatrix{<:Real};
          variable::AbstractString, start::Real, stop::Real, num_points::Integer,
          num_samples::Integer = 1, randomize::Bool = false,
          fixed_indices::AbstractVector{<:Integer} = Int[],
          uniform_indices::AbstractVector{<:Integer} = Int[]
         )::Vector{Matrix{Float64}}
```

Sweep values are `range(start, stop, length = num_points)`; the output order
is `(point, sample)` with `point` outer, matching the script's
`output_index*num_sample + i` numbering.

### Layer 2 orchestration (exported)

```julia
function sample_mfa_incar(incar_path::AbstractString;
        variable::AbstractString,
        start::Real, stop::Real, num_points::Integer,
        num_samples::Integer = 1,
        randomize::Bool = false,
        fix::AbstractString = "",
        uniform_atoms::AbstractString = "",
        outdir::AbstractString = ".",
        prefix::AbstractString = "sample",
    )::Vector{String}
```

Reads `MAGMOM` (else `M_CONSTR`) from the INCAR, builds the `3×N` matrix,
parses `fix` / `uniform_atoms` index specs (`"1-10,12,20-22"`), calls
`mfa_sweep`, and writes each config to `joinpath(outdir, "<prefix>-NN.INCAR")`.
File numbering follows the script exactly: for sweep point `p` (1-based over
`num_points`) and sample `s` (1-based over `num_samples`), the file index is
`(p − 1) * num_samples + s`, zero-padded to a width of
`ndigits(num_points * num_samples)`. The output INCAR copies all input keys,
sets both `MAGMOM` and `M_CONSTR` to the sampled vector, and uses
`wrap_vectors = true`. Returns the written paths. Standard Julia docstring
(`# Arguments` / `# Returns` / `# Examples`).

`prefix` is a Layer 2 keyword (default `"sample"`) for Julia callers; the CLI
leaf intentionally does **not** expose `--prefix` (keeps the command surface
minimal, matching the script's fixed `sample-NN.INCAR` naming).

### Layer 3 CLI

```julia
Comonicon.@cast function mfa(incar::String, variable::String;
        start::Float64, stop::Float64, num_points::Int,
        num_samples::Int = 1, randomize::Bool = false,
        fix::String = "", uniform_atoms::String = "", outdir::String = ".")
    variable in ("tau", "m") || error("variable must be \"tau\" or \"m\"")
    paths = Magesty.sample_mfa_incar(incar; variable, start, stop, num_points,
        num_samples, randomize, fix, uniform_atoms, outdir)
    println("Wrote $(length(paths)) INCAR file(s) to $(outdir).")
    return nothing
end
```

Comonicon flags: `--start --stop --num-points --num-samples --randomize --fix
--uniform-atoms --outdir`. Validation of `start ≤ stop` and `num_points ≥ 1`
lives in `sample_mfa_incar` (so Julia callers get it too).

## Types and conventions

- Spin matrix `3 × n_atoms`, columns are atoms; directions unit, magnitudes
  preserved. No convention change.
- vMF/MFA numerical conventions preserved exactly (see requirements
  Invariants), including the `MIN_TEMP`/`MAX_TEMP` clamps and the κ special
  branches.
- New invariant: `sample_vmf_direction` reproduces the vMF law with mean
  resultant length `A₃(κ) = coth κ − 1/κ`.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): none.
- [ ] SCE coefficient XML (`save` / `load`): none.
- [ ] `Fitting` <-> `SALCBasis`: none.
- [x] `.claude/agents/`: check module-name lists (new `MfaSampling`/`IncarIO`/
      `VaspSampling`) — sweep agent files if they enumerate modules.
- [x] `SPEC.md` / `docs/src/api.md`: add `sample_mfa_incar`; add the CLI
      narrative entry for `magesty vasp mfa`.

## Test strategy

- `test/component/test_MfaSampling.jl`:
  - vMF: with fixed `Random.seed!`, large-N sample mean resultant length
    matches `A₃(κ)` within tolerance; direction mean aligns with `mean_dir`;
    κ→0 isotropy.
  - `thermal_averaged_m` ↔ `tau_from_magnetization` round-trip; boundary
    clamps; `mfa_sample` magnitude preservation and zero-norm skip.
  - `fix` / `uniform` / `randomize` interaction; `parse_atom_index_spec` edge
    cases and error paths.
  - Expected values derived analytically — not pinned to current output.
- `test/component/test_IncarIO.jl`: parse → write → parse round-trip on an
  INCAR fixture; MAGMOM/M_CONSTR formatting and `wrap_vectors`.
- `cli/test/runtests.jl`: `command_main(["vasp","mfa", incar, "tau",
  "--start","0.1","--stop","0.3","--num-points","3","--outdir",tmp])`;
  assert exit 0 and `num_points * num_samples` files with the expected naming;
  assert per-atom magnitudes preserved. Use structural invariants rather than
  exact (stochastic) spin values.
- Fixture: small INCAR with `MAGMOM` + `M_CONSTR` under
  `test/component/fixtures/`.

## Risks and open items

- **Numerical**: the hand-rolled vMF replaces a library sampler. Mitigated by
  the analytic `A₃(κ)` test and the κ-branch preservation; `numerical-reviewer`
  to confirm equivalence.
- **Open**: export `mfa_sweep` for direct Julia reuse now, or keep it
  `Magesty.`-qualified until a second DFT adapter needs it? Default: keep
  non-exported initially (only `sample_mfa_incar` exported).
- **Open**: `IncarIO.jl` as a standalone module vs. folding into `VaspIO.jl`.
  Default: standalone `IncarIO.jl` (INCAR text format is unrelated to the
  `vasprun.xml` parsing in `VaspIO.jl`).
