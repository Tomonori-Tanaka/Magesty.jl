# Design: CLI foundation (Comonicon) with VASP-to-extxyz pilot

Status: draft (2026-05-20)

## Summary

Adopt Comonicon.jl as the CLI framework. Comonicon turns annotated
package functions (`@cast`) into a command tree with auto-generated help,
and installs a launcher into `~/.julia/bin/`. The launcher invokes
`julia` with `using Magesty` and the generated entrypoint — it resolves
the package *by name through an environment*, not by a baked absolute
path, so it follows `Pkg.update` instead of breaking (the failure mode of
the removed `install_tools`).

The pilot proves the full pattern each later tool spec will reuse:
**pure conversion logic moves into a `src/` submodule as plain functions;
a thin `@cast` subcommand wraps it.** For the VASP-to-extxyz converter the
logic is already cleanly separable — `VaspParser` is 100% library code and
`ExtXYZWriter` ~95% — so promotion is a move plus a thin high-level
wrapper, not a rewrite.

Subcommands are grouped per DFT code: the converter is reached as
`magesty vasp extxyz`, with `vasp` a Comonicon command group
(`@cast module Vasp`) and `extxyz` the leaf naming the target format.
Future codes follow the same shape — `magesty qe extxyz`, etc. — and the
code-independent `version` command stays at the top level.

Plain install only. Comonicon's `[sysimg]` / `[application]` modes bake a
frozen package snapshot and need a rebuild on every update; they are
deferred until CLI startup latency is shown to matter.

## Module layout

| Target | Change |
|---|---|
| `Project.toml` | Add `Comonicon` to `[deps]` + `[compat]`; add `[tool.comonicon]` (`name = "magesty"`). |
| `src/ExtXYZ.jl` | New. `module ExtXYZ` — **renamed from `ExtXYZWriter`** — `AtomFrame`, `write_extxyz`. Conversion logic unchanged from `tools/ExtXYZ.jl` (dep: `Printf`). |
| `src/VaspIO.jl` | New. `module VaspIO` — **renamed from `VaspParser`** — `VaspRunData`, `OszicarMagData`, `parse_vasprun`, `parse_oszicar_magdata`. Logic unchanged from `tools/vasp/VaspParser.jl` (deps: `EzXML`, `Printf` — both already in `[deps]`). |
| `src/CLI.jl` | Holds the high-level `vasp_to_extxyz` function, the `@cast module Vasp` command group, and the Comonicon `@cast` / `@main` definitions. Created as a skeleton in M1; the `Vasp` group is added in M3. |
| `src/Magesty.jl` | Insert `include` of `ExtXYZ.jl`, `VaspIO.jl` after `XMLIO.jl`; include `CLI.jl` last. Final order: `... → XMLIO.jl → ExtXYZ.jl → VaspIO.jl → FitCheckIO.jl → CLI.jl` (`FitCheckIO.jl` is the current last include). Export `vasp_to_extxyz`. |
| `deps/build.jl` | New. Calls `Magesty.comonicon_install()` to install the `magesty` launcher into `~/.julia/bin`. |
| `tools/ExtXYZ.jl`, `tools/vasp/VaspParser.jl` | Removed (logic now in `src/`). |
| `tools/vasp/vasp2extxyz.jl` | M2: the two `include` lines + `using .ExtXYZWriter` / `using .VaspParser` become `using Magesty.ExtXYZ` / `Magesty.VaspIO`, keeping it runnable. M3: removed once its logic is lifted into `vasp_to_extxyz` and the `magesty vasp extxyz` subcommand. |
| `tools/vasp/vasp2extxyz_recursive.jl` | Same M2 loading update as `vasp2extxyz.jl`. Its `using ArgParse` is **kept** — the script still needs `ArgParse` in its run environment until its own promotion spec. Stays a script; full promotion deferred. |
| `Makefile`, `CLAUDE.md` | The `make test-tools` target (and its `ci-local` reference) is dropped — `tools/test/` no longer exists; the CLAUDE.md test-target table is updated to match. |
| `test/component/test_ExtXYZ.jl`, `test/component/test_VaspIO.jl` | Moved from `tools/test/`, retargeted at `Magesty.ExtXYZ` / `Magesty.VaspIO`, and added to `test/runtests.jl`. `tools/test/fixtures/` is relocated to `test/component/fixtures/`. |
| `docs/src/installation.md`, `tools.md`, `api.md`, `SPEC.md` | Updated (see Impact). |

## API

```julia
# src/VaspIO.jl  (module VaspIO)
parse_vasprun(path::AbstractString)::VaspRunData
parse_oszicar_magdata(path::AbstractString, m_constr, num_atoms)::OszicarMagData

# src/ExtXYZ.jl  (module ExtXYZ)
write_extxyz(io::IO, frame::AtomFrame)
write_extxyz(path::AbstractString, frame::AtomFrame)

# src/CLI.jl — high-level API, exported from Magesty
"""
    vasp_to_extxyz(vasprun; oszicar=nothing, output=nothing) -> String

Convert a VASP run to extended XYZ. ...
"""
function vasp_to_extxyz(vasprun::AbstractString;
                        oszicar::Union{AbstractString,Nothing} = nothing,
                        output::Union{AbstractString,Nothing} = nothing)::String

# src/CLI.jl — Comonicon CLI, reached as `magesty vasp extxyz ...`
@cast module Vasp
    import Comonicon
    Comonicon.@cast function extxyz(vasprun::String;
                                    oszicar::String = "", output::String = "")
        ...
    end
    Comonicon.@main
end
```

The leaf subcommand is a thin shell: it normalizes empty-string options to
`nothing` and calls `Magesty.vasp_to_extxyz`. All conversion behavior
lives in the exported API function so it is testable without the CLI.
`@cast` / `@main` are qualified with `Comonicon.` because `@main` clashes
with `Base.@main` on Julia 1.12.

### Install path

`[tool.comonicon]` in `Project.toml` sets the command name (`magesty`).
`deps/build.jl` calls `Magesty.comonicon_install()`, so the launcher is
installed by `Pkg.build("Magesty")`; `installation.md` documents this.
The generated launcher resolves Magesty by name through a dedicated
environment — no absolute package path is baked in.

## Types and conventions

No physics-convention, unit, sign, or numerical impact. `VaspRunData` /
`OszicarMagData` / `AtomFrame` move unchanged; the conversion arithmetic
(`build_vasp_comment`, RWIGS per-atom expansion, magnetic-moment columns)
is copied verbatim from `tools/vasp/vasp2extxyz.jl` into `vasp_to_extxyz`.

New invariant: extxyz output is byte-identical to the pre-spec
`tools/vasp/vasp2extxyz.jl` for identical inputs (enforced by the
regression test).

## Impact on linked sites

- [x] Spherical-harmonics convention (`TesseralHarmonics`): N/A — not touched.
- [x] SCE coefficient XML (`save` / `load`): N/A — not touched.
- [x] `Fitting` <-> `SALCBasis`: N/A — not touched.
- [x] `.claude/agents/`: new `src/` modules (`ExtXYZ`, `VaspIO`, `CLI`)
      and a new Makefile/test layout — sweep agent files for module lists
      and `tools/test` references.
- [x] `SPEC.md` / `docs/src/api.md`: new modules in the directory layout;
      `vasp_to_extxyz` in the public API; CLI section.

## Test strategy

- Relocate `test_ExtXYZ.jl` / `test_VaspParser.jl` (renamed
  `test_VaspIO.jl`) to `test/component/`, retargeted at the `src/`
  modules; keep all existing assertions unchanged.
- New integration test: `vasp_to_extxyz` on a `tools/test/fixtures/`
  input produces output byte-identical to a stored reference file.
- CLI smoke test: `magesty vasp extxyz --help` exits 0, and a real
  conversion via the subcommand exits 0 with output matching the
  reference (`tools/test/` or a new `test/` target).
- `make test-all` / `test-aqua` / `test-jet`.

## Risks and open items

- **Comonicon is a non-trivial new dependency.** It is loaded by every
  `using Magesty`, including users who never touch the CLI, adding to
  precompile/load time. *Resolved in M1:* measured load times (warm,
  precompiled) — `using Comonicon` alone ~0.20 s; `using Magesty`
  ~1.05 s with Comonicon vs. ~0.85 s without (upper-bound delta ~0.20 s,
  the full Comonicon subtree). This is small in absolute terms against a
  workload whose real cost (SALC construction, fitting) is seconds to
  minutes, so Comonicon stays a plain `[deps]` entry; the
  package-extension path is **not** taken.
- Removing `tools/vasp/vasp2extxyz.jl` breaks anyone invoking that exact
  script path; note in `CHANGELOG.md` (as was done for `install_tools`).
- `sysimg` / `application` install modes deferred; revisit if startup
  latency is reported as a problem.
