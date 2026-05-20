# Extract the CLI into a separate `MagestyCLI` package

Status: not started — follow-up to spec `260520-cli-foundation`.

## Problem

The `magesty` command-line interface (`src/CLI.jl`) is built with
Comonicon.jl, whose `@main` macro runs command-tree code generation at
macro-expansion time. `JET.report_package(Magesty)` macro-expands every
top-level form in the package through a JuliaInterpreter-based loader, and
that code generation throws under the interpreter
(`Comonicon.AST.Entry has no field name`).

This is an external Comonicon/JET incompatibility, not a Magesty defect:
`using Magesty` and `make test-all` are unaffected. But while the CLI
lives inside the `Magesty` source tree, JET cannot analyze the package at
all — `test/jet.jl` is currently marked `@test_broken`.

No JET configuration avoids this: `target_modules` only filters reported
issues, and `concretization_patterns` acts after macro expansion. The
crash is during expansion, which JET performs unconditionally for all
package top-level code. The only fix is structural — move the Comonicon
code out of the `Magesty` source tree that `report_package` walks.

## Proposed approach: a separate `MagestyCLI` package

Extract the CLI into its own package, `MagestyCLI`, depending on `Magesty`
and `Comonicon`. The `magesty` launcher then runs
`using MagestyCLI; exit(MagestyCLI.command_main())`.

A Julia **package extension** was considered and rejected: an extension
module cannot be `using`-ed, so Comonicon's `comonicon_install` — which
generates a launcher naming the module `@main` ran in — does not work, and
the launcher would have to be hand-rolled. A separate package keeps
Comonicon's standard install machinery intact.

Scope sketch:

- `Magesty` core loses `src/CLI.jl`, the `Comonicon` dependency, and
  `deps/build.jl`; `vasp_to_extxyz` and the `VaspIO` / `ExtXYZ` modules
  stay in core (they carry no Comonicon dependency).
- `MagestyCLI` holds the `@cast` / `@main` definitions and the thin
  command wrappers over `Magesty`'s API.
- Installation docs gain a step to add `MagestyCLI`.
- `test/jet.jl` returns to a plain assertion (the `@test_broken` guard
  there reverts on its own once `report_package` stops throwing).

## Trade-off

A second package in the repository and a slightly longer install path,
versus restored JET coverage of the numerical core. Deferred to a
dedicated spec; not a release blocker.
