# `docs/specs/` index

Folders for mid-sized or larger development units. Operating rules live
in the "Managing development units" section of
[CLAUDE.md](../../CLAUDE.md). When starting a new spec, copy from
[`_template/`](_template/).

## Status table

| Spec | Status | One-line summary |
|---|---|---|
| [260513-estimator-dispatch](260513-estimator-dispatch/) | complete (2026-05-13) | Estimator-dispatch refactor for `fit(SCEFit, ...)` |
| [260513-write-xml-api](260513-write-xml-api/) | complete (2026-05-13) | Cleanup of `save` / `load` XML API |
| [260514-replace-sorted-container](260514-replace-sorted-container/) | complete (2026-05-14) | Introduce `SortedCounters` to drop a heavy dependency |
| [260514-sce-public-api](260514-sce-public-api/) | merged (2026-05-14, tasklist still says draft) | Public-API design: `SCEBasis` / `SCEDataset` / `SCEFit` / `SCEModel` |
| [260516-coupled-basis-atoms-svector](260516-coupled-basis-atoms-svector/) | abandoned (2026-05-16) | `SVector` for the atoms field; no performance win |
| [260516-coupled-basis-typeparam](260516-coupled-basis-typeparam/) | complete (2026-05-16) | Type parameterization `CoupledBasis{R}` for type stability |
| [260516-optimize-workspace](260516-optimize-workspace/) | complete (2026-05-16) | Scratch workspace for the Optimize hot path |
| [260516-src-layout-refactor](260516-src-layout-refactor/) | complete (2026-05-16) | Flatten `src/`, drop `common/` / `types/` / `utils/` |
| [260516-typed-input-spec](260516-typed-input-spec/) | merged (2026-05-16, tasklist still says draft) | Typed-value gate: `SystemSpec` / `InteractionSpec` / `SymmetryOptions` |
| [260517-elasticnet-estimator](260517-elasticnet-estimator/) | complete (2026-05-18) | `ElasticNet` (incl. `Lasso` convenience) via GLMNet.jl |
| [260518-energy-centered-design-matrix](260518-energy-centered-design-matrix/) | complete (2026-05-18) | Remove the bias column from `X_E`; replace with energy-only centering in `assemble_weighted_problem` |
| [260518-adaptive-lasso-oneshot](260518-adaptive-lasso-oneshot/) | complete (2026-05-19) | `AdaptiveLasso` oneshot estimator (Zou 2006 / ALAMODE recipe) via the existing GLMNet path |
| [260519-adaptive-lasso-precomputed-pilot](260519-adaptive-lasso-precomputed-pilot/) | complete (2026-05-19) | `PrecomputedPilot` adapter + `AdaptiveLasso(::SCEFit; ...)` / `AdaptiveLasso(::SCEModel; ...)` to reuse a fitted model as the pilot |
| [260519-pre-release-safety](260519-pre-release-safety/) | complete (2026-05-20) | `SpinConfig` unit-norm validation, in-memory basis fingerprint, integration tests for `SCEModel(fit)` / batched `predict_*` / save-load round trip (M3 `torque_weight` persistence deferred) |
| [260520-fitcheck-io-writers](260520-fitcheck-io-writers/) | complete (2026-05-20) | `write_energies` / `write_torques` text writers for the current API, feeding `FitCheck_*.py` visualization |
| [260520-cli-foundation](260520-cli-foundation/) | complete (2026-05-21) | Comonicon-based `magesty` CLI foundation; VASP-to-extxyz promoted to API + `magesty vasp extxyz` subcommand |
| [260521-cli-package-extraction](260521-cli-package-extraction/) | draft (2026-05-21) | Extract the Comonicon CLI into a `MagestyCLI` subdirectory package; de-Comonicon the core to restore JET coverage |

This table and the `Status:` line in each `tasklist.md` are duplicated
intentionally; update both when a spec lands. Reconcile in bulk when
drift is found.

## Completion criteria

Each `tasklist.md` ends with a shared exit checklist (see
[`_template/tasklist.md`](_template/tasklist.md)). In particular,
**`.claude/agents/` is easy to forget** — whenever module names or
Makefile targets change, sweep the agent files as well.
