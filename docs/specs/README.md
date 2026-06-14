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
| [260521-cli-package-extraction](260521-cli-package-extraction/) | complete (2026-05-21) | Extract the Comonicon CLI into a `MagestyCLI` subdirectory package; de-Comonicon the core to restore JET coverage |
| [260522-adaptive-ridge-iterative](260522-adaptive-ridge-iterative/) | complete (2026-05-22) | Iterative `AdaptiveRidge` estimator (Frommlet & Nuel 2016 L0 approximation) with an analytic weighted-ridge refit |
| [260522-cluster-generation-benchmark](260522-cluster-generation-benchmark/) | complete (2026-05-22) | `bench-cluster` benchmark: per-stage timing of `Cluster` construction to locate the three-body bottleneck |
| [260522-refit-estimator](260522-refit-estimator/) | complete (2026-05-23) | `refit` post-selection refit on the selected basis support, scaled-magnitude criterion |
| [260523-cluster-generation-perf](260523-cluster-generation-perf/) | complete (2026-05-23) | Replace `irreducible_clusters` O(N^2) linear scan with canonical-form `Dict` lookup; remove duplicate `set_mindist_pairs` call in `generate_clusters` |
| [260523-design-matrix-3body-perf](260523-design-matrix-3body-perf/) | complete (2026-05-23) | `Fitting.jl` hot-path rewrite for 3-body SALCs: replace `Set{UInt}` de-dup with `Vector{UInt}` linear scan, and replace `@views spin_directions[:, atom]` SubArrays with stack `SVector{3, Float64}` reads at three sites |
| [260524-design-matrix-restructuring](260524-design-matrix-restructuring/) | complete (2026-05-24) | Four-stage algorithmic restructuring of `build_design_matrix_*`: Mf folding, orbit pre-enumeration, per-spinconfig SH cache, torque as cluster-major reverse-mode Jacobian. ~12× speedup on FeGe 2x2x2 light. |
| [260526-solver-unification-and-memory-log](260526-solver-unification-and-memory-log/) | complete (2026-05-26) | Unify OLS / Ridge / AdaptiveRidge around Cholesky on `X'X`, with `PosDefException` rethrown as a Ridge-pointing `ArgumentError`; add pre / post memory log to the design-matrix builders. ~3× faster solver. |
| [260529-sce-sunny-export](260529-sce-sunny-export/) | merged (2026-05-29, PR #31, #32) | `sce_to_sunny` exporter: emit a runnable Sunny.jl LSWT script from a fitted `SCEModel` (bilinear pairs + single-ion anisotropy; round-trip energy gate; primitive-cell unfolding incl. multiplicity>1). |
| [260601-src-refactor](260601-src-refactor/) | merged (2026-06-01, PR #33) | Architecture + hot-path cleanup from an `src/` review panel: encapsulate the angular-momentum coupling cache (M1) and cut SALC / SH-cache allocations incl. a combined SH value+gradient entry point (M2, ~2.4× torque SH fill). Fitting energy/torque kernel dedup (M3) deferred — the shared skeleton does not cleanly exist. Behavior-preserving (bit-identical). |
| [260605-mfa-sampling-cli](260605-mfa-sampling-cli/) | merged (2026-06-06, PR #34) | Promote `tools/sampling_mfa.jl` to `sample_mfa_incar` + the `magesty vasp mfa` CLI. Code-agnostic MFA sampler (`MfaSampling`) with a hand-rolled exact p=3 von Mises-Fisher draw (drops `Distributions`; adds only `Roots`) and a VASP `IncarIO` adapter. Sweep specified by point count; corrects a latent boundary-sign bug in the MFA solver (no sampling-output change). |
| [260610-gcv-diagnostics](260610-gcv-diagnostics/) | merged (2026-06-12, PR #36) | Hat-matrix generalized cross-validation on the combined energy+torque weighted objective: `gcv(f)`, ridge `gcv_lambda` penalty path (one SVD), and `gcv_learning_curve` data-sufficiency sweep (random subsets, mean +/- std). Linear estimators only; FitCheck-style table writers + Python plots. |
| [260614-magesty-module-split](260614-magesty-module-split/) | complete (2026-06-14) | Split the oversized `src/Magesty.jl` via `include` into same-namespace files `Evaluation.jl` (prediction verbs + accuracy metrics) and `GCV.jl` (GCV diagnostics + result types); not submodules. Verbatim move — public API / signatures / docstrings unchanged; main file ~2088 -> ~1260 lines; full suite + JET + Aqua green. |

This table and the `Status:` line in each `tasklist.md` are duplicated
intentionally; update both when a spec lands. Reconcile in bulk when
drift is found.

## Completion criteria

Each `tasklist.md` ends with a shared exit checklist (see
[`_template/tasklist.md`](_template/tasklist.md)). In particular,
**`.claude/agents/` is easy to forget** — whenever module names or
Makefile targets change, sweep the agent files as well.
