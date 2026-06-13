# Whole-package review sweep — open findings backlog

**Status**: in progress (2026-06-13)

A four-axis review panel (numerical / maintainability / performance / API)
was run over the whole `src/` tree (27 modules, ~14k lines). This note
records the findings that have **not** yet been addressed, so they can be
picked up later without re-running the review. Each entry has a
`file:line`, the problem, a proposed fix, and any decision the change
needs. Sever­ity follows the panel schema (blocker / major / minor).

## Already landed

- **Blocker (numerical)** — XML lattice round-trip transposed
  `lattice_vectors` for low-symmetry cells. Fixed: writer now emits
  columns; regression test added for an asymmetric lattice.
  (commit `df5d941`)
- **Major (maintainability + performance)** — hidden global coupling
  cache. Replaced with a per-construction cache threaded through the
  build; no lock needed (population is single-threaded, before the SALC
  `@threads` loop). The performance reviewer's "concurrent write" race
  claim was investigated and found **not** to occur in the current call
  graph. (commit `f28a84d`)
- **Major (maintainability) — internal `_` prefixes + GCV dedup + minors.**
  Renamed 11 internal helpers in `SALCBases.jl` and 7 in `Clusters.jl` to
  the leading-`_` convention (call sites in `src/`, `test/`, and
  `tools/micromagnetics.jl` updated); removed dead `is_symmetric`, broke a
  long line, expanded the `SALCBases` module docstring; extracted
  `_gcv_core` to dedup the GCV setup; index loops `in`->`=` in `Fitting.jl`
  / `ExtXYZ.jl`; `DEFAULT_PERIODICITY` made an immutable tuple. Full unit
  suite green (22,772). **Still open below:** the OSZICAR parser dedup and
  the package-wide indentation normalization.
- **Major (numerical) — zero-moment atom handling.** `read_embset` no longer
  normalizes `0/0` for a non-magnetic site; it stores a placeholder unit
  direction (never consumed unless the site enters the basis) and the zero
  magnitude makes the torque vanish. The `SCEDataset` constructor now rejects
  a configuration that puts a (near-)zero moment on a SALC-referenced atom,
  catching both a non-magnetic site wrongly pulled into the basis and a
  magnetic site that collapsed to zero. Regression tests added for both.
- **API — validation + docstrings + convention.** `oszicar_to_embset` now
  validates `energy_kind in ("f","e0")` up front (was a misleading "energy
  not found" later); SPEC argument name aligned to `filename`;
  `Magesty.save` / `load` docstrings completed with `# Returns` / `# Examples`.
  The "inconsistent docstring sections" finding rested on a wrong premise —
  the `# Arguments` / `# Keyword arguments` split is the package majority (8
  modules), not three outliers — so instead of mass-editing, the split style
  was codified as canonical in `STYLE_GUIDE.md`. `Magesty.jl`'s unified
  docstrings may converge to it opportunistically; not blocking.
- **Performance — `_predict_energy` buffer removal.** Fused the per-call
  `design_vector` allocation plus `dot` into a single scalar accumulator in
  the SALC loop. Each term is unchanged; only the summation runs sequentially
  (ULP-level difference, `≈` tests unaffected). Per-config inference on the
  fept baseline: -5.8% time, -2 allocations/config. Recorded in
  `.claude/bench_log.md`; full unit suite green (22,780).
- **Performance — `assemble_weighted_problem` in-place assembly.** Replaced the
  six intermediate arrays + two `vcat`s (peak ~3x the design-matrix size) with
  fused broadcast writes into preallocated `X` / `y`. Output is bit-identical
  (verified `X == X_old && y == y_old` at weight 0.5 and 1.0), so all tests pass
  unchanged. fept fixture: -57% time, -49% peak memory; mitigates the OOM risk
  on large torque-rich datasets. Recorded in `.claude/bench_log.md`.

## Deferred (decision pending)

- **OSZICAR parser dedup** (`VaspIO.jl`). A shared `_scan_oszicar_records`
  scanner was prototyped and passed all tests byte-for-byte, but it added
  ~68 net lines and a four-predicate injection because the two parsers
  differ in header match, block-end rule, EOF handling, and column
  selection. Reverted pending a decision on whether the consolidation of
  the (genuinely triplicated) block-walking state machine is worth the
  added indirection. The `in`->`=` index-loop minor in this file went back
  with the revert and should ride along with whatever is decided.

## Numerical

- **Minor — overloaded constant.** `MfaSampling.jl:62-63`: the
  temperature guards `(MIN_TEMP, MAX_TEMP)` are reused as the bracket for
  the magnetization root `m`. Correct today only because `m in (0,1)`
  coincides with the guard interval. Introduce a named `(M_MIN, M_MAX)`
  bracket to decouple the two roles.
- **Minor — implicit unitarity assumption.**
  `SphericalHarmonicsTransforms.jl:10-14`: `r2c_sph_harm_matrix` returns
  `c2r_sph_harm_matrix(l)'`, which is the true inverse only because
  `c2r` is unitary. Add a comment stating the assumption (or an `@assert`
  in a test) so a future normalization change cannot silently break it.

## Maintainability

The internal `_`-prefix renames, GCV dedup, dead-code removal, long-line
break, module-docstring expansion, index-loop style, and the immutable
periodicity default are done (see "Already landed"). The OSZICAR parser
dedup is under "Deferred". Remaining:

- **Major — `ENV` lookup in a keyword default.** `SALCBases.jl` (the
  `check_irrep_unitary` keyword):
  `get(ENV, "MAGESTY_CHECK_IRREP_UNITARY", "0") == "1"` is a hidden,
  per-call side channel. Read the env var once at load time into a module
  constant, or make the knob an explicit, documented keyword.
- **Minor — indentation inconsistency.** Tabs vs 4-space mixed across the
  package (STYLE_GUIDE prescribes 4-space). Normalize in a single
  mechanical commit (re-baseline benchmarks if it touches hot-path files).
  Large and noisy; kept separate on purpose.
- **Minor (spec-level) — main module size.** `Magesty.jl` (~2000 lines)
  mixes type definitions, fitting, prediction, ~24 evaluation metrics,
  GCV, and I/O. Extracting `Metrics.jl` and a GCV-public file would each
  shed ~200 lines. Pursue via a spec if taken up.

## Performance

Hot-path items (high value):

- **Major — type erasure in `salc_list`.** `SALCBases.jl:75,262`: the
  containers hold `Vector{CoupledBasis_with_coefficient}` (the `{R,N}`
  parameter is dropped), so the inner design-matrix kernels dispatch
  dynamically per element instead of monomorphizing. Each key group is
  homogeneous in body count; split or dispatch on `first(key_group)` to a
  type-stable helper. **Needs numerical sign-off** (float summation order
  may change) and a before/after `@btime` in the bench log.

Setup-time items (lower value, but matter for large supercells):

- **Major — O(n^2) translational-equivalence dedup.**
  `SALCBases.jl:536-549`: linear scan over `orbit_basis_list`. Implement
  `Base.hash(::CoupledBasis)` over the integer fields (exclude the float
  tensor, matching `salc_fingerprint`) and dedup via a `Dict`.
- **Major — defensive buffered SH derivative variants.**
  `TesseralHarmonics.jl:549-`: the unbuffered `∂Zₗₘ_∂r̂{x,y,z}_unsafe` /
  `zzₗₘ_unsafe` allocate a `dnPl` vector per call. Not on the current hot
  path (the kernels use the buffered `∂ᵢZlm_unsafe`), but add buffered
  overloads to guard against a future hot-path caller.
- **Minor** — `√2` as a module constant
  (`TesseralHarmonics.jl`, several sites); `zzₗₘ_unsafe` recomputing the
  Legendre recursion three times (`TesseralHarmonics.jl:683-686`);
  `_canonicalize_eigenspace` forming an explicit n×n projector
  (`SALCBases.jl:188-206`) instead of `V * (V' * e_j)`.
- **Profiler recommended** — measure SALC-construction throughput
  (`construct_and_classify_coupled_basislist` /
  `projection_matrix_coupled_basis`) on a representative large supercell
  to confirm where the O(n^2) dedup and the projector dominate.

## API / UX

The `energy_kind` validation, the SPEC argument-name fix, the
`Magesty.save` / `load` docstring sections, and the docstring-section
convention are done (see "Already landed"). Remaining:

- **Minor** — `SpinConfig` docstring lacks `# Examples`
  (`SpinConfigs.jl:52`); `atol_unit_norm::Float64` should be `::Real` with
  internal conversion (`SpinConfigs.jl:110`); `saxis` has no length / norm
  validation, so `[0,0,0]` silently gives an identity rotation
  (`VaspConvert.jl:397`); `intercept` / `nobs` / `dof` docstrings lack the
  section layout (`Magesty.jl:935-949`); `AbstractEstimator` / `OLS` /
  `Ridge` docstrings lack `# Arguments` (`Fitting.jl:138-197`); `Lasso` is
  a factory function returning `ElasticNet`, not a type — worth a one-line
  SPEC annotation (`Magesty.jl:93`).

## Confirmed clean (for the record)

Physics conventions (3 x n_atoms layout, unit spin vectors, real
tesseral `Zₗₘ`, `Jφ` / `j0` separation), linked-site synchronization (the
spherical-harmonics convention `Z = C·Y` verified to ~2e-16 across
l = 0..4; SALC key-group ordering; jphi serialized in `salc_index`
order), and numerical stability (Cholesky SPD guarantees, ridge
conditioning, GCV denominator saturation, RotationMatrix gimbal-lock
branch, MfaSampling kappa / zero-norm guards) were all verified consistent
— the lattice transpose (now fixed) was the single synchronization break.
