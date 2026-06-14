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
- **API — remaining UX minors.** `SpinConfig` docstring gained `# Examples`;
  `atol_unit_norm` relaxed to `::Real` with internal `Float64` conversion;
  `_saxis_rotation` now rejects a zero / non-finite / wrong-length `saxis`
  (was a silent identity rotation) with a regression test;
  `intercept` / `nobs` / `dof` docstrings gained `# Returns`; SPEC notes
  that `Lasso` is a convenience function returning `ElasticNet`, not a type.
  The `AbstractEstimator` / `OLS` / `Ridge` "missing `# Arguments`" finding
  was declined: `OLS` is a zero-field singleton, `AbstractEstimator` is
  abstract, and `Ridge` documents `lambda` via `# Fields` (matching the
  `ElasticNet` estimator style) — forcing `# Arguments` would contradict the
  package's own convention.
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

All numerical minors are done.

- **Minor — overloaded constant.** *(done)* `MfaSampling.jl`: the
  magnetization root in `thermal_averaged_m` now brackets on named
  `(M_MIN, M_MAX)` constants instead of reusing the temperature guards
  `(MIN_TEMP, MAX_TEMP)`. The values are numerically coincident, so the
  result is bit-identical; only the role is decoupled.
- **Minor — implicit unitarity assumption.** *(done)*
  `SphericalHarmonicsTransforms.jl`: `r2c_sph_harm_matrix` now documents
  that its adjoint-equals-inverse identity holds only because `c2r` is
  unitary, and a named "Unitarity of c2r" testset asserts `c2r' * c2r ≈ I`
  (and `c2r * c2r' ≈ I`) for `l = 0..5`, so a normalization change that
  breaks unitarity fails here rather than silently downstream.

## Maintainability

The internal `_`-prefix renames, GCV dedup, dead-code removal, long-line
break, module-docstring expansion, index-loop style, and the immutable
periodicity default are done (see "Already landed"). The OSZICAR parser
dedup is under "Deferred". Remaining:

- **Major — `ENV` lookup in a keyword default.** *(done)* `SALCBases.jl`:
  the `check_irrep_unitary` keyword of `_projection_matrix_coupled_basis`
  no longer reads `ENV` on every call. The environment variable
  `MAGESTY_CHECK_IRREP_UNITARY` is read once at load time into the module
  constant `CHECK_IRREP_UNITARY_DEFAULT`, which is now the keyword default;
  the docstring documents this. An explicit keyword still overrides it.
- **Minor — indentation inconsistency.** *(done)* Leading tabs across
  `src/` and `test/` (50 files) normalized to 4 spaces per STYLE_GUIDE in a
  single mechanical commit. The change is whitespace-only (`git diff -w` is
  empty), so the parsed AST is unchanged and no benchmark re-baseline is
  needed despite touching hot-path files; the full suite (unit + integration)
  passes. `cli/` (a separate package), `tools/`, `examples/`, and `bench/`
  were left out of scope.
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

- **Major — O(n^2) translational-equivalence dedup.** *(prototyped,
  verified, reverted — no benefit at current sizes.)* The dedup uses
  `_is_translationally_equivalent_coupled_basis`, a symmetry-dependent
  predicate that is *not* structural equality (it has non-reflexive
  exclusions — same atom list or same first atom → not equivalent — and a
  float `coeff_tensor` comparison), so a plain `Base.hash`/`==` swap would
  change the merge result. The safe form keeps the predicate as the final
  arbiter and only prefilters candidates into buckets keyed by the
  integer-only necessary conditions `(Lf, Lseq, sorted (atom,l) pairs)`;
  this is provably result-identical (full suite passed, SALC counts and Jφ
  unchanged) and turns the scan into O(sum of bucket^2). But the profiler
  showed the dedup is not a bottleneck at the available fixtures
  (fept 16-atom: -2.5%; fege 64-atom: +1.0% — both within the ±2-4 ms noise
  floor; allocations -0.03..-0.4%). Reverted to avoid added complexity for
  no measured gain. The bucket-prefilter approach is the right one; revisit
  only if large-supercell profiling (e.g. 4x4x4) shows the orbit scan
  dominating.
- **Major — defensive buffered SH derivative variants.**
  `TesseralHarmonics.jl:549-`: the unbuffered `∂Zₗₘ_∂r̂{x,y,z}_unsafe` /
  `zzₗₘ_unsafe` allocate a `dnPl` vector per call. Not on the current hot
  path (the kernels use the buffered `∂ᵢZlm_unsafe`), but add buffered
  overloads to guard against a future hot-path caller.
- **Minor — `√2` as a module constant.** *(declined)* `@code_typed` shows
  `√2` is already constant-folded to the literal `1.4142135623730951` at
  compile time (`x * √2` lowers to `mul_float(x, 1.4142135623730951)`), so a
  module constant gives no runtime benefit. `√2` also mirrors the `√2` in the
  surrounding equation docstrings; replacing it with `SQRT2` would reduce
  readability. Not worth doing.
- **Minor — `zzₗₘ_unsafe` triple Legendre recompute.** *(declined / folded
  into the buffered-variants item)* The unbuffered `zzₗₘ_unsafe` and
  `∂Zₗₘ_∂r̂{x,y,z}_unsafe` are not called by the design-matrix kernels (the
  hot path uses the buffered `∂ᵢZlm_unsafe` / `Zₗₘ_grad_unsafe` in
  `Fitting.jl`); they are reached only through thin internal wrappers. The
  3× recompute is therefore off the hot path. If these gain a hot-path
  caller, address it together with the "defensive buffered SH derivative
  variants" item above.
- **Minor — `_canonicalize_eigenspace` explicit n×n projector.** *(deferred)*
  Replacing `P = V*V'` with on-demand columns `V*(V'eⱼ)` would avoid the n×n
  allocation and skip unused columns after the loop's early break, but it is
  setup-time only (once per degenerate eigenspace in SALC construction), the
  gain is limited unless n is large (large supercells), and it changes the
  BLAS reduction (gemm → gemv) so the canonical gauge shifts at the ULP
  level — a SALC linked site that needs numerical sign-off. Revisit only if
  large-supercell profiling shows the projector dominates.
- **Profiler recommended** — measure SALC-construction throughput
  (`construct_and_classify_coupled_basislist` /
  `projection_matrix_coupled_basis`) on a representative large supercell
  to confirm where the O(n^2) dedup and the projector dominate.

## API / UX

All API / UX findings are done (see "Already landed"): the `energy_kind`
validation, the SPEC argument-name fix, the `Magesty.save` / `load`
docstring sections, the docstring-section convention, and the remaining UX
minors (`SpinConfig` examples, `atol_unit_norm::Real`, `saxis` validation,
`intercept` / `nobs` / `dof` `# Returns`, the `Lasso` SPEC note). The
`AbstractEstimator` / `OLS` / `Ridge` "missing `# Arguments`" item was
declined with rationale recorded above.

## Confirmed clean (for the record)

Physics conventions (3 x n_atoms layout, unit spin vectors, real
tesseral `Zₗₘ`, `Jφ` / `j0` separation), linked-site synchronization (the
spherical-harmonics convention `Z = C·Y` verified to ~2e-16 across
l = 0..4; SALC key-group ordering; jphi serialized in `salc_index`
order), and numerical stability (Cholesky SPD guarantees, ridge
conditioning, GCV denominator saturation, RotationMatrix gimbal-lock
branch, MfaSampling kappa / zero-norm guards) were all verified consistent
— the lattice transpose (now fixed) was the single synchronization break.
