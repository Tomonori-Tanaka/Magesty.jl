# Requirements: SALC projection time-reversal dedup

Status: complete (2026-07-02)

## Goal

Remove the redundant recomputation across the `time_rev_sym in
(false, true)` loop in `_projection_matrix_coupled_basis`
(`src/SALCBases.jl`), roughly halving the allocation count and the
CPU-dominant tensor-contraction work of the historically most expensive
stage of `SALCBasis` construction, while keeping the output
bit-identical.

## Background

A whole-package performance review found that for each `(symop, cb1)`
pair the loop body recomputes work that does not depend on
`time_rev_sym` at all: the shifted atom list, the translation fold
(`_find_translation_atoms`), the atom reorder (`reorder_atoms`), and
the entire inner `cb2` loop (`_is_obviously_zero_coupled_basis_product`
+ `_tensor_inner_product`). Only the final sign
(`multiplier = time_rev_sym ? (-1)^sum(ls) : 1`) differs, and the sign
does not even need the reordered basis, because `sum(reordered_cb.ls) ==
sum(cb1.ls)` (summing a permuted vector is permutation-invariant).

Measured on `test/integration/fege_2x2x2` (12-basis `Lf = 0` key group,
`nsym = 96`): `reorder_atoms` alone costs ~1.09 KiB / 28 allocations per
call and is invoked `nsym × 2 × nbasis = 2304` times in that one group,
accounting for ~59% of `_projection_matrix_coupled_basis`'s 111,729
allocations / 4.15 MiB for the group. An `Lf = 2` group doubles the
total, consistent with the duplicated `cb2` loop scaling with
`submatrix_dim²`. This is the same redundancy class as the
`Δl`/`base_rot_mats` hoist recorded in `.claude/bench_log.md` entry #9,
which missed the reorder/fold/inner-product step.

## Scope

Includes:

- Restructuring the `(symop, time_rev_sym)` loop in
  `_projection_matrix_coupled_basis` so the `time_rev_sym`-independent
  work runs once per `(symop, cb1)`.
- A bit-identity regression test (before/after SALC coefficients and
  `salc_fingerprint`).
- A before/after entry in `.claude/bench_log.md` using
  `bench/benchmark_salcbasis_hotspots.jl`.

Excludes:

- A low-allocation specialized `reorder_atoms` variant
  (`src/CoupledBases.jl`). Deferred: re-measure after this spec lands
  and open a follow-up only if the remaining allocation floor still
  matters (it also trades against the single-generic-path preference).
- Any change to the SALC/CG conventions, the projection average
  (`/ (2*nsym)`), or the unitarity check semantics.
- Threading changes.

## Invariants

- **Bit-identical output on the development machine**: `projection_mat`
  (and every downstream SALC coefficient) must match the current
  implementation exactly when before/after snapshots are captured on
  the same machine. Floating-point accumulation order into
  `projection_mat` must be preserved: for each symmetry operation `n`,
  the `time_rev_sym = false` matrix is added before the `true` matrix,
  and block values are computed with the same product order
  (`base_rot_mat .* (multiplier * phase)`).
  This is a one-time development gate, not a committed test: the SALC
  pipeline goes through `eigen!` and trig-based Wigner matrices, which
  are not bit-stable across BLAS/libm implementations, so committed
  baselines follow the repository's `isapprox(atol = 1e-8)` convention
  (see `test/component/test_save_load.jl`).
  `salc_fingerprint` equality is additionally required, but it hashes
  only structural integers (it deliberately excludes the
  floating-point coefficients), so it is a structural sanity check,
  not evidence of numerical identity.
- The per-operation unitarity check (`check_irrep_unitary`) still runs
  on each of the `2*nsym` representation matrices `D(g)` before
  averaging.
- Physics conventions unchanged (tesseral `Zₗₘ`, SALC/CG, `3×n_atoms`
  layout). XML round-trips unchanged.

## Completion criteria

- [ ] One-time bit-identity gate passes on the development machine:
      `projection_mat` outputs and full SALC coefficients captured
      before the change compare `==` after the change, on at least the
      dimer and `fege_2x2x2` fixtures. Result recorded in the
      bench-log entry.
- [ ] Committed regression coverage: a projection-level baseline probe
      (extending `test/component/test_SALC_projection.jl`) compares
      against a committed baseline with `isapprox(atol = 1e-8)`,
      consistent with the existing fresh-build-vs-baseline test in
      `test_save_load.jl`; `salc_fingerprint` matches exactly.
- [ ] `make test-all` passes; `make test-jet` / `make test-aqua` clean.
- [ ] `.claude/bench_log.md` entry shows the expected reduction
      (target: roughly −50% allocations and ≥ 1.3× median speedup in
      the `_projection_matrix_coupled_basis` benchmark on
      `fege_2x2x2`; the `SALCBasis` constructor benchmark must not
      regress).

## References

- Related issues / PRs: none yet.
- Related specs / design notes: `.claude/bench_log.md` entry #9
  (`base_rot_mats` hoist); `docs/specs/260601-src-refactor/` (M2 SALC
  allocation work).
