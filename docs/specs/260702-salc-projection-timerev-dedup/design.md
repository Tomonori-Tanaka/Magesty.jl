# Design: SALC projection time-reversal dedup

Status: implemented (2026-07-02)

## Summary

Build the two per-operation representation matrices `D(g, trs=false)`
and `D(g, trs=true)` in a **single pass** over `(cb1, cb2)`, instead of
two full passes. The loop over `time_rev_sym in (false, true)`
disappears; for each symmetry operation `n` the body computes the
shift/fold/reorder and the `cb2` phase once, then writes each block
twice — once with `factor = phase` into the `false` buffer and once
with `factor = multiplier * phase` into the `true` buffer, where
`multiplier = (-1)^sum(cb1.ls)` (permutation-invariant, so the
reordered basis is not needed for the sign).

The alternative — keeping the two-pass loop and caching per-`(i, j)`
phases in a `nbasis × nbasis` matrix — was rejected: it saves the same
recomputation but keeps two passes over the block-write loop and needs
the same extra memory, with more bookkeeping.

Cost: one extra `full_matrix_dim²` `Matrix{Float64}` scratch buffer
(the `true`-variant twin of the existing `representation_mat`). For the
observed group sizes (`full_matrix_dim` up to a few hundred) this is
negligible next to the current per-call allocation traffic it removes.

## Module layout

| Target | Change |
|---|---|
| `src/SALCBases.jl` | Restructure the loop body of `_projection_matrix_coupled_basis` (currently lines ~862-900). No signature change. |
| `test/component/` | New bit-identity regression test (see Test strategy). |
| `.claude/bench_log.md` | Before/after benchmark entry. |

## API

No public or internal API change. `_projection_matrix_coupled_basis`
keeps its exact signature and return type:

```julia
function _projection_matrix_coupled_basis(
    coupled_basislist::SortedCounter{CoupledBases.CoupledBasis},
    symmetry::Symmetry;
    check_irrep_unitary::Bool = CHECK_IRREP_UNITARY_DEFAULT,
)::Matrix{Float64}
```

## Restructured loop (sketch)

```julia
representation_mat_no_time_rev = zeros(Float64, full_matrix_dim, full_matrix_dim)
representation_mat_time_rev   = zeros(Float64, full_matrix_dim, full_matrix_dim)
for (n, symop) in enumerate(symmetry.symdata)
    fill!(representation_mat_no_time_rev, 0.0)
    fill!(representation_mat_time_rev, 0.0)
    base_rot_mat = base_rot_mats[n]
    for (i, cb1) in enumerate(coupled_basislist)
        # -- time_rev_sym-independent work, now computed once --
        @inbounds for k = 1:N_atoms
            atoms_shifted_list[k] = symmetry.map_sym[cb1.atoms[k], n]
        end
        primitive_atoms = _find_translation_atoms(atoms_shifted_list,
                                                  cluster_atoms, symmetry)
        reordered_cb = reorder_atoms(cb1, primitive_atoms)
        multiplier = (-1)^sum(cb1.ls)   # == (-1)^sum(reordered_cb.ls)
        col_range = ((i-1)*submatrix_dim+1):(i*submatrix_dim)
        for (j, cb2) in enumerate(coupled_basislist)
            _is_obviously_zero_coupled_basis_product(reordered_cb, cb2) && continue
            phase = _tensor_inner_product(cb2.coeff_tensor,
                                          reordered_cb.coeff_tensor) / (2*Lf+1)
            row_range = ((j-1)*submatrix_dim+1):(j*submatrix_dim)
            @views representation_mat_no_time_rev[row_range, col_range] .=
                base_rot_mat .* phase
            @views representation_mat_time_rev[row_range, col_range] .=
                base_rot_mat .* (multiplier * phase)
        end
    end
    # unitarity check + accumulation, in the original order (false, then
    # true). The check is the existing inline pattern
    # `check_irrep_unitary && !_is_unitary(mat, tol = 1e-8) && error(...)`
    # (with the current D(g)'D(g) diagnostic), applied to each buffer
    # before it is accumulated; shown abbreviated here.
    check_irrep_unitary && _assert_unitary(representation_mat_no_time_rev, n)  # illustrative
    projection_mat .+= representation_mat_no_time_rev
    check_irrep_unitary && _assert_unitary(representation_mat_time_rev, n)   # illustrative
    projection_mat .+= representation_mat_time_rev
end
```

## Types and conventions

- No physics-convention change. The projector is still the average of
  the same `2*nsym` representation matrices divided by `2*nsym`.
- **Bit-identity argument.** Three properties guarantee an unchanged
  result:
  1. Block values: the current code computes
     `base_rot_mat .* (multiplier * phase)` with
     `multiplier = 1` (`false` case) or `(-1)^sum(ls)` (`true` case).
     The `false`-buffer write `base_rot_mat .* phase` is the identical
     expression with `multiplier` folded away
     (`1 * phase` and `phase` are the same Float64); the `true`-buffer
     write is verbatim (`multiplier` is an exactly representable ±1, so
     `multiplier * phase` is exact and order-insensitive).
  2. Block writes are assignments (`.=`), not accumulations, so
     within-`D(g)` construction order is irrelevant.
  3. Accumulation into `projection_mat` keeps the original order:
     for each `n`, `false` before `true`.
- `multiplier` via `sum(cb1.ls)` instead of `sum(reordered_cb.ls)`:
  `reorder_atoms` permutes `ls`; integer sums are permutation-exact.
- **Scope of the bit-identity claim.** The argument above covers
  `projection_mat` (pure elementwise arithmetic). Downstream SALC
  coefficients pass through `eigen!`, and the Wigner matrices through
  trig functions — neither is bit-stable across BLAS/libm
  implementations or architectures. Bit-identity is therefore verified
  as a one-time before/after comparison on a single machine; committed
  regression tests use `isapprox(atol = 1e-8)` per the existing
  `test_save_load.jl` precedent. `salc_fingerprint` hashes structural
  integers only (no floating-point payload), so its equality is a
  structural check, not numerical evidence.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): N/A.
- [ ] SCE coefficient XML (`save` / `load`): N/A (coefficients
      unchanged by construction; guarded by the regression coverage).
- [x] `Fitting` <-> `SALCBasis`: SALC ordering and coefficients must be
      unchanged — this is the central invariant of the spec.
- [ ] `.claude/agents/` references: N/A.
- [ ] `SPEC.md` / `docs/src/api.md` updates: N/A (internal change).

## Test strategy

- **One-time bit-identity gate (development machine, not committed)**:
  capture `projection_mat` outputs and full `SALCBasis` coefficients on
  the dimer and `fege_2x2x2` fixtures with the current implementation,
  re-capture after the change on the same machine, and compare with
  `==`. Any bit difference means stop and consult — do not proceed by
  loosening the comparison. The outcome is recorded in the
  `.claude/bench_log.md` entry.
- **Committed regression test**: extend
  `test/component/test_SALC_projection.jl` (which already exercises
  `_projection_matrix_coupled_basis` directly) with a baseline probe
  comparing a projection matrix against a committed baseline under
  `test/component/baselines/`, using `isapprox(atol = 1e-8)` — the
  tolerance convention of the existing fresh-build-vs-baseline test in
  `test_save_load.jl` (eigensolver/trig results are not bit-stable
  across BLAS/libm/architectures, so `==` would be flaky on CI).
  `salc_fingerprint` is compared with `==` (structural integers only).
- Existing suites: `make test-all`, `make test-jet`, `make test-aqua`.
- Benchmark: `make bench-salcbasis` /
  `bench/benchmark_salcbasis_hotspots.jl` (references verified against
  the current internal API; commit 1417c04) before and after; record
  in `.claude/bench_log.md`.

## Risks and open items

- **Numerical-review sign-off required before merge**: reordering
  floating-point work is exactly the class of change the review panel
  flags; the bit-identity argument above plus the one-time `==` gate
  and the committed regression probe are the evidence.
- The extra `full_matrix_dim²` buffer doubles the scratch memory of
  this function. If a pathological group makes this matter, fall back
  to the per-`(i, j)` phase-cache alternative (same dedup, one
  `nbasis × nbasis` cache instead of a second full matrix).
- Deferred follow-up (separate decision after re-measuring): a
  low-allocation `reorder_atoms` variant for this call site
  (`src/CoupledBases.jl`), which trades a second specialized path
  against the single-generic-path preference.
