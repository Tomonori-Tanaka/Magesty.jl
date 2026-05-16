# Post-Step 7 cleanup follow-ups

Status: open (created 2026-05-15 during the SCE public-API breaking commit).

Tracks `[B]` items from the Step 7 audit that were intentionally deferred
from the breaking commit (because they need benchmarks or because they
expand scope), plus follow-ups discovered while writing Step 7.

Related: `docs/specs/260514-sce-public-api/` (the parent refactor).

## Hot-path type stability (audit B1–B4)

Highest-priority follow-up — these affect every fit and every predict
call. Land as separate `perf(...)` commits with before/after benchmarks
in `.claude/bench_log.md`.

- ~~**B1**: `src/types/Basis.jl:32, :310` — `CoupledBasis.coeff_tensor` and
  `CoupledBasis_with_coefficient.coeff_tensor` are typed as bare
  `AbstractArray`. Concretize…~~ *(resolved 2026-05-16 on
  `refactor/coupled-basis-typeparam`. Structs parameterized as
  `CoupledBasis{R}`, `CoupledBasis_with_coefficient{R}`,
  `AngularMomentumCouplingResult{R}` with `coeff_tensor::Array{Float64,
  R}`. See `docs/specs/260516-coupled-basis-typeparam/`.)*
- ~~**B2**: `src/Optimize.jl:design_matrix_energy_element` and `calc_∇ₑu!`
  — per-translation allocation of `Vector{Vector{Float64}}` for
  spherical-harmonic values.~~ *(per-translation aspect resolved
  earlier by `260516-optimize-workspace`'s workspace pooling; the
  residual `Vector{Vector{Float64}}` layout was flattened on
  2026-05-16. `EnergyWorkspace` / `GradWorkspace` now hold a single
  contiguous `sh_values::Vector{Float64}` plus
  `sh_offsets::Vector{Int}`, where site `i`'s slice is
  `sh_values[sh_offsets[i] + 1 : sh_offsets[i + 1]]`. Removes outer
  pointer indirection, makes site values cache-adjacent, and replaces
  the per-site `resize!` loop in `_ensure_sh_buffer!` with one
  cumulative-offset pass + at most one `Vector{Float64}` resize.
  Numerics unchanged (full unit + integration suite passes); time
  noise-level, allocs ~-60 per call. See `.claude/bench_log.md` "B2".)*
- ~~**B3**: same files — `cbc.coeff_tensor[idx_buf..., mf_idx]` splat
  indexing is not statically resolvable.~~ *(resolved 2026-05-16
  together with B1: `idx_buf::MVector{R-1,Int}`,
  `dims_t::NTuple{R-1,Int}`, `CartesianIndices` now have statically
  known rank. Energy/torque design-matrix builders gained ×24 / ×8.6
  on the fept fixture.)*

## Minor follow-ups

- ~~**B5**: `src/utils/xml_io.jl` — `Vector{Any}` of `Colon()` for
  tensor slicing~~ *(resolved 2026-05-16: rewritten as
  `selectdim(coeff_tensor, ndims(coeff_tensor), mf_idx) .= tensor_slice`
  — the trailing dimension is the only one indexed by an integer, and
  every other dim takes the full range)*.

## Documentation cleanup (audit C2)

- **C2** *(interim fix landed 2026-05-15)*: the placeholders
  `Reference: Equation (***) in T. Tanaka and Y. Gohda, ***` in the `Zₗₘ`
  (`src/utils/MySphericalHarmonics.jl:488`) and `∂Zₗₘ_∂r̂x`
  (`:574`) docstrings now point to the technical notes
  (`https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/`). Swap
  in the published paper reference once it is out.

## SALC build determinism across platforms (resolved)

*Resolved 2026-05-15 in `src/SALCBases.jl::_canonicalize_eigenspace`.*

Degenerate eigenvalue subspaces of the SALC projection matrix used to
produce LAPACK-implementation-dependent eigenbases — fept and fege
both have `Lf=2` groups that span 2D degenerate eigenspaces, so the
macOS-generated XML baselines did not match Linux x86 builds. The fix
applies modified Gram-Schmidt on `P * eⱼ` (axis-ordered) where
`P = V Vᵀ` is the basis-invariant projector onto the eigenspace,
yielding a canonical orthonormal basis of the same subspace. SALC
coefficients are now bit-stable across BLAS/LAPACK implementations,
and the "fresh build" byte-diff regression test in
`test/component_test/test_save_load.jl` is restored for fept and
fege.

## Variable naming vs. mathematical content *(resolved 2026-05-16)*

A recurring smell exposed during the Step 7 audit: variables named after
one mathematical object actually hold a different one. Misleads readers
and silently invites incorrect changes (e.g. checking unitarity on a
"projection matrix" only makes sense if it's actually a representation).

Resolved mismatches:

- `src/SALCBases.jl:projection_matrix_coupled_basis::temp_projection_mat`
  → renamed to `representation_mat`.
- Same function, `check_unitary` kwarg → renamed to
  `check_irrep_unitary`; env var renamed to
  `MAGESTY_CHECK_IRREP_UNITARY`.

Broader sweep across `src/` (SALCBases, Optimize, Symmetries, Clusters,
Magesty, types/Basis, utils/xml_io) found no further mathematically
misleading `*_mat` / `*_matrix` / `*_list` / `*_dict` / `*_set` / `*_vec`
names — the remaining names accurately describe their runtime values.
Two non-naming items surfaced during the sweep and were fixed alongside:
the `interactiong_atoms` typo in `src/Clusters.jl`, and a redundant
`sort()` on an already-sorted vector in `irreducible_clusters`.

## Out-of-scope reminders

- ~~`tools/check_convergence_embset.jl`, `tools/convert2tensor.jl`,
  `tools/micromagnetics.jl`, `tools/plot_jphi_cluster_distance.jl`, and
  `tools/personal/*` still reference legacy symbols~~ *(resolved
  2026-05-16: `check_convergence_embset.jl` was deleted (recreate when
  needed); `micromagnetics.jl` was migrated to the `SCEBasis` API and
  now loads the basis from the same XML the coefficients live in.
  `convert2tensor.jl` and `plot_jphi_cluster_distance.jl` turned out
  to use only `Structure` / `Symmetry`, which survived Step 7 — they
  needed no migration. `tools/personal/*` is gitignored and migrated by
  the owning user when needed. `tools/CrossValidation.jl` (LOOCV
  module) was deleted in an earlier pass.)*
- ~~`docs/src/tools.md` still mentions the legacy `System` struct and
  `system.jld2` caching pattern~~ *(resolved 2026-05-16 alongside the
  tools migration.)*
- ~~`test/benchmark_*.jl` (except `benchmark_salcbasis_hotspots.jl`) and
  `test/profile_run.jl` still reference legacy symbols~~ *(resolved
  2026-05-16: deleted `benchmark_optimize.jl`,
  `benchmark_optimize_hotspots.jl`,
  `benchmark_optimize_sphericart.jl`, `benchmark_threads.jl`,
  `profile_run.jl`; the `bench-optimize-sphericart{,-fept}` Makefile
  targets and the `bench_optimize_sphericart` TEST_MODE branch were
  removed in the same pass. `benchmark_salcbasis_hotspots.jl`,
  `benchmark_sphericart.jl`, and `benchmark_spherical_harmonics.jl`
  remain as the modern templates; new perf work creates fresh benches
  from those.)*
