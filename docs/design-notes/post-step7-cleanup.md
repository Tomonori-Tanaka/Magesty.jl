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

- **B1**: `src/types/Basis.jl:32, :310` — `CoupledBasis.coeff_tensor` and
  `CoupledBasis_with_coefficient.coeff_tensor` are typed as bare
  `AbstractArray`. Concretize to at least `AbstractArray{Float64}` (or
  parameterize the structs by tensor rank so a concrete `Array{Float64,
  N}` is possible). Pervasive type-stability win across the whole
  design-matrix path.
- **B2**: `src/Optimize.jl:design_matrix_energy_element` (around L424) and
  `calc_∇ₑu!` (around L611) — per-translation allocation of
  `Vector{Vector{Float64}}` for spherical-harmonic values, and `@views`
  slices that don't go through `SVector`. Preallocate (`MMatrix` or
  fixed-size scratch) and lift slices into `SVector{3, Float64}`.
  Mirrors the pattern already used in `build_design_matrix_torque`
  (`dir_iatom_svec = SVector{3, Float64}(dir_iatom)`).
- **B3**: same files — `cbc.coeff_tensor[idx_buf..., mf_idx]` splat
  indexing is not statically resolvable. Try `CartesianIndex`-based
  indexing or `@generated` dispatch by rank.
- **B4**: `src/Magesty.jl:SCEFit.metrics` — `Dict{Symbol, Any}` →
  `Dict{Symbol, Float64}` (only four `Float64` values stored).

## I/O performance (audit B6)

- **B6**: `src/utils/xml_io.jl` — `read_model_components_from_xml` /
  `read_basis_components_from_xml` / `read_salcbasis_from_xml` /
  `_read_tolerance_sym` / `_read_isotropy` each call `readxml(path)`
  independently. A single large `load(SCEModel, path)` re-parses the XML
  up to four times. Split the internals into `(::EzXML.Document, ...)`
  overloads and parse the document once at the public-entry boundary.

## Minor follow-ups

- **B5**: `src/utils/xml_io.jl:394` — `Vector{Any}` of `Colon()` for
  tensor slicing. Cosmetic; the path is I/O-bound so perf impact is
  small. Replace with chained `selectdim`.
- **B14**: `src/utils/MySphericalHarmonics.jl:735–736` — `d_Zlm` /
  `d_Zlm_unsafe` are exported, mutable (no `const`), and typed implicitly
  as `Vector{Any}`. They are exported, so type-instability here leaks
  into user code. At minimum, `const ... ::Vector{Function} = [...]`;
  ideally an `NTuple` of concrete function types.

## Documentation cleanup (audit C2)

- **C2** *(interim fix landed 2026-05-15)*: the placeholders
  `Reference: Equation (***) in T. Tanaka and Y. Gohda, ***` in the `Zₗₘ`
  (`src/utils/MySphericalHarmonics.jl:488`) and `∂Zₗₘ_∂r̂x`
  (`:574`) docstrings now point to the technical notes
  (`https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/`). Swap
  in the published paper reference once it is out.

## SALC build determinism across platforms (degenerate eigenvalues)

The "fresh build" byte-diff regression test in
`test/component_test/test_save_load.jl` was dropped during the Step 7
push because both `fept` and `fege` have 2-body `Lf=2` SALC groups that
span 2D degenerate eigenvalue subspaces. LAPACK eigensolvers do not
promise a platform-stable basis within a degenerate subspace, so the
committed XML baseline (generated on macOS arm64) does not match the
output on Linux x86 — the values differ by a unitary rotation inside the
degenerate subspace (verified: the rotated entries satisfy
`a² + b² = const`).

Per-platform reproducibility is still verified by the
"basis identity (design matrices) after reload" testset (it compares
the energy/torque design matrices with `≈`, which is invariant under
rotations within a degenerate subspace).

Two possible follow-up fixes:

- **Subspace-span comparison**: group SALCs by `(atom-orbit, ls, Lf)` and
  compare the orthogonal projection onto each group's span (or the QR
  factorization's `R` block). Adds a helper to the test suite without
  touching `src/`. Lets the byte-diff "fresh build" test come back in a
  span-aware form.
- **Canonical gauge fix at SALC build time**: post-process eigenvectors
  within each degenerate subspace to a deterministic representative
  (e.g., sort by first-nonzero entry magnitude, sign-fix on a chosen
  pivot, then Gram-Schmidt). Makes the SALC output bit-stable across
  LAPACK implementations and reinstates simple byte-diff regression at
  the cost of changing a numerical convention (and regenerating
  baselines). Treat the gauge convention as a physics decision — flag
  it for review before landing.

The example that triggered this (fept):

```xml
<SALC index="6" num_basis="2" body="2" Lf="2">
  <basis multiplicity="2" atoms="9 11"  ls="1 1" ... >0.0 0.0 0.5076 0.0  0.4923</basis>
  <basis multiplicity="2" atoms="9 13"  ls="1 1" ... >0.0 0.0 0.5076 0.0 -0.4923</basis>
</SALC>
<SALC index="7" num_basis="2" body="2" Lf="2">
  <basis multiplicity="2" atoms="9 11"  ls="1 1" ... >0.0 0.0 0.4923 0.0 -0.5076</basis>
  <basis multiplicity="2" atoms="9 13"  ls="1 1" ... >0.0 0.0 0.4923 0.0  0.5076</basis>
</SALC>
```

SALC index 6 and 7 span a 4D subspace (2 groups × 2 basis each) that is
4-fold degenerate; any orthogonal rotation within it is a valid
eigenbasis.

## Variable naming vs. mathematical content

A recurring smell exposed during the Step 7 audit: variables named after
one mathematical object actually hold a different one. Misleads readers
and silently invites incorrect changes (e.g. checking unitarity on a
"projection matrix" only makes sense if it's actually a representation).

Known mismatches:

- `src/SALCBases.jl:projection_matrix_coupled_basis::temp_projection_mat`
  — holds the single-symmetry-operation representation matrix `D(g)`
  (Wigner D × `tensor_inner_product` phase), **not** the projector. The
  projector is the accumulated sum after `/2nsym`. Rename candidates:
  `representation_mat`, `D_g_mat`.
- Same function, `check_unitary` kwarg — the name does not say *what*
  must be unitary. The check tests "each `D(g)` is unitary as required
  by a unitary irrep", not "`P` is unitary" (which would imply `P = I`).
  Rename the kwarg or the helper (e.g. `check_irrep_unitary`) and update
  the docstring.

Broader sweep task: walk `src/` and flag `*_mat` / `*_matrix` / `*_list`
/ `*_dict` variables whose name implies one mathematical object but
stores another. Output a list of rename candidates plus the scope of
caller-site updates each rename would require.

## Out-of-scope reminders

- `tools/CrossValidation.jl`, `tools/check_convergence_embset.jl`,
  `tools/personal/*` still reference legacy symbols (`System`,
  `SpinCluster`, `fit_sce_model`, `_fit_sce_model_internal`,
  `get_j0_jphi`, `calc_energy`). These were intentionally left untouched
  during Step 7 — re-evaluate each tool's usefulness before migrating or
  deleting.
- `docs/src/tools.md` still mentions the legacy `System` struct and
  `system.jld2` caching pattern. Update alongside the tools migration
  above so the tool docs match the live scripts.
- `test/benchmark_*.jl` (except `benchmark_salcbasis_hotspots.jl`) and
  `test/profile_run.jl` still reference legacy symbols. Migrate when
  someone needs them next; the patterns are all `System(input)` →
  `SCEBasis(input)` and `SpinCluster(system, input)` →
  `SCEDataset(basis, ...)` + `fit(SCEFit, ...)`.
