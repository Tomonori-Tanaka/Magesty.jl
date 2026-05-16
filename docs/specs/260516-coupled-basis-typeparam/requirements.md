# Requirements: CoupledBasis type parameterization (B1)

Created 2026-05-16. Topic branch: `refactor/coupled-basis-typeparam`.

## Purpose

Concretize `coeff_tensor` storage in `Basis.CoupledBasis` and
`Basis.CoupledBasis_with_coefficient` to recover compile-time type
information across the entire design-matrix path. This is the highest
priority follow-up from the Step 7 audit (B1 in
`docs/design-notes/post-step7-cleanup.md`).

Current state (`src/types/Basis.jl:32, :310`):

```julia
coeff_tensor::AbstractArray
```

Element type is unknown, rank is unknown. Reads such as
`cbc.coeff_tensor[idx_buf..., mf_idx]` therefore return `Any`, which
prevents downstream type-stable arithmetic in
`design_matrix_energy_element` and `calc_∇ₑu!`.

## Scope

In scope:

- `Basis.CoupledBasis` — parameterize over tensor rank, concretize
  `coeff_tensor::Array{Float64, R}`.
- `Basis.CoupledBasis_with_coefficient` — same.
- `Basis.AngularMomentumCouplingResult` — same (currently
  `coeff_tensor::Array{Float64}` with unspecified rank).
- All construction sites in `src/` (SALCBases, Optimize, xml_io,
  AngularMomentumCoupling).
- Type annotations on function signatures that mention `CoupledBasis*`
  (use UnionAll `CoupledBasis` when rank-agnostic; constrain when
  rank-specific is appropriate).
- Tests under `test/component_test/`, `test/examples/`.

In scope (added 2026-05-16 after Step 2 measurement showed B1 alone
left wall-time unchanged because splat indexing still returns `Any`):

- B3: replace runtime-sized scratch buffers in
  `design_matrix_energy_element` and `calc_∇ₑu!` with `MVector{N,Int}`
  and `NTuple{N,Int}` (with `N = R - 1` derived from the type
  parameter). This lets `cbc.coeff_tensor[idx_buf..., mf_idx]` resolve
  statically once the field type is concrete.

Out of scope (separate follow-ups):

- B2: `sh_values::Vector{Vector{Float64}}` consolidation. The inner
  vector lengths vary per site (`2lᵢ+1`), so this needs a different
  layout (flattened storage + offsets, or `NTuple{N,Vector{Float64}}`).
- `atoms::Vector{Int}` → `SVector{N,Int}` (旧 #10).
- XML on-disk format changes (the format is rank-agnostic already).

## Invariants (must hold)

1. **Numerical equivalence**: every existing test passes with bit-stable
   outputs (no recomputation, no reordering of reductions). SCE
   coefficients and design matrices must match pre-refactor values to
   machine precision.
2. **XML round-trip byte stability**: `write_xml` /
   `build_sce_basis_from_xml` continue to round-trip identically. Use
   `selectdim` (already in place) so the XML schema is unchanged.
3. **Public API compatibility at the call-site level**: callers that
   write `::CoupledBasis` (without type parameters) continue to compile
   — Julia's UnionAll dispatch covers this. We do *not* require the
   struct signature to be non-breaking: this is a pre-public refactor
   and a struct-signature change is acceptable.
4. **Containers (`SortedCounter`, `IdDict`, `Vector{CoupledBasis}`)**
   continue to work with mixed-rank elements via the UnionAll element
   type.

## Completion criteria

- `make test-all` green.
- `make test-aqua` green (no new ambiguities).
- `make test-jet` shows no new type-instability flags inside
  `design_matrix_energy_element` / `calc_∇ₑu!` related to `coeff_tensor`
  access.
- Before/after `@time` (5+ trials) of `build_design_matrix_energy` on
  the fept example recorded in `.claude/bench_log.md`. Expected
  allocation reduction in the inner contraction loop.
- `docs/design-notes/post-step7-cleanup.md` B1 marked resolved.
- One breaking commit on `refactor/coupled-basis-typeparam`, merged via
  PR to `main` after review.

## Risks

- **Mixed-rank containers**: `SortedCounter{CoupledBasis}` will hold
  the UnionAll. Need to verify that `isless` still works correctly
  across ranks (it already compares by `length(ls)` first, so this
  should be safe).
- **XML reconstruction path** (`xml_io.jl:320`) creates the tensor with
  `zeros(Float64, full_shape...)` where `full_shape` is a `Vector{Int}`
  — this yields `Array{Float64,N}` where N is runtime-derived. The
  constructor needs to accept this without explicit rank narrowing
  (Julia infers).
- **`reorder_atoms` uses `permutedims`**: returns same rank, type
  preserved. Should be safe.
- **AngularMomentumCoupling internals**: `coeff_tensor_complex` and
  `build_all_real_bases` produce tensors of varying rank. Constructor
  must remain dispatch-flexible.
