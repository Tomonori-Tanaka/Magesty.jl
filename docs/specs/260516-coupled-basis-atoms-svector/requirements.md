# Requirements: CoupledBasis.atoms as SVector

Created 2026-05-16. Topic branch: `refactor/coupled-basis-atoms-svector`.

## Purpose

Convert the `atoms` field on the coupled-basis struct family from
`Vector{Int}` to `SVector{R-1, Int}`, eliminating heap-allocated
`Vector{Int}` traffic in symmetry-mapping and indexing hot paths.

`R` is the tensor rank parameter already present after the B1+B3
refactor (`coeff_tensor::Array{Float64, R}`); since
`length(atoms) == length(ls) == R - 1`, no new type parameter is
required.

## Scope

In scope (3 structs in `src/types/Basis.jl`):

- `CoupledBasis{R}.atoms::Vector{Int}` → `SVector{R-1, Int}`
- `CoupledBasis_with_coefficient{R}.atoms::Vector{Int}` → `SVector{R-1, Int}`
- `AngularMomentumCouplingResult{R}.atoms::Vector{Int}` → `SVector{R-1, Int}`

Inner constructors: accept `AbstractVector{<:Integer}` and convert to
`SVector{R-1, Int}` (R inferred from `coeff_tensor` ndims).

`reorder_atoms` helpers: accept `AbstractVector{<:Integer}`, produce a
new struct with `SVector{R-1, Int}` atoms.

Callers (verified by grep, transparent to the change):

- `src/SALCBases.jl`: `cb1.atoms != cb2.atoms`, indexing, zipping with
  `cb1.ls` — all work on `AbstractVector` API.
- `src/Optimize.jl`: `cbc.atoms[site_idx]` — pure indexing.
- `src/utils/xml_io.jl`: `length(cbc.atoms)`, `string.(cbc.atoms)` —
  read-only.

Hot spot the change targets:

- `projection_matrix_coupled_basis` (`SALCBases.jl:756`):
  `atoms_shifted_list[k] = symmetry.map_sym[cb1.atoms[k], n]` — the
  buffer `atoms_shifted_list` was hoisted earlier (former #10 quick fix).
  With `SVector` atoms, the natural rewrite is
  `atoms_shifted = SVector{R-1,Int}(symmetry.map_sym[a, n] for a in cb1.atoms)`,
  fully stack-allocated.

Out of scope:

- B2 (sh_values flat storage).
- Renaming `atoms` or changing semantics.
- Any public-API signature change beyond the breaking constructor /
  field type.

## Invariants

1. **Numerical equivalence**: every test passes bit-stable; XML output
   byte-identical to baseline.
2. **Type parameter unchanged externally**: structs remain
   `CoupledBasis{R}` etc.; `R = ndims(coeff_tensor)` still infers
   from the caller-supplied tensor.
3. **Constructor ergonomics**: callers passing `Vector{Int}` (the
   majority of in-tree call sites) continue to work — the inner
   constructor accepts `AbstractVector{<:Integer}` and converts.

## Completion criteria

- `make test-all`, `make test-aqua`, `make test-jet` green.
- XML byte-identical (XML roundtrip + fept fixture).
- Alloc reduction recorded for `projection_matrix_coupled_basis` (the
  primary target) and overall `build_sce_basis` / design-matrix build
  on the fept fixture, in `.claude/bench_log.md`.

## Risks

- **Constructor breakage**: any external caller building these structs
  directly will need to pass a vector that `SVector{R-1,Int}` can
  convert (a `Vector{Int}` of the right length suffices).
- **`reorder_atoms` length mismatch**: must validate that
  `length(new_atoms) == R - 1` (assert inside constructor).
- **Show / XML output formatting**: `join(string.(svec), " ")` works
  identically for `SVector`; verified by integration test.
