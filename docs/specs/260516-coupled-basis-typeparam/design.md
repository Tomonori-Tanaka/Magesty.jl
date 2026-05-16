# Design: CoupledBasis type parameterization

## Type parameter choice

Parameterize by **tensor rank** `R = ndims(coeff_tensor) = length(ls) +
1`. Single parameter, lets Julia store `coeff_tensor` as a fully
concrete `Array{Float64, R}`.

```julia
struct CoupledBasis{R}
    ls::Vector{Int}
    Lf::Int
    Lseq::Vector{Int}
    atoms::Vector{Int}
    coeff_tensor::Array{Float64, R}
end
```

**Why rank, not `Nsites`**: Julia does not allow arithmetic in field
type expressions (`Array{Float64, N+1}` is not a valid field type for
`struct CoupledBasis{N}`). Parameterizing by `R` directly avoids
needing two parameters or `@generated` constructors. Site count is
recoverable as `R - 1` and is checked at construct time via
`length(ls) == R - 1`.

Same shape for `CoupledBasis_with_coefficient{R}` and
`AngularMomentumCouplingResult{R}`.

## Constructor changes

User-facing constructor signature stays generic (accepts
`AbstractArray{<:Number}`) and converts to `Array{Float64,R}`
internally — `R` is inferred from `ndims(coeff_tensor)`:

```julia
function CoupledBasis(
    ls::AbstractVector{<:Integer},
    Lf::Integer,
    Lseq::AbstractVector{<:Integer},
    atoms::AbstractVector{<:Integer},
    coeff_tensor::AbstractArray{<:Number},
)
    R = ndims(coeff_tensor)
    tensor_concrete = convert(Array{Float64, R}, coeff_tensor)
    N = length(ls)
    # ...existing length / rank checks...
    return new{R}(
        collect(Int.(ls)),
        Int(Lf),
        collect(Int.(Lseq)),
        collect(Int.(atoms)),
        tensor_concrete,
    )
end
```

`convert(Array{Float64, R}, x)` is a no-op when `x` is already
`Array{Float64, R}` (common path from
`build_all_real_bases`/`coeff_tensor_complex`).

## Call-site annotation policy

- **Rank-agnostic containers** (`Vector{CoupledBasis}`,
  `SortedCounter{CoupledBasis}`, `IdDict{CoupledBasis, Int}`): use the
  UnionAll `Basis.CoupledBasis` without parameters. Julia treats this as
  `CoupledBasis where R`. Iteration uses runtime dispatch on `R`, which
  is unavoidable when ranks are mixed.
- **Hot-path inner-loop variables**: where a single rank dominates a
  call (e.g. inside `design_matrix_energy_element` after `cbc` is
  extracted), let Julia specialize via the loop body. Avoid
  `::CoupledBasis_with_coefficient` annotations that erase the
  parameter — prefer no annotation (so the parameter is inferred) or
  `::CoupledBasis_with_coefficient{R} where R`.
- **Function signatures**: where the function body uses
  `cbc.coeff_tensor` arithmetically, take the argument without
  rank-erasing annotation. Where the function just stores or forwards,
  use UnionAll.

## Container behavior verification

`Base.isless(::CoupledBasis, ::CoupledBasis)` already compares by
`length(ls)` first, then field-by-field. Mixed-rank instances in a
`SortedCounter` are well-ordered. No changes needed.

`Base.show` and `convert_to_coupled_basis` are unaffected by the type
parameter.

`reorder_atoms` calls `permutedims(coeff_tensor, dims_perm)` which
preserves rank → preserved `Array{Float64, R}`. Returns
`CoupledBasis{R}` with the same `R`. No changes needed beyond
threading the parameter through the constructor call.

## XML I/O

`xml_io.jl`:

- `write_xml` already uses `selectdim` over the trailing axis →
  rank-agnostic, no change.
- `build_sce_basis_from_xml` reconstructs via
  `zeros(Float64, full_shape...)` where `full_shape::Vector{Int}`. This
  returns `Array{Float64, N}` with `N` known at runtime, not compile
  time. Julia's `splat` on `Vector{Int}` yields `Array{Float64, N}`
  where `N` is inferred from the length at runtime — that's fine, the
  constructor accepts any `AbstractArray{<:Number}` and the type
  parameter `R` of the returned `CoupledBasis{R}` is correctly set
  from `ndims(coeff_tensor)`.

No XML schema changes.

## Linked sites that must stay in sync

Per `CLAUDE.md`'s "連動箇所" section, the SCE coefficient I/O loop
(write_xml / build_sce_basis_from_xml round-trip) must continue to
match integer-precision. Verified by
`test/component_test/test_save_load.jl` and
`test/examples/fept_tetragonal_2x2x2/test.jl`.

The Optimize ↔ SALCBasis (l, m, site) ordering is unaffected — we are
not touching basis order, only field type.

## Benchmarking

Before:

1. On `main` (cd72f68 or the commit just before the refactor),
   `make test-integration` once to populate the test cache.
2. Use a small script (or extend
   `test/benchmark_salcbasis_hotspots.jl`) to time
   `build_design_matrix_energy` on the fept example with `@time` (5
   trials each, report min and median).

After: same script post-refactor. Record in `.claude/bench_log.md`
under a "B1" section. Expectation: allocation count drops in the inner
contraction loop; wall time may or may not improve (depends on whether
`Any` arithmetic was being constant-folded).
