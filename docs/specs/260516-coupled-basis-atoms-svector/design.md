# Design: CoupledBasis.atoms as SVector

## Field change

```julia
# src/types/Basis.jl

struct CoupledBasis{R}
    atoms::SVector{R-1, Int}     # was Vector{Int}
    ls::Vector{Int}
    coeff_tensor::Array{Float64, R}
    # ...
end

struct CoupledBasis_with_coefficient{R}
    atoms::SVector{R-1, Int}     # was Vector{Int}
    # ...
end

struct AngularMomentumCouplingResult{R}
    atoms::SVector{R-1, Int}     # was Vector{Int}
    # ...
end
```

## Inner constructors

Accept `AbstractVector{<:Integer}`, convert to `SVector{R-1, Int}`:

```julia
function CoupledBasis(
    atoms::AbstractVector{<:Integer},
    ls::AbstractVector{<:Integer},
    coeff_tensor::AbstractArray,
    # ...
)
    R = ndims(coeff_tensor)
    @assert length(atoms) == R - 1
    @assert length(ls) == R - 1
    atoms_s = SVector{R-1, Int}(atoms)
    ls_v = Vector{Int}(ls)
    tensor = convert(Array{Float64, R}, coeff_tensor)
    return new{R}(atoms_s, ls_v, tensor, ...)
end
```

Same pattern for `CoupledBasis_with_coefficient` and
`AngularMomentumCouplingResult`.

## `reorder_atoms`

```julia
function reorder_atoms(cb::CoupledBasis{R}, new_atoms::AbstractVector{<:Integer}) where {R}
    @assert length(new_atoms) == R - 1
    # ... existing permutation logic ...
    return CoupledBasis(new_atoms, new_ls, new_tensor, ...)
end
```

The inner constructor handles `SVector` conversion; `reorder_atoms`
just forwards `new_atoms` as a vector.

## Hot-path rewrite (`SALCBases.jl:756`)

Before (current, post former #10 quick fix):

```julia
atoms_shifted_list = Vector{Int}(undef, length(cb1.atoms))  # hoisted outside loop
for n in 1:num_symops
    for k in 1:length(cb1.atoms)
        atoms_shifted_list[k] = symmetry.map_sym[cb1.atoms[k], n]
    end
    # ... use atoms_shifted_list
end
```

After:

```julia
for n in 1:num_symops
    atoms_shifted = SVector{R-1, Int}(symmetry.map_sym[a, n] for a in cb1.atoms)
    # ... use atoms_shifted (stack-allocated)
end
```

The exact local rewrite depends on the surrounding loop; the
hoisted buffer pattern goes away.

## Callers (verified transparent)

- `SALCBases.jl:553`: `cb1.atoms != cb2.atoms` — `SVector` defines `==`
  element-wise, works.
- `SALCBases.jl:838, :853`: `zip(cb1.atoms, cb1.ls)`, `cb1.atoms` as
  iterable — fine.
- `Optimize.jl:311, :504`: `cbc.atoms[site_idx]` — indexing fine.
- `types/Basis.jl:116`: `lc1.atoms < lc2.atoms` — `SVector` lex compare,
  fine.
- `types/Basis.jl:410`: `print(io, "atoms=$(cbc.atoms), ")` — `show` of
  `SVector` differs slightly from `Vector` (prefix `[…]` vs `Int[…]`).
  If a test asserts exact `Vector{Int}` show output, adjust there.
- `xml_io.jl:133, :149, :415`: `length`, `join(string.(...), " ")` —
  fine.

## Linked sites

- XML I/O: output format unchanged (`join` over either container is
  identical bytewise). Byte-identical roundtrip preserved.
- B1/B3 spec (`260516-coupled-basis-typeparam`): this change builds on
  the `{R}` parameter from that spec.

## Benchmarking plan

Same fixture / scripts as B1+B3 and the workspace refactor:

- `bench_b1_design_matrix.jl` on fept (energy + torque build).
- `bench_salcbasis_hotspots.jl` for `projection_matrix_coupled_basis`
  (primary alloc target).
- `profile_remaining_allocs.jl` for per-call alloc count on
  `design_matrix_energy_element` and `calc_∇ₑu!` (may already be at
  floor; recheck).

Expectation: modest. The hoisted-buffer (former #10 quick fix) already
absorbed most of the per-call alloc. Main win is in
`projection_matrix_coupled_basis` build phase and code clarity.
