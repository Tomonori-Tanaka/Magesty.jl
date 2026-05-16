# Design: Optimize hot-path scratch workspace

## Workspace structs

Two small structs in `src/Optimize.jl`. Mutable so `empty!` and resize
mutations are cheap.

```julia
mutable struct EnergyWorkspace
    searched_pairs::Set{UInt}
    sh_values::Vector{Vector{Float64}}

    EnergyWorkspace() = new(Set{UInt}(), Vector{Vector{Float64}}())
end

mutable struct GradWorkspace
    searched_pairs::Set{UInt}
    sh_values::Vector{Vector{Float64}}
    atom_grad_values::Vector{SVector{3, Float64}}

    GradWorkspace() = new(
        Set{UInt}(),
        Vector{Vector{Float64}}(),
        Vector{SVector{3, Float64}}(),
    )
end
```

Note: `∂ᵢZlm_unsafe` returns `SVector{3, Float64}`. Storing it in
`Vector{Vector{Float64}}` (as the pre-refactor code did) forced a
`convert` to `Vector{Float64}` on every assignment — that's where many
of the grad-side allocations originated. Switching to
`Vector{SVector{3, Float64}}` matches the producer's return type and
eliminates the conversion.

These are **not exported**. They are an internal optimization
contract between `build_design_matrix_*` and the per-element kernels.

## Resize helpers

```julia
function _ensure_sh_buffer!(buf::Vector{Vector{Float64}}, ls::Vector{Int})
    N = length(ls)
    # Grow outer if needed
    while length(buf) < N
        push!(buf, Vector{Float64}())
    end
    # Ensure each inner has capacity 2lᵢ+1
    @inbounds for i in 1:N
        needed = 2 * ls[i] + 1
        if length(buf[i]) != needed
            resize!(buf[i], needed)
        end
    end
    return buf
end
```

The outer never shrinks within a workspace's lifetime — once a thread
sees max-N cluster it keeps that capacity. Inner vectors resize as
needed; `resize!` on `Vector{Float64}` is amortized O(1) and
reuses the existing buffer when growing within capacity.

For `atom_grad_values` we need a separate helper sized by
`2 * l_atom + 1` (where l_atom is the differentiated site's l) —
since this is a single scalar size, it's a `resize!` of the outer
vector and each inner.

## Function signatures

```julia
function design_matrix_energy_element(
    cbc::Basis.CoupledBasis_with_coefficient{R},
    spin_directions::AbstractMatrix{<:Real},
    symmetry::Symmetry,
    ws::EnergyWorkspace,
)::Float64 where {R}
    empty!(ws.searched_pairs)
    _ensure_sh_buffer!(ws.sh_values, cbc.ls)
    # ... rest of body uses ws.searched_pairs and ws.sh_values
end

function calc_∇ₑu!(
    result::MVector{3, Float64},
    cbc::Basis.CoupledBasis_with_coefficient{R},
    atom::Integer,
    spin_directions::AbstractMatrix{<:Real},
    symmetry::Symmetry,
    ws::GradWorkspace,
) where {R}
    empty!(ws.searched_pairs)
    _ensure_sh_buffer!(ws.sh_values, cbc.ls)
    # ... atom_grad_values resized inside the loop after l_atom is known
end
```

## Caller changes

`build_design_matrix_energy`:

```julia
@threads for i = 1:num_salcs
    key_group = salc_list[i]
    ws = EnergyWorkspace()  # one per thread iteration; one alloc per thread/salc
    # ...
    for cbc in key_group
        group_value += design_matrix_energy_element(cbc, sd, sym, ws)
    end
end
```

Strictly speaking we want one workspace per *thread*, not per
salc-iteration. But each `@threads` iteration is one chunk
assigned to one thread, and Julia's default scheduler hands
contiguous chunks. Creating one workspace per `i` is slightly
wasteful (one extra allocation per salc) compared to per-thread,
but keeps the code simple. If profiling shows this matters, switch
to `Vector{EnergyWorkspace}(undef, nthreads())` allocated outside
the loop.

Initial implementation: per-`i` workspace. Switch to per-thread if
benchmarks justify it.

`build_design_matrix_torque`:

```julia
@threads for sc_idx in 1:num_spinconfigs
    ws = GradWorkspace()
    # ...
    for cbc in key_group
        calc_∇ₑu!(grad_u_buf, cbc, iatom, sd, sym, ws)
    end
end
```

`_predict_energy` is called once per inference, not threaded:

```julia
ws = EnergyWorkspace()
for i = 1:num_salcs
    # ...
end
```

## Non-bang wrapper

`calc_∇ₑu(cbc, atom, spin_directions, symmetry)` keeps its current
signature; internally allocates a `GradWorkspace`. Not hot path.

## Linked sites

None — workspace is an Optimize-internal contract. XML I/O, SALCBases,
type definitions all unchanged.

## Benchmarking plan

Pre/post on the same fept fixture:

- Per-call allocs for `design_matrix_energy_element` and `calc_∇ₑu!`
  (using the script in `test/develop_tmp/profile_remaining_allocs.jl`).
- Wall time for `build_design_matrix_energy` and
  `build_design_matrix_torque` (extend
  `test/develop_tmp/bench_b1_design_matrix.jl` or reuse).

Target: 50×+ alloc reduction per element call; wall-time speedup
modest (likely 1.2-2× on energy since contraction loop dominates
the remaining 2ms).
