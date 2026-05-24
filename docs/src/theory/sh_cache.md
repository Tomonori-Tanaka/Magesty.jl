# Spherical-harmonic cache

The two preceding hot-path pages — [Folded tensor](folded_tensor.md) and
[Orbit clusters](orbit_clusters.md) — moved work out of the
design-matrix kernel that depended on neither the spin configuration nor
the SALC channel. The remaining per-element work that *does* depend on
the spin configuration is the evaluation of the tesseral spherical
harmonics ``Z_{l,m}(\hat{\boldsymbol{e}}_{a})`` and, for the torque
kernel, their direction gradients
``\partial_i Z_{l,m}(\hat{\boldsymbol{e}}_{a})``. Even though these
*do* vary with the spin configuration, they depend only on
``(l, m, \hat{\boldsymbol{e}}_{a})`` — not on which SALC, which coupled
basis, or which cluster image is being processed. Inside the original
kernel they were recomputed once per
``(\text{spinconfig},\, \text{SALC},\, \text{coupled basis},\,
\text{cluster image},\, \text{site})`` tuple. Hoisting that work to
"once per ``(\text{spinconfig},\, l, m, \text{atom})``" is the job of
the *spherical-harmonic cache* (`SHCache`).

## What is stored

For one spin configuration, the cache stores the harmonic values for
every atom up to a basis-determined cutoff ``l_\text{max}``:

```julia
struct SHCache
    Z::Matrix{Float64}                # (l_max+1)² × num_atoms
    ∂Z::Matrix{SVector{3, Float64}}   # same shape, or 0×0 (energy path)
    l_max::Int
end
```

Both arrays are addressed by ``(\text{lm\_idx}, \text{atom})`` with the
standard flat tesseral ordering

```math
\text{lm\_idx}(l, m) = l^2 + l + m + 1,
\qquad 0 \le l \le l_\text{max},
\qquad -l \le m \le l.
```

The range ``[l^2 + 1,\, (l + 1)^2]`` therefore covers the ``2l + 1``
m-values for a fixed ``l``. Column-major storage (atom as the trailing
index) matches the access pattern — the inner contraction loops
``m`` while the atom is fixed by the cluster image.

`l_max` is the largest ``l`` value appearing in any `cbc.ls` across the
SALC list, computed once at the start of design-matrix construction by
`Magesty.Fitting._compute_l_max`. Entries with ``l`` smaller than every
``l`` in `cbc.ls` are still computed; the wasted work is bounded by
``(l_\text{max}+1)^2`` evaluations per atom and is dwarfed by the
savings on the hot path.

## Why precompute

The per-cluster contribution to ``\Phi_\nu`` (see
[Design matrix and fitting](design_matrix_and_fitting.md)) is

```math
\tilde\Phi_\nu^{(\text{cluster})}
  = \sum_{m_1, \ldots, m_N}
    \tilde T_{m_1, \ldots, m_N}
    \prod_{k=1}^{N}
      Z_{l_k, m_k}\!\bigl(\hat{\boldsymbol{e}}_{a_k}\bigr),
```

with ``a_k`` the ``k``-th atom of this cluster image. The values
``Z_{l_k, m_k}(\hat{\boldsymbol{e}}_{a_k})`` are functions of the
*spin configuration* and the *atom*; they are independent of which SALC
or which coupled basis is being summed. Yet the same
``Z_{l_k, m_k}(\hat{\boldsymbol{e}}_{a_k})`` was being recomputed for
every cluster image, every coupled basis, and every SALC that touched
atom ``a_k`` with that ``(l_k, m_k)``. The cache replaces this with a
single sweep at the start of each spin configuration:

```julia
sh_cache = build_sh_cache_energy(spinconfig.spin_directions, l_max)
# ... or build_sh_cache_torque for the gradient path.

# Inside the kernel, in place of
#   Zₗₘ_unsafe(l, m, dir, legendre_buf)
# we read
#   sh_cache.Z[l^2 + l + m + 1, atom]
```

The energy and torque paths use distinct builders because the gradient
kernel additionally needs ``\partial_i Z_{l, m}``:
`build_sh_cache_energy` populates only ``Z``, while
`build_sh_cache_torque` populates both ``Z`` and ``\partial Z``.

For the FeGe B20 2x2x2 light fixture (64 atoms, 100 spin configurations,
1232 SALC groups, body distribution
``[1 \!\to\! 2,\, 2 \!\to\! 330,\, 3 \!\to\! 900]``) the cache cut about
85 % off `build_design_matrix_energy` and about 60 % off
`build_design_matrix_torque`, on top of the
[Folded tensor](folded_tensor.md) and [Orbit clusters](orbit_clusters.md)
gains.

## Threading and lifetime

Energy and torque thread differently, and the cache is built to match.

For the energy design matrix, `build_design_matrix_energy` threads over
SALC columns (``@threads\ \text{for}\ i = 1{:}\text{num\_salcs}``), so a
given spin configuration is needed by every thread. The caches for all
spin configurations are therefore built up front in a separate
``@threads\ \text{for}\ j = 1{:}\text{num\_spinconfigs}`` loop and held
in a `Vector{SHCache}` shared read-only across the SALC threads:

```julia
sh_caches = Vector{SHCache}(undef, num_spinconfigs)
@threads for j = 1:num_spinconfigs
    sh_caches[j] = build_sh_cache_energy(spinconfig_list[j].spin_directions, l_max)
end
# Subsequent @threads over SALCs reads sh_caches[j] in its inner j-loop.
```

For the torque design matrix, `build_design_matrix_torque` threads over
spin configurations (``@threads\ \text{for}\ \text{sc\_idx} =
1{:}\text{num\_spinconfigs}``), so each thread builds and owns the
``Z + \partial Z`` cache for its own configuration only:

```julia
@threads for sc_idx = 1:num_spinconfigs
    sh_cache = build_sh_cache_torque(spinconfig.spin_directions, l_max)
    # ... use across all (iatom, salc_idx, cbc) inside this thread.
end
```

The cache lives for the duration of one spin configuration's work in
both layouts. It is never serialized: it is purely runtime memoization,
fully reconstructible from `spinconfig.spin_directions` and `l_max`.

## Hot-path simplification

With the cache in place, `design_matrix_energy_element` and `calc_∇ₑu!`
drop their `spin_directions`, `symmetry`, and per-thread workspace
arguments altogether: the only spin-configuration-dependent input they
need is the cache, and the only symmetry-dependent input
(`cbc.clusters`) lives on the coupled basis. The signatures collapse to

```julia
design_matrix_energy_element(cbc, sh_cache)::Float64
calc_∇ₑu!(result, cbc, atom, sh_cache)
```

and the per-call scratch buffers that used to live on `EnergyWorkspace`
/ `GradWorkspace` (`sh_values` / `sh_offsets` / `legendre_buf` /
`atom_grad_values`) are no longer needed and have been removed
together with the workspace structs themselves.

## Numerical note

The cache is *exact value memoization*: every entry holds the same
`Float64` that the previous per-call `Zₗₘ_unsafe` /
`∂ᵢZlm_unsafe` invocation would have returned, and no reduction order
changes. The design-matrix equivalence test consequently passes at
``\text{rtol} = 10^{-14},\, \text{atol} = 10^{-15}`` on both FeGe and
FePt regression fixtures — bit-identical to within a few ULPs.
