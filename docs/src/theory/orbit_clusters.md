# Orbit clusters

[Symmetry adaptation](symmetry_adaptation.md) explained that the SCE
expansion runs over *cluster orbits* — sets of clusters that the
crystal's symmetry maps onto one another. When the design feature
``\Phi_\nu`` is evaluated, every cluster in the orbit must contribute, so
the kernel has to walk through them in some order. The earlier
implementation re-derived that walk per matrix element, by sweeping
through every pure translation of the seed cluster and de-duplicating
the resulting atom tuples on the fly. The orbit itself, however, is a
*structural* object — it depends only on the seed atoms and the
crystal's symmetry, not on the spin configuration or the Mf channel —
so the walk is moved to SALC build time and the result is cached on
each `CoupledBasis_with_coefficient`. The hot path then iterates that
cached list directly.

## What is stored

The `clusters` field of `CoupledBasis_with_coefficient` is a
`Vector{Vector{Int}}`. Each inner vector is one cluster in the orbit
— an `N`-element list of atom indices, in the same site order as
`cbc.atoms`, so its `k`-th entry sits on the `l_k`-th
orbital-angular-momentum axis of the coupled tensor:

```julia
struct CoupledBasis_with_coefficient{R, N}
    ls::Vector{Int}
    Lf::Int
    Lseq::Vector{Int}
    atoms::Vector{Int}
    coeff_tensor::Array{Float64, R}
    coefficient::Vector{Float64}
    multiplicity::Int
    folded_tensor::Array{Float64, N}      # see Folded tensor
    clusters::Vector{Vector{Int}}         # all distinct images of `atoms`
end
```

For a seed cluster `cbc.atoms` and a supercell symmetry whose pure
translations are indexed by `symnum_translation`, the orbit is

```math
\bigl\{\,(\sigma\,\text{atoms}[1],\, \ldots,\, \sigma\,\text{atoms}[N])
  \,\bigm|\,
  \sigma \in \text{translations}\,\bigr\},
```

with images that produce the same multiset of atoms collapsed into a
single representative (the first occurrence in `symnum_translation`
order).

## How it is built

`CoupledBases.enumerate_orbit_clusters(atoms, map_sym, symnum_translation)`
runs the enumeration. It iterates the translations, applies
`map_sym[atoms[k], itrans]` site by site to get each image, hashes the
sorted image (so two translations whose images are the same multiset
are deduplicated), and returns the un-sorted images in site order:

```julia
clusters = enumerate_orbit_clusters(
    cbc.atoms,
    symmetry.map_sym,
    symmetry.symnum_translation,
)
```

This is called once per seed cluster during SALC construction
(`SALCBases._compute_salc_groups`) and once per
`CoupledBasis_with_coefficient` on XML load
(`XMLIO.read_salcbasis_from_xml`). At SALC build time it is hoisted
above the SALC's Mf channel loop, since every Mf channel of the same
orbit reuses the same `clusters`.

XML I/O does **not** serialize `clusters`: the field is fully determined
by `cbc.atoms` and the (independently loaded) `Symmetry`, so on-disk
files round-trip byte-for-byte across this change.

## Why pre-enumerate

The orbit walk used to live inside the inner kernel:

```julia
# Old design_matrix_energy_element (sketch)
for itrans in symmetry.symnum_translation
    for site_idx in 1:N
        translated_atoms[site_idx] = symmetry.map_sym[cbc.atoms[site_idx], itrans]
    end
    sorted = sort(copy(translated_atoms))
    key = hash(sorted)
    key in seen && continue
    push!(seen, key)
    # ... evaluate harmonics, contract ...
end
```

That block ran for every matrix element — every spin configuration,
every SALC key group, every coupled basis inside the group. None of
the work depends on the spin configuration or on which Mf channel of
the SALC is being processed; only the seed atoms and the symmetry's
`map_sym` enter. Hoisting it to construction time means the work runs
once per orbit, full stop.

For the FeGe B20 2×2×2 light fixture (64 atoms, 100 spin
configurations, 1232 SALC groups, body distribution
``[1 \!\to\! 2,\, 2 \!\to\! 330,\, 3 \!\to\! 900]``) this collapsed the
per-element overhead far enough to recover roughly 12 % on
`build_design_matrix_energy` and 26 % on `build_design_matrix_torque`,
on top of the [Folded tensor](folded_tensor.md) speedup. Per-element
allocations dropped accordingly (energy −50 %, torque −67 % vs the
pre-restructuring baseline) because the per-element sort buffer,
membership scan, and hash-key dedup are gone.

## Invariants worth knowing

- Atom indices inside a single `clusters[k]` are *distinct*. The
  cluster enumeration upstream (`Clusters.allunique_atoms` filter) only
  emits clusters with distinct atoms, and pure translations are
  bijections of the supercell, so the property propagates.
- `clusters` carries the unique members of the orbit only; the
  separately stored `multiplicity` field counts how many distinct
  *coupled bases* — not clusters — collapse to the same
  `(ls, Lf, Lseq, atoms)` tuple during basis construction. The two
  fields encode different aspects of the symmetry data and are not
  redundant.
