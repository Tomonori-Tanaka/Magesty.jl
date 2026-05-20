# Symmetry adaptation (SALC)

Crystal symmetry constrains the spin model: symmetry-equivalent clusters
must carry equal coefficients, and some basis functions are forbidden
outright. Imposing those constraints up front — rather than hoping the fit
discovers them — is the job of *symmetry-adapted linear combinations*
(SALCs). This page covers `src/Clusters.jl`, `src/Symmetries.jl`, and
`src/SALCBases.jl`.

## Clusters and orbits

The expansion runs over clusters: pairs, triplets, and higher groups of
atoms within distance cutoffs. `src/Clusters.jl` enumerates them, using
per-body-order, per-species-pair cutoff radii (`bodyn_cutoff` in
`InteractionSpec`); the body order `nbody` sets the largest cluster size.

The space group, obtained from Spglib through `src/Symmetries.jl`, sorts
these clusters into *orbits*: sets of clusters mapped onto one another by
symmetry operations. Every cluster in an orbit must contribute with the
same coefficient, because the crystal cannot tell its members apart. So
the free parameters of the model are indexed by orbits, not by individual
clusters.

## From coupled bases to SALCs

A coupled basis (see [Angular-momentum coupling](angular_momentum_coupling.md))
attached to one cluster is not yet symmetry-adapted: under a symmetry
operation it generally turns into a combination of coupled bases on other
clusters of the same orbit. Symmetry adaptation projects the coupled bases
onto the irreducible representations of the crystal symmetry group. The
projection produces linear combinations — SALCs — that are invariant under
the group, so each surviving SALC corresponds to exactly one independent
coefficient ``J_\nu``.

Concretely, each SALC is a set of coupled bases sharing an orbit, each
weighted by a SALC coefficient vector ``c_\nu^{M_f}`` (one entry per
``M_f``). Bases that the projection sends to zero are dropped: those
angular components are symmetry-forbidden. The net effect is a basis whose
size matches the true number of independent interactions, which both
shrinks the fit and guarantees the model respects the crystal symmetry
exactly.

The `isotropy` option (`SymmetryOptions.isotropy`) is a further
restriction: when set, only ``L_f = 0`` terms are kept, i.e. the purely
isotropic (Heisenberg-like) part of the model.

## Implementation

`src/SALCBases.jl` builds and stores the symmetry-adapted basis:

```julia
struct SALCBasis
    coupled_basislist::SortedCounter{CoupledBasis}
    salc_list::Vector{Vector{CoupledBasis_with_coefficient}}
    angular_momentum_couplings::Vector{AngularMomentumCouplingResult}
end
```

The central field is `salc_list`. Its outer index runs over *key groups*
(one per independent SALC / coefficient ``J_\nu``); each inner vector holds
the `CoupledBasis_with_coefficient` objects of that group — a `CoupledBasis`
extended with the SALC coefficient vector and a `multiplicity`:

```julia
struct CoupledBasis_with_coefficient{R}
    ls::Vector{Int}
    Lf::Int
    Lseq::Vector{Int}
    atoms::Vector{Int}
    coeff_tensor::Array{Float64,R}
    coefficient::Vector{Float64}   # SALC coefficient, length 2·Lf + 1
    multiplicity::Int              # translational-equivalence count
end
```

The `multiplicity` records how many times cluster enumeration produced this
exact `CoupledBasis`. In a periodic supercell, translationally equivalent
placements of a cluster — most visibly when its atoms lie on a cell
boundary — generate the identical coupled basis more than once; the basis
is kept once, with the repeat count stored as `multiplicity`. The
design-matrix construction weights each coupled basis by this count, so
every translationally equivalent copy contributes to ``\Phi_\nu``.

Constructing the SALC basis is the computationally expensive stage of the
pipeline. The result is wrapped in `SCEBasis`, which can be persisted with
`Magesty.save` and reloaded with `Magesty.load` to avoid recomputation.

!!! warning "The key-group order defines `Jφ`"
    The outer order of `salc_list` is the order of the fitted coefficient
    vector `jphi`: column ``\nu`` of the design matrix, the ``\nu``-th
    entry of `jphi`, and the ``\nu``-th key group all refer to the same
    interaction. The XML I/O serializes coefficients in this same order.
    Reordering SALCs therefore changes the physical meaning of every
    ``J_\nu``; the design-matrix construction, the fitted vector, and the
    on-disk format must stay synchronized.

With the symmetry-adapted basis in hand, the remaining step is to turn it
into a numerical design matrix and fit the coefficients:
[Design matrix and fitting](design_matrix_and_fitting.md).
