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

### Representable interaction range

Cluster distances are measured under the **minimum-image convention**. From
the 3×3×3 block of periodic images of the supercell (`calc_x_images` in
`src/Structures.jl` — the home cell plus its 26 neighbors), `set_mindist_pairs`
in `src/Clusters.jl` keeps, for each atom pair, the shortest image together
with every image tied with it to within tolerance. A pair displacement is
therefore represented uniquely only when it is the strictly shortest vector
among its periodic images, i.e. when it lies **inside the Wigner–Seitz cell of
the supercell lattice**. Beyond that region the shortest image folds onto a
closer lattice vector, so longer shells collapse onto shorter ones and the
model effectively treats interactions past the representable range as zero.

This region is anisotropic. For a cubic supercell of edge ``L`` (Wigner–Seitz
cell ``[-L/2,\,L/2]^3``) it reaches ``L/2`` along ``\langle 100\rangle`` but
``\sqrt{3}\,L/2`` along ``\langle 111\rangle``; in a bcc 2×2×2 cell
(``L = 2a``) the origin–body-corner pair at distance
``\sqrt{3}\,a = \sqrt{3}\,L/2`` sits exactly on the Wigner–Seitz corner
(eightfold degenerate) and is still kept. A convenient direction-independent
guarantee is **range ``< L/2``** (the inscribed-sphere radius): any
interaction shorter than this is always representable, which is the safe rule
when choosing a supercell. Pairs sitting on the Wigner–Seitz boundary are the
equal-distance degenerate case recorded by `multiplicity` (below).

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
struct CoupledBasis_with_coefficient{R, N}
    ls::Vector{Int}
    Lf::Int
    Lseq::Vector{Int}
    atoms::Vector{Int}
    coeff_tensor::Array{Float64, R}
    coefficient::Vector{Float64}      # SALC coefficient, length 2·Lf + 1
    multiplicity::Int                 # translational-equivalence count
    folded_tensor::Array{Float64, N}  # see Folded tensor; N = R - 1
    clusters::Vector{Vector{Int}}     # pre-enumerated symmetry orbit
end
```

The `multiplicity` records how many times cluster enumeration produced this
exact `CoupledBasis`. In a periodic supercell, translationally equivalent
placements of a cluster — most visibly when its atoms lie on a cell
boundary — generate the identical coupled basis more than once; the basis
is kept once, with the repeat count stored as `multiplicity`. The
design-matrix construction weights each coupled basis by this count, so
every translationally equivalent copy contributes to ``\Phi_\nu``.

### Cell-boundary pairs: odd ``L_f`` SALCs vanish

A symmetry cancellation is specific to those cell-boundary pairs. When a
pair's two atoms are translationally equivalent (they fold to the same
primitive atom) and the pair lies on a Wigner–Seitz *face*, both orderings
``(a, b)`` and ``(b, a)`` appear in the orbit (``multiplicity > 1``): a pure
supercell translation maps the ordered pair onto its own reverse. Under
exchange of the two sites the coupled pair tensor of folded angular momentum
``L_f`` picks up a factor ``(-1)^{L_f}`` (for equal site momenta; in general
``(-1)^{l_1 + l_2 - L_f}``, see
[Angular-momentum coupling](angular_momentum_coupling.md)). The projection
symmetrizes over the group, which now contains this site-swapping translation,
so for **odd ``L_f`` the tensor equals its own negative and is projected to
zero**; even ``L_f`` survives.

For an ``l_1 = l_2 = 1`` pair this removes the ``L_f = 1`` channel — the
antisymmetric, Dzyaloshinskii–Moriya part — while ``L_f = 0`` (Heisenberg) and
``L_f = 2`` (symmetric anisotropy) remain. It is the same mechanism that
forbids a Dzyaloshinskii–Moriya vector on an inversion-symmetric bond.
`SALCBases.filter_basisdict` drops these bases up front: a cluster is removed
when all its atoms fold to one primitive atom, ``L_f`` is odd, and the
multiplicity exceeds one. The FeGe 2×2×2 fixture exercises this.

The trailing `folded_tensor` is the contraction
``\tilde T = \sum_{M_f} c_\nu^{M_f}\, T^{(L_f, M_f)}``, precomputed in the
inner constructor and read by the design-matrix kernel in place of the
``(coeff\_tensor, coefficient)`` pair. The standalone
[Folded tensor](folded_tensor.md) page explains the linearity argument and
the kernel-side consequences; XML I/O persists only ``T`` and
``c_\nu^{M_f}`` so files round-trip unchanged.

`clusters` is the full symmetry-orbit list — the distinct atom tuples
produced by applying every pure translation to the seed `atoms` — built
once via `enumerate_orbit_clusters` during basis construction (and
recomputed lazily on XML load). The design-matrix kernel iterates this
list directly instead of redoing the translation walk and the
sorted-tuple dedup for each ``\Phi_\nu`` evaluation. The
[Orbit clusters](orbit_clusters.md) page covers the enumeration,
storage layout, and the cluster-distinctness invariant downstream
kernels rely on.

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
