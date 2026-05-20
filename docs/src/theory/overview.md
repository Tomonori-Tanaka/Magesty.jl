# Overview

This section explains the theory behind Magesty.jl — the spin-cluster
expansion (SCE) — and connects each concept to the module that implements
it. It is meant to be read in order; each page builds on the previous one.

## The problem

A noncollinear spin DFT calculation returns, for one fixed set of spin
directions ``\{\hat{e}_i\}``, a total energy and the torque on every
magnetic atom. We want a closed-form *effective spin model*: an energy
function

```math
E\bigl(\{\hat{e}_i\}\bigr)
```

that reproduces the DFT energies and torques over the whole space of
noncollinear configurations, not just the few that were computed. Here
``\hat{e}_i`` is the unit vector along the classical spin on atom ``i``
(the magnetic-moment magnitude is treated separately).

The SCE solves this by writing ``E`` as a linear model in a fixed,
physically motivated basis, then fitting the linear coefficients to a set
of DFT reference configurations.

## The spin-cluster expansion

The SCE expands the energy as a sum of contributions from *clusters* of
atoms — single sites, pairs, triplets, and so on:

```math
E\bigl(\{\hat{e}_i\}\bigr)
  = J_0 + \sum_{\nu} J_\nu \, \Phi_\nu\!\bigl(\{\hat{e}_i\}\bigr).
```

``J_0`` is a constant reference energy, and the sum index ``\nu`` labels
the terms of the expansion. Each ``\Phi_\nu`` is a known *basis function*
of the spin directions, and each ``J_\nu`` is a scalar coefficient to be
fitted. This is the angular-momentum analogue of the
cluster expansion used for substitutional alloys: instead of discrete
occupation variables, the degrees of freedom are continuous spin
directions on the unit sphere, so the per-site factors of ``\Phi_\nu`` are
spherical harmonics of ``\hat{e}_i`` rather than occupation numbers.

Concretely, each ``\Phi_\nu`` is built from per-site real spherical
harmonics ``Z_{l,m}``. For a cluster of ``n`` atoms the elementary building
block is a product of one harmonic per atom,

```math
(4\pi)^{n/2}\;\prod_{a=1}^{n} Z_{l_a, m_a}\bigl(\hat{e}_{i_a}\bigr),
```

where ``i_a`` is the ``a``-th atom of the cluster, ``(l_a, m_a)`` its
angular indices, and ``(4\pi)^{n/2}`` a fixed normalization factor (see
[Real spherical harmonics](spherical_harmonics.md)). A basis function
``\Phi_\nu`` is a specific linear combination of such products: the
per-site harmonics are *coupled* into a rotationally well-behaved object
and the result is *symmetrized* over the crystal, so that one coefficient
``J_\nu`` multiplies one independent interaction.

Building a usable model means answering three questions, each handled by a
later page:

1. **Which basis functions?** Per-site real spherical harmonics, multiplied
   together across the atoms of a cluster and coupled into rotationally
   well-behaved combinations.
2. **How many independent coefficients?** Crystal symmetry forces many
   ``J_\nu`` to be equal or zero; symmetry adaptation keeps only the
   independent ones.
3. **How are the coefficients found?** The model is linear in ``J_\nu``, so
   they follow from a (regularized) least-squares fit.

## Pipeline and modules

The implementation mirrors this structure. The construction stages, in
dependency order, are:

| Stage | Module | Theory page |
|---|---|---|
| Crystal structure and supercell | `src/Structures.jl` | — |
| Space-group symmetry | `src/Symmetries.jl` | [Symmetry adaptation](symmetry_adaptation.md) |
| Cluster enumeration and cutoffs | `src/Clusters.jl` | [Symmetry adaptation](symmetry_adaptation.md) |
| Single-site basis (``Z_{l,m}``) | `src/TesseralHarmonics.jl` | [Real spherical harmonics](spherical_harmonics.md) |
| Angular-momentum coupling | `src/AngularMomentumCoupling.jl`, `src/CoupledBases.jl` | [Angular-momentum coupling](angular_momentum_coupling.md) |
| Symmetry-adapted basis (SALC) | `src/SALCBases.jl` | [Symmetry adaptation](symmetry_adaptation.md) |
| Design matrix and regression | `src/Fitting.jl` | [Design matrix and fitting](design_matrix_and_fitting.md) |

The user-facing types follow the same path: `SCEBasis` carries the
symmetry-adapted basis, `SCEDataset` pairs it with reference
configurations and the design matrices, `SCEFit` holds the fitted
coefficients, and `SCEModel` is the lightweight predictor. See
[`SPEC.md`](https://github.com/Tomonori-Tanaka/Magesty.jl/blob/main/SPEC.md)
for the full type layout and public API.

Once the coefficients ``J_\nu`` are fitted, the
[Technical Notes](../technical_notes.md) page shows how the lowest-order
ones map onto conventional spin-model parameters (Heisenberg ``J``,
Dzyaloshinskii–Moriya vector ``\vec{D}``, anisotropic exchange
``\boldsymbol{\Gamma}``).

## References

1. R. Drautz and M. Fähnle, "Spin-cluster expansion: Parametrization of the
   general adiabatic magnetic energy surface with ab initio accuracy",
   *Phys. Rev. B* **69**, 104404 (2004).
   DOI: [10.1103/PhysRevB.69.104404](https://doi.org/10.1103/PhysRevB.69.104404)
2. T. Tanaka and Y. Gohda, "General spin models from noncollinear spin
   density functional theory and spin-cluster expansion",
   [arXiv:2512.04458](https://arxiv.org/abs/2512.04458)
