# Folded tensor

The previous two pages introduced two objects: the coupled tensor
``T^{(L_f, M_f)}`` (from [Angular-momentum coupling](angular_momentum_coupling.md))
and the SALC coefficient vector ``c_\nu^{M_f}`` (from
[Symmetry adaptation](symmetry_adaptation.md)). The design-matrix
construction uses them as a pair: every per-cluster contribution to
``\Phi_\nu`` is a linear combination, indexed by ``M_f``, of products of
site harmonics. Because that combination is linear in ``M_f``, the two
tensors can be contracted along ``M_f`` once, ahead of time. The
resulting rank-``(R-1)`` tensor — the *folded tensor* — is what the
design-matrix kernel actually reads.

## Definition

For one `CoupledBasis_with_coefficient` (the unit stored inside
`salc_list[ν]`), the folded tensor is

```math
\tilde T_{m_1, \ldots, m_N}
  = \sum_{M_f} c_\nu^{M_f}\; T^{(L_f, M_f)}_{m_1, \ldots, m_N},
```

where ``N`` is the number of sites in the cluster (so
``\tilde T`` has rank ``R - 1 = N``). It carries everything about this
SALC's angular dependence per cluster, with the ``M_f`` axis already
collapsed.

## Why fold

The per-cluster, per-`CoupledBasis_with_coefficient` contribution to the
design feature
``\Phi_\nu`` is (see [Design matrix and fitting](design_matrix_and_fitting.md))

```math
\sum_{M_f} c_\nu^{M_f} \sum_{m_1, \ldots, m_N}
  T^{(L_f, M_f)}_{m_1, \ldots, m_N}
  \prod_{k=1}^{N} Z_{l_k, m_k}\bigl(\hat{\boldsymbol{e}}_k\bigr).
```

The site harmonics do not depend on ``M_f``, and the SALC coefficient
``c_\nu^{M_f}`` is fixed once symmetry adaptation has run. Linearity
lets the ``M_f`` sum be carried out first, giving ``\tilde T`` and a
shorter per-cluster expression:

```math
\sum_{m_1, \ldots, m_N}
  \tilde T_{m_1, \ldots, m_N}\,
  \prod_{k=1}^{N} Z_{l_k, m_k}\bigl(\hat{\boldsymbol{e}}_k\bigr).
```

In `design_matrix_energy_element` and `calc_∇ₑu!`, this collapses what
used to be a nested ``M_f`` / site-``m`` contraction into a single
site-``m`` contraction. For the FeGe B20 2x2x2 fixture (with
``L_f \in \{0, 1, 2\}`` mixed across the body-2 and body-3 SALCs), the
measured wall-time improvement was about 30 % on
`build_design_matrix_energy` and about 15 % on
`build_design_matrix_torque`; allocations were unchanged because no
other loop structure changes.

## Storage

The folded tensor is computed once, in the inner constructor of
`CoupledBasis_with_coefficient`, and stored alongside the originals so
both representations stay available:

```julia
struct CoupledBasis_with_coefficient{R, N}
    ls::Vector{Int}
    Lf::Int
    Lseq::Vector{Int}
    atoms::Vector{Int}
    coeff_tensor::Array{Float64, R}    # T^(Lf, Mf); kept for reference / I/O
    coefficient::Vector{Float64}       # c_ν^Mf;    kept for reference / introspection
    multiplicity::Int
    folded_tensor::Array{Float64, N}   # T̃; consumed by the hot path. N = R - 1.
end
```

The second type parameter ``N = R - 1`` is the rank of ``\tilde T``; it
appears as its own parameter because Julia does not permit arithmetic on
the rank parameter ``R`` inside a struct field declaration.

XML I/O persists ``T`` and ``c_\nu^{M_f}``, not ``\tilde T``: the inner
constructor recomputes the folded tensor from those two fields on load,
so on-disk files round-trip byte-for-byte across this change.

## Numerical note

The fold is computed in increasing ``M_f`` order — the same order the
original in-kernel loop used. The change in floating-point reduction
order is therefore confined to distributing the ``c_\nu^{M_f}`` weight
outwards across the cluster-level summation. In practice the
post-folding design matrix matches the pre-folding one to within
``\sim 10^{-14}`` relative error (a few ULPs) on the FeGe B20 2x2x2 and
FePt L1_0 2x2x2 regression fixtures.
