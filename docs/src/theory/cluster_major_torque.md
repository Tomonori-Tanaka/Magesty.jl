# Cluster-major torque kernel

The torque design matrix ``X_T`` differs from the energy matrix ``X_E``
in one structural way: every row of ``X_T`` is keyed by an atom, not just
a configuration. For each spin configuration ``c`` and atom ``a``, the
three rows of ``X_T`` are the per-SALC torque components

```math
(\boldsymbol{\tau}_a)_\nu = \hat{\boldsymbol{e}}_a
  \times \nabla_{\hat{\boldsymbol{e}}_a}\!\Phi_\nu.
```

So the kernel that fills ``X_T`` has to produce ``N_\text{atoms}`` angular
gradients per ``(c, \nu)`` cell rather than the single scalar that ``X_E``
needs. How the loop is *nested* turns out to set the cost.

## The per-atom layout, and what it wastes

The natural-looking nesting is per-atom on the outside:

```julia
for sc_idx in spinconfigs
    for iatom in 1:num_atoms
        for ν in salcs
            for cbc in salc_list[ν]
                for cluster in cbc.clusters
                    if iatom ∉ cluster
                        continue
                    end
                    # contract folded_tensor with ∂Z at the site of iatom
                end
            end
        end
    end
end
```

Each cluster contains only ``N`` atoms (``N \le 3`` for the body-3 basis
this codebase ships), so for every ``iatom`` we visit, ``N_\text{atoms} -
N`` of the cluster scans hit the `continue` line and do no useful work.
For the FeGe B20 2x2x2 fixture (``N_\text{atoms} = 64``, ``N = 3`` for the
dominant 900 body-3 SALCs), about 95 % of the per-`(iatom, cluster)`
visits are membership checks that find nothing. Asymptotically the
wasted factor is ``N_\text{atoms} / N``.

Inside the kernel itself the wasted work is even more concrete: a
fresh per-call gradient buffer, a per-cluster `searched`/membership
scan, an `idx_buf` repopulation for *one* differentiated site, and a
final `cross(spin[iatom], grad) * scaling` per `(iatom, ν, cbc)`
tuple.

## Cluster-major reverse-mode

The fix is to make `iatom` an *output* of the loop nest, not an input.
A single sweep over the cluster orbit visits each cluster image once
and writes its gradient contribution to every site it contains:

```julia
for sc_idx in spinconfigs
    grad_buf .= 0          # shape (3, num_atoms, num_salcs)
    for ν in salcs
        for cbc in salc_list[ν]
            for cluster in cbc.clusters
                for m₁,…,m_N in multi-index of folded_tensor
                    for j in 1:N
                        a_j = cluster[j]
                        # leave-one-out: Π_{k ≠ j} Z[l_k, m_k][a_k]
                        coeff = folded_tensor[m₁,…,m_N] * multiplicity * Π_{k ≠ j} Z[…]
                        grad_buf[:, a_j, ν] .+= coeff .* ∂Z[l_j, m_j][a_j]
                    end
                end
            end
        end
    end
    # Single reduction: cross with each atom's spin direction.
    for ν in salcs
        for iatom in 1:num_atoms
            row = cross(spin[iatom], grad_buf[:, iatom, ν]) * scaling[ν]
            # write three rows of X_T
        end
    end
end
```

This is the standard reverse-mode-AD pattern: walk the forward graph
once, push contributions back to every input variable that touched it.
Here the "forward graph" of one cluster's contribution to ``\Phi_\nu``
is the product ``\prod_k Z_{l_k, m_k}(\hat{\boldsymbol{e}}_{a_k})``; the
gradient with respect to ``\hat{\boldsymbol{e}}_{a_j}`` is that product
with the ``j``-th factor replaced by ``\partial_i Z_{l_j, m_j}``. The
inner loop pays for that gradient at every site simultaneously.

Two invariants of [`cbc.clusters`](orbit_clusters.md) make this work:

- **Atoms in a cluster are pairwise distinct.** Otherwise two sites
  ``j_1, j_2`` of the same cluster would write to
  `grad_buf[:, a_{j_1}, ν]` *and* `grad_buf[:, a_{j_2}, ν]` with
  ``a_{j_1} = a_{j_2}``, and we would double-count the cluster's
  contribution. Pure translations are bijections of the supercell, so
  the distinctness of `cbc.atoms` (enforced upstream by
  `Clusters.allunique_atoms`) propagates to every translation image.
  A `@boundscheck @assert allunique(cluster_atoms)` guards the
  property at debug build time.
- **Translation symmetry only.** The hot-path symmetry group is the
  translation subgroup, so spin directions are not rotated between
  cluster images and the chain rule above stays linear in
  ``\partial_i Z_{l_j, m_j}(\hat{\boldsymbol{e}}_{a_j})``. Admitting
  non-translation operations in the hot-path symmetry set would
  invalidate this kernel.

## The accumulation buffer

The accumulation buffer `grad_buf` carries the *un-crossed* angular
gradient — the cross product with each atom's spin direction is taken
once at the end of every spinconfig, in the reduction phase. Pulling
the cross out is what makes the cluster-major layout efficient: a
cluster contributes to the gradient at each of its sites, so any cross
product would have to be deferred to "after the cluster sweep" anyway.

The shape is `(3, num_atoms, num_salcs)`, with the `xyz` axis innermost
so that reading the 3-vector at a single `(iatom, salc_idx)` cell costs
one cache line. For the FeGe 2x2x2 light fixture this is
`3 × 64 × 1232 × 8 ≈ 1.9 MB` per thread, comfortably inside L2.

## Threading

Threading is over spin configurations
(`@threads for sc_idx = 1:num_spinconfigs`), which matches the
[Spherical-harmonic cache](sh_cache.md) layout for the torque path:
each thread owns one `SHCache` and one `grad_buf` for its spinconfig,
and the design-matrix writes go to disjoint row blocks
(`row_offset = 3·num_atoms·(sc_idx - 1)`), so there is no shared
state to synchronize.

## Inference path

`_predict_torque` uses the same algorithm via a scalar-coefficient
variant `_accumulate_grad_torque_scaled!` that folds
``\text{scaling}(\nu) \cdot J_\nu`` into the per-cbc factor at
accumulation time. The buffer for inference is `(3, num_atoms)`
instead of `(3, num_atoms, num_salcs)`: the SALC dimension is summed
away as soon as it is encountered, since prediction does not need a
per-column Jacobian.

The kernel that fills ``X_T`` is `_accumulate_grad_torque_cluster!`;
its scalar-coefficient sibling `_accumulate_grad_torque_scaled!`
serves the inference path. Both are internal — nothing outside the
design-matrix and inference paths exposes a per-atom gradient API.

## Numerical note

The cluster-major restructure preserves the per-`(salc, iatom)`
accumulation order: an atom still receives contributions from
`(cbc, cluster_containing_atom, multi_idx)` tuples in the same outer
sequence the per-atom kernel walked, because the inner `for j` loop
writes only to the unique `a_j = cluster_atoms[j]`. The only
floating-point change is that `cbc.multiplicity` is now multiplied
into the per-multi-idx factor instead of applied once at the end of
each cluster — a redistribution of one rounding step. On both FeGe
and FePt regression fixtures the design-matrix equivalence test
therefore passes at ``\text{rtol} = 10^{-14}, \text{atol} = 10^{-15}``,
i.e. bit-identical to within a few ULPs.
