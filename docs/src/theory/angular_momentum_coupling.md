# Angular-momentum coupling

A single-site harmonic describes one spin. To describe a *cluster* — a
pair, triplet, or higher group of atoms — the per-site harmonics must be
combined. This page explains how, and which modules do it:
`src/AngularMomentumCoupling.jl` and `src/CoupledBases.jl`.

## Why couple

The naive cluster basis function is a plain product
``Z_{l_1,m_1}(\hat{\boldsymbol{e}}_1)\,Z_{l_2,m_2}(\hat{\boldsymbol{e}}_2)\cdots``. That product
works, but it does not transform simply under a global rotation of all
spins: the ``m`` indices mix among themselves in a way that obscures the
physics and makes symmetry adaptation awkward.

Coupling fixes this. By contracting the product with Clebsch–Gordan
coefficients, the site harmonics ``l_1, l_2, \ldots, l_n`` are combined
into an object of definite *total* angular momentum ``L_f``. Such a
coupled object transforms under a global spin rotation exactly like a
single harmonic of degree ``L_f`` — it is rotationally covariant. This is
what lets the later symmetry projection act cleanly, and it is what gives
the lowest coupled channels their direct physical reading: for a pair with
``l_1 = l_2 = 1``, the three channels ``L_f = 0, 1, 2`` are precisely
Heisenberg exchange, the Dzyaloshinskii–Moriya interaction, and
anisotropic symmetric exchange (see the
[Technical Notes](../technical_notes.md)).

This same covariance also explains the `isotropy` option. A coupled object
of total angular momentum ``L_f`` transforms like a degree-``L_f``
harmonic, so under a global rotation of all spins only the ``L_f = 0``
channel — the rotational scalar — stays invariant; every ``L_f \ge 1``
channel rotates into its other ``M_f`` components. Demanding that the
energy itself be invariant under a global spin rotation (an isotropic spin
model, with no preferred direction in spin space) therefore keeps only the
``L_f = 0`` coupled bases and discards all ``L_f \ge 1`` ones. Setting
`isotropy = true` imposes exactly this restriction, leaving the purely
isotropic, Heisenberg-like part of the model.

## Coupling scheme

Coupling ``n`` momenta is done one pair at a time, in a fixed *left-coupling*
order:

```math
l_1 \otimes l_2 \to L_{12},\qquad
L_{12} \otimes l_3 \to L_{123},\qquad \ldots,\qquad
L_{1\cdots(n-1)} \otimes l_n \to L_f.
```

Each step obeys the triangle rule ``|l_a - l_b| \le L \le l_a + l_b``, so a
given set ``\{l_i\}`` admits several admissible coupling sequences. A
sequence is recorded as a pair `(Lseq, Lf)`, where `Lseq` is the list of
intermediate momenta ``[L_{12}, L_{123}, \ldots]`` and `Lf` is the final
total ``L_f``. Inside a fixed tree, different ``(L_{12}, \ldots, L_f)``
labels produce orthogonal coupled states. Other coupling trees are related
to the left-coupling tree by Racah (6j) recoupling, so fixing one tree
loses no generality.

## Implementation

`src/AngularMomentumCoupling.jl` builds the coupled bases:

| Function | Role |
|---|---|
| `enumerate_paths_left_all(ls)` | Enumerate all admissible `(Lseq, Lf)` for site momenta `ls = [l₁,…,l_n]`. |
| `coeff_tensor_complex(ls, Lseq, Lf)` | Chain Clebsch–Gordan coefficients into a coupling tensor (complex basis). |
| `complex_to_real_tensor(tensor, ls, Lf)` | Transform that tensor from the complex ``Y_{l,m}`` basis to the real ``Z_{l,m}`` basis. |
| `build_all_real_bases(...)` | Build every real coupled basis up to the requested coupling order. |

The result of coupling is a real coefficient tensor ``T^{(L_f,M_f)}``
indexed by the per-site ``m`` values and the order ``M_f`` of the final
multiplet (``-L_f \le M_f \le L_f``). A cluster
basis function is the full contraction of that tensor with the site
harmonics:

```math
T^{(L_f, M_f)}_{m_1 m_2 \cdots m_n}\;
  Z_{l_1, m_1}(\hat{\boldsymbol{e}}_1)\, Z_{l_2, m_2}(\hat{\boldsymbol{e}}_2)\cdots
  Z_{l_n, m_n}(\hat{\boldsymbol{e}}_n),
```

summed over all site ``m`` indices.

`src/CoupledBases.jl` stores this as the type `CoupledBasis{R}` (rank
``R = n + 1``: one tensor mode per site plus one for ``M_f``):

```julia
struct CoupledBasis{R}
    ls::Vector{Int}                # per-site angular momenta
    Lf::Int                        # final total L
    Lseq::Vector{Int}              # intermediate L sequence
    atoms::Vector{Int}             # atom indices of the cluster
    coeff_tensor::Array{Float64,R} # T^{(Lf,Mf)}, real (tesseral) basis
end
```

A `CoupledBasis` is still tied to specific atom indices but is not yet
symmetry-adapted. The next stage — projecting onto crystal symmetry — is
covered in [Symmetry adaptation](symmetry_adaptation.md).
