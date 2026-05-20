# Real spherical harmonics

The single-site building block of every SCE basis function is a spherical
harmonic of the spin direction. This page describes that block and the
module that computes it, `src/TesseralHarmonics.jl`.

## Why a spherical-harmonic basis

The energy depends on each spin only through its direction ``\hat{e}_i``, a
point on the unit sphere. Any well-behaved function on the sphere expands
in spherical harmonics, so they are the natural per-site basis: a finite
set ``\{Z_{l,m}\}`` with ``l \le l_{\max}`` spans all angular dependence up
to angular resolution ``l_{\max}``. Here ``l = 0, 1, 2, \ldots`` is the
degree (the angular momentum) and ``m`` the order, with ``-l \le m \le l``,
so the number of ``(l, m)`` pairs for a single site is ``(l_{\max}+1)^2``.

Truncating at ``l_{\max}`` is the central modeling choice for the per-site
factor. ``l = 1`` already captures all bilinear (two-spin) physics —
Heisenberg exchange, the Dzyaloshinskii–Moriya interaction, and
single-pair anisotropy — as shown on the
[Technical Notes](../technical_notes.md) page. Higher ``l_{\max}`` adds
higher-order angular terms. The per-species value is set by `body1_lmax`
in `InteractionSpec`.

## Real (tesseral) harmonics

Magesty.jl uses *real* spherical harmonics ``Z_{l,m}`` (tesseral
harmonics), not the complex ``Y_{l,m}``. The two are related by a fixed
unitary transformation (`src/SphericalHarmonicsTransforms.jl`,
`r2c_sph_harm_matrix` / `c2r_sph_harm_matrix`). Working in the real basis
keeps every coefficient tensor, design-matrix entry, and fitted ``J_\nu``
real, which is both physically transparent and faster.

The ``Z_{l,m}`` are orthonormal over the unit sphere; each carries a
normalization factor ``1/\sqrt{4\pi}``. For ``l = 1`` on a unit spin
vector ``\hat{e} = (e_x, e_y, e_z)``,

```math
Z_{1,m}(\hat{e}) = \sqrt{\tfrac{3}{4\pi}}\;\hat{e}_{\mu(m)},
\qquad
\mu(m) = \begin{cases} y & m = -1 \\ z & m = 0 \\ x & m = +1. \end{cases}
```

The per-site ``1/\sqrt{4\pi}`` is exact bookkeeping: the product of ``n``
single-site harmonics for an ``n``-site cluster carries ``(4\pi)^{-n/2}``.
The design-matrix construction cancels it with a ``(4\pi)^{n/2}`` factor.
That factor is dimensionless — the fitted ``J_\nu`` carry the input energy
unit (typically eV) either way — but it fixes the *normalization* of
``J_\nu``, so their magnitudes map cleanly onto conventional spin-model
parameters without a stray ``4\pi``. That step is described on the
[Design matrix and fitting](design_matrix_and_fitting.md) page.

## Implementation

`src/TesseralHarmonics.jl` follows the formalism of Drautz,
*Phys. Rev. B* **102**, 024104 (2020). The relevant entry points:

| Function | Role |
|---|---|
| `Zₗₘ(l, m, uvec)` | Tesseral harmonic ``Z_{l,m}``; validates inputs. |
| `Zₗₘ_unsafe(l, m, uvec)` | Same value, no validation — for hot paths. |
| `∂ᵢZlm(l, m, uvec)` | Cartesian gradient ``\partial Z_{l,m}/\partial\hat{e}`` as a `Vector{Float64}`; validates inputs. |
| `∂ᵢZlm_unsafe(l, m, uvec)` | Cartesian gradient as an `SVector{3,Float64}`, no validation — for hot paths. |

The unsafe variants are semantically identical to the checked ones but
skip bounds and argument validation; the basis-evaluation loop in
`Fitting.jl` uses them, together with pre-allocated scratch buffers, to
avoid per-call heap allocation. The gradient ``\partial Z_{l,m}/\partial\hat{e}`` is what
makes torque prediction possible — see
[Design matrix and fitting](design_matrix_and_fitting.md).

!!! note "Convention-linked code"
    The normalization and signs of the ``Z_{l,m}`` family are a fixed
    convention shared with the angular-momentum coupling, the
    symmetry-adapted basis, and the SpheriCart agreement test. They must
    not be changed in isolation. See the
    [technical notes](https://Tomonori-Tanaka.github.io/Magesty.jl/technical_notes/)
    for the precise definitions.

Next: [Angular-momentum coupling](angular_momentum_coupling.md).
