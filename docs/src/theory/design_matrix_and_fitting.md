# Design matrix and fitting

The symmetry-adapted basis fixes *which* functions ``\Phi_\nu`` appear in
the expansion. This final page turns those functions into numbers — a
design matrix — and fits the coefficients ``J_\nu``. The relevant module
is `src/Fitting.jl`.

## The design feature

Recall the SCE energy of a configuration ``\{\hat{e}_i\}``:

```math
E = J_0 + \sum_{\nu} J_\nu \, \Phi_\nu\!\bigl(\{\hat{e}_i\}\bigr).
```

Each ``\Phi_\nu`` is the basis function of one SALC key group ``\nu``. For
an ``n``-site cluster it is evaluated by contracting the coupled tensor
``T^{(L_f,M_f)}`` with the site harmonics, weighting by the SALC
coefficient ``c_\nu^{M_f}``, and summing over every cluster in the
symmetry orbit. For a 2-body SALC this reads

```math
\Phi_\nu = (4\pi)\sum_{M_f} c_\nu^{M_f}
  \!\!\sum_{\text{clusters in orbit}}\;\sum_{m_1,m_2}
  T^{(L_f,M_f)}_{m_1 m_2}\;
  Z_{1,m_1}(\hat{e}_i)\,Z_{1,m_2}(\hat{e}_j).
```

The prefactor generalizes to ``(4\pi)^{n/2}`` for an ``n``-site cluster.
It cancels the ``(4\pi)^{-n/2}`` carried by the ``n`` per-site tesseral
harmonics (see [Real spherical harmonics](spherical_harmonics.md)). Being
dimensionless it does not change the unit of ``J_\nu`` — always the input
energy unit, typically eV — but it fixes their normalization, so the
low-order coefficients map onto conventional spin-model parameters through
the factors derived in the [Technical Notes](../technical_notes.md). In the
code this is the helper

```julia
@inline _cluster_scaling(n_sites::Integer)::Float64 = (4π)^(n_sites / 2)
```

## Building the design matrix

A linear model in ``J_\nu`` becomes a matrix equation once ``\Phi_\nu`` is
tabulated over the reference configurations. `build_design_matrix_energy`
produces the energy design matrix ``X_E``: one row per spin configuration,
one column per SALC key group, entry ``(c, \nu) = \Phi_\nu`` for
configuration ``c``. Each column entry is the sum of
`design_matrix_energy_element` over the coupled bases of the group, times
`_cluster_scaling`.

DFT also provides per-atom torques, and these constrain the same
coefficients. `build_design_matrix_torque` produces ``X_T`` with
``3 \times n_\text{atoms}`` rows per configuration. The torque is the
angular gradient of the energy,

```math
\boldsymbol{\tau}_i = \hat{e}_i \times \nabla_{\hat{e}_i} E,
```

which is why the spherical-harmonic gradient ``\partial Z_{l,m}/\partial\hat{e}``
(`∂ᵢZlm`) is needed alongside the values.

The constant ``J_0`` is *not* a column of ``X_E``. The energy matrix is
mean-centered instead, and ``J_0`` is recovered analytically after the
solve by `extract_j0_jphi` as ``J_0 = \operatorname{mean}(y_E - X_E\,J)``,
where ``y_E`` is the vector of observed (DFT) energies and ``J = (J_\nu)``
the fitted coefficient vector.

These matrices, together with the reference data, are held by `SCEDataset`
(fields `X_E`, `X_T`, `y_E`, `y_T`). They are stored *unweighted*, so the
energy/torque balance can be re-tuned at fit time without rebuilding them.

## Fitting the coefficients

`fit(SCEFit, dataset, estimator; torque_weight)` solves the regression.
The energy and torque residuals are combined into a single loss by a
convex weight ``w =`` `torque_weight` ``\in [0, 1]``:

```math
\text{loss}(J, J_0) = (1 - w)\,\text{MSE}_\text{energy}
                    + w\,\text{MSE}_\text{torque},
```

with each term a *per-sample* mean squared error
(``\text{MSE}_\text{energy}`` averaged over ``n_E`` configurations,
``\text{MSE}_\text{torque}`` over ``n_T = 3\,n_\text{atoms}\,n_E`` torque
components). The per-sample normalization keeps the two terms commensurate
regardless of system size.

The `estimator` selects how the least-squares problem is regularized:

| Estimator | Penalty |
|---|---|
| `OLS` | none — ordinary least squares |
| `Ridge(lambda)` | ``L_2`` |
| `Lasso(lambda)` | ``L_1`` (sparse coefficients) |
| `ElasticNet(alpha, lambda)` | mix of ``L_1`` and ``L_2`` |
| `AdaptiveLasso(pilot, lambda, gamma)` | weighted ``L_1`` (Zou, 2006) |
| `PrecomputedPilot(beta)` | adapter to reuse prior coefficients as an `AdaptiveLasso` pilot |

In these constructors ``lambda`` is the overall regularization strength;
``alpha`` (in ``[0, 1]``) is the ``L_1`` fraction of the `ElasticNet`
penalty, with ``alpha = 1`` recovering `Lasso`; ``gamma`` is the exponent
applied to the pilot coefficients when forming the adaptive-Lasso penalty
weights (Zou, 2006); ``pilot`` is the estimator that supplies those pilot
coefficients; and ``beta`` is a precomputed coefficient vector reused
directly as a pilot.

Regularization matters here because the SALC basis can be large relative
to the number of DFT configurations; `Lasso` / `AdaptiveLasso` additionally
drive small coefficients to zero, yielding a sparse, interpretable model.

## Result types

The fit returns an `SCEFit` (fitted coefficients `j0` / `jphi`, the
estimator, `torque_weight`, residuals; it implements the StatsAPI
`RegressionModel` interface — `coef`, `intercept`, `nobs`, `dof`). For
deployment, `SCEModel(f::SCEFit)` extracts the lightweight predictor —
just the basis and the coefficients — used by `predict_energy` and
`predict_torque`, and persisted via `Magesty.save` / `Magesty.load`.

Once `jphi` is fitted, the lowest-order coefficients can be read as
conventional spin-model parameters. That mapping — SCE coefficients into
Heisenberg ``J``, the Dzyaloshinskii–Moriya vector ``\vec{D}``, and the
anisotropic exchange tensor ``\boldsymbol{\Gamma}`` — is the subject of the
[Technical Notes](../technical_notes.md).
