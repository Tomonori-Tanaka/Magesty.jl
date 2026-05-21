# Design matrix and fitting

The symmetry-adapted basis fixes *which* functions ``\Phi_\nu`` appear in
the expansion. This final page turns those functions into numbers — a
design matrix — and fits the coefficients ``J_\nu``. The relevant module
is `src/Fitting.jl`.

## The design feature

Recall the SCE energy of a configuration ``\{\hat{\boldsymbol{e}}_i\}``:

```math
E = J_0 + \sum_{\nu} J_\nu \, \Phi_\nu\!\bigl(\{\hat{\boldsymbol{e}}_i\}\bigr).
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
  Z_{1,m_1}(\hat{\boldsymbol{e}}_i)\,Z_{1,m_2}(\hat{\boldsymbol{e}}_j).
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
\boldsymbol{\tau}_i = \hat{\boldsymbol{e}}_i \times \nabla_{\hat{\boldsymbol{e}}_i} E,
```

which is why the spherical-harmonic gradient ``\partial Z_{l,m}/\partial\hat{\boldsymbol{e}}``
(`∂ᵢZlm`) is needed alongside the values.

The constant ``J_0`` is *not* a column of ``X_E``. The energy matrix is
mean-centered instead, and ``J_0`` is recovered analytically after the
solve by `extract_j0_jphi` as ``J_0 = \operatorname{mean}(y_E - X_E\,J)``,
where ``y_E`` is the vector of observed (DFT) energies and ``J = (J_\nu)``
the fitted coefficient vector.

In other words: once the SCE coefficients ``J`` have been fixed by the fit,
``J_0`` is the value that minimizes the energy misfit for that ``J``.
Minimizing the sum of squares
``\sum_c \bigl(y_E^{(c)} - J_0 - (X_E\,J)_c\bigr)^2`` over the single
parameter ``J_0`` (set its derivative to zero) gives exactly the mean
residual ``\operatorname{mean}(y_E - X_E\,J)`` — an ordinary least-squares
estimate of ``J_0`` with ``J`` held fixed. Because the loss is quadratic in
``J_0``, this is not a two-stage approximation: substituting
``J_0 = J_0(J)`` back and fitting ``J`` on the mean-centered system yields
the same ``(J, J_0)`` as minimizing the loss jointly over both.

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

where each term is a *per-sample* mean squared error,

```math
\text{MSE}_\text{energy}
  = \frac{1}{n_E}\sum_{c=1}^{n_E}
    \bigl(y_E^{(c)} - J_0 - (X_E\,J)_c\bigr)^2,
\qquad
\text{MSE}_\text{torque}
  = \frac{1}{n_T}\sum_{r=1}^{n_T}
    \bigl(y_T^{(r)} - (X_T\,J)_r\bigr)^2 .
```

Here ``y_E^{(c)}`` is the DFT energy of configuration ``c``, ``y_T^{(r)}``
the ``r``-th observed torque component, ``J = (J_\nu)`` the coefficient
vector, and ``(X_E\,J)_c`` / ``(X_T\,J)_r`` the model predictions.
Averaging over ``n_E`` configurations and over
``n_T = 3\,n_\text{atoms}\,n_E`` torque components — rather than summing —
keeps the two terms commensurate regardless of system size.

The intercept ``J_0`` is eliminated analytically (the mean-centering of
the previous section), so the solve is over ``J`` alone; write
``\text{loss}(J)`` for the loss with
``J_0 = \operatorname{mean}(y_E - X_E\,J)`` substituted back. Every
estimator first reduces the problem to a single *augmented*
least-squares system — the mean-centered energy rows and the torque rows
are stacked and row-weighted so that

```math
\lVert X\,J - y\rVert^2 = \text{loss}(J),
\qquad N = n_E + n_T,
```

with ``X`` the augmented design matrix, ``y`` the augmented observation
vector, and ``N`` its total row count.

The `estimator` then selects how that loss is regularized:

| Estimator | Penalty |
|---|---|
| `OLS` | none — ordinary least squares |
| `Ridge(lambda)` | ``L_2`` |
| `Lasso(lambda)` | ``L_1`` (sparse coefficients) |
| `ElasticNet(alpha, lambda)` | mix of ``L_1`` and ``L_2`` |
| `AdaptiveLasso(pilot, lambda, gamma)` | weighted ``L_1`` (Zou, 2006) |
| `PrecomputedPilot(beta)` | adapter to reuse prior coefficients as an `AdaptiveLasso` pilot |

The fitted coefficient vector ``\hat{J}`` minimizes the loss plus the
estimator's penalty. With ``\lVert J\rVert_1 = \sum_\nu |J_\nu|`` and
``\lVert J\rVert_2^2 = \sum_\nu J_\nu^2``, the objective of each estimator
is:

```math
\textbf{OLS:}\qquad
\hat{J} = \arg\min_{J}\ \text{loss}(J),
```

```math
\textbf{Ridge:}\qquad
\hat{J} = \arg\min_{J}\ \text{loss}(J) + \lambda\,\lVert J\rVert_2^2,
```

```math
\textbf{ElasticNet:}\qquad
\hat{J} = \arg\min_{J}\ \frac{1}{2N}\,\text{loss}(J)
  + \lambda\sum_{\nu}\Bigl[\tfrac{1-\alpha}{2}\,J_\nu^2 + \alpha\,|J_\nu|\Bigr],
```

```math
\textbf{Lasso}\ (\alpha = 1):\qquad
\hat{J} = \arg\min_{J}\ \frac{1}{2N}\,\text{loss}(J)
  + \lambda\,\lVert J\rVert_1,
```

```math
\textbf{AdaptiveLasso:}\qquad
\hat{J} = \arg\min_{J}\ \frac{1}{2N}\,\text{loss}(J)
  + \lambda\sum_{\nu} w_\nu\,|J_\nu|,
\qquad
w_\nu = \frac{1}{\max\bigl(|J_\nu^\text{pilot}|,\ \epsilon\bigr)^{\gamma}}.
```

`OLS` is `Ridge` at ``\lambda = 0``, and `Lasso` is `ElasticNet` at
``\alpha = 1``. In `AdaptiveLasso`, ``J^\text{pilot}`` is the coefficient
vector of a preliminary *pilot* fit: a large pilot coefficient gets a small
weight ``w_\nu`` and is penalized lightly, while a small one gets a large
weight and is pushed harder toward zero — the mechanism behind the oracle
property of Zou (2006). `PrecomputedPilot` has no objective of its own; it
simply supplies a stored ``J^\text{pilot}`` to an `AdaptiveLasso`.

In these constructors ``lambda`` (``\lambda``) is the overall
regularization strength; ``alpha`` (``\alpha \in [0, 1]``) is the ``L_1``
fraction of the `ElasticNet` penalty, with ``\alpha = 1`` recovering
`Lasso`; ``gamma`` (``\gamma \ge 0``, default ``1``) is the exponent in the
adaptive weight ``w_\nu``; ``epsilon`` (``\epsilon``, default machine
epsilon) floors the pilot magnitude so ``w_\nu`` stays finite when a pilot
coefficient is zero; ``pilot`` is the estimator that supplies
``J^\text{pilot}``; and ``beta`` is a precomputed coefficient vector reused
directly as that pilot.

!!! note "Comparing `lambda` across estimators"
    `Ridge` adds its penalty to the un-normalized loss, whereas
    `ElasticNet` / `Lasso` / `AdaptiveLasso` are solved by GLMNet, whose
    objective divides the data term by ``2N`` and standardizes the columns
    of ``X`` internally. GLMNet also rescales the `AdaptiveLasso` weights
    so that they sum to the number of SALC key groups. Because of these
    differing conventions, a given numerical value of `lambda` does *not*
    carry the same meaning from one estimator to another; it must be tuned
    per estimator.

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

## References

1. H. Zou, "The Adaptive Lasso and Its Oracle Properties",
   *J. Am. Stat. Assoc.* **101**, 1418 (2006).
   DOI: [10.1198/016214506000000735](https://doi.org/10.1198/016214506000000735)
