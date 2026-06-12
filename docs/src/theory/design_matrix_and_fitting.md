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

In the kernel, the inner ``M_f`` sum over
``c_\nu^{M_f}\, T^{(L_f, M_f)}`` is not iterated per element: the two
factors are precontracted into a single rank-``(N)`` tensor at SALC build
time, so `design_matrix_energy_element` and
`_accumulate_grad_torque_cluster!` read that folded tensor directly. The math is the same — only the order of operations
changes. See [Folded tensor](folded_tensor.md) for the derivation and
storage layout.

The outer sum over "clusters in the orbit" is treated the same way: the
list of symmetry-equivalent clusters is enumerated once during SALC
construction and cached on each `CoupledBasis_with_coefficient`, so the
kernel iterates it directly instead of re-running the translation walk
and the sorted-tuple dedup per matrix element. See
[Orbit clusters](orbit_clusters.md) for the layout and the build-time
invariants.

The per-site tesseral harmonics
``Z_{l, m}(\hat{\boldsymbol{e}}_a)``, which *do* depend on the spin
configuration, are still constant across the SALCs, coupled bases, and
cluster images that share an atom and an ``(l, m)``. They are
precomputed once per spin configuration into an `SHCache` keyed by
``(l^2 + l + m + 1,\, \text{atom})`` and read by the kernel in place
of a per-call `Zₗₘ_unsafe` invocation. The torque path additionally
caches the direction gradients
``\partial_i Z_{l, m}(\hat{\boldsymbol{e}}_a)``. See
[Spherical-harmonic cache](sh_cache.md) for the storage layout, the
threading rationale, and the hot-path simplification it enables.

DFT also provides per-atom torques, and these constrain the same
coefficients. `build_design_matrix_torque` produces ``X_T`` with
``3 \times n_\text{atoms}`` rows per configuration. The torque is the
angular gradient of the energy,

```math
\boldsymbol{\tau}_i = \hat{\boldsymbol{e}}_i \times \nabla_{\hat{\boldsymbol{e}}_i} E,
```

which is why the spherical-harmonic gradient ``\partial Z_{l,m}/\partial\hat{\boldsymbol{e}}``
(`∂ᵢZlm`) is needed alongside the values.

The torque kernel is *cluster-major*: a single sweep over each cluster
image writes the gradient contribution to every site it contains, with
the cross product against each atom's spin direction deferred to a
final reduction. The per-atom kernel that the energy block naively
suggests would scan ``N_\text{atoms} / N`` clusters that contribute
nothing, which the cluster-major layout avoids. See
[Cluster-major torque kernel](cluster_major_torque.md) for the
derivation, the accumulation buffer's shape, and the cluster-distinctness
invariant the chain rule relies on.

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
| `AdaptiveRidge(lambda, epsilon)` | iteratively reweighted ``L_2`` — an ``L_0`` approximation (Frommlet & Nuel, 2016) |

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

```math
\textbf{AdaptiveRidge:}\qquad
\hat{J} = \arg\min_{J}\ \text{loss}(J) + \lambda\sum_{\nu} w_\nu\,J_\nu^2,
\qquad
w_\nu = \frac{1}{J_\nu^2 + \epsilon}.
```

`OLS` is `Ridge` at ``\lambda = 0``, and `Lasso` is `ElasticNet` at
``\alpha = 1``. In `AdaptiveLasso`, ``J^\text{pilot}`` is the coefficient
vector of a preliminary *pilot* fit: a large pilot coefficient gets a small
weight ``w_\nu`` and is penalized lightly, while a small one gets a large
weight and is pushed harder toward zero — the mechanism behind the oracle
property of Zou (2006). `PrecomputedPilot` has no objective of its own; it
simply supplies a stored ``J^\text{pilot}`` to an `AdaptiveLasso`.

`AdaptiveRidge` is iterative: its weights ``w_\nu`` depend on the very
coefficients being solved for, so they are found self-consistently. A plain
ridge fit initializes ``J``; then ``w_\nu`` and ``J`` are refit alternately
until ``J`` stops changing. Because the penalty term
``w_\nu J_\nu^2 = J_\nu^2 / (J_\nu^2 + \epsilon)`` saturates at ``1`` for
``J_\nu^2 \gg \epsilon`` and vanishes for ``J_\nu^2 \ll \epsilon``, the
converged penalty approximates ``\lambda`` times the number of nonzero
coefficients — an ``L_0`` penalty (Frommlet & Nuel, 2016). `AdaptiveRidge`
is `OLS` at ``\lambda = 0``.

In these constructors ``lambda`` (``\lambda``) is the overall
regularization strength; ``alpha`` (``\alpha \in [0, 1]``) is the ``L_1``
fraction of the `ElasticNet` penalty, with ``\alpha = 1`` recovering
`Lasso`; ``gamma`` (``\gamma \ge 0``, default ``1``) is the exponent in the
adaptive weight ``w_\nu``; ``pilot`` is the estimator that supplies
``J^\text{pilot}``; and ``beta`` is a precomputed coefficient vector reused
directly as that pilot. The parameter ``epsilon`` (``\epsilon``) appears in
both adaptive estimators but plays distinct roles: in `AdaptiveLasso`
(default machine epsilon) it floors the pilot magnitude so ``w_\nu`` stays
finite when a pilot coefficient is zero; in `AdaptiveRidge` (default
``10^{-8}``) it sets the magnitude scale ``\sqrt{\epsilon}`` below which a
coefficient is treated as negligible. `AdaptiveRidge` additionally takes
``max_iter`` (iteration cap, default ``50``) and ``tol`` (relative
coefficient-change convergence threshold, default ``10^{-6}``).

!!! note "Comparing `lambda` across estimators"
    `Ridge` and `AdaptiveRidge` add their penalty to the un-normalized
    loss and are solved analytically, whereas
    `ElasticNet` / `Lasso` / `AdaptiveLasso` are solved by GLMNet, whose
    objective divides the data term by ``2N`` and standardizes the columns
    of ``X`` internally. GLMNet also rescales the `AdaptiveLasso` weights
    so that they sum to the number of SALC key groups. Because of these
    differing conventions, a given numerical value of `lambda` does *not*
    carry the same meaning from one estimator to another; it must be tuned
    per estimator.

Regularization matters here because the SALC basis can be large relative
to the number of DFT configurations; `Lasso` / `AdaptiveLasso` produce
exactly sparse solutions, and `AdaptiveRidge` drives small coefficients
toward zero as an ``L_0`` approximation, yielding a sparse, interpretable
model.

### Post-selection refit

`Lasso` and `AdaptiveLasso` shrink the surviving coefficients toward
zero (the bias inherent in any ``L_1`` penalty), and `AdaptiveRidge`
likewise leaves them attenuated by its ``L_0`` approximation. Once the
support has been chosen, that shrinkage is no longer useful: a plain
`OLS` (or low-``\lambda`` `Ridge`) restricted to the selected columns
recovers an unbiased estimate of the coefficients on that support.
`refit` performs this post-selection refit:

```math
\text{support} = \{\,\nu : |J_\nu|\,\lVert X_{\cdot,\nu}\rVert > \tau\,\},
\qquad
J^\text{refit} = \arg\min_{J}\,
  \tfrac{1-w}{n_E}\,\lVert y_E - X_E J\rVert^2
+ \tfrac{w}{n_T}\,\lVert y_T - X_T J\rVert^2
\quad\text{s.t.}\ J_\nu = 0\ \forall\, \nu \notin \text{support},
```

with ``X`` the same energy-centered, row-whitened design matrix that the
original fit assembled (the input fit's `torque_weight` is reused). The
selection criterion ``|J_\nu|\,\lVert X_{\cdot,\nu}\rVert`` weights the
raw coefficient by its column norm to cancel the per-cluster
``(4\pi)^{N/2}`` factor that otherwise makes ``|J_\nu|`` magnitudes
incomparable across cluster sizes. At the default ``\tau = 0`` an
``L_1``-selected support survives intact (exact zeros are dropped, all
other coefficients kept), while an `AdaptiveRidge` user passes a
positive ``\tau`` to drop near-negligible bases. Coefficients of dropped
SALC key groups remain ``0`` in the returned `SCEFit`, so the
coefficient ordering is preserved.

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

## Cross-validation diagnostics

Generalized cross-validation (GCV) estimates the out-of-sample prediction
error of a fit *from the fit itself*, without a held-out set. It is
computed on the **same weighted, energy-centered augmented system** that
the fit minimizes — energy rows stacked above torque rows, row-whitened by
``\sqrt{(1-w)/n_E}`` and ``\sqrt{w/n_T}`` — so the score reflects the
combined energy+torque objective and the chosen `torque_weight`.

For a linear estimator the fitted values are a linear image of the
observations, ``\hat{y} = H y`` with a *hat matrix* ``H`` independent of
``y``, and

```math
\mathrm{GCV} = \frac{\lVert y - \hat{y}\rVert^2 / N}
                    {\bigl(1 - \mathrm{tr}(H)/N\bigr)^2},
```

where ``\mathrm{tr}(H)`` is the *effective degrees of freedom* and ``N`` is
the number of **live** rows. A block whose whitening scale is zero
(``\text{torque\_weight} = 1`` zeroes the energy rows;
``\text{torque\_weight} = 0`` zeroes the torque rows) contributes only dead
all-zero rows, which add nothing to the residual or to ``\mathrm{tr}(H)``;
they are excluded from ``N`` (so ``N = n_E + n_T`` for ``0 < w < 1``,
``N = n_T`` for ``w = 1``, ``N = n_E`` for ``w = 0``). The reference energy
``J_0`` was eliminated by the energy-block centering, so it contributes one
degree of freedom to ``\mathrm{tr}(H)`` **only when the energy block is
live** (``w < 1``). The score is in the weighted-objective unit, **not**
``\mathrm{eV}^2`` — it is meant for *comparing* configurations (across
penalties or data sizes), not for reading an absolute error.

To make a single value interpretable, normalize the GCV score against the
**null-model mean square** ``\mathrm{MSY} = \lVert y\rVert^2 / N`` — the
prediction error of the ``\beta = 0`` baseline on the centered, whitened system
(energy predicted at its mean, torque predicted as zero). The **GCV-based
predictive R²** is

```math
R^2_{\mathrm{gcv}} = 1 - \frac{\mathrm{GCV}}{\mathrm{MSY}},
```

which lives on a fixed scale: ``1`` is a perfect fit, ``0`` matches the null
model, and a negative value means the fit predicts worse than the null
(over-parameterized or too little data). Both terms share the weighted-objective
unit, so the ratio is dimensionless. `gcv_r2` returns it for a fit; the sweep
results carry it per-point (`GCVLambdaPath.gcv_r2`, `GCVSizeCurve.gcv_r2_mean`).
Because ``\mathrm{MSY}`` is fixed along a ``\lambda`` path, the GCV minimizer is
exactly the ``R^2_{\mathrm{gcv}}`` maximizer.

The effective dof depends on the estimator (the ``+1`` below is the
conditional intercept term, present only when the energy block is live):

- **`OLS`**: ``\mathrm{tr}(H) = 1 + \mathrm{rank}(X)``.
- **`Ridge`** with penalty ``\lambda``: from the singular values
  ``\sigma_k`` of ``X``,
  ``\mathrm{tr}(H) = 1 + \sum_k \sigma_k^2/(\sigma_k^2 + \lambda)``.
- **`AdaptiveRidge`**: conditional dof
  ``1 + \mathrm{tr}\bigl((X^\top X + \lambda D)^{-1} X^\top X\bigr)`` with
  the converged diagonal reweighting ``D``. The weights are recovered
  exactly from the fitted coefficients
  (``D_{\nu\nu} = 1/(J_\nu^2 + \varepsilon)``); only treating ``D`` as
  fixed in ``H`` is an approximation.

Non-linear estimators (`ElasticNet`, the `Lasso` alias, `AdaptiveLasso`)
have no exact hat matrix, so GCV is undefined for them and the entry
points raise an `ArgumentError`.

Two uses share this machinery. `gcv_lambda` sweeps the ridge penalty: a
single SVD of ``X`` yields ``\mathrm{tr}(H)`` and the residual norm for
every ``\lambda`` from the closed forms above, and the GCV minimizer is
returned as `lambda_best`. `gcv_learning_curve` sweeps the training-set size: at
each size it draws several random subsets (reusing the prebuilt design
matrix by row slicing, with a seeded RNG), fits each, and averages their
GCV scores. A curve that flattens with size indicates that enough training
data is present. `gcv` returns the single score for an existing fit, and
`gcv_r2` its predictive-R² companion.

The sweep results carry to text via `write_gcv_lambda` /
`write_gcv_learning_curve` for the `FitCheck_gcv_lambda.py` /
`FitCheck_gcv_learning_curve.py` plotters under `tools/`.

## References

1. H. Zou, "The Adaptive Lasso and Its Oracle Properties",
   *J. Am. Stat. Assoc.* **101**, 1418 (2006).
   DOI: [10.1198/016214506000000735](https://doi.org/10.1198/016214506000000735)
2. F. Frommlet and G. Nuel, "An Adaptive Ridge Procedure for L0
   Regularization", *PLoS ONE* **11**, e0148620 (2016).
   DOI: [10.1371/journal.pone.0148620](https://doi.org/10.1371/journal.pone.0148620)
3. G. H. Golub, M. Heath, and G. Wahba, "Generalized Cross-Validation as a
   Method for Choosing a Good Ridge Parameter", *Technometrics* **21**, 215
   (1979). DOI: [10.1080/00401706.1979.10489751](https://doi.org/10.1080/00401706.1979.10489751)
