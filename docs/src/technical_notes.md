# Technical notes

## Conversion between the lowest SCE model and conventional spin model

### Setup

The SCE model energy per supercell is

$$E = J_0 + \sum_\nu J_\nu \,\Phi_\nu\!\left(\{\hat{e}_i\}\right),$$

where $J_0$ = `optimizer.reference_energy` and $J_\nu$ = `optimizer.SCE[ν]`.

The design feature $\Phi_\nu$ for a 2-body SALC is computed as

$$\Phi_\nu = (4\pi) \sum_{M_f} c_\nu^{M_f} \!\!\!\sum_{\text{clusters in orbit}} \;\sum_{m_1,m_2} T^{(L_f,M_f)}_{m_1 m_2}\; Z_{1,m_1}(\hat{e}_i)\, Z_{1,m_2}(\hat{e}_j),$$

where $T^{(L_f,M_f)}$ is the real (tesseral) coupled angular momentum tensor, $c_\nu^{M_f}$ is the SALC coefficient vector, and the sum runs over all clusters in the symmetry orbit.

### Tesseral harmonics for $l = 1$

The $l = 1$ tesseral harmonics evaluated on a unit spin vector $\hat{e} = (e_x, e_y, e_z)$ are

$$Z_{1,m}(\hat{e}) = \sqrt{\frac{3}{4\pi}}\,\hat{e}_{\mu(m)}, \qquad
\mu(m) = \begin{cases} y & m = -1 \\ z & m = 0 \\ x & m = +1 \end{cases}$$

### Design features for $l_1 = l_2 = 1$

The three channels $L_f = 0, 1, 2$ arising from coupling $l_1 = l_2 = 1$ yield the following design features for a single pair cluster $(i, j)$:

#### $L_f = 0$ — Heisenberg exchange

$$\boxed{\Phi^{(0)} = \sqrt{3}\;\hat{e}_i \cdot \hat{e}_j}$$

#### $L_f = 1$ — Dzyaloshinskii–Moriya interaction

The three components ($M_f = -1,0,+1$ correspond to $\mu = y, z, x$, respectively):

$$\boxed{\Phi^{(1)}_\mu = \frac{3}{\sqrt{2}}\,(\hat{e}_i \times \hat{e}_j)_\mu}$$

#### $L_f = 2$ — Anisotropic symmetric exchange

The five components are

| $M_f$ | $\Phi^{(2)}_{M_f}$ |
|:-----:|:-------------------|
| $-2$ | $\dfrac{3}{\sqrt{2}}(e_{ix}\,e_{jy}+e_{iy}\,e_{jx})$ |
| $-1$ | $\dfrac{3}{\sqrt{2}}(e_{iy}\,e_{jz}+e_{iz}\,e_{jy})$ |
|  $0$ | $\sqrt{6}\!\left[e_{iz}\,e_{jz}-\dfrac{1}{2}(e_{ix}\,e_{jx}+e_{iy}\,e_{jy})\right]$ |
| $+1$ | $\dfrac{3}{\sqrt{2}}(e_{ix}\,e_{jz}+e_{iz}\,e_{jx})$ |
| $+2$ | $\dfrac{3}{\sqrt{2}}(e_{ix}\,e_{jx}-e_{iy}\,e_{jy})$ |

> **Verification** (from `test/examples/dimer/test.jl`):
> - $L_f=0$: `SCE[1]`$\times\sqrt{3} = -1$ for FM/AFM energies $\mp 1$ ✓  
> - $L_f=1$: `SCE[1]`$\times 3/\sqrt{2} = -1$ for planar DMI with $D_z=-1$ ✓

### Conversion to conventional spin model parameters

The bilinear exchange Hamiltonian for a pair $(i,j)$ is often written as

$$E_{ij} = J\,\hat{e}_i\cdot\hat{e}_j \;+\; \vec{D}_{ij}\cdot(\hat{e}_i\times\hat{e}_j) \;+\; \hat{e}_i^{\,\top}\boldsymbol{\Gamma}_{ij}\,\hat{e}_j,$$

where $\boldsymbol{\Gamma}_{ij}$ is symmetric and traceless.

The correspondence with SCE coefficients is:

$$\boxed{J \;=\; \sqrt{3}\;J^{\mathrm{SCE}}_{L_f=0}}$$

$$\boxed{\vec{D} \;=\; \frac{3}{\sqrt{2}}\;\vec{J}^{\mathrm{SCE}}_{L_f=1}}$$

where $\vec{J}^{\mathrm{SCE}}_{L_f=1} = \bigl(J^{\mathrm{SCE}}_{M_f=+1},\, J^{\mathrm{SCE}}_{M_f=-1},\, J^{\mathrm{SCE}}_{M_f=0}\bigr)$ maps to $(D_x, D_y, D_z)$.

The symmetric traceless tensor $\boldsymbol{\Gamma}$ is assembled from the five $L_f=2$ coefficients (abbreviating $J_m \equiv J^{\mathrm{SCE}}_{L_f=2,\,M_f=m}$):

$$\boldsymbol{\Gamma} = \begin{pmatrix}
-\dfrac{\sqrt{6}}{2}\,J_0 + \dfrac{3}{\sqrt{2}}\,J_{+2} & \dfrac{3}{\sqrt{2}}\,J_{-2} & \dfrac{3}{\sqrt{2}}\,J_{+1} \\[6pt]
\dfrac{3}{\sqrt{2}}\,J_{-2} & -\dfrac{\sqrt{6}}{2}\,J_0 - \dfrac{3}{\sqrt{2}}\,J_{+2} & \dfrac{3}{\sqrt{2}}\,J_{-1} \\[6pt]
\dfrac{3}{\sqrt{2}}\,J_{+1} & \dfrac{3}{\sqrt{2}}\,J_{-1} & \sqrt{6}\,J_0
\end{pmatrix}$$

One can verify that $\mathrm{Tr}(\boldsymbol{\Gamma}) = 0$ and $\boldsymbol{\Gamma}^{\!\top} = \boldsymbol{\Gamma}$.

### Full bilinear exchange tensor

Combining all three channels, the full exchange tensor $\mathbf{J}_{ij}$ (so that $E_{ij} = \hat{e}_i^{\,\top}\mathbf{J}_{ij}\,\hat{e}_j$) is

$$\mathbf{J}_{ij} = J\,\mathbf{I} + \begin{pmatrix} 0 & D_z & -D_y \\ -D_z & 0 & D_x \\ D_y & -D_x & 0 \end{pmatrix} + \boldsymbol{\Gamma}_{ij},$$

where $J = \sqrt{3}\,J^{\mathrm{SCE}}_0$, $\vec{D} = \frac{3}{\sqrt{2}}\vec{J}^{\mathrm{SCE}}_1$, and $\boldsymbol{\Gamma}$ is given above.

> **Note on sign convention**: In the literature the Heisenberg term is often written as $-J\hat{e}_i\cdot\hat{e}_j$. With that convention replace $J \to -J$ in the formula above, giving $J_{\mathrm{Heis}} = -\sqrt{3}\,J^{\mathrm{SCE}}_0$.

### Summary table

| Physical parameter | Conversion from SCE |
|:---|:---|
| Isotropic exchange $J$ ($E = J\,\hat{e}_i\cdot\hat{e}_j$) | $J = \sqrt{3}\;J^{\mathrm{SCE}}_{L_f=0}$ |
| DM vector component $D_\mu$ | $D_\mu = \dfrac{3}{\sqrt{2}}\;J^{\mathrm{SCE}}_{L_f=1,\,\mu}$ |
| Anisotropic exchange $\Gamma_{zz}$ | $\Gamma_{zz} = \sqrt{6}\;J^{\mathrm{SCE}}_{L_f=2,\,M_f=0}$ |
| Anisotropic exchange $\Gamma_{xz}$ | $\Gamma_{xz} = \dfrac{3}{\sqrt{2}}\;J^{\mathrm{SCE}}_{L_f=2,\,M_f=+1}$ |
| Anisotropic exchange $\Gamma_{yz}$ | $\Gamma_{yz} = \dfrac{3}{\sqrt{2}}\;J^{\mathrm{SCE}}_{L_f=2,\,M_f=-1}$ |
| Anisotropic exchange $\Gamma_{xx}-\Gamma_{yy}$ | $\Gamma_{xx}-\Gamma_{yy} = \dfrac{3}{\sqrt{2}}\;J^{\mathrm{SCE}}_{L_f=2,\,M_f=+2}$ |
| Anisotropic exchange $\Gamma_{xy}$ | $\Gamma_{xy} = \dfrac{3}{\sqrt{2}}\;J^{\mathrm{SCE}}_{L_f=2,\,M_f=-2}$ |

The $M_f$-to-Cartesian mapping for the tesseral basis follows: $M_f = +1 \leftrightarrow x$, $M_f = -1 \leftrightarrow y$, $M_f = 0 \leftrightarrow z$.
