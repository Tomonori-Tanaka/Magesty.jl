# Mean-Field Analysis: Finding the Ordering Wave Vector

Theoretical background for `tools/mfa_analysis.jl`.

---

## 1. Spin Hamiltonian

When the Magesty SCE model is restricted to pairwise exchange interactions (i.e., `[interaction] nbody = 2` and `[interaction.body2] lsum = 2`), the Hamiltonian is written in the following full tensorial form:

$$H = \sum_{i < j} \boldsymbol{S}_i^\top \mathbf{J}_{ij} \, \boldsymbol{S}_j$$

where $\mathbf{J}_{ij}$ is the $3 \times 3$ exchange tensor (unit: meV) returned by `convert2tensor`.
It can be decomposed into three physically distinct contributions:

$$\mathbf{J}_{ij} = \underbrace{J_{ij}^\text{iso} \mathbf{I}}_{\text{isotropic}} + \underbrace{\mathbf{W}_{ij}}_{\text{antisymmetric (DMI)}} + \underbrace{\boldsymbol{\Gamma}_{ij}}_{\text{anisotropic symmetric}}$$

| Component | Definition | Physical origin |
|---|---|---|
| $J_{ij}^\text{iso} = \tfrac{1}{3}\operatorname{tr}(\mathbf{J}_{ij})$ | scalar | Heisenberg exchange |
| $\mathbf{W}_{ij} = \tfrac{1}{2}(\mathbf{J}_{ij} - \mathbf{J}_{ij}^\top)$ | antisymmetric | Dzyaloshinskii–Moriya interaction (DMI), $\boldsymbol{D}_{ij} \cdot (\boldsymbol{S}_i \times \boldsymbol{S}_j)$ |
| $\boldsymbol{\Gamma}_{ij} = \tfrac{1}{2}(\mathbf{J}_{ij} + \mathbf{J}_{ij}^\top) - J_{ij}^\text{iso}\mathbf{I}$ | symmetric traceless | Anisotropic symmetric exchange |

The sign convention is $J_{ij}^\text{iso} < 0$ for ferromagnetic coupling.

> **Note on double counting**  
> `convert2tensor` uses the no-double-counting convention (sum over $i < j$).
> The with-double-counting convention uses $\mathbf{J}_{ij}^\text{(DC)} = \mathbf{J}_{ij}/2$.
> Both conventions yield the same MFA transition temperature.

---

## 2. Fourier Transform and the $J(q)$ Matrix

For a system with $n_\text{sub}$ sublattices, $J(q)$ is a $3n_\text{sub} \times 3n_\text{sub}$ block matrix.
Each block is a $3\times 3$ matrix that retains the **full tensorial** structure of $\mathbf{J}_{ij}$,
including isotropic, DMI, and anisotropic symmetric contributions:

$$\mathbf{J}_{\mu\nu}(\boldsymbol{q}) = \sum_{\substack{j \in \text{sublattice } \nu \\ j \neq \text{prim}_\mu}} \mathbf{J}(\text{prim}_\mu,\, j) \, e^{2\pi i \, \boldsymbol{q} \cdot \boldsymbol{R}_{\mu j}}$$

- $\mu, \nu$: sublattice indices ($1, \ldots, n_\text{sub}$)
- $\text{prim}_\mu$: representative atom of sublattice $\mu$ in the primitive cell (`symmetry.atoms_in_prim[μ]`)
- $\boldsymbol{R}_{\mu j}$: minimum-image displacement in supercell fractional coordinates

$$\boldsymbol{R}_{\mu j} = \boldsymbol{x}_j^\text{frac} - \boldsymbol{x}_{\text{prim}_\mu}^\text{frac} - \operatorname{round}\!\left(\boldsymbol{x}_j^\text{frac} - \boldsymbol{x}_{\text{prim}_\mu}^\text{frac}\right)$$

- $\boldsymbol{q}$: wave vector in supercell fractional reciprocal coordinates

The sublattice index of atom $j$ is identified via the symmetry mapping:

$$\nu = \texttt{symmetry.map\_s2p[j].atom} \quad (1 \leq \nu \leq n_\text{sub})$$

---

## 3. MFA Ground State and Spin Spirals

Assuming a spin-spiral state with ordering vector $\boldsymbol{q}$, the mean-field effective field on sublattice $\mu$ is

$$\boldsymbol{h}_{\text{eff},\mu} = -\sum_{\nu} \mathbf{J}_{\mu\nu}(\boldsymbol{q}) \langle \boldsymbol{S}_\nu \rangle$$

Because $\mathbf{J}_{\mu\nu}(\boldsymbol{q})$ is a full $3\times 3$ tensor, the MFA eigenvalue problem is formulated on the $3n_\text{sub} \times 3n_\text{sub}$ Hermitian matrix

$$\mathcal{J}(\boldsymbol{q}) = \frac{\mathbf{J}(\boldsymbol{q}) + \mathbf{J}(\boldsymbol{q})^\dagger}{2}$$

The ground-state energy density is determined by the minimum eigenvalue of $\mathcal{J}(\boldsymbol{q})$:

$$E(\boldsymbol{q}) \propto S^2 \cdot \lambda_\text{min}(\mathcal{J}(\boldsymbol{q}))$$

Each eigenvalue corresponds to a spin-polarization mode (mixing isotropic, DMI-driven, and anisotropic channels).
The optimal ordering wave vector $\boldsymbol{q}_\text{opt}$ minimizes $\lambda_\text{min}$ (makes it most negative):

$$\boldsymbol{q}_\text{opt} = \underset{\boldsymbol{q}}{\operatorname{argmin}} \left[ \lambda_\text{min}(\mathcal{J}(\boldsymbol{q})) \right]$$

| $\boldsymbol{q}_\text{opt}$ | Ordering type |
|---|---|
| $(0,0,0)$ | Ferromagnetic (FM) |
| $(1/2, 1/2, 1/2)$, etc. | Antiferromagnetic (AFM) |
| Other rational or irrational values | Spin spiral |

---

## 4. Classical-Spin MFA Transition Temperature

Linearizing the classical-spin (Langevin) self-consistency equation gives

$$k_B T_C = \frac{S^2}{3} \left| \lambda_\text{min}(\mathcal{J}(\boldsymbol{q}_\text{opt})) \right|$$

- $S$: spin magnitude (classical normalization: $|\boldsymbol{S}| = S$)
- $k_B = 8.617333 \times 10^{-2}$ meV/K

> For quantum spins (quantum Heisenberg model), replace $S^2 \to S(S+1)$.

---

## 5. Wavevector Units and Wavelength

Let $\{\boldsymbol{a}_i^\text{prim}\}$ be primitive lattice vectors (Å), and $\boldsymbol{b}_i$ be reciprocal vectors with
$\boldsymbol{b}_i \cdot \boldsymbol{a}_j = \delta_{ij}$ (unit Å$^{-1}$).

First compute

$$\boldsymbol{q}_\text{cart} = \sum_i q_i^\text{prim} \, \boldsymbol{b}_i \quad (\text{unit: } 1/\text{Å})$$

Then define the angular wavevector

$$\boldsymbol{k}_\text{cart} = 2\pi \, \boldsymbol{q}_\text{cart} \quad (\text{unit: rad/Å})$$

The script prints:

- `Ordering vector |q| (Cartesian)`: $|\boldsymbol{q}_\text{cart}|$ reported with unit $2\pi/\text{Å}$ (numerically $|\boldsymbol{k}_\text{cart}|/(2\pi)$),
- `Angular wavevector |k|`: $|\boldsymbol{k}_\text{cart}|$ in rad/Å.

The wavelength is

$$\lambda = \frac{2\pi}{|\boldsymbol{k}_\text{cart}|} = \frac{1}{|\boldsymbol{q}_\text{cart}|}$$

For $|\boldsymbol{q}|=0$ (ferromagnetic), $\lambda = \infty$.

---

## 6. Brillouin Zone Scan

The primitive BZ is sampled on an $n_k^3$ uniform grid. With $\boldsymbol{q}^\text{prim} \in [0,1)^3$ as fractional coordinates of the primitive reciprocal lattice, the conversion to supercell fractional coordinates is

$$\boldsymbol{q}_\text{sc} = A_\text{prim} \, \boldsymbol{q}_\text{prim}$$

where $A_\text{prim}$ is the $3 \times 3$ matrix whose columns are the primitive lattice vectors in supercell fractional coordinates, constructed from the pure-translation symmetry operations.

The Fourier phase used in $J_{\mu\nu}(\boldsymbol{q})$ is

$$\text{phase} = \exp(2\pi i \, \boldsymbol{q}_\text{sc} \cdot \boldsymbol{R}_{\mu j})$$

---

## 7. Usage

```bash
julia --project=. tools/mfa_analysis.jl \
    --xml path/to/jphi.xml \
    --nk 20 \
    --spin 1.0 \
    --eigvec 1
```

| Option | Description | Default |
|---|---|---|
| `--xml` / `-x` | Path to an XML file containing SCE model information | (required) |
| `--nk` / `-n` | Number of $q$-points per reciprocal direction (positive integer) | 20 |
| `--spin` / `-s` | Spin magnitude $S$ (if omitted: classical-spin estimate; if specified: quantum-spin estimate with $S(S+1)$) | omitted |
| `--eigvec` / `-e` | Print the $N$-th eigenvector of $\mathcal{J}(\boldsymbol{q}_\text{opt})$ as a spin configuration ($N=1$: minimum eigenvalue, $N=2$: second smallest, …) | omitted |
| `--no-refine` | Skip gradient-descent refinement after the grid scan | false |

---

## 8. Reading the Eigenvectors: Spin Configuration at $\boldsymbol{q}_\text{opt}$

### What is the eigenvector $\boldsymbol{v}$?

After the optimal ordering wave vector $\boldsymbol{q}_\text{opt}$ is found, `mfa_analysis` can compute the eigenvectors of $\mathcal{J}(\boldsymbol{q}_\text{opt})$.
Each eigenvector $\boldsymbol{v}$ is a $3n_\text{sub}$-dimensional complex vector.
Its components are grouped by sublattice:

$$\boldsymbol{v} = \bigl(\underbrace{v_x^{(1)},\, v_y^{(1)},\, v_z^{(1)}}_{\text{sublattice 1}},\;\underbrace{v_x^{(2)},\, v_y^{(2)},\, v_z^{(2)}}_{\text{sublattice 2}},\;\ldots\bigr)$$

The 3-component slice $\boldsymbol{v}_\mu = (v_x^{(\mu)}, v_y^{(\mu)}, v_z^{(\mu)})$ describes the spin state of sublattice $\mu$.

### How does $\boldsymbol{v}_\mu$ map to real-space spins?

MFA assumes a spin-spiral (plane-wave) ansatz.
The spin on sublattice $\mu$ at unit-cell position $\boldsymbol{R}$ is

$$\boldsymbol{S}_\mu(\boldsymbol{R}) = \mathrm{Re}\!\left[\boldsymbol{v}_\mu\, e^{2\pi i\,\boldsymbol{q}\cdot\boldsymbol{R}}\right]$$

Expanding $\boldsymbol{v}_\mu = \boldsymbol{a}_\mu + i\,\boldsymbol{b}_\mu$ and using Euler's formula:

$$\boldsymbol{S}_\mu(\boldsymbol{R}) = \boldsymbol{a}_\mu\cos(2\pi\boldsymbol{q}\cdot\boldsymbol{R}) - \boldsymbol{b}_\mu\sin(2\pi\boldsymbol{q}\cdot\boldsymbol{R})$$

where $\boldsymbol{a}_\mu = \mathrm{Re}[\boldsymbol{v}_\mu]$ and $\boldsymbol{b}_\mu = \mathrm{Im}[\boldsymbol{v}_\mu]$.

The physical meaning of each quantity is summarised below.

| Quantity | Symbol | Physical meaning |
|---|---|---|
| Real part | $\boldsymbol{a}_\mu = \mathrm{Re}[\boldsymbol{v}_\mu]$ | Spin direction at $\boldsymbol{q}\cdot\boldsymbol{R}=0$ (the cosine amplitude) |
| Imaginary part | $\boldsymbol{b}_\mu = \mathrm{Im}[\boldsymbol{v}_\mu]$ | Spin direction at $\boldsymbol{q}\cdot\boldsymbol{R}=\tfrac{1}{4}$ (the sine amplitude, $\tfrac{1}{4}$ period away) |
| Magnitude | $|\boldsymbol{v}_\mu|$ | Overall spin amplitude on sublattice $\mu$ (proportional to the ordered moment) |

**Intuition.**
Think of the phase $\varphi = 2\pi\boldsymbol{q}\cdot\boldsymbol{R}$ as a clock hand that advances as you move along the $\boldsymbol{q}$ direction.
At each unit cell the spin is the projection of the complex vector $\boldsymbol{v}_\mu$ onto the real axis after rotating by $\varphi$:

| Phase $\varphi$ | Spin direction |
|---|---|
| $0$ | $+\boldsymbol{a}_\mu$ |
| $\pi/2$ | $-\boldsymbol{b}_\mu$ |
| $\pi$ | $-\boldsymbol{a}_\mu$ |
| $3\pi/2$ | $+\boldsymbol{b}_\mu$ |

As $\varphi$ sweeps from $0$ to $2\pi$, the spin tip traces an **ellipse** in the plane spanned by $\boldsymbol{a}_\mu$ and $\boldsymbol{b}_\mu$.

---

## 9. Identifying the Magnetic Structure from $\boldsymbol{v}_\mu$

The shape and orientation of the ellipse traced by $\boldsymbol{S}_\mu(\boldsymbol{R})$ reveals the type of magnetic order.

### Diagnostic quantities

For each sublattice $\mu$, define $\boldsymbol{a}_\mu = \mathrm{Re}[\boldsymbol{v}_\mu]$ and $\boldsymbol{b}_\mu = \mathrm{Im}[\boldsymbol{v}_\mu]$, and compute:

$$\text{ellipticity} = \frac{\min(|\boldsymbol{a}_\mu|,\,|\boldsymbol{b}_\mu|)}{\max(|\boldsymbol{a}_\mu|,\,|\boldsymbol{b}_\mu|)}, \qquad \text{non-orthogonality} = \frac{|\boldsymbol{a}_\mu \cdot \boldsymbol{b}_\mu|}{|\boldsymbol{a}_\mu||\boldsymbol{b}_\mu|}$$

For non-collinear structures where $\boldsymbol{a}_\mu \perp \boldsymbol{b}_\mu$, the spin rotation plane has a well-defined normal

$$\hat{n}_\mu = \frac{\boldsymbol{a}_\mu \times \boldsymbol{b}_\mu}{|\boldsymbol{a}_\mu \times \boldsymbol{b}_\mu|}$$

The relationship between $\hat{n}_\mu$ and the propagation direction $\hat{q}$ distinguishes the two canonical helical types:

| Structure type | $|\boldsymbol{b}_\mu|$ | $\boldsymbol{a}_\mu \cdot \boldsymbol{b}_\mu$ | Ellipticity | $\hat{n}_\mu$ vs $\hat{q}$ | Description |
|---|---|---|---|---|---|
| **Collinear** (FM/AFM) | $\approx 0$ | — | $0$ | — | Spins point along a fixed axis; no transverse oscillation |
| **Proper screw** (Bloch-type) | $\approx |\boldsymbol{a}_\mu|$ | $\approx 0$ | $\approx 1$ | $\hat{n}_\mu \parallel \hat{q}$ | Rotation plane ⊥ propagation direction; typical of B20 compounds along $\langle 111 \rangle$ |
| **Cycloid** (Néel-type) | $\approx |\boldsymbol{a}_\mu|$ | $\approx 0$ | $\approx 1$ | $\hat{n}_\mu \perp \hat{q}$ | Propagation direction lies within the rotation plane; driven by bond-direction DM interaction |
| **Elliptical helix** | $0 < |\boldsymbol{b}_\mu| < |\boldsymbol{a}_\mu|$ | $\approx 0$ | $0 < \varepsilon < 1$ | either | Spins precess on an ellipse; arises when easy-axis anisotropy competes with DMI |
| **Conical helix** | $> 0$ | $\neq 0$ | — | — | $\boldsymbol{a}_\mu$ and $\boldsymbol{b}_\mu$ are not orthogonal; spins lie on a cone |

### Example: B20 compound

The output below is from a B20-type material with $q \approx [0.073,\,0,\,0]$ (propagation along $x$).
Sublattices $\mu = 1$–$4$ are the magnetically active sites.

```
sublattice μ |  Re[vx]   Re[vy]   Re[vz]  |  Im[vx]   Im[vy]   Im[vz]  | |v_μ|
μ =  1       | -0.0010  +0.0031  +0.3480  | +0.3500  -0.0030  +0.0007  | 0.4936
```

Extracting the key quantities for $\mu = 1$:

$$\boldsymbol{a}_1 \approx (0,\;0,\;+0.348) \quad\text{(spin along }z\text{)}$$
$$\boldsymbol{b}_1 \approx (+0.350,\;0,\;0) \quad\text{(spin along }x\text{)}$$
$$|\boldsymbol{a}_1| \approx 0.348,\quad |\boldsymbol{b}_1| \approx 0.350,\quad \boldsymbol{a}_1\cdot\boldsymbol{b}_1 \approx 0$$

- **Ellipticity** $= 0.348/0.350 \approx 0.994 \approx 1$ → nearly circular
- **Non-orthogonality** $\approx 0$ → $\boldsymbol{a}_1 \perp \boldsymbol{b}_1$
- **Rotation plane normal**: $\hat{n}_1 = \boldsymbol{a}_1 \times \boldsymbol{b}_1 / |\cdots| \approx \hat{z} \times \hat{x} = \hat{y}$
- $\hat{n}_1 \approx \hat{y} \perp \hat{q} = \hat{x}$ → the propagation direction lies **within** the rotation plane

**Conclusion:** the spin traces a circle in the $xz$-plane (which contains $\boldsymbol{q}$) — a **cycloid (Néel-type helix)**.
The rotation plane normal $\hat{n} \approx \hat{y}$ is perpendicular to $\hat{q} = \hat{x}$, the defining signature of a cycloid.

The near three-fold degeneracy of the lowest eigenvalues ($\lambda_1 \approx \lambda_2 \approx \lambda_3$) reflects the nearly isotropic energy landscape in cubic B20: the minimum lies on a sphere of radius $|\boldsymbol{q}_\text{opt}|$ in reciprocal space, with the type of spiral (proper screw vs.\ cycloid) and its orientation determined by weak magnetic anisotropy beyond the isotropic MFA.
