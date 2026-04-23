# Mean-Field Analysis: Finding the Ordering Wave Vector

Theoretical background for `tools/mfa_analysis.jl`.

---

## 1. Spin Hamiltonian

The Magesty SCE model describes pairwise exchange interactions via the full tensorial Hamiltonian

$$H = \sum_{i > j} \boldsymbol{S}_i^\top \mathbf{J}_{ij} \, \boldsymbol{S}_j$$

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
> `convert2tensor` uses the no-double-counting convention (sum over $i > j$).
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

$$k_B T_C = \frac{S^2}{3} \left| \lambda_\text{min}(J(\boldsymbol{q}_\text{opt})) \right|$$

- $S$: spin magnitude (classical normalization: $|\boldsymbol{S}| = S$)
- $k_B = 8.617333 \times 10^{-2}$ meV/K

> For quantum spins (quantum Heisenberg model), replace $S^2 \to S(S+1)$.

---

## 5. Wavelength

Let $\{\boldsymbol{a}_i^\text{prim}\}$ be the primitive lattice vectors (Å) and $\boldsymbol{b}_i$ the corresponding reciprocal vectors ($\boldsymbol{b}_i \cdot \boldsymbol{a}_j = \delta_{ij}$, unit Å$^{-1}$).

The Cartesian representation of $\boldsymbol{q}_\text{opt}$ (unit Å$^{-1}$) is

$$\boldsymbol{q}_\text{cart} = \sum_i q_i^\text{prim} \, \boldsymbol{b}_i$$

and the spiral wavelength (Å) is

$$\lambda = \frac{1}{|\boldsymbol{q}_\text{cart}|}$$

For $|\boldsymbol{q}| = 0$ (ferromagnetic), the wavelength is infinite.

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
    --spin 1.0
```

| Option | Description | Default |
|---|---|---|
| `--xml` / `-x` | Path to `jphi.xml` or `scecoeffs.xml` | (required) |
| `--nk` / `-n` | Number of $q$-points per reciprocal direction | 20 |
| `--spin` / `-s` | Spin magnitude $S$ (classical spin) | 1.0 |
