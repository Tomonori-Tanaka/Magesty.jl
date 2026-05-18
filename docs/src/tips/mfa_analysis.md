# Mean-Field Analysis: Finding the Ordering Wave Vector

Theoretical background for `tools/mfa_analysis.jl`. The derivation follows
the formulation of Mendive-Tapia (2021), adapted to Magesty's per-bond
storage convention for exchange parameters.

---

## 1. Classical Spin Hamiltonian

The classical spin model analyzed here is written in the ordered double-sum
form

```math
E = -\sum_{\boldsymbol{R},\boldsymbol{R}'}\sum_{\mu\nu} \boldsymbol{e}_{\boldsymbol{R}\mu}^{\top}\, \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{R}'-\boldsymbol{R})\, \boldsymbol{e}_{\boldsymbol{R}'\nu},
```

where $\boldsymbol{e}_{\boldsymbol{R}\mu}$ is the unit spin direction at
sublattice $\mu$ in Bravais cell $\boldsymbol{R}$, and
$\boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{\Delta})$ is the $3\times 3$
exchange-coupling tensor. The reality of $E$ enforces the symmetry condition

```math
\boldsymbol{\mathcal{J}}_{\mu\nu}(-\boldsymbol{\Delta}) = \boldsymbol{\mathcal{J}}_{\nu\mu}^{\top}(\boldsymbol{\Delta}), \qquad \boldsymbol{\Delta} = \boldsymbol{R}'-\boldsymbol{R}.
```

### Relation to the per-bond tensor stored in `jphi.xml`

Magesty's SCE fit produces the per-bond tensor
$\mathbf{J}_{ij}^{\text{per-bond}}$ (unit: meV), returned by
`convert2tensor`, defined for the **unordered**-pair Hamiltonian summed
over $i<j$:

```math
H = +\sum_{i<j} \boldsymbol{S}_i^{\top}\, \mathbf{J}_{ij}^{\text{per-bond}}\, \boldsymbol{S}_j,
```

with the convention $J_{ij}^\text{iso} < 0$ for ferromagnetic coupling.
Matching the same physical energy with the ordered double-sum form above
gives the conversion

```math
\boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{\Delta}) = -\tfrac{1}{2}\, \mathbf{J}_{\mu\nu}^{\text{per-bond}}(\boldsymbol{\Delta}).
```

Under this convention $\boldsymbol{\mathcal{J}} > 0$ corresponds to
ferromagnetic coupling. `precompute_pairs` in `mfa_analysis.jl` applies
the $-\tfrac{1}{2}$ factor internally, so all downstream quantities
(eigenvalues, $T_C$) are already in the paper convention.

### Decomposition of $\mathbf{J}_{ij}^{\text{per-bond}}$

The per-bond tensor splits into three physically distinct contributions:

$$\mathbf{J}_{ij}^{\text{per-bond}} = \underbrace{J_{ij}^\text{iso} \mathbf{I}}_{\text{isotropic}} + \underbrace{\mathbf{W}_{ij}}_{\text{antisymmetric (DMI)}} + \underbrace{\boldsymbol{\Gamma}_{ij}}_{\text{anisotropic symmetric}}$$

| Component | Definition | Physical origin |
|---|---|---|
| $J_{ij}^\text{iso} = \tfrac{1}{3}\operatorname{tr}(\mathbf{J}_{ij}^{\text{per-bond}})$ | scalar | Heisenberg exchange |
| $\mathbf{W}_{ij} = \tfrac{1}{2}(\mathbf{J}_{ij}^{\text{per-bond}} - \mathbf{J}_{ij}^{\text{per-bond}\,\top})$ | antisymmetric | Dzyaloshinskii–Moriya interaction (DMI), $\boldsymbol{D}_{ij} \cdot (\boldsymbol{S}_i \times \boldsymbol{S}_j)$ |
| $\boldsymbol{\Gamma}_{ij} = \tfrac{1}{2}(\mathbf{J}_{ij}^{\text{per-bond}} + \mathbf{J}_{ij}^{\text{per-bond}\,\top}) - J_{ij}^\text{iso}\mathbf{I}$ | symmetric traceless | Anisotropic symmetric exchange |

The $3\times 3$ block structure of $\boldsymbol{\mathcal{J}}(\boldsymbol{q})$
built below retains all three channels.

---

## 2. Mean-Field Decomposition

Introduce the local order parameter

```math
\boldsymbol{m}_{\boldsymbol{R}\mu} = \langle \boldsymbol{e}_{\boldsymbol{R}\mu} \rangle,
```

and split the spin direction into mean value plus fluctuation,
$\boldsymbol{e}_{\boldsymbol{R}\mu} = \boldsymbol{m}_{\boldsymbol{R}\mu} + \delta\boldsymbol{e}_{\boldsymbol{R}\mu}$.
Substituting into a single pair term and dropping the fluctuation-fluctuation
correlation gives the mean-field approximation

```math
\boldsymbol{e}_{\boldsymbol{R}\mu}^{\top} \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{\Delta})\, \boldsymbol{e}_{\boldsymbol{R}'\nu} \simeq \boldsymbol{e}_{\boldsymbol{R}\mu}^{\top} \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{\Delta})\, \boldsymbol{m}_{\boldsymbol{R}'\nu} + \boldsymbol{m}_{\boldsymbol{R}\mu}^{\top} \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{\Delta})\, \boldsymbol{e}_{\boldsymbol{R}'\nu} - \boldsymbol{m}_{\boldsymbol{R}\mu}^{\top} \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{\Delta})\, \boldsymbol{m}_{\boldsymbol{R}'\nu}.
```

Summing over all ordered pairs and using the reality condition, the two
linear terms give identical contributions. The mean-field energy becomes

```math
E^{\text{MF}} = -2\sum_{\boldsymbol{R}\mu} \boldsymbol{h}^{\text{MF}}_{\boldsymbol{R}\mu} \cdot \boldsymbol{e}_{\boldsymbol{R}\mu} + \sum_{\boldsymbol{R},\boldsymbol{R}'}\sum_{\mu\nu} \boldsymbol{m}_{\boldsymbol{R}\mu}^{\top}\, \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{R}'-\boldsymbol{R})\, \boldsymbol{m}_{\boldsymbol{R}'\nu},
```

with the **mean-field effective field**

```math
\boldsymbol{h}^{\text{MF}}_{\boldsymbol{R}\mu} = \sum_{\boldsymbol{R}'\nu} \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{R}'-\boldsymbol{R})\, \boldsymbol{m}_{\boldsymbol{R}'\nu}.
```

---

## 3. Self-Consistent Equation and Linearization Near $T_C$

For classical (unit-vector) spins on the unit sphere, the single-site
partition function in the effective field is

```math
Z_{\boldsymbol{R}\mu} = \int \mathrm{d}\Omega_{\boldsymbol{R}\mu}\, \exp\!\left(2\beta\, \boldsymbol{h}^{\text{MF}}_{\boldsymbol{R}\mu} \cdot \boldsymbol{e}_{\boldsymbol{R}\mu}\right) = \frac{2\pi}{\beta h^{\text{MF}}_{\boldsymbol{R}\mu}}\, \sinh\!\bigl(2\beta h^{\text{MF}}_{\boldsymbol{R}\mu}\bigr).
```

The order parameter then satisfies the self-consistent equation

```math
\boldsymbol{m}_{\boldsymbol{R}\mu} = L\!\bigl(2\beta h^{\text{MF}}_{\boldsymbol{R}\mu}\bigr)\, \hat{\boldsymbol{h}}^{\text{MF}}_{\boldsymbol{R}\mu}, \qquad L(x) = \coth x - \frac{1}{x},
```

where $L$ is the Langevin function. Near the magnetic transition the order
parameter is small, and the Langevin expansion $L(x) = x/3 - x^3/45 + \mathcal{O}(x^5)$
gives, to leading order,

```math
\boldsymbol{m}_{\boldsymbol{R}\mu} = \frac{2\beta}{3}\, \boldsymbol{h}^{\text{MF}}_{\boldsymbol{R}\mu} + \mathcal{O}(|\boldsymbol{m}|^3).
```

Substituting the definition of $\boldsymbol{h}^{\text{MF}}$ yields the
linearized real-space mean-field equation

```math
\boldsymbol{m}_{\boldsymbol{R}\mu} = \frac{2\beta}{3} \sum_{\boldsymbol{R}'\nu} \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{R}'-\boldsymbol{R})\, \boldsymbol{m}_{\boldsymbol{R}'\nu}.
```

---

## 4. Fourier Representation and the Matrix Equation

Define the lattice Fourier transforms

```math
\boldsymbol{m}_{\mu}(\boldsymbol{q}) = \frac{1}{\sqrt{N}} \sum_{\boldsymbol{R}} e^{-i\boldsymbol{q}\cdot\boldsymbol{R}}\, \boldsymbol{m}_{\boldsymbol{R}\mu}, \qquad \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{q}) = \sum_{\boldsymbol{\Delta}} e^{+i\boldsymbol{q}\cdot\boldsymbol{\Delta}}\, \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{\Delta}),
```

where $N$ is the number of Bravais cells. The linearized equation in
Fourier space reads

```math
\boldsymbol{m}_{\mu}(\boldsymbol{q}) = \frac{2\beta}{3} \sum_{\nu} \boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{q})\, \boldsymbol{m}_{\nu}(\boldsymbol{q}).
```

Collecting the spin and sublattice components into the
$3n_\text{sub}$-dimensional vector

```math
\boldsymbol{M}(\boldsymbol{q}) = \bigl(\boldsymbol{m}_1(\boldsymbol{q}),\, \boldsymbol{m}_2(\boldsymbol{q}),\, \ldots,\, \boldsymbol{m}_{n_\text{sub}}(\boldsymbol{q})\bigr)^{\top},
```

and assembling the $3n_\text{sub} \times 3n_\text{sub}$ block matrix
$\boldsymbol{\mathcal{J}}(\boldsymbol{q})$ whose $(\mu,\nu)$ block is the
$3\times 3$ tensor $\boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{q})$, the
linearized self-consistency equation becomes the eigenvalue problem

```math
\boldsymbol{M}(\boldsymbol{q}) = \frac{2\beta}{3}\, \boldsymbol{\mathcal{J}}(\boldsymbol{q})\, \boldsymbol{M}(\boldsymbol{q}).
```

The reality condition makes $\boldsymbol{\mathcal{J}}(\boldsymbol{q})$
Hermitian with real eigenvalues $\eta_n(\boldsymbol{\mathcal{J}}(\boldsymbol{q}))$.
Each eigenvector mixes the isotropic, DMI, and anisotropic-symmetric
channels of the underlying $3\times 3$ blocks.

---

## 5. Ordering Wavevector and Transition Temperature

A non-trivial ordered solution appears whenever
$\tfrac{2\beta}{3}\eta_n(\boldsymbol{\mathcal{J}}(\boldsymbol{q})) = 1$ for
some band $n$ and wavevector $\boldsymbol{q}$. The **first instability**
from the paramagnetic state therefore picks the mode that maximizes the
largest eigenvalue:

```math
\boldsymbol{q}_\text{max} = \underset{\boldsymbol{q}}{\operatorname{argmax}}\, \eta_\text{max}\!\left(\boldsymbol{\mathcal{J}}(\boldsymbol{q})\right).
```

| $\boldsymbol{q}_\text{max}$ | Ordering type |
|---|---|
| $(0,0,0)$ | Ferromagnetic (FM) |
| $(1/2, 1/2, 1/2)$, etc. | Antiferromagnetic (AFM) |
| Other rational or irrational values | Spin spiral (helimagnetic) |

The corresponding mean-field transition temperature is

```math
k_B T_C^{\text{MF}} = \frac{2}{3}\, \eta_\text{max}\!\left(\boldsymbol{\mathcal{J}}(\boldsymbol{q}_\text{max})\right) \quad (k_B = 8.617333 \times 10^{-2}\,\text{meV/K}).
```

For quantum Heisenberg spins of magnitude $S$ the prefactor becomes
$\tfrac{2}{3}\, S(S+1)$, which is what the `--spin` option uses.

The spin-spiral period (wavelength) in real space is

```math
\lambda = \frac{2\pi}{|\boldsymbol{k}_\text{cart}|}, \qquad \boldsymbol{k}_\text{cart} = 2\pi\, \boldsymbol{q}_\text{cart},
```

where $\boldsymbol{q}_\text{cart} = \sum_i q_i^\text{prim}\, \boldsymbol{b}_i$ is
$\boldsymbol{q}_\text{max}$ in Cartesian reciprocal space (unit Å$^{-1}$),
expressed via the primitive reciprocal basis $\{\boldsymbol{b}_i\}$ that
satisfies $\boldsymbol{b}_i \cdot \boldsymbol{a}_j^\text{prim} = \delta_{ij}$.
For $|\boldsymbol{q}|=0$ (ferromagnetic), $\lambda = \infty$.

> **Equivalent per-bond form.** Using
> $\boldsymbol{\mathcal{J}} = -\tfrac{1}{2}\mathbf{J}^{\text{per-bond}}$,
> the same $T_C$ can be written as
> $k_B T_C = \tfrac{1}{3}\bigl(-\lambda_\text{min}(\mathbf{J}^{\text{per-bond}}(\boldsymbol{q}))\bigr)$.
> Both forms are numerically identical.

---

## 6. Implementation in `mfa_analysis.jl`

### Materializing $\boldsymbol{\mathcal{J}}(\boldsymbol{q})$ from `jphi.xml`

The Fourier sum is built by looping over all neighbor atoms $j$ of
$\text{prim}_\mu$ that fall in sublattice $\nu$:

```math
\boldsymbol{\mathcal{J}}_{\mu\nu}(\boldsymbol{q}) = \sum_{\substack{j \in \text{sublattice } \nu \\ j \neq \text{prim}_\mu}} \left[-\tfrac{1}{2}\, \mathbf{J}^{\text{per-bond}}(\text{prim}_\mu,\, j)\right] e^{2\pi i\, \boldsymbol{q}_\text{sc} \cdot \boldsymbol{R}_{\mu j}}.
```

- $\text{prim}_\mu$: representative atom of sublattice $\mu$ in the primitive cell (`symmetry.atoms_in_prim[μ]`)
- $\nu = \texttt{symmetry.map\_s2p[j].atom}$: sublattice index of neighbor $j$
- $\boldsymbol{R}_{\mu j}$: minimum-image displacement in supercell fractional coordinates,

$$\boldsymbol{R}_{\mu j} = \boldsymbol{x}_j^\text{frac} - \boldsymbol{x}_{\text{prim}_\mu}^\text{frac} - \operatorname{round}\!\left(\boldsymbol{x}_j^\text{frac} - \boldsymbol{x}_{\text{prim}_\mu}^\text{frac}\right).$$

The Hermitian symmetrization
$\boldsymbol{\mathcal{J}}(\boldsymbol{q}) \leftarrow \tfrac{1}{2}\bigl(\boldsymbol{\mathcal{J}}(\boldsymbol{q}) + \boldsymbol{\mathcal{J}}(\boldsymbol{q})^\dagger\bigr)$
is applied before diagonalization to suppress floating-point asymmetry.

### Wavevector coordinates

The primitive BZ is sampled on an $n_k^3$ uniform grid with
$\boldsymbol{q}_\text{prim} \in [0,1)^3$ (fractional coordinates in the
primitive reciprocal basis). The supercell-fractional wavevector used in
the Fourier phase is

$$\boldsymbol{q}_\text{sc} = A_\text{prim}\, \boldsymbol{q}_\text{prim},$$

where $A_\text{prim}$ is the $3\times 3$ matrix whose columns are the
primitive lattice vectors in supercell fractional coordinates (constructed
from the pure-translation symmetry operations in `prim_lattice_in_sc_frac`).
The script prints both

- `Ordering vector |q| (Cartesian)`: $|\boldsymbol{q}_\text{cart}|$ in units of $2\pi/\text{Å}$ (numerically $|\boldsymbol{k}_\text{cart}|/(2\pi)$), and
- `Angular wavevector |k|`: $|\boldsymbol{k}_\text{cart}|$ in rad/Å.

### Grid scan and gradient-ascent refinement

After the $n_k^3$ grid scan locates the coarse maximum, the script refines
$\boldsymbol{q}_\text{max}$ by steepest ascent on
$\eta_\text{max}(\boldsymbol{\mathcal{J}}(\boldsymbol{q}))$ with Armijo
backtracking. The gradient uses the Hellmann–Feynman theorem

```math
\frac{\partial \eta_\text{max}}{\partial q_k} = \mathrm{Re}\!\left[\boldsymbol{v}^\dagger \frac{\partial \boldsymbol{\mathcal{J}}(\boldsymbol{q})}{\partial q_k} \boldsymbol{v}\right], \qquad \frac{\partial \boldsymbol{\mathcal{J}}(\boldsymbol{q})}{\partial q_k^{\text{sc}}} = 2\pi i\, R_k\, \boldsymbol{\mathcal{J}}(\boldsymbol{q}),
```

with $\boldsymbol{v}$ the dominant eigenvector and the chain rule
$\partial \eta / \partial \boldsymbol{q}_\text{prim} = A_\text{prim}^{\top}\, \partial \eta / \partial \boldsymbol{q}_\text{sc}$.
Refinement can be disabled with `--no-refine`.

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
| `--eigvec` / `-e` | Print the $N$-th eigenvector of $\boldsymbol{\mathcal{J}}(\boldsymbol{q}_\text{max})$ ($N=1$: largest eigenvalue $\eta_\text{max}$, $N=2$: second largest, …) | omitted |
| `--no-refine` | Skip gradient-ascent refinement after the grid scan | false |
