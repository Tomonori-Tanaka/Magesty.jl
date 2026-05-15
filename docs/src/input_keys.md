# Input Keys Reference

Magesty.jl reads structure and interaction settings from a TOML file
(consumed by `SCEBasis(toml_path)` / `SCEBasis(input_dict)`). Fit
parameters (estimator, regularization, torque weight) are **not** part
of the TOML — they are passed in Julia at `fit` time, e.g.
`fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)`.

## Annotated Example

```toml:input.toml
[general]
name = "bccfe"
kd   = ["Fe"]           # list of element names
nat  = 16               # total number of atoms in the supercell
periodicity = [true, true, true]  # apply periodic boundary to all (x, y, z) directions

[symmetry]
tolerance = 1e-5        # symmetry detection tolerance (optional, default 1e-3)
isotropy = true         # restrict to Lf = 0 (isotropic exchange) terms

[interaction]
nbody = 2               # maximum interaction body
[interaction.body1]
lmax.Fe = 0             # 1-body maximum angular momentum per element, i.e. this represents on-site anisotropy
[interaction.body2]
lsum = 2                # cutoff summation of l values for basis functions
cutoff."Fe-Fe" = -1     # pairwise cutoff radius in Å (-1 uses all possible pairs)

[structure]
kd_list  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # element index per atom
lattice  = [            # lattice parameters
  [5.66, 0.0, 0.0],
  [0.0, 5.66, 0.0],
  [0.0, 0.0, 5.66],
]
position = [            # fractional coordinates
  [0.00, 0.00, 0.00],
  [0.25, 0.25, 0.25],
  [0.00, 0.00, 0.50],
  [0.50, 0.00, 0.00],
  [0.00, 0.50, 0.00],
  [0.25, 0.25, 0.75],
  [0.75, 0.25, 0.25],
  [0.25, 0.75, 0.25],
  [0.00, 0.50, 0.50],
  [0.50, 0.00, 0.50],
  [0.50, 0.50, 0.00],
  [0.25, 0.75, 0.75],
  [0.75, 0.25, 0.75],
  [0.75, 0.75, 0.25],
  [0.50, 0.50, 0.50],
  [0.75, 0.75, 0.75],
]
```

## Key Reference

### `[general]`

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `name` | String | yes | System name |
| `kd` | Vector{String} | yes | Element names |
| `nat` | Int | yes | Number of atoms in supercell |
| `periodicity` | Vector{Bool} | no | Periodic boundary conditions (default: `[true,true,true]`) |

### `[symmetry]`

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `tolerance` | Float64 | no | Symmetry detection tolerance (default: `1e-3`) |
| `isotropy` | Bool | no | Only include isotropic (`Lf = 0`) terms (default: `false`) |

### `[interaction]`

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `nbody` | Int | yes | Maximum interaction order |
| `body1.lmax.<elem>` | Int | no | On-site max angular momentum per element |
| `body<n>.lsum` | Int | yes (n≥2) | Max L sum for n-body basis |
| `body<n>.cutoff."<e1>-<e2>"` | Float64 | yes (n≥2) | Pairwise cutoff radius in Å (`-1` = include all pairs) |

### `[structure]`

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `kd_list` | Vector{Int} | yes | Element index (1-based, into `kd`) per atom |
| `lattice` | 3×3 Float64 | yes | Lattice vectors in Å; each of the three rows in the TOML array defines one lattice vector (a₁, a₂, a₃) |
| `position` | Vector of 3-vectors | yes | Fractional atomic coordinates |

## Fit parameters

The TOML covers the *material* — structure, symmetry, interaction. The
*fit* — estimator, regularization, torque weight, training data path —
is configured in Julia at `fit` time. The training data path is also
passed in Julia (to `SCEDataset`), not in the TOML.

```julia
basis   = SCEBasis("input.toml")
dataset = SCEDataset(basis, "EMBSET.dat")
f = fit(
    SCEFit, dataset,
    Ridge(lambda = 1e-4);   # estimator: OLS() or Ridge(; lambda)
    torque_weight = 0.5,    # ∈ [0, 1]: 0 = energy only, 1 = torque only,
                            # 0.5 = balanced (per-sample MSE convex combination)
)
```

- **`torque_weight`** ∈ `[0, 1]` — the convex weight applied to the
  per-sample torque MSE; the energy MSE gets `1 - torque_weight`.
  `0` = fit energies only; `1` = fit torques only. Default: `0.5`.
- **Estimator** — `OLS()` for no regularization, or `Ridge(lambda = λ)`
  for L2 with strength `λ ≥ 0`. The bias column (`j0`) is excluded
  from the penalty.
