# Input Keys Reference

Magesty.jl uses TOML configuration files. Below is an annotated example for a BCC Fe supercell,
followed by a full reference for every supported key.

## Annotated Example

```toml:input.toml
[general]
name = "bccfe"
kd   = ["Fe"]           # list of element names
nat  = 16               # total number of atoms in the supercell
periodicity = [true, true, true]  # apply periodic boundary to all (x, y, z) directions

[symmetry]
tolerance = 1e-5        # symmetry detection tolerance (optional, default 1e-3)
isotropy = true

[interaction]
nbody = 2               # maximum interaction body
[interaction.body1]
lmax.Fe = 0             # 1-body maximum angular momentum per element, i.e. this represents on-site anisotropy
[interaction.body2]
lsum = 2                # cutoff summation of l values for basis functions
cutoff."Fe-Fe" = -1     # pairwise cutoff radius in Å (-1 uses all possible pairs)

[regression]
datafile = "EMBSET" # path to training data
weight   = 0.5          # 0 = torque only, 1 = energy only, 0.5 = balanced
alpha    = 0.0          # elastic-net mixing (0 = ridge)
lambda   = 0.0          # regularization strength (0 = no regularization)

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
| `isotropy` | Bool | no | Only include isotropic (L=0) terms (default: `false`) |

### `[interaction]`

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `nbody` | Int | yes | Maximum interaction order |
| `body1.lmax.<elem>` | Int | no | On-site max angular momentum per element |
| `body<n>.lsum` | Int | yes (n≥2) | Max L sum for n-body basis |
| `body<n>.cutoff."<e1>-<e2>"` | Float64 | yes (n≥2) | Pairwise cutoff radius in Å |

### `[regression]`

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `datafile` | String | yes | Path to EMBSET training data |
| `ndata` | Int | no | Number of data points to use (default: `-1` = all) |
| `weight` | Float64 | no | Energy/torque balance: 0=torque only, 1=energy only (default: `0.0`) |
| `alpha` | Float64 | no | Elastic-net mixing parameter (default: `0.0`) |
| `lambda` | Float64 | no | Regularization strength (default: `0.0`) |

### `[structure]`

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `kd_list` | Vector{Int} | yes | Element index (1-based, into `kd`) per atom |
| `lattice` | 3×3 Float64 | yes | Lattice vectors in Å; each of the three rows in the TOML array defines one lattice vector (a₁, a₂, a₃) |
| `position` | Vector of 3-vectors | yes | Fractional atomic coordinates |
