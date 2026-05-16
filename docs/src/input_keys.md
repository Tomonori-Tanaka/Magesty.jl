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
| `body1.lmax.<elem>` | Int | no | On-site max angular momentum per element (wildcard `*` accepted, see below) |
| `body<n>.lsum` | Int | yes (n≥2) | Max L sum for n-body basis |
| `body<n>.cutoff."<e1>-<e2>"` | Float64 | yes (n≥2) | Pairwise cutoff radius in Å (`-1` = include all pairs; wildcards accepted, see below) |

#### Wildcard species notation

For `body1.lmax.<species>` and `body<n>.cutoff."<e1>-<e2>"` the species
name may be replaced with `*`, matching every species. When several keys
cover the same species or pair, the **most specific** key wins (it is
not the textual order of the keys). For pair keys the specificity order
is:

| Tier | Pattern | Example |
|------|---------|---------|
| 0 (most specific) | both concrete | `"Co-Ni"` |
| 1 | one wildcard | `"Fe-*"`, `"*-Fe"` |
| 2 (least specific) | both wildcards | `"*-*"` |

Pair keys are **unordered**: `"Co-Ni"` and `"Ni-Co"` are equivalent (the
parser rejects them as duplicates if both are supplied). The species
notation for `body1.lmax` has only two tiers: concrete species
(`"Fe"`) wins over the wildcard (`"*"`).

If two keys at the same tier cover the same pair with different values
(for instance `"Fe-*"` and `"*-Ni"` both covering `Fe-Ni`), the parser
reports an ambiguity and requires the user to add a tier-0 entry for
that pair. Unknown species names in a non-wildcard key are an error.

```toml
[interaction]
nbody = 2

[interaction.body1.lmax]
"*"  = 2          # default for every species
"Fe" = 4          # Fe overrides the default

[interaction.body2]
lsum = 4

[interaction.body2.cutoff]
"*-*"  = 8.0      # default for every pair
"Fe-*" = 10.0     # every pair involving Fe (incl. Fe-Fe)
"Co-Ni" = 12.0    # tier-0 entry — wins over any wildcard match
```

### Interaction cluster definition

Given the cutoff table for body `n`, an `n`-body atom set is accepted
as an interaction cluster if **every one of its `C(n, 2)` pairwise
cartesian distances** is less than or equal to the corresponding
pair-specific cutoff `body<n>.cutoff."<kᵢ>-<kⱼ>"`. The same rule applies
uniformly to every `n ≥ 2`: pair, triplet, quartet, and so on. No body
order treats the criterion differently.

Concretely, for body `n` and an `n`-body candidate `(a₁, …, aₙ)` of
species `(k₁, …, kₙ)`:

> the candidate is accepted iff `‖xᵢ − xⱼ‖ ≤ body<n>.cutoff."<kᵢ>-<kⱼ>"` for every pair `(i, j)` with `1 ≤ i < j ≤ n`.

If even one pair distance exceeds its cutoff the candidate is rejected
in full, regardless of how many other pairs fit. Setting a pair cutoff
to `-1` disables the constraint for that pair only — the rule is still
evaluated for every other pair in the cluster.

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
