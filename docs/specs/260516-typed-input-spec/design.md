# Design: Typed input spec + wildcard species + cluster definition

Status: draft (2026-05-16)
See requirements.md for goal / invariants / completion criteria.

## 1. Module layout

New file: `src/InputSpecs.jl` (module `InputSpecs`). Contains:

- Three typed value structs: `SystemSpec`, `InteractionSpec`, `SymmetryOptions`.
- Public entry point `parse_toml_inputs(dict)` that returns the
  three-tuple `(SystemSpec, InteractionSpec, SymmetryOptions)`.
- Internal helpers: `expand_pair_table`, `expand_species_table` (wildcard
  expansion), `_validate_*` (per-spec validation routines moved from
  `ConfigParser`).

`ConfigParser.jl` is deleted once all references are gone.
`AtomsBaseAdapter.jl` keeps its module name and public function names
(`atomsbase_to_specs`, `kwargs_to_specs` — renamed from
`system_to_input_dict` / `kwargs_to_input_dict`).

## 2. Typed value definitions

### 2.1 `SystemSpec`

```julia
struct SystemSpec
    name::String
    num_atoms::Int
    kd_name::Vector{String}              # unique species labels, in order
    kd_int_list::Vector{Int}             # per-atom index into kd_name (1-based)
    lattice_vectors::Matrix{Float64}     # 3 × 3, column = a / b / c (Å)
    x_fractional::Matrix{Float64}        # 3 × num_atoms
    is_periodic::Vector{Bool}            # length 3, per-direction
end
```

Inner constructor validation (matches the existing `Config4System`
behaviour to preserve numerical regression):

- `!isempty(name)`.
- `num_atoms > 0`.
- `!isempty(kd_name)`; `allunique(kd_name)`.
- `size(lattice_vectors) == (3, 3)`.
- `size(x_fractional) == (3, num_atoms)`.
- `length(kd_int_list) == num_atoms`.
- `1 ≤ kd_int_list[i] ≤ length(kd_name)` for all i.
- `length(is_periodic) == 3`.

The current code does not range-check `x_fractional`; we do not add new
checks beyond what `Config4System` enforced. No cross-spec dependencies.

### 2.2 `InteractionSpec`

```julia
struct InteractionSpec
    nbody::Int
    body1_lmax::Vector{Int}                                   # length = nkd
    bodyn_lsum::OffsetArray{Int, 1}                           # axis: body 2..nbody
    bodyn_cutoff::OffsetArray{Float64, 3}                     # axes: (body 2..nbody, species i, species j); -1.0 = include all
end
```

(Axis layout matches the current `Config4System` exactly, so downstream
code that indexed `bodyn_cutoff[body, ki, kj]` and `bodyn_lsum[body]`
keeps working with a field-name change only.)

Inner constructor signature:

```julia
function InteractionSpec(
    ; nbody::Int,
      body1_lmax::AbstractVector,
      bodyn_lsum::AbstractArray,
      bodyn_cutoff::AbstractArray,
      kd_name::Vector{String},          # for cross-validation against species
)
```

Validation (matches current `Config4System` behaviour):

- `nbody ≥ 1`.
- `length(body1_lmax) == length(kd_name)`.
- For `nbody == 1`: `bodyn_lsum` and `bodyn_cutoff` are empty
  OffsetArrays over range `2:1` (mirrors current code).
- For `nbody ≥ 2`: `bodyn_lsum` has axis `2:nbody` and `bodyn_cutoff`
  has axes `(2:nbody, 1:nkd, 1:nkd)`; `bodyn_cutoff[body, i, j] ==
  bodyn_cutoff[body, j, i]` for all body, i, j; every pair entry
  filled (no `0.0`-as-uninitialised sentinel — the parser fills with
  user values, with `-1.0` meaning "all pairs").
- All pair entries are concrete numbers; no wildcards survive (the
  parser handles expansion).

`kd_name` is consumed for validation but not stored.

### 2.3 `SymmetryOptions`

```julia
struct SymmetryOptions
    tolerance_sym::Float64
    isotropy::Bool
end
```

Defaults (when not supplied, matching `DEFAULT_VALUES_SYSTEM` in current
`ConfigParser`): `tolerance_sym = 1e-3`, `isotropy = false`.

Validation: `tolerance_sym > 0`.

## 3. `parse_toml_inputs(dict)` — dict / TOML parser

Single TOML-shape dict entry point. Pseudocode:

```julia
function parse_toml_inputs(dict::AbstractDict)
    _require_sections(dict, ["general", "symmetry", "interaction", "structure"])

    general = dict["general"]
    sym     = dict["symmetry"]
    inter   = dict["interaction"]
    str_    = dict["structure"]

    # 1. SystemSpec.
    kd_name = general["kd"] :: Vector{String}
    lattice_vectors = hcat(str_["lattice"]...)
    positions = [Float64.(vec) for vec in str_["position"]]
    x_fractional = parse_position(positions, general["nat"])  # reuse current helper
    system = SystemSpec(
        name = general["name"],
        num_atoms = general["nat"],
        kd_name = kd_name,
        kd_int_list = str_["kd_list"],
        lattice_vectors = lattice_vectors,
        x_fractional = x_fractional,
        is_periodic = get(general, "periodicity", [true, true, true]),
    )

    # 2. InteractionSpec with wildcard expansion.
    nbody = inter["nbody"]::Int
    body1_lmax = expand_species_table(
        kd_name,
        haskey(inter, "body1") ? inter["body1"]["lmax"] : Dict{String, Int}(),
    )                                                       # nbody == 1 with empty body1 → all zeros (current behavior)
    bodyn_lsum, bodyn_cutoff = expand_bodyn_tables(kd_name, nbody, inter)

    interaction = InteractionSpec(
        nbody = nbody,
        body1_lmax = body1_lmax,
        bodyn_lsum = bodyn_lsum,
        bodyn_cutoff = bodyn_cutoff,
        kd_name = kd_name,
    )

    # 3. SymmetryOptions.
    options = SymmetryOptions(
        tolerance_sym = get(sym, "tolerance", 1e-3),
        isotropy = get(sym, "isotropy", false),
    )

    return (system, interaction, options)
end
```

Key paths match the current `ConfigParser` schema verbatim
(`general.{name,nat,kd,periodicity}`, `symmetry.{tolerance,isotropy}`,
`interaction.nbody` / `interaction.body1.lmax.<species>` /
`interaction.body<n>.{lsum,cutoff."<a>-<b>"}`,
`structure.{lattice,kd_list,position}`). No schema breaking changes;
we only newly accept `*` in pair / species keys.

## 4. Wildcard expansion

### 4.1 Specificity rule (pair table)

For a pair key string, normalize species to alphabetic order
(`"Ni-Fe"` → `"Fe-Ni"`), then classify:

| Tier | Pattern | Example |
|---|---|---|
| 0 (highest) | both concrete | `"Fe-Ni"` |
| 1 | one wildcard | `"Fe-*"` (same as `"*-Fe"` after normalization to "left wildcard if present") |
| 2 (lowest) | both wildcards | `"*-*"` |

For each concrete pair `(i, j)` with `i ≤ j` in `kd_name`:
1. Collect every key that matches it: a tier-0 key with the same pair,
   tier-1 keys where the concrete side matches `i` or `j`, and the
   tier-2 key.
2. Pick the highest-tier match (smallest tier number).
3. If multiple keys exist at that same tier with different values
   covering this pair → raise `WildcardConflictError` with the offending
   keys named.
4. Same key listed twice (e.g. both `"Fe-Ni"` and `"Ni-Fe"`) → raise
   `DuplicatePairError`.

Missing all tiers → raise `MissingPairError(pair)`.

Unknown species in a non-wildcard key (e.g. `"Fe-Xx"` with `Xx ∉ kd_name`)
→ raise `UnknownSpeciesError`.

### 4.2 Specificity rule (single-species table, for `body1_lmax`)

Same idea with two tiers: concrete species (tier 0), `"*"` (tier 1).
Conflicts and unknowns handled identically.

### 4.3 Algorithm shape

```julia
function expand_pair_table(kd_name::Vector{String}, entries::Dict{String, T}) where {T}
    pairs_by_tier = classify(entries, kd_name)   # validates unknown species
    out = OffsetArray(...)
    for i in 1:length(kd_name), j in i:length(kd_name)
        winners = pick_highest_tier(pairs_by_tier, i, j)
        length(winners) == 0 && error(MissingPairError(...))
        length(winners) > 1 && check_or_error(winners)   # tier collision
        out[i, j] = out[j, i] = entries[winners[1]]
    end
    return out
end
```

Helpers live as internal functions in `InputSpecs`. Per-body wrapping
(`expand_pair_table_per_body`) loops over `body in 2:nbody` and calls
`expand_pair_table` per body.

## 5. AtomsBase / kwargs paths

Adapter functions move from "dict builder" to "typed value builder":

```julia
# Before:  system_to_input_dict(system; kwargs...) :: Dict
# After:
function atomsbase_to_specs(
    system::AbstractSystem;
    nbody::Int,
    body1_lmax,           # Dict{String,Int} or Int (broadcast)
    bodyn_lsum,
    bodyn_cutoff,
    symmetry_options = SymmetryOptions(),
)::Tuple{SystemSpec, InteractionSpec, SymmetryOptions}
```

The implementation:
1. Pull lattice / positions / species from `system` via existing
   AtomsBase boundary code (preserve unit handling: `ustrip` to Å).
2. Build `SystemSpec` directly.
3. Build `InteractionSpec` directly. If `bodyn_cutoff` is provided as a
   `Dict{String, ...}` with wildcard keys, route it through
   `expand_pair_table_per_body` (shared helper). If provided as a fully
   filled `OffsetArray`, pass through.
4. Return the tuple.

`kwargs_to_specs` is structurally the same but takes raw `lattice`,
`positions`, `kd`, etc. kwargs.

The `SCEBasis(toml_path)` / `SCEBasis(dict)` / `SCEBasis(system; ...)` /
`SCEBasis(; lattice, kd, ...)` public surface in `Magesty.jl` keeps its
signature; only the inner glue switches from `Config4System` to the
three-tuple.

## 6. Downstream signature changes

| Function | Old | New |
|---|---|---|
| `Structure(config)` | `Structure(config::Config4System)` | `Structure(system::SystemSpec)` |
| `Symmetry(structure, config)` | `Symmetry(structure, config::Config4System)` | `Symmetry(structure::Structure, options::SymmetryOptions)` |
| `Cluster(structure, symmetry, config)` | `Cluster(structure, symmetry, config::Config4System; verbosity)` | `Cluster(structure::Structure, symmetry::Symmetry, interaction::InteractionSpec; verbosity)` |
| `SALCBasis(structure, symmetry, cluster, config)` | uses `config.body1_lmax`, `config.bodyn_lsum`, `config.nbody`, `config.isotropy` | `SALCBasis(structure, symmetry, cluster, interaction::InteractionSpec, options::SymmetryOptions)` |

`Fitting.jl` and `XMLIO.jl` keep their public signatures (they consume
`SCEBasis`, not `Config4System`).

`SCEBasis` internal constructor gathers the three specs and calls
the chain above.

## 7. Migration order (matches tasklist.md milestones)

1. **M1**: Add `InputSpecs.jl` with the three structs and validation.
   Keep `Config4System` untouched. New module is unreferenced.
2. **M2**: Add `parse_toml_inputs` + wildcard helpers. Rewire
   `SCEBasis(toml_path)` and `SCEBasis(dict)` to call
   `parse_toml_inputs` → then build `Config4System` from the tuple
   (temporary bridge). Tests for wildcard expansion go in here.
3. **M3**: Rewrite `atomsbase_to_specs` / `kwargs_to_specs` to typed
   construction. Delete dict builders. AtomsBase / kwargs paths now
   produce typed values; still bridged to `Config4System` for downstream.
4. **M4**: Migrate downstream one at a time:
   - `Structure(system::SystemSpec)` first.
   - `Symmetry(structure, options::SymmetryOptions)`.
   - `Cluster(structure, symmetry, interaction::InteractionSpec)`.
   - `SALCBasis(...)` last (largest surface).
5. **M5**: Remove the `Config4System` bridge from `SCEBasis`. Delete
   `Config4System`, `ConfigParser.jl`, and the old `test_ConfigParser.jl`.
6. **M6**: Docs (`docs/src/input_keys.md`, `SPEC.md`, `DESIGN_NOTES.md`).

Each milestone leaves the test suite green.

## 8. Cluster definition (documentation only — no code change)

Existing `Clusters.is_within_cutoff` (`src/Clusters.jl`) iterates over
every 2-combination of atoms in an n-body cluster candidate and accepts
the cluster only if every pair distance is within its pair-specific
cutoff. The rule is uniform across body orders: for every n ≥ 2, all
`C(n, 2)` pairwise distances must satisfy the criterion.

In M6, `docs/src/input_keys.md` gains a new section "Interaction cluster
definition" stating this rule with one general formula (the
all-pairs-within-cutoff criterion) and, optionally, a brief example for
a specific n; the example must not imply that other body orders use a
different rule.

`SPEC.md` cross-references this section from the `Cluster` type
description.

## 9. Linked invariants checklist

- TOML schema unchanged except: pair / species keys may use `*`.
- `(l, m, site)` ordering in `SALCBasis`: untouched.
- `Clusters.is_within_cutoff` source: untouched (documentation only).
- XML I/O fingerprint: identical for existing inputs.
- `cutoff[body, ki, kj]` axis layout: unchanged.

## 10. Future work (out of this spec)

- Wildcards for `bodyn_lsum` (when use cases arise).
- Multi-error aggregation in `parse_toml_inputs`.
- Migrating `SystemSpec.lattice_vectors` to `SMatrix{3, 3, Float64}`
  (defer until hot-path effect measured).
- Additional input formats (CIF, JSON) — now straightforward because
  typed value construction is decoupled from dict parsing.
