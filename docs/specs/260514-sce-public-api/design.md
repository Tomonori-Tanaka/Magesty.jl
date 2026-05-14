# Design: SCE public API refactor

Status: draft (2026-05-14)
Requirements: `./requirements.md`

## Module / file layout

New / changed files under `src/`:

- `src/SCEBasis.jl` — new module `SCEBases`(?) or keep types in
  `Magesty.jl`. **Decision pending** (see "Open questions"): the four
  public types are small wrappers; they may live directly in
  `Magesty.jl` rather than each getting a module. Current leaning:
  define `SCEBasis`, `SCEDataset`, `SCEFit` in `Magesty.jl` (top-level,
  like `System`/`SpinCluster` today), keep `SCEModel` where it is.
- `src/utils/atomsbase_adapter.jl` — new. Isolates all `AtomsBase` /
  `Unitful` calls: `AbstractSystem` → `Structure` conversion, species
  sublabel extraction, unit stripping.
- `src/Optimize.jl` — `fit(::Type{SCEFit}, ...)` and the StatsAPI verb
  methods; `assemble_weighted_problem` already does per-sample MSE
  normalization (PR #3), only the kwarg name `weight` → `torque_weight`
  threads through.
- `src/utils/xml_io.jl` — `save`/`load` backends; `SCEBasis` XML
  (basis-only) writer/reader; `SCEModel` XML stays byte-identical.
- `src/Magesty.jl` — new exports, old API removed in the final step.

New dependencies: `AtomsBase`, `StatsAPI`, `Unitful` (the latter
transitively via AtomsBase but used directly in the adapter). `JLD2`
is intentionally *not* a dependency — see "save / load".

## Type definitions

```julia
struct SCEBasis
    structure::Structure          # Magesty internal representation
    symmetry::Symmetry            # derived from structure + symmetry.tol; cached
    salcbasis::SALCBasis          # heavy SALC result
end
```

`SCEBasis` stores only what the finished basis *is*, not the
construction scaffolding. `Cluster` is **not** a field: it is a
construction step (`structure` + `symmetry` + interaction params →
`Cluster` → `SALCBasis`), computed inside the constructor and
discarded once `SALCBasis` is built. `Cluster` stays an internal
struct — only its role changes from stored state to construction
intermediate. `Symmetry` *is* stored (it is used on every
`SCEDataset` build and every `predict_*` call) but is conceptually a
cached derivation of `structure` + `symmetry.tol`.

The AtomsBase `AbstractSystem` is *not* stored: it is converted to
`Structure` at construction and discarded (keeps Unitful out of the
stored type, per requirements "Unitful boundary"). Species sublabels
(`Fe_4a`, `Fe_8e`) are carried into `Structure.kd_name`.

```julia
struct SCEDataset
    basis::SCEBasis
    spinconfigs::Vector{SpinConfig}
    X_E::Matrix{Float64}   # energy design matrix, bias column at col 1, UNWEIGHTED
    X_T::Matrix{Float64}   # torque design matrix, no bias column, UNWEIGHTED
    y_E::Vector{Float64}   # observed energies, length n_configs
    y_T::Vector{Float64}   # observed torques flattened, length 3*n_atoms*n_configs
end
```

`X_E` / `X_T` are stored **unweighted**. `torque_weight` is applied at
`fit` time (via `assemble_weighted_problem`), so a weight sweep reuses
the same `SCEDataset` without rebuilding design matrices.

```julia
struct SCEFit <: StatsAPI.RegressionModel
    dataset::SCEDataset
    j0::Float64
    jphi::Vector{Float64}
    estimator::AbstractEstimator
    torque_weight::Float64
    residuals::Vector{Float64}    # over the augmented (weighted) system
    metrics::Dict{Symbol, Any}
end
```

`SCEFit <: StatsAPI.RegressionModel` (concrete, no intermediate
abstract type — see requirements "Out of scope", Bayesian).

```julia
struct SCEModel
    basis::SCEBasis
    j0::Float64                   # was: reference_energy
    jphi::Vector{Float64}         # was: SCE
end
```

`SCEModel` is reshaped to hold a `SCEBasis` plus the fitted bias and
coefficients. Holding the whole `basis` (rather than flat `salcbasis`
/ `symmetry` / `num_atoms` fields) gives `SCEModel` the `structure` it
needs to write a complete XML, and lets the XML `save` path share the
`SCEBasis` writer. `num_atoms` becomes
`model.basis.structure.supercell.num_atoms`; `SCEModel(f::SCEFit) =
SCEModel(f.dataset.basis, f.j0, f.jphi)`. The historical field renames
(`reference_energy` → `j0`, `SCE` → `jphi`) still apply to the names
kept.

## Interaction parameters

The interaction specification is hierarchical — 1-body and N-body
terms take different parameters — so it is passed as a nested
`NamedTuple` keyed `body1`, `body2`, ..., structurally identical to
the TOML `[interaction]` table.

| Parameter | Body level | Meaning | Type |
|---|---|---|---|
| `lmax` | `body1` only | max angular momentum `l` per atomic species (on-site anisotropy order) | `Dict{Symbol, Int}` (species-keyed) |
| `lsum` | `body2`+ | upper bound on the **sum** of the angular momenta of the N atoms in the term | `Int` (one scalar per body) |
| `cutoff` | `body2`+ | distance cutoff per species pair (`-1.0` = no cutoff) | `Dict{Tuple{Symbol,Symbol}, V}`, `V <: Union{Real, Unitful.Length}` |

```julia
interaction = (
    body1 = (lmax = Dict(:Fe => 2, :Pt => 2),),
    body2 = (lsum = 2,
             cutoff = Dict((:Fe, :Fe) => -1.0,
                           (:Fe, :Pt) => -1.0,
                           (:Pt, :Pt) => -1.0)),
    # body3 = (lsum = ..., cutoff = ...),   # add for nbody >= 3
)
```

`nbody` is inferred from the `body1`, `body2`, ... keys present.
`lmax` (1-body) caps a single atom's `l`; `lsum` (N-body) caps the
sum over the N atoms — they live at different body levels and are not
interchangeable. Species keys may use sublabels (`:Fe_4a`); an
element-only key fans out to all sublabels of that element.

`cutoff` values may be written as bare reals (taken as angstrom) or as
`Unitful.Length` quantities (auto-converted to angstrom); the `-1.0`
no-cutoff sentinel stays a bare real. See "Unitful boundary".

**Time-reversal symmetry**: SALC construction enforces time-reversal
symmetry, which removes all odd-`l` terms. An odd `lmax` or `lsum`
therefore produces exactly the same basis as `value - 1`. The
`SCEBasis` constructor emits a `@warn` when it sees an odd `lmax` or
`lsum`, e.g.
`lmax = 3 for species :Fe is odd; time-reversal symmetry drops odd-l terms, so this is equivalent to lmax = 2.`
It is a warning, not an error — the run still proceeds. The check is
in the `SCEBasis` constructor, so both the AtomsBase and the TOML
input paths get it.

## Constructors

```julia
# SCEBasis — material + basis. Four input paths.

# (1) AtomsBase system — Unitful required only if the user hand-writes
#     the system; the common case is the return value of AtomsIO
#     (CIF/POSCAR/extxyz), where the user never touches Unitful.
SCEBasis(system::AtomsBase.AbstractSystem;
         interaction::NamedTuple,           # see "Interaction parameters"
         tolerance_sym::Real = 1e-5,
         isotropy::Bool = false,
         verbosity::Bool = true) -> SCEBasis

# (2) TOML template — no Unitful
SCEBasis(toml_path::AbstractString; verbosity = true) -> SCEBasis

# (3) parsed TOML dict — no Unitful
SCEBasis(input_dict::AbstractDict; verbosity = true) -> SCEBasis

# (4) Julia-native keyword constructor — no Unitful. Numbers are taken
#     as angstrom (lattice) and fractional (positions); Magesty applies
#     the units internally. This is the dict of path (3) spread into
#     keyword arguments — the recommended way to build a basis in pure
#     Julia without importing Unitful.
SCEBasis(; lattice::AbstractMatrix,         # 3x3; bare reals = angstrom, Unitful lengths also OK
           kd::AbstractVector{Symbol},      # unique species/sublabel names
           kd_list::AbstractVector{<:Integer},  # per-atom index into `kd`
           positions::AbstractVector,       # fractional coordinates, one per atom
           periodicity::NTuple{3,Bool} = (true, true, true),
           interaction::NamedTuple,
           tolerance_sym::Real = 1e-5,
           isotropy::Bool = false,
           verbosity::Bool = true) -> SCEBasis

# SCEDataset — basis + data + design matrices
SCEDataset(basis::SCEBasis, spinconfigs::AbstractVector{SpinConfig}) -> SCEDataset
SCEDataset(basis::SCEBasis, embset_path::AbstractString) -> SCEDataset
# sugar: builds SCEBasis internally
SCEDataset(system::AtomsBase.AbstractSystem,
           spinconfigs::AbstractVector{SpinConfig};
           interaction, tolerance_sym, isotropy, ...) -> SCEDataset
SCEDataset(toml_path::AbstractString,
           spinconfigs::AbstractVector{SpinConfig}) -> SCEDataset
```

The sugar constructors build a throwaway `SCEBasis` internally and
embed it in `dataset.basis`. For multi-dataset workflows that share a
basis (FM + AFM), the user constructs `SCEBasis` explicitly once and
passes it to several `SCEDataset(basis, ...)` calls.

## Verbs

Energy and torque are two separate response blocks. **Any quantity
that differs between them is always suffixed `_energy` / `_torque`** —
there is no bare `r2` / `rss` / `residuals` / `predict` that silently
means "energy". Only response-block-independent verbs are taken from
StatsAPI unmodified.

### StatsAPI verbs (response-block independent)

```julia
# Magesty imports `fit` / `coef` / `nobs` / `dof` from StatsAPI, adds
# methods, and re-exports them — `using Magesty` is enough; the user
# never imports StatsAPI. Because GLM.jl etc. re-export the same
# StatsAPI.fit, `using GLM, Magesty` produces no name clash.
# `intercept` is NOT a StatsAPI verb — StatsAPI has no intercept concept
# (GLM folds the intercept into `coef`). Magesty defines `intercept`
# natively because the SCE model keeps `j0` separate from `jphi`.
import StatsAPI: fit, coef, nobs, dof
export fit, coef, intercept, nobs, dof   # intercept is Magesty-native

fit(::Type{SCEFit}, dataset::SCEDataset, estimator::AbstractEstimator;
    torque_weight::Real = 0.5) -> SCEFit

coef(f::SCEFit)        -> jphi::Vector{Float64}   # SCE coefficients, one shared set
coef(m::SCEModel)      -> jphi::Vector{Float64}
intercept(f::SCEFit)   -> j0::Float64             # Magesty-native verb
intercept(m::SCEModel) -> j0::Float64
nobs(f::SCEFit)        -> Int                     # n_configs (energy-block observations)
dof(f::SCEFit)         -> Int                     # length(jphi) + 1
```

`coef` / `intercept` are bare (no suffix): the SCE coefficients are a
single set shared by the energy and torque models. `nobs` returns the
number of spin configurations (`n_configs`), i.e. the energy-block
observation count — torque components are derived per-config quantities,
not independent observations. The bare StatsAPI verbs `r2` / `rss` /
`residuals` / `predict` are **not** implemented — they assume a single
response and would force an implicit "energy" default.

### Evaluation & prediction verbs (Magesty, `(predictor, data)` form)

Every evaluation verb takes a predictor (`SCEModel` or `SCEFit`) and
evaluation data, following one regular overload pattern (shown for
`r2_energy`; identical for the whole family):

```julia
const SCEEvalData = Union{SCEDataset, AbstractVector{SpinConfig}, AbstractString}

r2_energy(predictor::Union{SCEModel, SCEFit}, data::SCEEvalData) -> Float64
r2_energy(f::SCEFit)                                            -> Float64  # in-sample: f's own training data
```

One `(predictor, data)` method handles `SCEModel` / `SCEFit` predictors
and `SCEDataset` / `Vector{SpinConfig}` / EMBSET-path data. `SCEFit`
predictors delegate through `SCEModel(f)`; the in-sample form is
`r2_energy(f) = r2_energy(f, f.dataset)`. The `SCEModel(f::SCEFit)`
conversion constructor (a few lines, add-only) is implemented in this
step rather than step 5, since the delegation depends on it.

The family, each with the overload set above:

| Verb | Returns |
|---|---|
| `r2_energy` / `r2_torque` | `Float64` — coefficient of determination |
| `rss_energy` / `rss_torque` | `Float64` — residual sum of squares |
| `residuals_energy` / `residuals_torque` | `Vector{Float64}` — residual vector |
| `rmse_energy` / `rmse_torque` | `Float64` — `√(rss / nobs)` convenience |
| `predict_energy` | `Float64` (single config) / `Vector{Float64}` (dataset) |
| `predict_torque` | `Matrix{Float64}` 3×n_atoms (single) / `Vector{Matrix{Float64}}` (dataset) |

`predict_energy` / `predict_torque` additionally accept the lowest-
level form — a single spin-direction matrix or a `SpinConfig`:

```julia
predict_energy(m::SCEModel, spin_directions::AbstractMatrix{<:Real}) -> Float64
predict_torque(m::SCEModel, spin_directions::AbstractMatrix{<:Real}) -> Matrix{Float64}
predict_energy(m::SCEModel, sc::SpinConfig) -> Float64
predict_torque(m::SCEModel, sc::SpinConfig) -> Matrix{Float64}
```

There is no `metrics` aggregator: callers ask for exactly the quantity
they want, for exactly the block they mean.

**Constraint**: the `data` argument must be built from the *same*
`SCEBasis` as the predictor — the design-matrix `(l, m, site)`
ordering must match. Enforced by a runtime identity check on the SALC
basis (`predictor`'s `salcbasis === data.basis.salcbasis`) when `data`
is a `SCEDataset`. `Vector{SpinConfig}` / EMBSET-path data carry no
design matrices and are evaluated through the predictor's own basis, so
they are compatible by construction and need no check.

## AtomsBase integration & Unitful boundary

All `AtomsBase` / `Unitful` use is confined to
`src/utils/atomsbase_adapter.jl`:

- `position(atom)` and `cell(system)` are `ustrip`-ed to angstrom
  `Float64` (`ustrip.(u"Å", ...)`). Non-length units error out via
  Unitful; nm / bohr auto-convert.
- Species / kind name: `AtomsBase.species(atom)` returns a
  `ChemicalSpecies`, but its `atom_name` is capped at 4 chars and is
  ignored by `==`, so it cannot carry a Wyckoff-style sublabel. The
  adapter instead reads the per-atom `atom.data[:atom_name]` entry
  (free-form `Symbol`, no length limit): the Magesty kind name is
  `get(atom.data, :atom_name, Symbol(atomic_symbol(species(atom))))`.
  Atoms sharing an element but carrying different `:atom_name` data
  get distinct kind indices; atoms with no `:atom_name` fall back to
  the element symbol.
- `lmax` / `lsum` / `cutoff` dicts are keyed by species `Symbol`. An
  element-only key (`:Fe`) fans out to every sublabel of that element
  as sugar; mixed element/species keys are an error.
- Distance-valued arguments the user writes directly — `cutoff` in
  `interaction`, and `lattice` in the keyword constructor — accept
  either a bare `Real` (taken as angstrom) or a `Unitful.Length`
  (auto-converted). A two-method helper normalizes them:
  `_to_angstrom(x::Unitful.Length) = ustrip(u"Å", x)` and
  `_to_angstrom(x::Real) = Float64(x)`. The `-1.0` no-cutoff sentinel
  passes through `_to_angstrom(::Real)` unchanged.
- Internal representation stays `Float64`; outputs stay `Float64`.
  Unitful never enters the stored fields of the four types.

**User perspective — when is `using Unitful` needed?**

| Input path | `using Unitful` |
|---|---|
| `SCEBasis("input.toml")` / `SCEBasis(dict)` | not needed |
| `SCEBasis(; lattice, kd, positions, ...)` (kwarg) | not needed |
| `SCEBasis(load_system("x.cif"))` (AtomsIO) | not needed |
| `SCEBasis(system)` with a hand-written AtomsBase system | needed — AtomsBase's `Atom` constructor requires Unitful positions |

The Magesty API never requires or returns a Unitful `Quantity`: inputs
are `ustrip`-ed at the boundary, outputs are `Float64`. Unitful surfaces
only if the user chooses to hand-write an AtomsBase system; the TOML,
keyword, and AtomsIO paths never expose it.

## save / load

```julia
save(obj, path::AbstractString)         # XML; obj in {SCEBasis, SCEModel}
load(::Type{T}, path::AbstractString)   # XML; T in {SCEBasis, SCEModel}
```

| Type | `.xml` |
|---|---|
| `SCEBasis` | structure + symmetry params + SALC basis (human-readable, full round-trip) |
| `SCEModel` | `SCEBasis` XML + `<JPhi>` block (full round-trip) |
| `SCEFit` | — (not covered; see below) |

Magesty's `save` / `load` handle XML only. The path must end in
`.xml`; any other extension is a clear error. `save` / `load` are kept
as the (ecosystem-aligned) names even though there is one format, so a
second format could be added later without an API change.

**`SCEFit` persistence is left to the user.** `SCEFit` / `SCEBasis` /
`SCEModel` are plain structs, so a user who wants to persist a fit
result writes `using JLD2; jldsave("fit.jld2"; fit = f)` directly —
JLD2 serializes the whole struct graph generically. Magesty does not
depend on JLD2 or wrap it; the recommended pattern is shown in
`examples/` (step 6). This keeps Magesty's `save` / `load` a single
XML code path and adds no dependency.

**XML writer.** The `xml_io.jl` writer is refactored to take
`SCEBasis` / `SCEModel` directly (it currently takes loose
`structure` / `symmetry` / `salc_basis` / `Optimizer` arguments). The
basis-only path (`SCEBasis`) and the basis+coeff path (`SCEModel`)
share the basis-writing code; `SCEModel` appends the `<JPhi>` block.

The schema gains two attributes so the file is self-describing:
- `tolerance_sym` on `<Symmetry>` — consumed on load.
- `isotropy` on `<SCEBasis>` — provenance metadata (a human cannot
  otherwise tell `isotropy = true` from `false` without inspecting
  every basis function for `Lf = 0`); not used by `load`.

Adding these attributes means the XML is no longer byte-identical to
the pre-refactor `write_xml` output — see requirements invariant 2 (it
is reframed to "stable schema, byte-regression against a committed
baseline").

**XML reader** reconstructs rather than blindly deserializing:
- `structure` ← `Structure(xml)` (constructor already exists).
- `symmetry` ← **recomputed** via `Symmetry(structure, tolerance_sym)`.
  `Symmetry` is a deterministic function of `structure` + `tol` (via
  Spglib), so it is not serialized in full — only `tolerance_sym` is
  stored. The `<Symmetry>` subtree the writer still emits (translation
  map etc.) is for human / external-tool consumption; `load` does not
  parse it.
- `salcbasis` ← `read_salcbasis_from_xml` (already exists).
- `SCEModel` additionally reads `<ReferenceEnergy>` / `<jphi>`.
- `load(SCEBasis, "x.xml")` accepts an `SCEModel` XML too (drops the
  `<JPhi>` block).

`write_xml` / `read_*_from_xml` stay during steps 1-6 (old API
coexistence) and are removed from the public surface in step 7; their
internals are reused by the `save` / `load` backends.

## Slicing semantics

```julia
length(d::SCEDataset)            -> Int               # n_configs
getindex(d::SCEDataset, i::Int)  -> SCEDataset         # 1 config, copy
getindex(d::SCEDataset, r)       -> SCEDataset         # copy (r: range, Vector, BitVector)
vcat(d1::SCEDataset, d2::SCEDataset) -> SCEDataset     # requires d1.basis === d2.basis
```

`getindex` copies (Julia convention). Row slicing of `X_E` (one row per
config) is mirrored on `X_T` (a `3*n_atoms` block per config) so the
two stay synchronized. `vcat` runtime-checks `basis` identity.
`view`-backed slicing is deferred to the CV follow-up spec (see
"Resolved decisions" #3).

## Usage examples

Representative scripts in the new API. These drive the design; the
`examples/` directory is rebuilt around them (tasklist Step 6).

### Basic flow

```julia
basis   = SCEBasis("input.toml")                                # heavy: SALC
dataset = SCEDataset(basis, "EMBSET.dat")                       # build design matrices
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)
println("RMSE energy: ", rmse_energy(f), "  R2 energy: ", r2_energy(f))
save(SCEModel(f), "model.xml")
```

### Building from a CIF file (AtomsIO)

```julia
using AtomsIO                                   # reads CIF / POSCAR / extxyz
system  = load_system("FePt.cif")               # -> AtomsBase.AbstractSystem
basis   = SCEBasis(system;
                   interaction = (body1 = (lmax = Dict(:Fe => 2, :Pt => 2),),
                                  body2 = (lsum = 2,
                                           cutoff = Dict((:Fe, :Fe) => -1.0,
                                                         (:Fe, :Pt) => -1.0,
                                                         (:Pt, :Pt) => -1.0))))
dataset = SCEDataset(basis, "EMBSET.dat")
f       = fit(SCEFit, dataset, Ridge(lambda = 1e-4); torque_weight = 0.5)
```

`AtomsIO` is not a dependency of Magesty itself; it is the recommended
reader and is used in `examples/`. The user installs it alongside
Magesty when reading structure files.

### Estimator comparison (reuse one dataset)

```julia
dataset = SCEDataset(SCEBasis("input.toml"), "EMBSET.dat")
for est in (OLS(), Ridge(lambda = 1e-4), Ridge(lambda = 1e-2))
    f = fit(SCEFit, dataset, est; torque_weight = 0.5)
    println(est, ": ", rmse_energy(f))
end
```

### Train / test split (slice + out-of-sample evaluation)

```julia
dataset = SCEDataset(SCEBasis("input.toml"), "EMBSET.dat")
n_train = round(Int, 0.8 * length(dataset))
train, test = dataset[1:n_train], dataset[n_train+1:end]
f = fit(SCEFit, train, Ridge(lambda = 1e-4); torque_weight = 0.5)
println("in-sample:  ", rmse_energy(f))
println("out-sample: ", rmse_energy(f, test))
```

### Cache the SALC basis, reuse in a later session

```julia
save(SCEBasis("input.toml"), "basis.jld2")    # cache the heavy SALC result
# --- later session ---
basis = load(SCEBasis, "basis.jld2")
f = fit(SCEFit, SCEDataset(basis, "EMBSET.dat"), Ridge(lambda = 1e-4))
save(SCEModel(f), "model.xml")                # human-readable, byte-identical to old write_xml
```

## TOML as a template

TOML is not just a back-compat shim — it is the recommended
*declarative template* for the material + interaction definition. A
new user copies an example `input.toml`, edits the structure and
interaction sections, and loads it with `SCEBasis("input.toml")`. The
schema loses `[regression]`:

```toml
[general]      # kd, name, nat, periodicity            -> SCEBasis
[symmetry]     # tolerance, isotropy                   -> SCEBasis
[interaction]  # nbody, body1.lmax, bodyN.lsum/cutoff  -> SCEBasis
[structure]    # kd_list, lattice, position            -> SCEBasis
# [regression] REMOVED — lambda / torque_weight / datafile are Julia-side
```

The TOML `[interaction]` table maps 1:1 onto the nested-`NamedTuple`
`interaction` argument of the AtomsBase constructor
(`[interaction.body1]` ↔ `body1 = (...)`, etc.), so the two input
paths are structurally the same. `SCEBasis("input.toml")` parses the
four sections; fit parameters are passed in Julia at `fit` time. TOML
`lmax` / `cutoff` keys stay element-only (no sublabels in the TOML
path).

## Migration map (old API -> new API)

| Old | New |
|---|---|
| `System(input_dict)` | `SCEBasis(input_dict)` |
| `build_sce_basis(input_dict)` | `SCEBasis(input_dict)` |
| `build_sce_basis_from_xml(input_dict, xml)` | `load(SCEBasis, xml)` |
| `SpinCluster(input_dict)` | `SCEDataset(...)` then `fit(SCEFit, ...)` |
| `fit_sce_model(system, configs, est, weight)` | `fit(SCEFit, dataset, est; torque_weight)` |
| `write_xml(sc, path)` | `save(model, path)` |
| `Optimizer` (struct/export) | internal; `SCEFit` is the user-facing result |
| `predict_energy` | `predict(model, ...)` |
| `get_j0(sc)` / `get_jphi(sc)` | `intercept(fit)` / `coef(fit)` |
| `get_j0_jphi(sc)` | `(intercept(fit), coef(fit))` |
| `calc_energy(sc, dirs)` | `predict(model, dirs)` |
| `calc_torque(sc, dirs)` | `predict_torque(model, dirs)` |

## Affected coupling points

- **`xml_io.jl`**: `save`/`load` backends. The writer is refactored to
  take `SCEBasis` / `SCEModel` (not loose args + `Optimizer`); the
  schema gains `tolerance_sym` / `isotropy` attributes. New readers:
  the `Symmetry` recompute path and the `<JPhi>` block reader. Byte-
  regression against a regenerated `fept` / `fege` baseline.
- **`Optimize.jl`**: `assemble_weighted_problem` already does
  per-sample MSE normalization; `fit(SCEFit, ...)` calls it with
  `torque_weight`. `extract_j0_jphi` unchanged.
- **Design matrix construction**: `build_design_matrix_energy` /
  `build_design_matrix_torque` move from being called inside
  `Optimizer`/`fit_sce_model` to being called inside the `SCEDataset`
  constructor (so `X_E`/`X_T` are built once and reused).
- **`SCEModel` reshape** to `{basis::SCEBasis, j0, jphi}`: update
  `predict_energy` / `predict_torque` (read `model.basis.salcbasis` /
  `model.basis.symmetry`), `SCEModel(f::SCEFit)`, and every
  `SCEModel(...)` construction site.
- **`SCEBasis` drops `cluster`**: the constructor still computes a
  `Cluster` locally to build `SALCBasis`, then discards it. Nothing
  reads `SCEBasis.cluster` post-construction (verified across `src/`).

## Test strategy

- Each step (1-6) keeps `make test-all` green; step 7 (old API
  removal) is the breaking commit.
- New unit tests per type: `test_SCEBasis.jl`, `test_SCEDataset.jl`
  (incl. slicing / `vcat`), `test_SCEFit.jl` (StatsAPI verbs),
  `test_save_load.jl` (extension dispatch, round-trip).
- `SCEModel` XML byte-identical: a round-trip diff test against a
  committed baseline for `fept` and `fege`.
- `test/examples/*` drivers rewritten in step 6; their assertions
  (where they have any) must still pass.

## Resolved decisions

1. **Type placement**: `SCEBasis` / `SCEDataset` / `SCEFit` / `SCEModel`
   all live at `Magesty.jl` top level (like `System` / `SpinCluster`
   today). They are integration types depending on multiple submodules,
   so a submodule home would create a dependency cycle. `SCEModel` was
   moved out of `Optimize` in step 5a: reshape #5 makes it hold a
   `SCEBasis`, which `Optimize` (included before `SCEBasis` is defined)
   cannot reference.
2. **No `metrics` aggregator**: dropped in favor of individual
   StatsAPI-style verbs (`r2_energy`, `rmse_torque`, ...). See "Verbs".
3. **`SCEDataset` slicing**: `getindex` (copy) only in this spec.
   `view`-backed slicing is deferred to the CV follow-up spec, which
   will decide the type design (parametrized `SCEDataset` vs a
   dedicated view type) against real k-fold requirements.
4. **`Cluster` is not a stored field of `SCEBasis`** (decided in step
   5, after reading `xml_io.jl`). `Cluster` is a construction step,
   not persistent state — nothing reads `SCEBasis.cluster` after
   construction (verified across `src/`). It stays an internal struct;
   only its role changes. This also makes `SCEBasis` exactly the
   serializable triple `structure` / `symmetry` / `salcbasis`.
5. **`SCEModel` holds a `SCEBasis`** rather than flat `salcbasis` /
   `symmetry` / `num_atoms` fields. It needs `structure` for a
   complete XML, and holding `basis` lets `save` reuse the `SCEBasis`
   XML writer.
6. **`Symmetry` is recomputed, not serialized, on XML load.**
   `Symmetry` is a deterministic function of `structure` +
   `tolerance_sym`. Storing `tolerance_sym` (one attribute) is enough;
   serializing the full `Symmetry` struct would be redundant. This
   drove the invariant 2 reframe: the schema gains attributes, so
   strict byte-identity with the pre-refactor output is dropped in
   favor of byte-regression against a new committed baseline.
7. **JLD2 is not a Magesty dependency.** Magesty's `save` / `load` are
   XML-only. `SCEFit` / `SCEBasis` / `SCEModel` are plain structs, so
   users who want struct-level persistence (the only way to persist a
   `SCEFit`) use `JLD2` directly. This keeps `save` / `load` a single
   code path and avoids a heavy dependency for a one-line user-side
   call.
