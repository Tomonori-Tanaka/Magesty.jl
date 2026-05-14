# Tasklist: SCE public API refactor

Status: draft (2026-05-14)
Requirements: `./requirements.md` — Design: `./design.md`

Coarse milestones. Day-to-day work tracked with `TaskCreate`.
Branch: `refactor/sce-public-api`. Through step 5 the old API keeps
working alongside the new one (coexistence); step 5 also reshapes the
*new* types (`SCEBasis` / `SCEModel`) as the design firms up; step 6
rewrites tests; step 7 is the breaking removal of the old API.
`make test-all` must stay green at every commit.

## Phase 0 — sign-off

- [x] User agrees on the four-type shape (`SCEBasis` / `SCEDataset` /
      `SCEFit` / `SCEModel`).
- [x] User agrees on `SCEBasis` storing `Structure` (not the raw
      `AbstractSystem`) — Unitful kept out of stored fields.
- [x] User agrees on the resolved decisions in design.md (see
      "Resolved decisions" 1-7 — type placement, no `metrics`
      aggregator, `view` deferred, `cluster` not stored, `SCEModel`
      holds `SCEBasis`, `Symmetry` recomputed on load, JLD2 not a dep).

## Prerequisites (done)

- [x] `SALCBases` rename (PR #2).
- [x] per-sample MSE normalization in `assemble_weighted_problem`
      (PR #3).

## Step 1 — `SCEBasis` type + AtomsBase input

- [x] Add deps `AtomsBase`, `Unitful` to `Project.toml` with compat
      bounds.
- [x] New `src/utils/atomsbase_adapter.jl`: `AbstractSystem` →
      `Structure` conversion, `ustrip` to angstrom, `ChemicalSpecies`
      sublabel → kind name, element-only `lmax`/`cutoff` fan-out.
- [x] Define `SCEBasis` (structure + symmetry + cluster + salcbasis).
      Note: `cluster` is dropped from the fields in step 5a.
- [x] Constructors: `SCEBasis(::AbstractSystem; nbody, lmax, ...)`,
      `SCEBasis(::AbstractString)` (TOML), `SCEBasis(::AbstractDict)`.
- [x] `SCEBasis` constructor: `@warn` on odd `lmax` / `lsum` values
      (time-reversal symmetry note); covers both input paths.
- [x] `SCEModel` field renames: `reference_energy` → `j0`, `SCE` →
      `jphi`, `basisset` → `salcbasis`. Update `EnergyTorque.jl` /
      `Optimize.jl` read & construction sites. XML tag names unchanged.
- [x] Add `test_SCEBasis.jl`: construct from AtomsBase system, from
      TOML; species sublabel handling; element-only fan-out.
- [x] `make test-all` green (old API still exported and working).

## Step 2 — `SCEDataset` type + slicing

- [x] Define `SCEDataset` (basis + spinconfigs + unweighted
      `X_E`/`X_T`/`y_E`/`y_T`).
- [x] Move `build_design_matrix_energy` / `build_design_matrix_torque`
      calls into the `SCEDataset` constructor.
- [x] Constructors: `SCEDataset(basis, spinconfigs)`,
      `SCEDataset(basis, embset_path)`, sugar
      `SCEDataset(system, spinconfigs; nbody, ...)`,
      `SCEDataset(toml_path, spinconfigs)`.
- [x] Slicing: `length`, `getindex` (copy), `vcat` (basis-identity
      check). `X_E` row slice mirrored to `X_T` block slice. `view`
      deferred to the CV follow-up spec.
- [x] Add `test_SCEDataset.jl`: construction, slicing copy semantics,
      `vcat`, row/block synchronization.
- [x] `make test-all` green.

## Step 3 — `SCEFit` + StatsAPI verbs

- [x] Add dep `StatsAPI` to `Project.toml`.
- [x] Define `SCEFit <: StatsAPI.RegressionModel`.
- [x] `fit(::Type{SCEFit}, dataset, estimator; torque_weight = 0.5)`:
      `assemble_weighted_problem` + `solve_coefficients` +
      `extract_j0_jphi`, packaging residuals + metrics.
- [x] Rename the `weight` kwarg to `torque_weight` along the call
      path. (New `fit` uses `torque_weight`; the internal
      `assemble_weighted_problem` keeps its positional `weight` arg,
      removed/refactored in Step 7.)
- [x] Implement the response-independent verbs `coef`, `intercept`,
      `nobs`, `dof` for `SCEFit`; `coef` / `intercept` for `SCEModel`.
      `fit` / `coef` / `nobs` / `dof` come from StatsAPI; `intercept`
      is Magesty-native (StatsAPI has no intercept verb). Bare `r2` /
      `rss` / `residuals` / `predict` are not implemented — see
      design.md "Verbs".
- [x] Add `test_SCEFit.jl`: fit an example, check `coef` / `intercept`
      / `nobs` / `dof`; golden `(j0, jphi)` vs `fit_sce_model` on the
      same inputs.
- [x] `make test-all` green.

## Step 4 — evaluation & prediction verbs

- [x] `SCEModel(f::SCEFit)` lightweight conversion. Pulled forward from
      step 5: the `(fit, data)` overloads delegate through it.
- [x] Implement the `(predictor, data)` verb family: `r2_energy` /
      `r2_torque`, `rss_energy` / `rss_torque`, `residuals_energy` /
      `residuals_torque`, `rmse_energy` / `rmse_torque`.
- [x] One `(predictor::Union{SCEModel,SCEFit}, data::SCEEvalData)`
      method per verb, plus `(f::SCEFit)` in-sample. `SCEEvalData =
      Union{SCEDataset, Vector{SpinConfig}, AbstractString}`. SALC
      basis-identity check when `data` is a `SCEDataset`.
- [x] `predict_torque(::SCEModel, spin_directions)` added to
      `Optimize.jl` (mirrors `predict_energy`). `predict_energy` /
      `predict_torque` overloads: `(model, sc::SpinConfig)`,
      `(fit, ...)`, and dataset-batch forms.
- [x] Extend `test_SCEFit.jl`: SCEModel conversion, predict_* agreement
      across forms, in-sample verbs vs cached metrics, configs/dataset
      paths, basis-identity check error path.
- [x] `make test-all` green (unit 6259, jet 0 issues, aqua 10,
      integration 155).

## Step 5 — type reshape + `save` / `load`

Driven by reading `xml_io.jl`: `Cluster` is a construction step, not
persistent state, and `Symmetry` is a derived quantity. See design.md
"Resolved decisions" 4-7. JLD2 is *not* a dependency — `SCEFit`
persistence is left to the user (`jldsave` on a plain struct).

### Step 5a — type reshape (no new deps)

- [x] `SCEModel(f::SCEFit)` lightweight conversion — done in step 4
      (re-pointed at the reshaped `SCEModel` in this step).
- [x] Drop `cluster` from `SCEBasis`; the constructor computes it
      locally for `SALCBasis` construction and discards it. `Cluster`
      stays an internal struct (a construction intermediate).
- [x] Reshape `SCEModel` to `{basis::SCEBasis, j0, jphi}`; `num_atoms`
      becomes `basis.structure.supercell.num_atoms`. `SCEModel` moves
      from `Optimize` to `Magesty.jl` top level (it now references
      `SCEBasis`); a checked inner constructor enforces
      `length(jphi) == n_salc`.
- [x] Update `predict_energy` / `predict_torque` / `SCEModel(f)` and
      `test_SCEFit.jl` for the reshaped types. `predict_*` base methods
      delegate to `Optimize._predict_energy` / `_predict_torque`
      (raw-argument internal cores).
- [x] `make test-all` green (old API untouched; coexistence holds).

### Step 5b — XML `save` / `load`

- [ ] Extend the XML schema: `tolerance_sym` attribute on `<Symmetry>`,
      `isotropy` attribute on `<SCEBasis>`. `isotropy` is provenance
      metadata; `tolerance_sym` is consumed on load.
- [ ] Refactor the `xml_io.jl` writer to take `SCEBasis` / `SCEModel`
      (not loose args + `Optimizer`); basis-only and basis+coeff paths
      share code.
- [ ] XML readers: reconstruct `SCEBasis` via `Structure(xml)` +
      `Symmetry(structure, tolerance_sym)` (recompute, not
      deserialize) + `read_salcbasis_from_xml`. `SCEModel` adds the
      `<JPhi>` reader. `load(SCEBasis, xml)` accepts an `SCEModel` XML.
- [ ] `save(obj, path)` / `load(::Type{T}, path)` — XML only, `.xml`
      extension required, clear error otherwise.
- [ ] `make test-all` green.

### Step 5c — tests + baselines

- [ ] Add `test_save_load.jl`: `.xml` round-trip for `SCEBasis` /
      `SCEModel`, basis-identity after reload, non-`.xml` extension
      error path.
- [ ] Regenerate the committed baseline XML for `fept` / `fege` under
      the new schema; byte-diff regression test against it.
- [ ] `make test-all` green.

## Step 6 — migrate `test/examples/*`, rebuild `examples/`, docs

- [ ] Rewrite each `test/examples/*/test.jl` (and `run.jl`) to the
      new API: `SCEBasis` / `SCEDataset` / `fit` / `save`.
- [ ] Remove `[regression]` from every `test/examples/*/input.toml`.
- [ ] Delete `examples/fept/` (stale data dump); add runnable new-API
      example scripts under `examples/` — a light system with EMBSET
      data, covering basic flow / building from a CIF file (AtomsIO) /
      estimator comparison / train-test split / save-load (mirroring
      design.md "Usage examples"). `AtomsIO` is an `examples/`-only
      dependency, not a package dependency.
- [ ] Update `docs/src/api.md`, `docs/src/examples.md`, `SPEC.md`,
      `README.md`.
- [ ] `make test-all` green; example assertions still pass;
      `SCEModel` / `SCEBasis` XML still matches the committed baseline.

## Step 7 — remove old API (breaking commit)

- [ ] Remove from `Magesty`: `System`, `SpinCluster`,
      `build_sce_basis`, `build_sce_basis_from_xml`, `fit_sce_model`,
      `write_xml`, `predict_energy`, `Optimizer` export, `get_j0` /
      `get_jphi` / `get_j0_jphi`, `calc_energy` / `calc_torque`.
- [ ] Delete now-dead code paths; `Optimizer` struct becomes internal
      or is removed if fully superseded by `SCEFit`.
- [ ] `grep` for old symbol names across `src/` `test/` `docs/`
      returns no hits (except historical mentions in `docs/specs/`
      and `docs/design-notes/`).
- [ ] `make test-all`, `make test-jet`, `make test-aqua` green.
- [ ] Mark `docs/design-notes/sce-public-api.md` and `DESIGN_NOTES.md`
      complete.
- [ ] `code-reviewer` agent pass on the cumulative diff.
- [ ] Single breaking commit via `git-helper`; merge to `main` after
      explicit user confirmation.

## After merge

- [ ] Follow-up specs (out of scope here): molecular / aperiodic
      support, `CV` wrapper estimator + `kfold`, Bayesian
      `AbstractSCEFit`.
- [ ] Add an `examples/` script for a system where `atom_name`
      sublabels matter — same element on inequivalent Wyckoff sites
      that the user wants treated as distinct kinds (distinct `lmax`
      etc.). Deferred until a concrete such system is on hand.

## Risks tracked

- AtomsBase API churn (pre-1.0): all calls isolated in
  `atomsbase_adapter.jl`; compat bound pinned.
- XML schema drift: caught by the byte-diff round-trip test in step 5c
  against the regenerated `fept` / `fege` baselines — run it before
  every commit that touches `xml_io.jl`.
- Intermediate commits non-importable: old / new APIs coexist through
  step 5;
  `make test-all` at every commit.
