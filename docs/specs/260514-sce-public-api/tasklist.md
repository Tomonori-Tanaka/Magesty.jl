# Tasklist: SCE public API refactor

Status: draft (2026-05-14)
Requirements: `./requirements.md` — Design: `./design.md`

Coarse milestones. Day-to-day work tracked with `TaskCreate`.
Branch: `refactor/sce-public-api`. Steps 1-5 are add-only (old and new
APIs coexist); step 6 rewrites tests; step 7 is the breaking removal.
`make test-all` must stay green at every commit.

## Phase 0 — sign-off

- [ ] User agrees on the four-type shape (`SCEBasis` / `SCEDataset` /
      `SCEFit` / `SCEModel`).
- [ ] User agrees on `SCEBasis` storing `Structure` (not the raw
      `AbstractSystem`) — Unitful kept out of stored fields.
- [ ] User agrees on the resolved decisions in design.md (top-level
      type placement, no `metrics` aggregator, `view` deferred to the
      CV spec).

## Prerequisites (done)

- [x] `SALCBases` rename (PR #2).
- [x] per-sample MSE normalization in `assemble_weighted_problem`
      (PR #3).

## Step 1 — `SCEBasis` type + AtomsBase input

- [ ] Add deps `AtomsBase`, `Unitful` to `Project.toml` with compat
      bounds.
- [ ] New `src/utils/atomsbase_adapter.jl`: `AbstractSystem` →
      `Structure` conversion, `ustrip` to angstrom, `ChemicalSpecies`
      sublabel → kind name, element-only `lmax`/`cutoff` fan-out.
- [ ] Define `SCEBasis` (structure + symmetry + cluster + salcbasis).
- [ ] Constructors: `SCEBasis(::AbstractSystem; nbody, lmax, ...)`,
      `SCEBasis(::AbstractString)` (TOML), `SCEBasis(::AbstractDict)`.
- [ ] `SCEBasis` constructor: `@warn` on odd `lmax` / `lsum` values
      (time-reversal symmetry note); covers both input paths.
- [ ] `SCEModel` field renames: `reference_energy` → `j0`, `SCE` →
      `jphi`, `basisset` → `salcbasis`. Update `EnergyTorque.jl` /
      `Optimize.jl` read & construction sites. XML tag names unchanged.
- [ ] Add `test_SCEBasis.jl`: construct from AtomsBase system, from
      TOML; species sublabel handling; element-only fan-out.
- [ ] `make test-all` green (old API still exported and working).

## Step 2 — `SCEDataset` type + slicing

- [ ] Define `SCEDataset` (basis + spinconfigs + unweighted
      `X_E`/`X_T`/`y_E`/`y_T`).
- [ ] Move `build_design_matrix_energy` / `build_design_matrix_torque`
      calls into the `SCEDataset` constructor.
- [ ] Constructors: `SCEDataset(basis, spinconfigs)`,
      `SCEDataset(basis, embset_path)`, sugar
      `SCEDataset(system, spinconfigs; nbody, ...)`,
      `SCEDataset(toml_path, spinconfigs)`.
- [ ] Slicing: `length`, `getindex` (copy), `vcat` (basis-identity
      check). `X_E` row slice mirrored to `X_T` block slice. `view`
      deferred to the CV follow-up spec.
- [ ] Add `test_SCEDataset.jl`: construction, slicing copy semantics,
      `vcat`, row/block synchronization.
- [ ] `make test-all` green.

## Step 3 — `SCEFit` + StatsAPI verbs

- [ ] Add dep `StatsAPI` to `Project.toml`.
- [ ] Define `SCEFit <: StatsAPI.RegressionModel`.
- [ ] `fit(::Type{SCEFit}, dataset, estimator; torque_weight = 0.5)`:
      `assemble_weighted_problem` + `solve_coefficients` +
      `extract_j0_jphi`, packaging residuals + metrics.
- [ ] Rename the `weight` kwarg to `torque_weight` along the call
      path.
- [ ] Implement the StatsAPI response-independent verbs `coef`,
      `intercept`, `nobs`, `dof` for `SCEFit`; `coef` / `intercept`
      for `SCEModel`. (Bare `r2` / `rss` / `residuals` / `predict` are
      not implemented — see design.md "Verbs".)
- [ ] Add `test_SCEFit.jl`: fit an example, check `coef` / `intercept`
      / `nobs` / `dof`; golden `(j0, jphi)` vs `fit_sce_model` on the
      same inputs.
- [ ] `make test-all` green.

## Step 4 — evaluation & prediction verbs

- [ ] Implement the `(predictor, data)` verb family: `r2_energy` /
      `r2_torque`, `rss_energy` / `rss_torque`, `residuals_energy` /
      `residuals_torque`, `rmse_energy` / `rmse_torque`.
- [ ] Overload set for each: `(model, dataset)`, `(fit, dataset)`,
      `(fit)` (in-sample), `(target, embset_path)`,
      `(target, configs)`. Normalize data to `SCEDataset`; runtime
      `basis`-identity check.
- [ ] `predict_energy` / `predict_torque`: `(model, spin_directions)`,
      `(model, sc::SpinConfig)`, and dataset-batch forms.
- [ ] Extend `test_SCEFit.jl`: in-sample vs out-of-sample evaluation,
      `predict_torque` shape/values vs `calc_torque`, basis-identity
      check error path.
- [ ] `make test-all` green.

## Step 5 — `SCEModel(fit)` + `save` / `load`

- [ ] `SCEModel(f::SCEFit)` lightweight conversion.
- [ ] `save(obj, path)` / `load(::Type{T}, path)` with extension
      dispatch.
- [ ] XML backend: refactor `xml_io.jl` so basis-only (`SCEBasis`) and
      basis+coeff (`SCEModel`) share code; `SCEModel` XML
      byte-identical.
- [ ] JLD2 backend (add `JLD2` dep): `SCEBasis` / `SCEFit` /
      `SCEModel`.
- [ ] Add `test_save_load.jl`: extension dispatch, round-trip,
      `SCEModel` XML byte-identical vs committed baseline (`fept`,
      `fege`).
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
      `SCEModel` XML still byte-identical.

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
- `SCEModel` XML byte-identical regression: caught only by the
  round-trip diff test in Step 5 — run it before every commit that
  touches `xml_io.jl`.
- Intermediate commits non-importable: Steps 1-5 are add-only;
  `make test-all` at every commit.
