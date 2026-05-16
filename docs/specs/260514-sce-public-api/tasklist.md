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

- [x] Extend the XML schema: `tolerance_sym` attribute on `<Symmetry>`,
      `isotropy` attribute on `<SCEBasis>`. `tolerance_sym` is consumed
      on load; `isotropy` is stored as the new 4th `SCEBasis` field
      (decided in step 5b — provenance, not derived from `salcbasis`).
- [x] Refactor the `xml_io.jl` writer to drop the `Optimizer`
      argument; `XMLIO` writers take loose components (it is included
      before the `SCEBasis` / `SCEModel` types), `SCEBasis` /
      `SCEModel` unpacking is in the `Magesty.jl` `save` / `load`
      methods. Basis-only and basis+coeff paths share
      `_write_system_subtree!`.
- [x] XML readers: reconstruct `SCEBasis` via `Structure(xml)` +
      `Symmetry(structure, tolerance_sym)` (recompute, not
      deserialize) + `read_salcbasis_from_xml`. `SCEModel` adds the
      `<JPhi>` reader. `load(SCEBasis, xml)` accepts an `SCEModel` XML.
- [x] `save(obj, path)` / `load(::Type{T}, path)` — XML only, `.xml`
      extension required, clear error otherwise.
- [x] `make test-all` green (unit 6257, jet 0 issues, aqua 10,
      integration 155).

### Step 5c — tests + baselines

- [x] Add `test_save_load.jl`: `.xml` round-trip for `SCEBasis` /
      `SCEModel`, basis-identity (design matrices) after reload,
      non-`.xml` extension error path, `load(SCEBasis)` from an
      `SCEModel` XML.
- [x] Regenerate the committed baseline XML for `fept` / `fege` under
      the new schema (`test/component_test/baselines/`); byte-diff
      regression test against it. Two checks: fept fresh-build
      (`SCEBasis` + `SCEModel`, byte-reproducible — small,
      degeneracy-free) and a load -> re-save round-trip for fept +
      fege (deterministic, no SALC eigensolve — robust on the large
      `fege` system). Baselines generated under `--check-bounds=yes`
      to match the `Pkg.test` environment. `read_salcbasis_from_xml`
      now preserves the `<AngularMomentumCouplings>` file order
      (was `collect(values(::Dict))`) so the round-trip is byte-exact.
- [x] `make test-all` green (unit 6314, jet 0 issues, aqua 10,
      integration 155).

## Step 6 — migrate `test/examples/*`, rebuild `examples/`, docs

- [x] Rewrite each `test/examples/*/test.jl` to the new API
      (`SCEBasis` / `SCEDataset` / `fit` / `save`). Six examples migrated;
      `test/examples/2d_fcc_2x2x2/` deleted (was orphaned — not wired
      into `runtests.jl`, content duplicated by `chain`). The legacy
      `chain` design-matrix factor assertion (`4√3`) was replaced by the
      physics-level `coef * √3 ≈ -1.0` check — the legacy expectation
      was already stale relative to the current `build_design_matrix_*`
      output.
- [x] Remove `[regression]` from every `test/examples/*/input.toml`
      (5 files). `Config4System` does not consume the `[regression]`
      table, so the new API ignores it; stripping it makes the schema
      match what the new TOML template advertises.
- [x] Delete `examples/fept/` (stale data dump); add 3 runnable new-API
      example scripts: `examples/01_basic_flow.jl`,
      `examples/02_cif_input.jl`, `examples/03_save_load.jl`. The
      basic-flow and save-load scripts reuse the FePt L1_0 fixture
      under `test/examples/`; the CIF script writes a minimal CIF
      inline and uses `AtomsIO.load_system`. `AtomsIO` is an
      `examples/`-only dependency (not in `Project.toml`).
      Estimator-comparison and train-test split are deferred — covered
      by `docs/src/examples.md` but not as standalone scripts.
- [x] Update `docs/src/api.md`, `docs/src/examples.md`, `SPEC.md`.
      `README.md` left to Step 7 (folded into the breaking commit).
- [x] `make test-all` green (unit 6314, integration 813); JET 0
      issues; Aqua 10.

## Step 7 — remove old API (breaking commit)

- [x] Remove from `Magesty`: `System`, `SpinCluster`,
      `build_sce_basis`, `build_sce_basis_from_xml`, `fit_sce_model`,
      `write_xml`, `Optimizer` export, `get_j0` /
      `get_jphi` / `get_j0_jphi`, `calc_energy` / `calc_torque`,
      `write_energies` / `write_torques`. `predict_energy` /
      `predict_torque` survive as the new-API verbs (the legacy
      `Optimizer`-typed methods are gone with `Optimizer`).
- [x] Delete now-dead code paths. `Optimizer` struct removed (fully
      superseded by `SCEFit`); `_fit_sce_model_internal`,
      `fit_sce_model_ols`, `fit_sce_model_ridge`, `_default_estimator`,
      `print_sce_coeffs`, `print_metrics` deleted. The entire
      `src/utils/EnergyTorque.jl` module deleted (only the legacy
      `calc_energy`/`calc_torque` wrappers used it). Legacy
      `write_xml(::System)` / `write_xml(::SpinCluster)` / `write_xml(...,
      ::Optimizer, ...)` overloads removed from `xml_io.jl`. Removed
      `_isotropy_from_salcbasis` (legacy-only) and the `using
      ..Optimize` dependency in `xml_io.jl`.
- [x] Audit fold-ins (see Step 7 audit triage, 2026-05-15):
      - B7+B8 SALCBases debug prints (`@show` / `display`) removed.
      - B9 dead helpers (`push_unique_coupled_basis!`,
        `flip_vector_if_negative_sum`) deleted.
      - B10+B11 export trim (`detect_num_atoms`, `cluster_orbits`
        no longer exported).
      - B12 `classify_coupled_basislist_test` moved out of
        `src/SALCBases.jl` into `test/benchmark_salcbasis_hotspots.jl`
        as a local helper.
      - B13 `read_embset` re-export pattern fixed
        (`const = …; export` → `using .SpinConfigs: read_embset; export
        read_embset`).
      - C1 module docstring updated to the new API.
      - C3 `x_frac` shape docstring tightened to `[3 × num_atoms]`.
      - C4 docstring style normalized (`**Arguments:**` → `# Arguments`).
      - C5 SALCBases unitary mismatch — message corrected (the
        variable is the single-op representation matrix `D(g)`, not the
        projector; check is for unitary irrep, not `P` unitarity);
        `@warn "Critical error"` replaced with `error(...)`.
      - Critical [B] items (B1–B6, B14) deferred to follow-ups —
        tracked in `docs/design-notes/post-step7-cleanup.md`.
      - `Config4Optimize` (and `REQUIRED_SECTIONS_OPTIMIZE` /
        `DEFAULT_VALUES_OPTIMIZE` / `validate_optimize_parameters`)
        deleted as dead code: the only consumer was the legacy
        `Optimizer` constructor, gone with this commit. The
        `Config4Optimize` test set in `test_ConfigParser.jl` was
        removed in turn.
      - `docs/src/input_keys.md` rewritten: the `[regression]`
        section is no longer part of the TOML schema (the new API
        does not consume it). The `weight` description had been
        inverted (it now correctly says `0 = energy only`,
        `1 = torque only` and is moved to a Julia-side
        `torque_weight` example with default `0.5`).
- [x] `grep` for old symbol names across `src/`, `test/component_test/`,
      `test/examples/`, `docs/src/`, `SPEC.md`, `README.md` returns no
      hits. `test/develop_tmp/` (gitignored, slated for deletion),
      `test/benchmark_*.jl`, `test/profile_run.jl`, and `tools/` are
      out of scope.
- [x] `make test-all` (unit 6314, integration 813), `make test-jet`
      (0 issues), `make test-aqua` (10) all green.
- [x] `docs/design-notes/sce-public-api.md` and `DESIGN_NOTES.md`
      updated. Post-Step-7 follow-up note:
      `docs/design-notes/post-step7-cleanup.md`.
- [x] `code-reviewer` agent pass on the cumulative diff.
- [x] Single breaking commit via `git-helper` (`091e9f8
      refactor(sce-api)!: remove legacy API, finalize four-type public
      surface`); merged to `main` via PR #4 (`f4ae8cf`).

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
