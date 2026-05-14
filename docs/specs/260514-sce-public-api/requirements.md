# Requirements: SCE public API refactor (four types + StatsAPI + AtomsBase)

Status: draft (2026-05-14)
Owner: T. Tanaka
Branch: `refactor/sce-public-api` (cut when this spec is agreed)
Design notes: `docs/design-notes/sce-public-api.md` (original proposal, superseded by this spec on naming and a few decisions)

## Goal

Replace the organically-grown public API (`System` / `SpinCluster` /
`build_sce_basis` / `build_sce_basis_from_xml` / `fit_sce_model` /
`write_xml`) with four clearly-separated types and verbs aligned to the
Julia statistics ecosystem (`StatsAPI.jl`). Make material input
first-class via `AtomsBase.jl` so TOML becomes one input path among
several rather than the only one. The package is pre-release
(`0.1.0-DEV`); breaking changes are accepted now in favor of the ideal
shape.

### The four types

| Type | Replaces | Content (concept) | Weight |
|---|---|---|---|
| `SCEBasis` | `System` (basis side), `build_sce_basis`, `build_sce_basis_from_xml` | material (structure) + symmetry + SALC basis | heavy (SALC), independently persistable |
| `SCEDataset` | (new abstraction) | `basis::SCEBasis` + unweighted `X_E`/`X_T`/`y_E`/`y_T` + spinconfigs | heavy, reusable, sliceable |
| `SCEFit` | `SpinCluster`, `fit_sce_model`, `Optimizer` (export) | fit result: coefficients + metrics + residuals + estimator/weight used | medium; `<: StatsAPI.RegressionModel` |
| `SCEModel` | (kept, already public) | lightweight predictor: `basis::SCEBasis` + `j0` + `jphi` | lightweight, serializable |

`SCEBasis` was called `SCEProblem` in the original design note; renamed
because it holds no data and cannot be "solved" — "Problem" carries the
wrong connotation (cf. `ODEProblem`). `SCEBasis` follows the package
convention that the central type of a module/concept is the singular
noun.

## Prerequisites (already merged)

These were completed as separate PRs before this spec and are assumed
done:

- **`SALCBases` rename** (PR #2): internal `BasisSets`/`BasisSet` →
  `SALCBases`/`SALCBasis`, XML tag `SCEBasisSet` → `SCEBasis`. This
  freed the name `SCEBasis` for the public type introduced here.
- **per-sample MSE normalization** (PR #3):
  `Optimize.assemble_weighted_problem` now scales energy rows by
  `√((1-w)/n_E)` and torque rows by `√(w/n_T)`. The `torque_weight`
  semantics this spec exposes are already implemented at the solver
  level; this spec only renames the kwarg (`weight` → `torque_weight`)
  and threads it through the new `fit` signature.

## Scope

In scope:

- **New types**: `SCEBasis`, `SCEDataset`, `SCEFit`; `SCEModel` kept
  but reshaped to hold a `SCEBasis` plus the fitted `j0` / `jphi`
  (was: flat `j0` / `jphi` / `salcbasis` / `symmetry` / `num_atoms`).
  `SCEBasis` stores `structure` / `symmetry` / `salcbasis` only —
  `Cluster` is a construction step, not a stored field.
- **AtomsBase.jl** as a first-class dependency. `SCEBasis` accepts any
  `AtomsBase.AbstractSystem`. Boundary conversion only: `ustrip` units
  at the constructor, internal representation stays `Float64` (see
  design.md "Unitful boundary").
- **Species-based** parameter specification: `lmax` / `cutoff` keyed by
  kind-name `Symbol`. The kind name comes from the per-atom
  `atom.data[:atom_name]` (AtomsBase), or the element symbol when that
  entry is absent — so `:Fe_4a` / `:Fe_8e` are distinct kinds.
  (`ChemicalSpecies.atom_name` is unusable: 4-char cap, ignored by
  `==`.) An element-only key (`Dict(:Fe => 2)`) is sugar that fans out
  to all sublabels of that element.
- **Constructor-level validation**: the `SCEBasis` constructor emits a
  `@warn` on odd `lmax` / `lsum` (time-reversal symmetry removes
  odd-`l` terms, so an odd value is equivalent to `value - 1`).
- **StatsAPI.jl** as a dependency. Implement the response-block-
  independent verbs `fit`, `coef`, `intercept`, `nobs`, `dof`. The
  bare StatsAPI `r2` / `rss` / `residuals` / `predict` are NOT
  implemented — they assume a single response and would force an
  implicit "energy" default.
- **Magesty evaluation & prediction verbs** in `(predictor, data)`
  form, energy and torque always explicit: `r2_energy` / `r2_torque`,
  `rss_energy` / `rss_torque`, `residuals_energy` / `residuals_torque`,
  `rmse_energy` / `rmse_torque`, `predict_energy` / `predict_torque`.
  No `metrics` aggregator.
- **`save` / `load`** (`.xml` only). `SCEBasis` and `SCEModel`
  round-trip through XML; the schema records `tolerance_sym` and
  `isotropy` so a saved basis reloads without the original input
  (`Symmetry` is recomputed from `structure` + `tolerance_sym`, not
  deserialized; `isotropy` is provenance metadata). `SCEFit` is **not**
  covered by Magesty's `save` / `load`: it is a plain struct, so users
  who want to persist a fit use `JLD2` directly. JLD2 is intentionally
  not a Magesty dependency.
- **TOML** kept as a thin convenience layer (`SCEBasis("input.toml")`),
  but the `[regression]` section is removed — fit parameters
  (`lambda`, `torque_weight`, datafile) move to Julia code.
- **Old API removal**: `System`, `SpinCluster`, `build_sce_basis`,
  `build_sce_basis_from_xml`, `fit_sce_model`, `write_xml`,
  `predict_energy`, `Optimizer` (export), `get_j0` / `get_jphi` /
  `get_j0_jphi`, `calc_energy` / `calc_torque` are removed from the
  `Magesty` module.
- **Slicing**: `dataset[idx]` returns a copy `SCEDataset`;
  `view(dataset, idx)` returns a view. `vcat(d1, d2)` requires the
  same `basis` (runtime check).
- **Migration**: `test/examples/*` Julia drivers rewritten to the new
  API; their TOML files lose `[regression]`; the `examples/` directory
  (currently a stale data dump under `examples/fept/`) rebuilt as
  runnable new-API example scripts; `docs/src/`, `SPEC.md`,
  `README.md` updated.

Out of scope (deferred to follow-up specs):

- **Molecule / aperiodic systems** (design-note "molecular support"):
  `SCEBasis` will accept `AtomsBase.isolated_system(...)` in API shape,
  but the SALC construction's reliance on translational symmetry is a
  separate effort.
- **`AtomsIO.jl` as a hard dependency**: documented as the recommended
  way to read CIF / POSCAR / extxyz and used in `examples/`, but not a
  `Project.toml` dependency of the package itself.
- **Cross-validation wrapper estimator** (`CV(...)`) and the `kfold`
  helper: the type hierarchy supports it, implementation is later.
- **`view`-backed `SCEDataset` slicing**: `getindex` returns a copy in
  this spec; `view` (for fold-heavy CV workflows) is deferred to the
  CV follow-up spec, which will settle the view type design.
- **Bayesian estimators / `BayesianSCEFit`**: `SCEFit` stays a concrete
  `struct <: StatsAPI.RegressionModel`. An `AbstractSCEFit` abstract
  type can be introduced non-breakingly when Bayesian regression is
  added (it is on the roadmap, but not this spec).
- **MLJ.jl / StatsModels.jl** integration.
- **Unitful output**: outputs stay `Float64`; users who want units
  wrap them themselves.

## Invariants — what must NOT change

1. **SALC / coupled-basis numerics**: basis function values, signs,
   normalization, `(l, m, site)` ordering. Guarded by `test_Basis.jl`,
   `test_SALCBases_l13_regression.jl`, `test_sphericart_agreement.jl`.
2. **`SCEModel` / `SCEBasis` XML schema**: the XML schema is stable and
   content-faithful, and `save(obj, "x.xml")` is regression-tested
   byte-for-byte against a committed baseline. The baseline is
   regenerated once for the schema introduced in this spec (which adds
   `tolerance_sym` / `isotropy` attributes and is therefore *not*
   byte-identical to the pre-refactor `write_xml` output); thereafter
   any drift from the committed baseline is a regression.
3. **Physics conventions**: unit-vector spin directions, real tesseral
   spherical harmonics, `3 × n_atoms` layout, energy in eV, torque in
   eV.
4. **Numerical fit results**: for the existing `test/examples/*`
   (all `torque_weight ∈ {0, 1}`, `lambda = 0`), `j0` / `jphi` stay
   within rounding error of the current baseline. (The per-sample MSE
   normalization that could change intermediate-weight results already
   landed in PR #3.)
5. **`make test-all` green at every commit** on the topic branch.

## Completion criteria

- [ ] The four types, their constructors, accessors, and `save`/`load`
      cover examples 1-9 in `docs/design-notes/sce-public-api.md`.
- [ ] Old public surface (listed under Scope) removed from `Magesty`.
- [ ] All `test/examples/*/input.toml` have `[regression]` removed;
      their Julia drivers build via `SCEBasis` / `SCEDataset` / `fit`.
- [ ] `make test-all`, `make test-jet`, `make test-aqua` green.
- [ ] `SCEModel` / `SCEBasis` XML output byte-identical to the
      committed (new-schema) baseline, and round-trips through `load`.
- [ ] `docs/src/api.md`, `SPEC.md`, `README.md`, `docs/src/examples.md`
      reflect the new API.
- [ ] `docs/design-notes/sce-public-api.md` and `DESIGN_NOTES.md`
      marked complete.

## Risk

| Risk | Mitigation |
|------|------------|
| AtomsBase API churn (it is pre-1.0) | Pin a compat bound; isolate all AtomsBase calls in the `SCEBasis` constructor and a small adapter file so churn is localized. |
| `SCEModel` / `SCEBasis` XML schema drift | `xml_io.jl` is the single source of the schema; a byte-diff round-trip test on `fept` / `fege` against committed baselines catches unintended drift. Baselines are regenerated deliberately when the schema changes. |
| Intermediate commits leave the package in a non-importable state | Steps 1-5 are add-only; old and new APIs coexist until step 7. `make test-all` runs at every commit. |
| `SCEDataset` design matrix memory blows up for large EMBSET | `X_E`/`X_T` stored unweighted and once; slicing returns views where asked. Document the memory model in design.md. |
| Species-based `lmax` keying confuses users coming from element-only TOML | Provide the element-only fan-out sugar and document both forms; TOML path keeps element-only keys. |
