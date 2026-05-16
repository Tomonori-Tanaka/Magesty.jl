# Tasklist: Typed input spec + wildcard species + cluster definition

Status: draft (2026-05-16)
Branch: `refactor/typed-input-spec`

Milestones are coarse. Day-to-day sub-tasks are tracked via in-session
TaskCreate, not committed here. Each milestone ends with `make test-all`
green and a `code-reviewer` agent pass.

## M0. Spec agreement (this folder)

- [ ] User reviews requirements.md / design.md / tasklist.md.
- [ ] Open questions in design.md resolved.
- [ ] Branch `refactor/typed-input-spec` cut from `main`.

## M1. Typed value structs

- [ ] Create `src/InputSpecs.jl` (module + three structs +
      inner-constructor validation, no wildcard logic yet).
- [ ] Move species-set / lattice-shape / x_fractional validation logic
      from `ConfigParser` into the appropriate inner constructor.
- [ ] New tests: `test/component_test/test_input_specs.jl` covering
      typed-value validation rules (error cases derived from spec, not
      from impl).
- [ ] `make test-unit` green.

## M2. Dict parser + wildcard expansion

- [ ] Add `parse_toml_inputs(dict)` to `InputSpecs`.
- [ ] Implement `expand_pair_table` / `expand_species_table` with
      specificity-based resolution. Errors: `WildcardConflictError`,
      `DuplicatePairError`, `MissingPairError`, `UnknownSpeciesError`.
- [ ] Wire `SCEBasis(toml_path)` / `SCEBasis(dict)` through
      `parse_toml_inputs` (still bridged to `Config4System` downstream).
- [ ] Wildcard expansion tests (spec-driven; see test catalogue below).
- [ ] Equivalence test: every TOML in `test/examples/` produces an XML
      output identical to the pre-refactor baseline.
- [ ] `make test-all` green.

## M3. AtomsBase / kwargs adapter rewrite

- [ ] Rename `system_to_input_dict` → `atomsbase_to_specs` (now returns
      the three-tuple).
- [ ] Rename `kwargs_to_input_dict` → `kwargs_to_specs` likewise.
- [ ] Delete `_interaction_section` and other dict-building helpers.
- [ ] AtomsBase / kwargs paths now reuse `expand_pair_table_per_body`
      when wildcard dicts are provided.
- [ ] Adapter tests adjusted to assert on typed values (not dicts).
- [ ] `make test-all` green.

## M4. Downstream migration

For each step: switch the constructor signature, update every caller,
keep tests green.

- [ ] `Structure(system::SystemSpec)`.
- [ ] `Symmetry(structure, options::SymmetryOptions)`.
- [ ] `Cluster(structure, symmetry, interaction::InteractionSpec; verbosity)`.
- [ ] `SALCBasis(structure, symmetry, cluster, interaction, options)`.
- [ ] `Fitting` / `XMLIO`: only update internal field accesses if any.
- [ ] `SCEBasis` internal constructor calls the new chain directly
      (no `Config4System` bridge).
- [ ] `make test-all`, `make test-jet`, `make test-aqua` all green.

## M5. Delete `Config4System`

- [ ] Confirm via `grep -rn Config4System src/ test/ tools/` that no
      references remain (except possibly `tools/personal/` — migrate or
      flag).
- [ ] Delete `src/ConfigParser.jl` and `test/component_test/test_ConfigParser.jl`.
- [ ] Remove `Config4System` export from any remaining module.
- [ ] `make test-all` green.

## M6. Documentation

- [ ] `docs/src/input_keys.md`:
  - Add wildcard syntax for `cutoff.<pair>` and `body1.lmax.<species>`.
  - Add "Interaction cluster definition" section. Wording must apply
    uniformly to every n ≥ 2; do not phrase as if only one specific body
    order uses the all-pairs criterion.
- [ ] `SPEC.md`: replace `Config4System` references with the three
      typed values; cross-reference the cluster definition section.
- [ ] `DESIGN_NOTES.md` index: update
      `config4system-role.md` Status to "spec 化済み (260516-typed-input-spec)".
- [ ] `docs/design-notes/config4system-role.md`: add a note pointing to
      this spec as the realization (already done in Step 0 of the plan
      — re-verify).

## Test catalogue (M1 + M2)

`test/component_test/test_input_specs.jl` covers, at minimum:

**SystemSpec validation**
- lattice not 3×3 → error.
- `length(kd_int_list) != num_atoms` → error.
- `kd_int_list[i] > length(kd_name)` → error.
- `x_fractional` outside `[0,1)` under periodic → error.
- duplicate species in `kd_name` → error.

**InteractionSpec validation**
- `nbody = 0` → error.
- `body1_lmax` length != `length(kd_name)` → error.
- negative `lmax` or `lsum` → error.
- asymmetric `bodyn_cutoff` (i,j ≠ j,i) → error.
- missing pair entry → error.

**SymmetryOptions**
- `tolerance ≤ 0` → error.

**Wildcard expansion (`parse_toml_inputs`)**
- 3-species system with only `"*-*" = 8.0` for body 2 → every pair
  expands to 8.0.
- `"*-*" = 8.0` + `"Fe-Fe" = 12.0` → Fe-Fe gets 12.0, others 8.0.
- `"*-*" = 8.0` + `"Fe-*" = 10.0` → Fe-X (X ≠ Fe) and Fe-Fe get 10.0, X-Y get 8.0.
- `"Fe-*" = 10.0` + `"*-Ni" = 11.0` covering Fe-Ni → `WildcardConflictError`.
- `"Fe-Ni" = 10.0` + `"Ni-Fe" = 12.0` → `DuplicatePairError`.
- `"Fe-Xx" = 10.0` with Xx unknown → `UnknownSpeciesError`.
- `body1.lmax`: `"*" = 2` + `"Fe" = 4` → Fe = 4, others 2.

**Cluster cutoff criterion (reaffirmation, not new behaviour)**
- 3-atom candidate, 1 pair distance > pair cutoff → cluster rejected.
- 3-atom candidate, all 3 pair distances ≤ pair cutoffs → cluster accepted.
- Same pattern at n = 4 if a fixture allows: 1 pair out → rejected;
  all pairs in → accepted. Important so the test suite documents that
  the rule is uniform across n, not specific to n = 3.

**Cross-equivalence**
- For each TOML in `test/examples/`: build via
  `parse_toml_inputs` and via the pre-refactor pipeline (during M2 the
  bridge still allows both); resulting XML output is bit-identical.

## Verification gate before PR

- [ ] `make test-all` green.
- [ ] `make test-jet` clean.
- [ ] `make test-aqua` clean.
- [ ] `make test-integration` against full `test/examples/` set.
- [ ] `code-reviewer` agent run on the full diff (commit-by-commit if large).
- [ ] No `Config4System` left in `src/`, `test/`, `docs/src/`, or `SPEC.md`.

## Out of this spec

See design.md §10. Tracked in `docs/design-notes/` only if and when
demand appears.
