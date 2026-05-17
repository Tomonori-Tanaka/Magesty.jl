# Requirements: Typed input spec + wildcard species + cluster definition

Status: draft (2026-05-16)
Owner: T. Tanaka
Branch: `refactor/typed-input-spec` (cut when this spec is agreed)
Design note: `docs/design-notes/config4system-role.md` (case B deep-dive)

## Goal

Replace the god-config `Config4System` with three focused typed value
objects (`SystemSpec`, `InteractionSpec`, `SymmetryOptions`) that share a
single dict-parsing layer for TOML / Dict inputs, while AtomsBase and
kwargs paths construct typed values directly without a dict round-trip.

At the same time, introduce ALAMODE-style wildcard species notation
(`"*-*"`, `"Fe-*"`, `"Fe-Fe"`) for `bodyn_cutoff` and `body1_lmax`,
resolved by specificity (more specific overrides less specific).

Finally, document the existing interaction-cluster definition (all
pairwise distances within their pair-specific cutoffs, applied uniformly
for every body order n ≥ 2) in the user-facing docs alongside the new
cutoff syntax.

## Motivation

1. **god-config**: `Config4System` packs four orthogonal concerns
   (system / symmetry / interaction / structure) into one struct.
   Downstream consumers receive fields they do not use (`Symmetry` only
   reads `tolerance_sym`; `Cluster` only reads `nbody` / `bodyn_cutoff`).
2. **dict round-trip**: AtomsBase and kwargs paths build a TOML-shape
   dict purely to feed `Config4System`. ~150 lines of dict serializers
   that produce nothing of value.
3. **Wildcards need parser logic**: the new `*-*` / `Fe-*` override
   resolution must live in the dict parser. Adding it on top of the
   current god-config grows the dict layer further; doing it together
   with the case-B move keeps the parser slim and well-scoped.

## Scope

In scope:

- New module `src/InputSpecs.jl` with three typed value objects.
- `parse_toml_inputs(dict::AbstractDict) -> (SystemSpec, InteractionSpec, SymmetryOptions)`
  as the single TOML / Dict entry point, with wildcard expansion built in.
- AtomsBase adapter rewritten to construct typed values directly.
- kwargs `SCEBasis(; lattice, kd, ...)` path rewritten to construct
  typed values directly.
- Downstream modules (`Structure`, `Symmetry`, `Cluster`, `SALCBasis`,
  `Fitting`, `XMLIO`) updated to consume the relevant spec(s) instead of
  `Config4System`.
- `Config4System` removed from the package (BREAKING for any caller that
  named it explicitly; it was never in `Magesty`'s export list).
- Tests for `InputSpecs` (typed-value validation, wildcard expansion,
  parse equivalence with current behaviour, cluster cutoff criterion).
- User-facing docs updated: `docs/src/input_keys.md` gains wildcard
  syntax description and a section on interaction cluster definition;
  `SPEC.md` updated to reflect the new types.

Out of scope:

- Wildcards for `bodyn_lsum` (deferred; covered by Future Work section
  in design.md).
- New input formats (CIF, JSON, ...). The architecture supports adding
  them post-refactor.
- Performance tuning of hot paths (`Fitting`, `SALCBases`).
- TOML schema changes beyond accepting `*` in existing pair-keyed
  tables. Existing TOML inputs must keep working unchanged.

## Invariants (must be preserved)

- **Numerical regression**: every TOML in `test/examples/` must produce
  an `SCEBasis` whose XML fingerprint matches the pre-refactor output
  bit-for-bit (XML I/O round-trip determinism is the cross-check).
- **`InteractionSpec` post-construction**: all `(species_i, species_j)`
  pairs (with `i ≤ j`) are concrete and symmetric; all body1 species
  have concrete `lmax`. No wildcards survive into the typed value.
- **Cluster definition**: existing `is_within_cutoff` behaviour
  (all pairwise distances ≤ pair cutoff for every n-body cluster, with
  the same rule for every n ≥ 2) is unchanged. Only documented, not
  modified.
- **SALCBasis (l, m, site) ordering**: untouched. Fitting / XML I/O rely
  on it (see CLAUDE.md "Linked sites").

## Out-of-scope invariants (explicitly may change)

- `Config4System` itself: removed. External callers naming it directly
  break (no shim).
- Internal validation entry points (`parse_interaction_body1` etc.):
  moved into the inner constructors of the new specs. Their names go
  away.

## Completion criteria

- All `make test-unit`, `make test-integration`, `make test-jet`,
  `make test-aqua` targets pass on the topic branch.
- New `test/component_test/test_input_specs.jl` covers each spec rule
  (typed-value validation, wildcard expansion, cluster cutoff criterion)
  with tests derived from this spec, not from the implementation.
- `docs/src/input_keys.md` and `SPEC.md` reflect the new types and
  wildcard syntax. `DESIGN_NOTES.md` index updated.
- `Config4System` no longer appears in `src/`, `test/`, or `docs/src/`.
- Existing TOML examples in `test/examples/` produce identical XML
  output (fingerprint check, integration test).

## Non-goals

- Backward compatibility shim for `Config4System` (the type goes away).
- Multi-error aggregation in `parse_toml_inputs` (first-error throw is
  acceptable; future work).
