# Requirements: write_xml public API cleanup

Status: draft (2026-05-13)
Owner: T. Tanaka

## Goal

Reduce the public surface of `write_xml` from four overlapping
signatures to two, dispatched on a single top-level argument
(`System` or `SpinCluster`). The four-argument
`write_xml(structure, symmetry, basis_set, optimize, ...)` form is
removed from the public API; it stays as an internal helper that the
`SpinCluster` overload delegates to.

## Scope

In scope:

- `src/Magesty.jl` — top-level `write_xml` methods (delete the public
  four-argument method; keep the two type-dispatched methods).
- `src/utils/xml_io.jl` — confirm the underlying `XMLIO.write_xml`
  implementation methods stay as internal helpers (no public
  semantics change).
- `docs/src/examples.md` line 186 — the single in-repo caller of the
  public four-argument form. Migrate to the `SpinCluster` form.
- The exported names in `src/Magesty.jl` `export write_xml` — no
  rename, only the signature surface narrows.

Out of scope:

- Renaming `write_xml` to `save` (that lives in the user-facing API
  spec).
- Filename-as-kwarg style (`write_xml(target; filename=...)`) — Julia
  ecosystem convention keeps filename positional; this spec follows
  that convention.
- Round-trip semantics or XML schema changes (already isolated in
  `xml_io.jl` per the R7 refactor).
- Documenter docs auto-generation other than the one example edit.

## Invariants

1. **XML output is byte-identical** before and after for any input
   that goes through one of the surviving overloads. The four
   non-`Optimizer` writer methods produce the same bytes; the new
   internal helper that backs them must too.
2. **`SpinCluster` write semantics unchanged**: full
   structure/symmetry/basis/Optimizer + JPhi (gated by
   `write_jphi=true`).
3. **`System` write semantics unchanged**: structure/symmetry/basis,
   no JPhi.
4. **Default filename unchanged** at `"jphi.xml"` for both surviving
   public methods (revisit only as part of a separate spec).
5. **Round-trip** with `build_sce_basis_from_xml` keeps passing for
   the integration examples.

## Non-goals

- No new positional or kwarg added (e.g., `compress=`, `format=`).
- No performance optimization.
- No deprecation period for the four-argument form. The single
  in-repo caller is updated in the same PR; any out-of-repo caller
  will see a `MethodError` and can re-route via `SpinCluster(...)`
  or by calling the (internal, unsupported) `XMLIO.write_xml`
  directly.

## Public API after

```julia
write_xml(system::System, filename::AbstractString = "jphi.xml")
write_xml(sc::SpinCluster, filename::AbstractString = "jphi.xml";
          write_jphi::Bool = true)
```

Both exported via `Magesty`. No `write_xml(structure, symmetry,
basis_set, optimize, ...)` exposed.

## Completion criteria

- [ ] `src/Magesty.jl` exports only the two surviving methods; the
      four-argument public overload is removed.
- [ ] `src/utils/xml_io.jl` internal helpers are unchanged in
      behavior; the `SpinCluster` overload still calls them with the
      same arguments.
- [ ] `docs/src/examples.md` L186 uses the `SpinCluster` form.
- [ ] `make test-all`, `make test-jet`, `make test-aqua` green.
- [ ] No XML byte-level diff vs the pre-refactor output on the
      integration examples (`dimer`, `chain`, `fept`, `fege`,
      `square_lattice`, `febcc_pm`, `2d_fcc`).
- [ ] `DESIGN_NOTES.md` R2 entry marked complete.

## Risk

| Risk | Mitigation |
|------|------------|
| External callers depending on the four-arg form break silently | Removal is explicit and the only in-repo caller is rerouted in the same PR; external callers see a `MethodError`. Document the migration in the commit message and the DESIGN_NOTES R2 entry. |
| Behavior drift between System and SpinCluster overloads after refactor | Run `make test-integration`; the example tests round-trip the XML for both code paths. |
