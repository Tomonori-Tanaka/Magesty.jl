# Design: write_xml public API cleanup

Status: draft (2026-05-13)
Spec: [`requirements.md`](requirements.md)

## Before

Today's public API exposes four `write_xml` methods:

```julia
# src/Magesty.jl
write_xml(sc::SpinCluster, filename = "jphi.xml"; write_jphi = true)              # (1)
write_xml(system::System, filename = "jphi.xml")                                  # (2)
write_xml(structure, symmetry, basis_set, optimize, filename = "jphi.xml";        # (3)
          write_jphi = true)

# src/utils/xml_io.jl  (XMLIO module, not exported but accessible)
XMLIO.write_xml(structure, symmetry, basis_set, optimize, filename;               # (3-impl)
               write_jphi = true)
XMLIO.write_xml(structure, symmetry, basis_set, filename)                         # (4)
```

In `Magesty.jl`, method (1) delegates to (3); method (2) delegates to
the four-arg `XMLIO.write_xml` (signature 4); method (3) delegates to
the five-arg `XMLIO.write_xml` (signature 3-impl).

Method (3) is the noisy one — it exposes the internal
(`structure`, `symmetry`, `basis_set`, `optimize`) tuple to users and
gives `write_xml` four public arities at once.

## After

```julia
# src/Magesty.jl
write_xml(system::System, filename = "jphi.xml")                                  # (2) unchanged
write_xml(sc::SpinCluster, filename = "jphi.xml"; write_jphi = true)              # (1) unchanged

# src/utils/xml_io.jl  (still XMLIO, still internal)
XMLIO.write_xml(structure, symmetry, basis_set, optimize, filename;               # (3-impl) unchanged
               write_jphi = true)
XMLIO.write_xml(structure, symmetry, basis_set, filename)                         # (4) unchanged
```

The public four-arg method `write_xml(structure, symmetry,
basis_set, optimize, filename; ...)` is **deleted from
`src/Magesty.jl`**. The XMLIO-side implementations remain — they are
internal helpers and `Magesty.write_xml(sc::SpinCluster, ...)` keeps
calling them.

## Implementation steps

1. **Delete** the four-arg public method in `src/Magesty.jl`. The
   exact block is the function definition that currently sits after
   the two type-dispatched methods (around L540).
2. **Keep** `src/Magesty.jl` `write_xml(sc::SpinCluster, ...)`
   forwarding to `XMLIO.write_xml(sc.structure, sc.symmetry,
   sc.basisset, sc.optimize, filename; write_jphi = write_jphi)`.
   This is the only path that needs the five-arg XMLIO method.
3. **Keep** `src/Magesty.jl` `write_xml(system::System, filename)`
   forwarding to `XMLIO.write_xml(system.structure, system.symmetry,
   system.basisset, filename)`. This uses the four-arg XMLIO method
   (no `Optimizer`).
4. **Update** the `# Examples` block in the surviving docstring(s) so
   no example shows the deleted form.
5. **Update** `docs/src/examples.md` L186 from
   `write_xml(system.structure, system.symmetry, system.basisset,
   optimizer, "new_results.xml")` to
   ```julia
   sclus = SpinCluster(system.structure, system.symmetry,
                       system.cluster, system.basisset, optimizer)
   write_xml(sclus, "new_results.xml")
   ```
   The example is illustrating "I have already fitted an Optimizer
   from a System loaded via `build_sce_basis_from_xml`" — wrapping
   into a `SpinCluster` is the natural way to express that.

## Docstring policy

After this change there are exactly two documented public arities:

- `write_xml(system::System, filename::AbstractString="jphi.xml")`
- `write_xml(sc::SpinCluster, filename::AbstractString="jphi.xml";
            write_jphi::Bool=true)`

The docstring on the `SpinCluster` method should explain both forms
under one heading (similar to `Base.write` having one docstring with
multiple signatures). The internal `XMLIO.write_xml` helpers do not
need user-facing docstrings; they keep their existing inline docs.

## Open questions

- ~~Should the four-arg public form get a deprecation warning instead
  of straight removal?~~ — No. Only one in-repo caller, no published
  release yet (`0.1.0-DEV`), and a `MethodError` is a clearer signal
  than a warning during early development.
- ~~Should default filename change?~~ — Out of scope for this spec
  (would require coordinating with the user-facing API spec that
  might rename `write_xml` to `save` and pick a fresh default).

## Linked sections

- `SPEC.md` line 87 — already lists `write_xml(sc, filename)` only;
  no change needed.
- `docs/src/api.md` line 48 — references `write_xml` for Documenter
  to render its docstrings; will pick up the updated docstring
  automatically.
- `docs/src/tutorial.md` and `docs/src/examples.md` use the
  `SpinCluster` and `System` forms already; the only edit is
  `examples.md` L186.

## Test plan

- `make test-unit` + `make test-integration` cover both surviving
  overloads via the existing example tests.
- `make test-jet` confirms no new method-error inference holes.
- `make test-aqua` confirms ambiguity / export cleanliness.
- A manual byte-level diff of one XML artifact before vs after
  serves as a quick numerical-equivalence check; the integration
  examples already check XML round-trip via
  `build_sce_basis_from_xml` for the `fept` case.
