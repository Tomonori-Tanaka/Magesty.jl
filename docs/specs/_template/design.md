# Design: <title>

Status: draft (YYYY-MM-DD)

## Summary

<!-- 1-3 paragraphs describing the chosen design. List alternatives only
     when you need to explain why they were rejected. -->

## Module layout

<!-- Affected files / modules / types and their responsibilities. -->

| Target | Change |
|---|---|
| `src/Foo.jl` | <!-- e.g., change the signature of `bar` --> |

## API

<!-- Public / internal API additions, changes, deletions. Note type
     annotations and docstring style. -->

```julia
# Example
function new_fn(x::Foo; opt::Bool = false)::Bar
```

## Types and conventions

<!-- Impact on physics conventions, units, numerical conventions. List
     any new invariants. -->

## Impact on linked sites

<!-- Which "Linked sites" in CLAUDE.md does this touch?
     Spherical-harmonics convention / SCE XML / SALC <-> Fitting / other. -->

- [ ] Spherical-harmonics convention (`TesseralHarmonics`):
- [ ] SCE coefficient XML (`save` / `load`):
- [ ] `Fitting` <-> `SALCBasis`:
- [ ] `.claude/agents/` references:
- [ ] `SPEC.md` / `docs/src/api.md` updates:

## Test strategy

<!-- New tests, modified tests, benchmarks. -->

## Risks and open items

<!-- Anything that may change numerical results, unresolved decisions,
     deferred alternatives. -->
