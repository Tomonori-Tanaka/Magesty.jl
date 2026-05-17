# Style guide

Follows the official Julia style guide and the DFTK style guide.
Guiding principle: **readability beats consistency**.

## Naming conventions

| Target | Rule | Example |
|---|---|---|
| Modules and types | `UpperCamelCase` | `SpinCluster`, `SALCBasis` |
| Functions and variables | `lowercase` / `snake_case` | `calc_energy`, `num_atoms` |
| Constants | `UPPER_SNAKE_CASE` | `MAX_ITER` |
| Internal helpers | leading `_` | `_compute_aux` |
| Mutating functions | trailing `!` | `update_spin!` |

## `for` loops

- Use `=` for numeric **indices**.
- Use `in` for **collection elements**.

```julia
# Index loops -> =
for i = 1:num_atoms
for isym = 1:n_operations

# Element loops -> in
for atom in atoms
for (salc_idx, key_group) in enumerate(salc_list)
```

The distinction makes it obvious at read time whether the loop variable
is an index or a value.

## Argument order

Follow the Julia-official ordering:

```
function -> IO -> mutated args -> types -> non-mutated args -> keyword args
```

## Type annotations

- Annotate arguments of public APIs.
- Annotate the return type as well (`::Float64`, etc.).
- Always brace the `where` clause.

```julia
# Good
function calc_energy(sc::SpinCluster, spin_config::AbstractMatrix{<:Real})::Float64
function foo(x::T) where {T <: AbstractFloat}

# Bad
function foo(x::T) where T <: AbstractFloat
```

## Named tuples

Use the explicit syntax `(; var=val)`:

```julia
# Good
return (; energy=E, forces=F)
(; energy, forces) = result

# Bad
return (energy=E, forces=F)
```

## Miscellaneous

- Indent with 4 spaces.
- Line length up to 92 characters.
- Do not parenthesize `if` / `while` conditions.
- Use `isa` / `<:` for type checks (not `==`).
- Pass `identity` for empty callbacks.
- Avoid type piracy (adding `Base` methods on types defined by other
  packages).

## Hot paths and physics conventions

Guidance for `SVector` / `MVector` usage in hot paths, `@views` /
`@inbounds`, and the physics conventions (unit-vector spin directions,
`3 × n_atoms` layout, real spherical harmonics `Zₗₘ`, energy units) lives
in [`CLAUDE.md`](CLAUDE.md) under "Physics conventions" and
"Performance guidelines". These directly affect numerical results — read
before editing.
