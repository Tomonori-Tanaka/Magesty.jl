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

### Public-API constructors and verbs

On top of the language-level ordering above, Magesty's public
constructors and verbs follow a consistent shape. Read it as "what the
result is built from" first, "how it should be built" next, and
"behavior knobs" last:

```
(data inputs) -> interaction spec -> symmetry / tolerance options
              -> cosmetic options -> verbosity
```

Examples (`SCEBasis` has four input paths; all four keep this order):

- `SCEBasis(system_spec, interaction, options; verbosity=true)`
- `SCEBasis(system; interaction, name="system", tolerance_sym=1e-5,
  isotropy=false, verbosity=true)` (AtomsBase path)
- `SCEBasis(; lattice, kd, kd_list, positions, periodicity, interaction,
  name="system", tolerance_sym=1e-5, isotropy=false, verbosity=true)`
  (kwargs path)
- `SCEDataset(basis, spinconfigs)` — basis already carries the
  symmetry / interaction settings.
- `SCEDataset(system, spinconfigs; interaction, name, tolerance_sym,
  isotropy, verbosity)` — same trailing kwargs as `SCEBasis`.
- `SCEModel(basis, j0, jphi)` — basis first, then intercept, then
  coefficients (scalar before vector).

Verbs follow the StatsBase / StatsAPI shape:

- `fit(SCEFit, dataset, estimator; torque_weight=1.0)` — the only
  tunable lives in kwargs.
- `predict_energy(predictor, data)`, `r2_energy(predictor, data)`, …
  always take the predictor first and the evaluation data second.

When adding a new input path or overload, keep this order. Behavior
flags (booleans, `verbosity`) belong at the end as keyword arguments
with sensible defaults so that the happy path stays positional.

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

## Docstrings

Public-API docstrings follow the standard Julia format. List positional
arguments under `# Arguments` and keyword arguments under a separate
`# Keyword arguments` section, mirroring the `;` in the signature; then
`# Returns`, and `# Examples`. A `# Throws` section is optional. Omit
`# Keyword arguments` when the function takes no keyword arguments.

```julia
"""
	f(x; scale = 1.0) -> Float64

# Arguments
- `x::Real`: the input value.

# Keyword arguments
- `scale::Real = 1.0`: multiplier applied to `x`.

# Returns
- `Float64`: `scale * x`.
"""
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
