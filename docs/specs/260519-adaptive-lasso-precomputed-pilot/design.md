# Design: AdaptiveLasso precomputed-pilot path

Status: draft (2026-05-19)

## Summary

Introduce a tiny adapter estimator `PrecomputedPilot <: AbstractEstimator`
that carries a fixed coefficient vector and returns it from
`solve_coefficients` unchanged. `AdaptiveLasso` already accepts any
`AbstractEstimator` as its `pilot` field, so the existing kernel needs
zero changes -- the precomputed path is just a different choice of
`pilot`.

Two convenience constructors,
`AdaptiveLasso(model::SCEModel; ...)` and
`AdaptiveLasso(fit::SCEFit; ...)`, wrap `coef(model)` / `coef(fit)` in a
`PrecomputedPilot` and forward all remaining kwargs to the existing
keyword constructor. This keeps the user-facing call sites symmetric
with the current `AdaptiveLasso(; pilot = OLS(), lambda = ...)` form.

Alternatives considered and rejected:

- Widening `AdaptiveLasso.pilot` to `Union{AbstractEstimator,
  AbstractVector{<:Real}}`. This breaks the "pilot is an estimator
  recipe" abstraction and forces a runtime branch inside
  `solve_coefficients`. The `PrecomputedPilot` wrapper keeps the
  `AbstractEstimator` invariant.
- Adding a new field `pilot_coef::Union{Nothing, Vector{Float64}}` to
  `AdaptiveLasso`. Same downsides plus a wider struct on the hot path.

## Module layout

`SCEFit` and `SCEModel` are defined in the parent `Magesty` module
*after* `include("Fitting.jl")`, and there is no forward-declaration in
Julia. Therefore:

- `PrecomputedPilot` (struct, inner ctor, `solve_coefficients` method)
  lives entirely inside `src/Fitting.jl`, alongside the other
  estimators. It is exported from `Fitting` and re-exported from
  `Magesty`.
- The two convenience constructors `AdaptiveLasso(::SCEFit; ...)` and
  `AdaptiveLasso(::SCEModel; ...)` live in `src/Magesty.jl`, placed
  after `SCEFit` / `SCEModel` and after `coef(::SCEFit)` /
  `coef(::SCEModel)` are defined (currently `Magesty.jl:624-625`). They
  add methods to the existing `AdaptiveLasso` constructor that
  `Fitting` exports.

| Target | Change |
|---|---|
| `src/Fitting.jl` | Add `PrecomputedPilot` struct (with copying inner ctor), `solve_coefficients(::PrecomputedPilot, X, y)` method; export `PrecomputedPilot`. No new `AdaptiveLasso` constructors here. |
| `src/Magesty.jl` | Re-export `PrecomputedPilot`. Add the two convenience `AdaptiveLasso` constructors after the `coef(::SCEFit)` / `coef(::SCEModel)` methods. |
| `test/component/test_fitting_estimators.jl` | New `@testset "PrecomputedPilot"` (parity, length-mismatch, ctor smoke). |
| `test/integration/febcc_2x2x2_pm/test.jl` | Extend the existing AdaptiveLasso smoke testset with an `AdaptiveLasso(fit_ols; lambda = ...)` reuse case. |
| `docs/src/api.md` | Add `PrecomputedPilot` to the `## Estimators` `@docs` block. |
| `SPEC.md` | Update the `# Estimators` API line. |
| `CHANGELOG.md` | New `[Unreleased] ### Added` entry. |
| `docs/design-notes/<lasso-family note>` | Note that `PrecomputedPilot` is the iterative-variant primitive. |

## API

```julia
# In src/Fitting.jl
"""
    PrecomputedPilot(beta::AbstractVector{<:Real})

Estimator adapter that returns a fixed coefficient vector from
`solve_coefficients`, ignoring the supplied `(X, y)` except for a length
check against `size(X, 2)`. Designed for `AdaptiveLasso.pilot`: lets the
adaptive call reuse coefficients from a previous fit instead of running
a fresh pilot.

The vector is copied at construction so later mutation of the caller's
storage does not leak into `PrecomputedPilot.beta`.
"""
struct PrecomputedPilot <: AbstractEstimator
    beta::Vector{Float64}
    PrecomputedPilot(b::AbstractVector{<:Real}) = new(Vector{Float64}(b))
end

function solve_coefficients(
    e::PrecomputedPilot,
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
)::Vector{Float64}
    length(e.beta) == size(X, 2) || throw(DimensionMismatch(
        "PrecomputedPilot coefficient length $(length(e.beta)) does not " *
        "match design-matrix column count $(size(X, 2)); the pilot was " *
        "likely fit on a different SCEBasis."))
    return e.beta
end
```

```julia
# In src/Magesty.jl, AFTER coef(::SCEFit) / coef(::SCEModel) are defined.
# Forward kwargs unchanged.
AdaptiveLasso(model::SCEModel; kwargs...) =
    AdaptiveLasso(; pilot = PrecomputedPilot(coef(model)), kwargs...)

AdaptiveLasso(fit::SCEFit; kwargs...) =
    AdaptiveLasso(; pilot = PrecomputedPilot(coef(fit)), kwargs...)
```

Notes on the convenience constructors:

- They forward `kwargs...` (including `lambda`, which is required).
  Passing `pilot` again via these kwargs raises Julia's standard
  `ErrorException("duplicate keyword argument: pilot")` because the
  underlying keyword ctor would receive `pilot` twice; that is the
  intended error (the precomputed pilot is the whole point).
- They do not infer `lambda` from the prior fit -- `lambda` belongs to
  the *adaptive* call, not the pilot.

## Types and conventions

- `PrecomputedPilot.beta` is stored as `Vector{Float64}` (eager copy of
  the input, enforced by the inner constructor). This matches the type
  that `solve_coefficients(::OLS, ...)` and friends return and avoids
  holding a reference to caller-mutable storage. The field is named
  `beta` rather than `coef` to avoid visual collision with the
  `StatsAPI.coef` function that Magesty extends.
- No physics-convention surface area is touched. SALC ordering is
  unchanged.
- The length check in `solve_coefficients` is the only safety net for
  basis mismatch; the docstring on `AdaptiveLasso(::SCEModel; ...)`
  must say "the model must have been fit on the same `SCEBasis`."

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): unaffected.
- [ ] SCE coefficient XML (`save` / `load`): unaffected. `PrecomputedPilot`
      is an estimator recipe, never serialized.
- [ ] `Fitting` <-> `SALCBasis`: unaffected at the design-matrix level;
      the new code path runs strictly inside `solve_coefficients`.
- [ ] `.claude/agents/` references: no module / Makefile change, but
      sweep for stale "AdaptiveLasso has one pilot field" claims.
- [ ] `SPEC.md` / `docs/src/api.md` updates: add `PrecomputedPilot` and
      the convenience constructors.

## Test strategy

Component tests in `test/component/test_fitting_estimators.jl`:

1. **Pilot-reuse parity (bit-exact).** Compute `beta_ols =
   solve_coefficients(OLS(), X, y)` directly on the fixture. Then:
   - `b_a = solve_coefficients(AdaptiveLasso(pilot = OLS(), lambda = L,
     gamma = 1), X, y)`
   - `b_b = solve_coefficients(AdaptiveLasso(pilot =
     PrecomputedPilot(beta_ols), lambda = L, gamma = 1), X, y)`
   - Assert `b_a == b_b` (no `atol` -- they are computing the same `pf`
     and calling the same GLMNet path with the same inputs).
2. **Length-mismatch error.** Call `solve_coefficients(PrecomputedPilot(zeros(k)),
   X_with_p_cols, y)` where `k != p`; expect `DimensionMismatch`.
3. **Constructor smoke.** Build a tiny `SCEModel` (or skip and fit a real
   one in the integration test below) and check that the constructor
   returns an `AdaptiveLasso` whose `pilot isa PrecomputedPilot` and
   carries the model's coefficients.

Integration test (`test/integration/febcc_2x2x2_pm/test.jl`): inside the
existing AdaptiveLasso testset, also exercise the reuse path:
`AdaptiveLasso(fitted_ols; lambda = 1e-4)` runs end-to-end and the
resulting `coef` is finite (no parity claim against the from-scratch
path; the from-scratch pilot inside the existing test uses `OLS()` too,
so equality holds but we keep the assertion at the component level).

## Risks and open items

- **Basis-mismatch silent corruption.** If the user passes an `SCEFit`
  fit on a different `SCEBasis` (different SALC ordering or count) the
  length check catches the most common case (different `num_salcs`),
  but a same-length-different-meaning case can still slip through. We
  rely on documentation here; a tighter check would require carrying
  the `SCEBasis` identity around, which is out of scope.
- **Method-ambiguity risk on the convenience constructors.**
  `AdaptiveLasso(model::SCEModel; ...)` and
  `AdaptiveLasso(fit::SCEFit; ...)` are positional-arg constructors; the
  existing keyword-only ctor `AdaptiveLasso(; pilot = OLS(), lambda, ...)`
  has no positional argument and so cannot collide. Verified by reading
  the current Fitting.jl ctor; no ambiguity warning expected.
- **`Vector{Float64}` eager copy.** Deliberate, enforced by the inner
  constructor so the default no-copy path is unreachable. The iterative
  variant will call this constructor inside its inner loop, so the
  per-call cost is one `Vector{Float64}` of length `num_salcs` --
  negligible compared to a GLMNet call.
