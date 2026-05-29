# Design: Solver unification (Cholesky) and design-matrix memory logging

Status: draft (2026-05-26)

## Summary

The OLS, Ridge, and AdaptiveRidge solvers all reduce to the same algebraic
problem: solve `(X'X + Λ) β = X'y` with `Λ = 0` (OLS), `Λ = λI` (Ridge), or
`Λ = λ · Diagonal(w)` (AdaptiveRidge). We route all three through
`cholesky(Symmetric(A)) \ b`, which is the fastest stable solver for
symmetric positive-definite systems and the most memory-frugal (the only
matrix that needs to live is the `num_salcs × num_salcs` Gram matrix).

For OLS, the new failure mode is `PosDefException` when `X'X` is not
strictly positive definite — i.e. when `X` is rank-deficient or numerically
collinear. We catch this and rethrow `ArgumentError` with a message that
explicitly names `Ridge(lambda = ...)` as the remedy. This is by design:
the previous QR + pivoted-QR fallback gave a numerical answer even on
rank-deficient input, but in the SCE context that answer is physically
meaningless (the SCE coefficients are not identifiable). Forcing the user
to declare a regularizer is the correct behavior.

The memory log is two `println`s under the existing `verbosity = true`
gate inside the design-matrix builders. Pre-construction the size is
computed analytically; post-construction we report `Base.summarysize`.

## Module layout

| Target | Change |
|---|---|
| `src/Fitting.jl` (`solve_coefficients(::OLS, ...)`) | Replace QR path with Cholesky on `X'X`; rethrow `PosDefException` as `ArgumentError`. |
| `src/Fitting.jl` (`solve_coefficients(::Ridge, ...)`) | Replace `MultivariateStats.ridge(...)` (and the `λ ≈ 0` `X \ y` fast path) with explicit `cholesky(Symmetric(X'X + λI))`; for `λ ≈ 0` delegate to the OLS Cholesky path so `PosDefException` is surfaced consistently. |
| `src/Fitting.jl` (`solve_coefficients(::AdaptiveRidge, ...)`) | Use `cholesky` for the iteration-zero ridge and the reweighted step. The reweighted-step matrix is SPD by construction (`λ > 0`, `w > 0`), so no `PosDefException` catch is needed there; iteration zero also has `λ > 0` (since we check `λ ≈ 0` upfront and delegate to OLS), so the catch is only on the OLS path. The `λ ≈ 0` fast path delegates to the OLS Cholesky path. |
| `src/Fitting.jl` (`build_design_matrix_energy`) | Add pre/post memory log under `verbosity = true`. |
| `src/Fitting.jl` (`build_design_matrix_torque`) | Same. |
| `src/Fitting.jl` (`OLS` / `Ridge` docstrings) | Reflect the Cholesky solver and the new error behavior. |
| `src/Fitting.jl` (top of file) | Drop `using MultivariateStats`. |
| `Project.toml` | Remove the `MultivariateStats` dependency (`[deps]` + `[compat]`). |
| `test/component/test_fitting_estimators.jl` | Add cross-check and PosDef tests alongside the existing OLS / Ridge / GLMNet tests. (If the file structure suggests a separate `test_fitting_solver.jl` makes the new tests easier to find, create it and register in `test/runtests.jl`.) |

## API

```julia
# Public API surface unchanged. The internal solve becomes:

function solve_coefficients(::OLS, X, y)::Vector{Float64}
    Xf = Matrix{Float64}(X)
    yf = Vector{Float64}(y)
    XtX = Symmetric(Xf' * Xf)
    Xty = Xf' * yf
    try
        return cholesky(XtX) \ Xty
    catch e
        e isa PosDefException || rethrow()
        throw(ArgumentError(
            "OLS solve failed: the normal-equation matrix X'X is not " *
            "positive definite, which means the design matrix has " *
            "linearly dependent (or numerically near-dependent) columns. " *
            "Typical causes: (1) num_salcs >= num_spinconfigs at " *
            "torque_weight = 0 (no torque rows to add rank); (2) redundant " *
            "SALCs in the basis — this can fire even with torque_weight > 0 " *
            "if the basis itself contains duplicate or symmetry-degenerate " *
            "columns. Use Ridge(lambda = ε) with a small ε to regularize, " *
            "or reduce the basis size."
        ))
    end
end

function solve_coefficients(e::Ridge, X, y)::Vector{Float64}
    if e.lambda ≈ 0.0
        return solve_coefficients(OLS(), X, y)
    end
    Xf = Matrix{Float64}(X)
    yf = Vector{Float64}(y)
    XtX = Xf' * Xf
    Xty = Xf' * yf
    # X'X + λI is strictly SPD for λ > 0 — Cholesky cannot fail here.
    return cholesky(Symmetric(XtX + e.lambda * I)) \ Xty
end

function solve_coefficients(e::AdaptiveRidge, X, y)::Vector{Float64}
    if e.lambda ≈ 0.0
        return solve_coefficients(OLS(), X, y)
    end
    Xf = Matrix{Float64}(X)
    yf = Vector{Float64}(y)
    XtX = Xf' * Xf
    Xty = Xf' * yf
    p = size(XtX, 1)
    # Iteration zero (plain ridge) and each reweighted step are SPD by
    # construction (λ > 0, w > 0), so no PosDef catch is needed.
    beta = cholesky(Symmetric(XtX + e.lambda * I)) \ Xty
    w = Vector{Float64}(undef, p)
    A = Matrix{Float64}(undef, p, p)
    iters = 0
    rel = Inf
    for outer iters in 1:e.max_iter
        @. w = 1.0 / (beta^2 + e.epsilon)
        copyto!(A, XtX)
        @inbounds for j in 1:p
            A[j, j] += e.lambda * w[j]
        end
        beta_new = cholesky(Symmetric(A)) \ Xty
        delta = mapreduce((a, b) -> abs(a - b), max, beta_new, beta)
        rel = delta / max(maximum(abs, beta_new), eps(Float64))
        beta = beta_new
        rel < e.tol && break
    end
    @debug "AdaptiveRidge solve" iterations = iters relative_change = rel
    return beta
end
```

```julia
# Memory log emitted by build_design_matrix_energy:

if verbosity
    bytes      = num_spinconfigs * num_salcs * sizeof(Float64)
    gram_bytes = num_salcs * num_salcs * sizeof(Float64)
    println(@sprintf(
        "Memory estimate: energy design matrix %s (%d × %d Float64), " *
        "Cholesky Gram %s (%d × %d). Total ~%s.",
        _format_bytes(bytes), num_spinconfigs, num_salcs,
        _format_bytes(gram_bytes), num_salcs, num_salcs,
        _format_bytes(bytes + gram_bytes),
    ))
end
# ... after construction ...
if verbosity
    println(@sprintf(
        "Energy design matrix built: %s actual.",
        _format_bytes(Base.summarysize(design_matrix)),
    ))
end

# Torque builder uses the same pattern but with the torque row count:
# rows = num_spinconfigs * 3 * num_atoms; the Gram matrix size
# (num_salcs²) is the same in both builders, so the torque pre-log
# typically dominates total memory when torque_weight > 0.
if verbosity
    rows       = num_spinconfigs * 3 * num_atoms
    bytes      = rows * num_salcs * sizeof(Float64)
    gram_bytes = num_salcs * num_salcs * sizeof(Float64)
    println(@sprintf(
        "Memory estimate: torque design matrix %s (%d × %d Float64), " *
        "Cholesky Gram %s (%d × %d). Total ~%s.",
        _format_bytes(bytes), rows, num_salcs,
        _format_bytes(gram_bytes), num_salcs, num_salcs,
        _format_bytes(bytes + gram_bytes),
    ))
end
```

A small `format_bytes(::Integer) -> String` helper (private,
`_format_bytes`) renders bytes as `KB` / `MB` / `GB` with one decimal.

## Types and conventions

- No physics-convention changes. `Jφ` unit and sign unchanged.
- `num_salcs²` Gram matrix is the new memory hot spot. Document this in the
  OLS docstring so users sizing problems know to budget it.
- Cholesky is called via `LinearAlgebra.cholesky` on `Symmetric(...)`; we
  rely on Julia stdlib, no new dependency.

## Impact on linked sites

- [ ] Spherical-harmonics convention (`TesseralHarmonics`): no.
- [ ] SCE coefficient XML (`save` / `load`): no.
- [x] `Fitting` <-> `SALCBasis`: design-matrix shape unchanged, only the
      memory log around its construction is new.
- [ ] `.claude/agents/` references: none, the agent file list still aligns.
- [x] `SPEC.md` / `docs/src/api.md` updates: refresh the OLS / Ridge
      paragraphs to mention Cholesky on the normal equations and the new
      error behavior. `docs/src/theory/` describes the regression model
      (objective and weighting), not the factorization method, so no
      update needed there.
- [x] `CHANGELOG.md` `[Unreleased]`: note the OLS solver change as
      "may surface previously-silent rank deficiency as an explicit
      error" — a behavioral change worth flagging even though API and
      numerics are otherwise preserved.

## Test strategy

- **Cross-check test**: on a known well-conditioned `(X, y)` (e.g. a
  random `Float64` design with `cond < 1e6`), assert
  `cholesky-OLS ≈ qr-OLS` to `1e-10`.
- **PosDef test**: construct a rank-deficient `X` (duplicate a column),
  call `solve_coefficients(OLS(), X, y)`, assert it throws `ArgumentError`
  whose message contains `"Ridge"` and `"linearly dependent"`.
- **Ridge cross-check**: confirm the new Cholesky-Ridge agrees with
  `MultivariateStats.ridge` to `1e-10` on a representative case.
- **Memory log**: a smoke test that captures stdout during
  `build_design_matrix_energy` under `verbosity = true` and asserts the
  estimate line is present and contains the expected byte count.

## Risks and open items

- **Behavioral regression for users who relied on the rank-deficient
  fallback**: anyone fitting at `torque_weight = 0` with
  `num_salcs >= num_spinconfigs` will now see an `ArgumentError` where they
  previously got coefficients (often garbage anyway). The error message
  tells them to switch to `Ridge`, which is the correct fix. Document in
  `CHANGELOG.md` `[Unreleased]`.
- **Cost of forming `X'X`**: for very tall thin `X`, the `O(np²)` Gram
  formation can compete with QR's blocked-Householder cost. In practice
  Cholesky still wins for SCE-sized problems, but this should be confirmed
  in the bench log entry.
- **Numerical agreement tolerance**: `1e-10` is loose enough to cover
  Cholesky vs QR rounding differences but tight enough to catch a sign or
  scale bug. If a real test case fails the tolerance, the right response is
  to investigate, not to loosen.
- **`Base.summarysize` cost**: O(matrix size), called once per build under
  the verbosity gate. Negligible compared to the construction itself.
