# Julia 1.10 (LTS) compatibility for the General-registry floor

**Status**: complete (2026-06-15)

Conclusion: the `[compat]` `julia` floor is kept at **1.12** for the initial
General-registry release. Lowering to 1.10 (LTS) is attractive for adoption and is
unblocked at the dependency-resolution level, but CI on Julia 1.10 surfaced two
1.10-specific runtime failures that need separate, reviewed fixes. Revisit 1.10/1.11
support as a follow-up.

## Background

The package developed against Julia 1.12. The only thing pinning the floor at 1.12 was
the test-only `JET` dependency: `JET` `0.10 – 0.11.x` requires `julia = "1.12"`, while
`JET` `0.9.x` supports `julia >= 1.10`. Widening `JET` to `"0.9, 0.10, 0.11"` lets Pkg
pick a per-version-compatible JET, and a registry check confirmed every production
dependency supports 1.10 (ProgressMeter and Roots newest versions require `>= 1.10`,
fixing 1.10 as the lowest feasible floor). The package source has no 1.12-only syntax
(no `public` statements, no `VERSION`/`@static` gates).

So lowering looked safe at the metadata level. CI on 1.10 (ubuntu + macos) proved
otherwise: the floor change was attempted and reverted in the same registry-prep branch.

## Findings — two 1.10-specific failures

Both pass on 1.12 and latest; both fail only on 1.10 (ubuntu and macos identically).

### 1. Solver rank-deficiency detection depends on LAPACK behavior

`test/component/test_fitting_solver.jl` — testsets "OLS PosDefException → ArgumentError
pointing to Ridge" (line ~76) and "Ridge(λ ≈ 0) delegates to OLS Cholesky" (line ~119).

The OLS/Ridge path solves via `cholesky(Symmetric(X'X + Λ)) \ (X'y)` (see
`src/Fitting.jl`). For a rank-deficient design it relies on `cholesky` throwing
`PosDefException`, which is caught and rethrown as an `ArgumentError` whose message names
`Ridge` and "linearly dependent". On Julia 1.10's LAPACK the factorization of the
rank-deficient fixture (`_rank_deficient_system(rng, 100, 8)`) does **not** throw, so the
solver returns a (min-norm-ish) result and the error path never runs — the test captures
`err === nothing`.

Possible fix, to make detection version-independent: replace the exception-driven check
with an explicit one, e.g. `F = cholesky(A; check = false); issuccess(F) || throw(...)`.
Caveat: for a *numerically borderline* (not exactly) singular Gram, `potrf`'s `info`
return can still differ across LAPACK versions, so `issuccess` alone may not be enough —
a rank/conditioning-based guard may be needed. This touches the solver's numerical error
semantics and must go through the numerical-correctness review before landing.

### 2. `Structure(::SystemSpec; verbosity)` MethodError on 1.10

`test/component/test_cluster_cutoff_criterion.jl:57` calls
`Magesty.Structure(system; verbosity = false)` where `system::SystemSpec`
(`Magesty.InputSpecs.SystemSpec`). On 1.10 this raises
`MethodError: no method matching Structure(::SystemSpec; verbosity::Bool)` even though the
method is defined at `src/Structures.jl:189` (`Structure(system::SystemSpec; verbosity)`).
It resolves fine on 1.12. The signature being present yet not matched points at a
version-dependent dispatch / type-identity issue (candidate: the runtests.jl pattern of
`include`-ing internal module files alongside `using Magesty` producing distinct type
bindings, surfacing differently across Julia versions). Needs a focused repro on 1.10.

## Verdict / next steps

- Initial General registration: floor stays `julia = "1.12"`, `JET = "0.11"`, CI matrix
  `1.12` (minimum) + `1` (latest).
- A future 1.10/1.11 effort should: (a) make solver rank-deficiency detection robust and
  version-independent (numerical review required); (b) root-cause and fix the `Structure`
  dispatch failure on 1.10; (c) add 1.10 and 1.11 to the CI matrix and confirm green
  before lowering the floor. Probing 1.11 alone first is cheap and would show whether
  either issue is strictly 1.10-specific.
