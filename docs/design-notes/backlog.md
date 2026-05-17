# Performance backlog

Lightweight notes and on-hold candidates from the second wave of
profiling. When picking one up, promote it to its own design note or to
a spec.

## DONE — Candidate F: fused `P̄ₗₘ` / `dP̄ₗₘ` in `∂ᵢZlm_unsafe`

**Status**: **complete** (commits `153e354` / `131d9a4`, around 2026-05).
`src/TesseralHarmonics.jl` gained a private helper
`_legendre_pair_unsafe!(buf, x, l, am) -> (P_am_l, P_am1_l)`: it calls
`_unsafednPl!` once to populate the `am`-degree cache, then advances
in-place by one step via `dPl_recursion` to recover `(am+1)`. The
buffered `∂ᵢZlm_unsafe(l, m, uvec, buf)` now goes through this helper,
collapsing two Legendre cache builds into one.

**Implementation notes (retained)**:
- `LegendrePolynomials._unsafednPl!` and `dPl_recursion` are both
  non-exported internals; upstream API changes need to be tracked.
- Buffer-size requirement is centralized in the module docstring
  ("Buffer requirements"): `length(buf) >= l - |m| + 1`. Guarded by
  `@boundscheck checkbounds(buf, _required_buf_size(l, m))`.
- Numerical results unchanged (covered by equivalence tests).

## DONE — Candidate E: pre-allocated buffer for `Zₗₘ_unsafe` (Magesty side complete)

**Status**: **complete on the Magesty.jl side** (commit `8a4a17d`,
2026-05-11). `P̄ₗₘ` / `dP̄ₗₘ_unsafe` / `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe`
gained 4-argument buffered overloads that go through
`LegendrePolynomials.dnPl(x, l, n, A)`, removing the per-call
`zeros(l-|m|+1)` heap allocation. The original 3-argument API stays.
Buffer-size requirement is centralized in the module docstring
(`length(buf) >= l - |m| + 1`); guarded since `131d9a4` by
`@boundscheck checkbounds(buf, _required_buf_size(l, m))`.

**Bench**: `Zₗₘ_unsafe` (l=4, m=2) standalone: 35.2 ns -> 20.8 ns
(1.70x, 80 B -> 0). Numerically equivalent (l=0..8 x all m x 9
directions equivalence test plus `@allocated == 0`).

**Remaining work (outside this repo)**: callers
(`SpinClusterMC` / `JPhiMagestyCarlo.jl`) need to pre-allocate
per-thread `Vector{Float64}(undef, max_l + 1)` and feed it to the
buffered overload (per-thread buffer is required in parallel sections).

## ON HOLD — Cache argument for SH buffers (formerly #1/#2)

**Target**: `src/Fitting.jl` `design_matrix_energy_element`, `calc_∇ₑu!`.

Idea: introduce `mutable struct SHCache` holding `sh_values` /
`atom_grad_values`, allocate per-thread in `build_design_matrix_*`, and
pass it through. Tried once; verdict:

- `calc_∇ₑu!` returns early on `atom_site_idx == 0` in many calls, so
  eager allocation at function entry hurt.
- SVector flavor was (-4.4% time / +14% memory); Vector flavor was
  (+6% time / -3% memory) — opposing trade-offs.
- BenchmarkTools noise is ~±5%, so significance was unclear.

Reasons it could come back:
- Lazy allocation matched to the early-return pattern (`Ref{Bool}` flag
  plus allocate on first work iteration).
- Parameterize `cb1.atoms` as `SVector{N, Int}` with `N` in the type and
  stack-allocate via `MVector{N, Int}`.
- Identify why the `Vector{SVector{3, Float64}}` accounting rose
  (`@code_warntype` / `@allocated` per line).

History: `.claude/bench_log.md` "#1/#2" section.

## RULED OUT — `cb1.atoms` as SVector (formerly #10)

**Status**: tried on 2026-05-16 and abandoned. Spec:
`docs/specs/260516-coupled-basis-atoms-svector/`.

Parameterized `CoupledBasis{R, N}.atoms` as `SVector{N, Int}` and re-ran
benchmarks. SALC build, including `projection_matrix_coupled_basis`,
regressed by **+12% time / +6 MB / +140K allocs**; the design-matrix
path was unchanged.

**Root cause**: `CoupledBasis` is stored inside
`SortedCounter{CoupledBasis}` / `Vector{CoupledBasis}` (UnionAll
containers). On iteration the type parameter `N` of `cb.atoms` is
erased, so `SVector{N, Int}` indexing is not specialized and becomes
slower than `Vector{Int}` (fixed element type, dynamic length). The
design-matrix side has a `where {R, N}` function barrier that
specializes the hot path, so it is unaffected.

Prerequisites for a re-attempt:
- Concrete containers (e.g., `Vector{CoupledBasis{R, N}}` per `(R, N)`,
  or grouped iteration with a function barrier).
- Avoid closure capture of `map(closure, atoms)` via `let` bindings.

Without those it is net negative. See the spec's "Outcome" section.

## CANDIDATES (verification needed; memory / alloc focus)

The second-wave investigation adopted only B and F. The items below
helped memory / allocations but the time difference was within noise
(±5%), so they were not adopted. Re-evaluate under GC-pressured
scenarios.

- **In-place broadcast for `representation_mat[row_range, col_range] = rot_mat * phase`**:
  SALCBasis memory -3.9%, allocs -224K. One-line change, no time effect.
  (Renamed 2026-05-16 from `temp_projection_mat`.)
- ~~**Replace `sum(conj.(t1) .* t2)` in `tensor_inner_product` with a
  fused `@simd` loop**: SALCBasis memory -8.8%, allocs -224K. Time in
  noise.~~ *(2026-05-16: already replaced with a deterministic scalar
  loop. Motivation was cross-platform reproducibility (reduction order
  independent of BLAS/LAPACK SIMD lane width). Memory / alloc savings
  are a side effect. Plain `@inbounds for`, not `@simd`.)*

Numbers: `.claude/bench_log.md` "C" / "D" sections (preserved in
perf-2-branch commit history).

## OUT OF SCOPE

- Former #7 (SALCBases `round.()` / masking fusion): not a bottleneck.
- Former #8 (`eigenvec` slice `@views`): not a bottleneck.

## LESSONS — approaches we tried and rejected

These regressed under the second-wave benchmarks. Julia's compiler
already optimizes many of these patterns; naive hand replacement makes
things worse.

- **Linearize `coeff_tensor[idx_buf..., mf_idx]` splatting via manual
  strides**: torque time +34%, energy time +1%. Julia's vararg dispatch
  is already fast at small N (=2); manual strides are slower.
- **Hand-unroll the three components of
  `mf_grad_contribution .+= coeff_val * product_other .* grad_atom`**:
  time +38%, allocs +88%, memory +38%. Broadcast `.+=` is fused
  in-place by the compiler; hand-rolled `MVector[i] += x` incurs
  `setindex!` overhead that can outweigh the benefit.

Takeaway: profiler sample counts are not proportional to time. Always
measure before / after with `@time` (5+ trials).
