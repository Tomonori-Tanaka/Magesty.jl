# Tasklist: CoupledBasis type parameterization

## Step 1 — baseline benchmark

- [x] Record pre-refactor `@time` for `build_design_matrix_energy` on
      the fept example (5 trials, min + median). Saved to
      `.claude/bench_log.md` "B1 baseline" section.
- [x] Capture `@code_warntype` excerpt of
      `design_matrix_energy_element` showing the `Any` return on
      `coeff_tensor[...]`. (line `%46 ... ::ABSTRACTARRAY` recorded.)

## Step 2 — struct definitions (done)

- [x] `src/types/Basis.jl`: parameterize `CoupledBasis{R}`,
      `CoupledBasis_with_coefficient{R}`,
      `AngularMomentumCouplingResult{R}` with
      `coeff_tensor::Array{Float64, R}`.
- [x] Constructors derive `R` from `ndims(coeff_tensor)` and call
      `new{R}(...)`. Input tensors `convert`-ed to
      `Array{Float64, R}`.
- [x] `make test-unit` (20607/20607) and `make test-integration`
      (813/813) green.
- [x] Mid-step warntype check: `cbc.coeff_tensor::Array{Float64, 3}`,
      `Mf_size::Int64`, `mf_idx::Int64` now type-stable.
- [x] Mid-step bench: allocs −1.3% (energy) / −4.6% (torque), wall
      time in noise. Confirms Step 3 (splat fix) is needed to convert
      the type-stability win into a wall-time win.

## Step 2 — struct definitions

- [ ] `src/types/Basis.jl`: parameterize `CoupledBasis{R}`,
      `CoupledBasis_with_coefficient{R}`, with
      `coeff_tensor::Array{Float64, R}`.
- [ ] Update both constructors to derive `R` from `ndims(coeff_tensor)`
      and call `new{R}(...)`.
- [ ] Update `AngularMomentumCouplingResult{R}` similarly.
- [ ] Adjust `reorder_atoms`, `convert_to_coupled_basis`,
      `CoupledBasis_with_coefficient(cb::CoupledBasis, ...)` to
      preserve / forward `R`.

## Step 3 — splat indexing (B3 merged into spec 2026-05-16)

Step 2 alone left wall-time unchanged because `coeff_tensor[idx_buf...,
mf_idx]` still returns `Any` (splat over a `Vector{Int}` defeats static
dispatch). Step 3 fixes this by giving the indexing path
compile-time-sized buffers, exploiting `R` from the type parameter.

- [x] `src/Optimize.jl::design_matrix_energy_element` —
      parameterize the function on `R`, use
      `dims_t = ntuple(i -> 2*cbc.ls[i]+1, Val(R-1))`,
      `idx_buf::MVector{R-1, Int}`,
      `translated_atoms / atoms_sorted_buf::MVector{R-1, Int}`,
      `other_dims_t::NTuple{R-2, Int}`,
      `CartesianIndices(other_dims_t)`.
- [x] Same pattern in `src/Optimize.jl::calc_∇ₑu!` (`other_sites_buf`,
      `other_dims_buf` → `MVector{R-2, Int}`, with
      `CartesianIndices(Tuple(other_dims_buf))` giving
      `CartesianIndices{R-2}`).
- [x] `src/SALCBases.jl`, `src/utils/xml_io.jl`: no changes needed;
      UnionAll usages compile unchanged, reconstruction path infers `R`
      at construct time. Verified by passing test suites.
- [x] Re-run `@code_warntype`: `ANY count: 0`. `tensor_result`,
      `mf_contribution`, `product_other` all `Float64`.

## Step 4 — tests

- [x] `make test-unit` green (20607/20607, after Step 2 and Step 3).
- [x] `make test-integration` green (813/813); XML byte-stability test
      passes on `fept_tetragonal_2x2x2`. FeGe end-to-end test dropped
      26.9s → 4.2s as a side effect.
- [x] `make test-aqua` green (10/10).
- [x] `make test-jet`: no issues found.

## Step 5 — post benchmark (done in Step 3)

- [x] Final benchmark recorded under "After B1 + B3 (Step 3)" in
      `.claude/bench_log.md`. Wall time improved ×24 (energy) / ×8.6
      (torque) on the fept fixture — well above the noise floor.

## Step 6 — close out

- [x] `code-reviewer` agent on cumulative diff (2 minor issues raised
      and resolved: AMCResult inner constructor parity, rank-erasing
      annotations in `build_design_matrix_energy` removed).
- [x] `docs/design-notes/post-step7-cleanup.md` B1/B3 marked resolved.
- [x] Single breaking commit via `git-helper`: `d3c2e8f
      refactor(basis)!: parameterize CoupledBasis on tensor rank`.
- [x] PR #7 opened and squash-merged to `main` as `6a062b5` on
      2026-05-16. Topic branch deleted.
