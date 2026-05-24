# Tasklist: Design-matrix algorithmic restructuring

Status: draft (2026-05-24)

Branch: `refactor/design-matrix-restructuring`.

Coarse-grained, commit-sized milestones. Day-to-day tracking goes
through `TaskCreate` in-session.

Each implementation milestone (M1–M4) is gated by a bench entry in
`.claude/bench_log.md` and a per-stage design-matrix equivalence
test. Do not start the next milestone until the prior gate is met.

## Milestones

### M0 — Regression and benchmark harness

- [ ] Create `test/regression_fixtures/` directory. Capture
      pre-restructuring reference design matrices for FeGe 2x2x2
      and FePt tetragonal 2x2x2 fixtures via `writedlm`
      (`DelimitedFiles` stdlib; no new dependency). Each file
      carries a header line with the generating commit hash.
      Commit the fixtures to the repo (expected size: hundreds of
      KB total).
- [ ] Add a `design-matrix equivalence` testset in
      `test/component/test_fitting.jl` parameterized by `(rtol,
      atol)` so each stage can call it with the appropriate
      tolerance from the design.md table.
- [ ] Record baseline (current main, commit `4809a5a`) bench
      numbers in `.claude/bench_log.md` under this spec's header,
      so subsequent entries can reference them.

### M1 — A: Mf folding into `folded_tensor` (strategy a)

- [ ] `src/CoupledBases.jl`: add `folded_tensor::Array{Float64,
      R-1}` field; constructor computes it from `(coeff_tensor,
      coefficient)` in the existing Mf-iteration order (preserves
      bit-identity).
- [ ] `src/Fitting.jl` `design_matrix_energy_element` /
      `calc_∇ₑu!`: drop the `mf_idx` loop, read from
      `cbc.folded_tensor` instead of contracting
      `coeff_tensor[..., mf_idx] * coefficient[mf_idx]`.
- [ ] `src/XMLIO.jl`: serialize / deserialize `folded_tensor`
      (additive change; bump version field but accept old files
      and recompute `folded_tensor` lazily on load).
- [ ] Regression: equivalence test passes at `rtol = 1e-14, atol =
      1e-15` (expect bit-identical).
- [ ] Bench: `.claude/bench_log.md` entry, FeGe 2x2x2, 5-trial
      median. No minimum speedup gate (clean-up step).
- [ ] Commit + push.

### M2 — C: Pre-enumerate orbit clusters

- [ ] `src/CoupledBases.jl`: add `clusters::Vector{NTuple{N,
      Int}}` field; document the distinctness + canonical-sort
      invariants in the docstring.
- [ ] `src/SALCBases.jl`: during basis build, enumerate translated
      cluster atom-tuples via `symmetry.map_sym` and
      `symnum_translation`, sort canonically, dedup, and store on
      each `cbc`.
- [ ] `src/Fitting.jl` `design_matrix_energy_element` /
      `calc_∇ₑu!`: replace `for itrans in symnum_translation` + the
      `searched_pairs` / `_atoms_hash_key` / sort machinery with
      `for atoms in cbc.clusters`. Drop the now-dead workspace
      fields (`searched_pairs`).
- [ ] `src/XMLIO.jl`: recompute `clusters` lazily on load
      (default; flip to serialize-and-load if load time becomes a
      bottleneck).
- [ ] Regression: equivalence test passes at `rtol = 1e-13, atol =
      1e-14`.
- [ ] Bench: gate is **torque ≥ 1.15× over baseline**. If gate
      fails, profile and decide before M3.
- [ ] Commit + push.

### M3 — B: Per-spinconfig SH cache

- [ ] Decide energy threading axis (design.md "Risks and open
      items"). Default: keep `@threads for i = 1:num_salcs` and
      rebuild cache per spinconfig inside the inner `j` loop.
      Benchmark both before committing.
- [ ] Add an internal cache type (e.g., `SHCache`) with `Z::Matrix`
      / `∂Z::Matrix` and the layout-offset map. Not exported.
- [ ] `src/Fitting.jl` `build_design_matrix_*`: build the cache
      once per spinconfig from the union of `cbc.ls` values
      present in `salc_list`.
- [ ] `src/Fitting.jl` `design_matrix_energy_element` /
      `calc_∇ₑu!`: accept the cache, read from it instead of
      calling `Zₗₘ_unsafe` / `∂ᵢZlm_unsafe`. Workspace
      `sh_values` / `sh_offsets` / `legendre_buf` no longer needed
      in the hot path (delete from `EnergyWorkspace` /
      `GradWorkspace`).
- [ ] Regression: equivalence test passes at `rtol = 1e-14, atol =
      1e-15` (caching is exact).
- [ ] Bench: gate is **energy ≥ 2× over post-C**. If gate fails,
      profile and decide before M4.
- [ ] Commit + push.

### M4 — D: Torque as cluster-major reverse-mode Jacobian

- [ ] `src/Fitting.jl` `build_design_matrix_torque`: rewrite the
      loop to `(sc_idx, salc_idx, cbc, cluster) → site_in_cluster
      → iatom`. Allocate the per-thread per-spinconfig
      accumulation buffer of shape `(3, num_atoms, num_salcs)`.
- [ ] Replace `calc_∇ₑu!` with a cluster-major kernel: forward
      `Φ_cbc(cluster)` and the `N` partial-product slices assembled
      together; accumulate `∂Φ/∂spin[a_k]` for each
      site-in-cluster into the buffer. Add `@assert
      allunique(translated_atoms)` in `@debug` build.
- [ ] At the end of each spinconfig, reduce the buffer into
      `design_matrix` rows via one `cross(spin[a_k], ·)` per
      `(atom, salc_idx)` cell.
- [ ] Drop the dead `calc_∇ₑu!` API (or keep as a thin wrapper if
      called externally — check `tools/personal/` first).
- [ ] Regression: equivalence test passes at energy `rtol = 1e-12,
      atol = 1e-13` and torque `rtol = 1e-10, atol = 1e-12`.
- [ ] Bench: gate is **torque ≥ 5× over post-B**.
- [ ] Commit + push.

### M5 — Tier 2 review + merge

- [ ] Run the four-axis Tier 2 review panel (numerical /
      maintainability / performance / API) in parallel on the
      cumulative branch diff. Resolve blockers and major findings;
      apply numerical findings unconditionally.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `SPEC.md` updated only if needed (no public-API change
      planned).
- [ ] Open PR; resolve review; merge to main.
- [ ] Update `Status:` in this file and the row in
      `docs/specs/README.md` to `complete (YYYY-MM-DD)` together.
- [ ] Delete
      `docs/design-notes/design-matrix-algorithmic-restructuring.md`
      and remove its row from `DESIGN_NOTES.md` (per the
      design-note operating rule: completed design notes are
      folded into their spec and removed from the index).

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~
items that do not apply.

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] Regression tests added per stage (`test/component/test_fitting.jl`
      equivalence harness; reference fixtures in
      `test/regression_fixtures/`).
- [ ] ~~If public API changed: `SPEC.md` and `docs/src/api.md`
      updated.~~ (No public API change.)
- [ ] Before / after recorded in `.claude/bench_log.md` for each of
      M1, M2, M3, M4.
- [ ] Tier 2 review panel run after M4 lands; findings resolved.
- [ ] ~~If module names or Makefile targets changed:
      `.claude/agents/` swept and updated.~~ (No rename.)
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hashes appended below.

## Implementation commits

(Append commit hash + one-line summary as each milestone lands.)

- M0: ...
- M1: ...
- M2: ...
- M3: ...
- M4: ...
