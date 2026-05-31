# Tasklist: src architecture and hot-path cleanup

Status: draft (2026-06-01)

This file holds coarse-grained, commit-sized milestones. Day-to-day
tracking goes through `TaskCreate` in-session.

## Milestones

### M1 — A1: cache encapsulation

- [ ] Add `CoupledBases.angular_momentum_coupling_results(ls, isotropy)`.
- [ ] Repoint the `SALCBases` lookup to the accessor; drop the direct
      `_angular_momentum_cache` / `:none` reference.
- [ ] `make test-all` green; commit.

### M2 — C: hot-path allocations

- [ ] C5 `(false, true)` literal (trivial).
- [ ] C3 / C4 SALCBases translation-equivalence / `find_translation_atoms`
      allocation removal.
- [ ] C2 `projection_matrix_coupled_basis` in-place matrix write.
- [ ] C1 `build_sh_cache_torque` single Legendre recursion.
- [ ] `bench-salcbasis` / `bench-sphericart` before/after to
      `.claude/bench_log.md`; `make test-all` green; commit.

### M3 — A2: Fitting kernel de-duplication (conditional)

- [ ] Extract the shared folded-tensor contraction skeleton (`@inline`).
- [ ] Micro-benchmark energy + torque kernels; land only if neutral-or-better.
- [ ] If regressed: revert and record the decision here.
- Exit: either a commit landing the helper, or an explicit decision note
  (no commit). Mark `[x]` in either case.

## Exit checklist

- [ ] `make test-all` passes.
- [ ] `make test-aqua` / `make test-jet` clean (no new warnings).
- [ ] ~~If results changed: regression or validation test added.~~ (results
      must not change; guarded by the existing suite)
- [ ] ~~If public API changed: `SPEC.md` and `docs/src/api.md` updated.~~
      (internal-only)
- [ ] If a hot path was touched: before / after recorded in
      `.claude/bench_log.md`.
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [ ] ~~If module names or Makefile targets changed: `.claude/agents/`
      swept.~~ (none)
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in
      `docs/specs/README.md` updated in sync.
- [ ] Implementation commit hash appended below.
