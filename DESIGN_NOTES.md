# Design notes index

Index of design discussions, investigations, and on-hold ideas. Follow the
links for detail. Active development units live under `docs/specs/`;
historical benchmarks live in `.claude/bench_log.md` and `git log`.

Operating rules: [`docs/design-notes/README.md`](docs/design-notes/README.md).

## Design proposals

| Topic | Status | Last update |
|---|---|---|

Completed proposals are folded into their corresponding spec under
`docs/specs/` and removed from this index. The spec folder is the
canonical history; the original design-note body is dropped (git history
preserves it).

## Investigations

- [`src/` refactoring review sweep](docs/design-notes/investigations/src-refactor-review.md)
  — multi-agent review retrospective; behavior-preserving cleanups landed,
  A2 kernel dedup and `kd_name` `<SpeciesOrder>` declined with reasoning
  (2026-06-01)
- [SpheriCart.jl adoption assessment](docs/design-notes/investigations/sphericart-adoption.md)
  — conclusion: not adopted (2–3x slower) (2026-05-11)

## Performance backlog

- [`Zₗₘ` / Legendre / SH buffer improvement ideas](docs/design-notes/backlog.md)
  — lightweight notes and on-hold candidates
