# Design notes index

Index of design discussions, investigations, and on-hold ideas. Follow the
links for detail. Active development units live under `docs/specs/`;
historical benchmarks live in `.claude/bench_log.md` and `git log`.

Operating rules: [`docs/design-notes/README.md`](docs/design-notes/README.md).

## Design proposals

| Topic | Status | Last update |
|---|---|---|
| [Design-matrix algorithmic restructuring](docs/specs/260524-design-matrix-restructuring/) — Mf folding, per-spinconfig SH cache, orbit pre-enumeration, torque as Jacobian (body still in `docs/design-notes/design-matrix-algorithmic-restructuring.md`) | in progress (spec) | 2026-05-24 |

Completed proposals are folded into their corresponding spec under
`docs/specs/` and removed from this index. The spec folder is the
canonical history; the original design-note body is dropped (git history
preserves it).

## Investigations

- [SpheriCart.jl adoption assessment](docs/design-notes/investigations/sphericart-adoption.md)
  — conclusion: not adopted (2–3x slower) (2026-05-11)

## Performance backlog

- [`Zₗₘ` / Legendre / SH buffer improvement ideas](docs/design-notes/backlog.md)
  — lightweight notes and on-hold candidates
