# Design notes index

Index of design discussions, investigations, and on-hold ideas. Follow the
links for detail. Active development units live under `docs/specs/`;
historical benchmarks live in `.claude/bench_log.md` and `git log`.

Operating rules: [`docs/design-notes/README.md`](docs/design-notes/README.md).

## Design proposals

| Topic | Status | Last update |
|---|---|---|
| [LASSO / Adaptive LASSO / Adaptive Ridge estimators](docs/design-notes/lasso-adaptive-estimators.md) | not started (GLMNet.jl chosen; pre-spec draft) | 2026-05-16 |
| [Pre-release cleanup (v0.1.0)](docs/design-notes/pre-release-cleanup.md) | not started | 2026-05-17 |
| [Julia compat policy](docs/design-notes/julia-compat-policy.md) | operating rule | 2026-05-17 |

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
