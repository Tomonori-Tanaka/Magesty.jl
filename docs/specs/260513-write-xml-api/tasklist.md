# Tasklist: write_xml public API cleanup

Status: draft (2026-05-13)
Spec: [`requirements.md`](requirements.md) ·
[`design.md`](design.md)

Coarse-grained milestones. Daily todos use `TaskCreate`.

## M0 — Baseline capture

- [ ] `make test-all` green on the branch tip before any edit; record
      pass/fail counts for diff comparison after M2.
- [ ] Capture a byte-level snapshot of one example XML artifact (e.g.
      the file produced by the `dimer` integration test) so the
      post-refactor XML can be `diff`-ed against it.

## M1 — Delete four-arg public `write_xml`; reroute the doc example

- [ ] Remove the four-arg `write_xml(structure, symmetry, basis_set,
      optimize, filename; write_jphi = true)` method from
      `src/Magesty.jl`.
- [ ] Keep `XMLIO.write_xml` five-arg helper as-is so the surviving
      `SpinCluster` overload can still delegate to it.
- [ ] Update the surviving docstring(s) so the `# Examples` block no
      longer shows the deleted four-arg form. Consolidate into one
      docstring on the `SpinCluster` method (the `System` method
      keeps its existing concise docstring).
- [ ] Update `docs/src/examples.md` line 186 to wrap the optimizer
      into a `SpinCluster` and call `write_xml(sclus, filename)`.

## M2 — Validation

- [ ] `make test-unit` / `make test-integration` / `make test-jet` /
      `make test-aqua` green.
- [ ] Byte-level XML diff vs the M0 snapshot is empty.
- [ ] Confirm there are no other in-repo callers of the four-arg
      form (one final grep `write_xml(.*,.*,.*,.*` to be sure).

## M3 — Wrap-up

- [ ] Mark `DESIGN_NOTES.md` section R2 with a Status line pointing
      to the merge commit.
- [ ] Commit, fast-forward into `main` (or open a PR if user
      prefers), push. Branch deleted afterwards per the recent
      R-series workflow.

## Out-of-scope reminders

- Do NOT rename `write_xml` to `save` — that is the user-facing API
  spec.
- Do NOT change `xml_io.jl` schema constants — that landed in R7.
- Do NOT change the default filename `"jphi.xml"` — out of scope.
- Do NOT add a deprecation warning for the four-arg form — the
  decision is straight removal (only one in-repo caller, project at
  `0.1.0-DEV`).
