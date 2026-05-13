# Tasklist: write_xml public API cleanup

Status: **complete (2026-05-13)** — implementation commit `d9ae9ba`,
fast-forwarded into main and pushed; branch deleted.
Spec: [`requirements.md`](requirements.md) ·
[`design.md`](design.md)

Coarse-grained milestones. Daily todos use `TaskCreate`. Items below
are kept as-is for historical reference; the actual per-step progress
was tracked in the session-scoped task list.

## M0 — Baseline capture  ✅

- [x] `make test-all` green on the branch tip before any edit; record
      pass/fail counts for diff comparison after M2.
- [x] Capture a byte-level snapshot of one example XML artifact (e.g.
      the file produced by the `dimer` integration test) so the
      post-refactor XML can be `diff`-ed against it. Snapshots saved
      to `/tmp/dimer_baseline.xml` and `/tmp/fept_baseline.xml`.

## M1 — Delete four-arg public `write_xml`; reroute the doc example  ✅

- [x] Remove the four-arg `write_xml(structure, symmetry, basis_set,
      optimize, filename; write_jphi = true)` method from
      `src/Magesty.jl`.
- [x] Keep `XMLIO.write_xml` five-arg helper as-is so the surviving
      `SpinCluster` overload can still delegate to it. The
      `SpinCluster` overload now calls `XMLIO.write_xml(...)`
      directly (previously bounced through the deleted public method).
- [x] Update the surviving docstring(s) so the `# Examples` block no
      longer shows the deleted four-arg form. Consolidated into one
      docstring on the `SpinCluster` method that also explains the
      migration path for code that used to call the four-arg form.
      The `System` method keeps its existing concise docstring.
- [x] Update `docs/src/examples.md` line 186 to wrap the optimizer
      into a `SpinCluster` and call `write_xml(sclus, filename)`.

## M2 — Validation  ✅

- [x] `make test-unit` / `make test-integration` / `make test-jet` /
      `make test-aqua` green.
- [x] Byte-level XML diff vs the M0 snapshot is empty.
- [x] Confirm there are no other in-repo callers of the four-arg
      form. Grep `write_xml(.*,.*,.*,.*` returned only the two
      internal `XMLIO.write_xml` definitions plus the surviving
      `System` overload (`XMLIO.write_xml(structure, symmetry,
      basisset, filename)` — note this is the 4-arg *without*
      Optimizer, which is the internal helper for the `System` path).

## M3 — Wrap-up  ✅

- [x] Mark `DESIGN_NOTES.md` section R2 with a Status line pointing
      to the merge commit (`d9ae9ba`).
- [x] Commit, fast-forward into `main`, push. Branch
      `refactor/write-xml-public-api` deleted afterwards.

## Out-of-scope reminders

- Do NOT rename `write_xml` to `save` — that is the user-facing API
  spec.
- Do NOT change `xml_io.jl` schema constants — that landed in R7.
- Do NOT change the default filename `"jphi.xml"` — out of scope.
- Do NOT add a deprecation warning for the four-arg form — the
  decision is straight removal (only one in-repo caller, project at
  `0.1.0-DEV`).
