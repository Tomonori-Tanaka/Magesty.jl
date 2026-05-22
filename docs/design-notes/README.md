# Design notes

Home for design discussions, investigations, and the performance
backlog. `DESIGN_NOTES.md` at the repository root is the index; each
topic body lives here as one file.

## When to use which directory

| Location | Contents | Example |
|---|---|---|
| `docs/specs/` | Agreed-on, in-progress, or completed development units (`requirements.md` / `design.md` / `tasklist.md`) | `docs/specs/260513-write-xml-api/` |
| `docs/design-notes/` | Pre-spec design discussions, investigations, backlog | this directory |
| `docs/src/` | User-facing docs built by Documenter | `docs/src/api.md` |

## File naming

- English kebab-case (`<topic>.md` or `<topic>-<aspect>.md`).
- No date prefix (topic-first; the timestamp lives in the in-file
  `Status:` line).
- One topic per file (a single proposal, or a tight group of
  investigations).
- Post-mortem materials (investigations, benchmarks, adoption
  assessments) go under `investigations/`.

## Operating rules

1. **New topic**: create `docs/design-notes/<topic>.md` and add a row
   to `DESIGN_NOTES.md`.
2. **When the topic becomes a spec**: point the `DESIGN_NOTES.md` row
   at the new `docs/specs/...` path. The design-note body stays in
   place only while the spec is active.
3. **When the spec completes**: delete the design-note file and remove
   its row from `DESIGN_NOTES.md`. The spec folder is the canonical
   history; git preserves the original design-note body.
4. **Length**: if a single file exceeds ~500 lines, consider splitting.

## File header template

```markdown
# <Title>

**Status**: not started / in progress / complete (YYYY-MM-DD)

[1-2 line conclusion summary, if useful]

## Background

...
```

## Internal references

Do not reference files in this directory from `.jl` source (CLAUDE.md
rule). This directory is Claude-collaboration scaffolding; source code
must read independently.
