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
2. **When the topic becomes a spec**: replace the link in
   `DESIGN_NOTES.md` with the new `docs/specs/...` path. Keep the
   design-note body as a historical record (do not delete it).
3. **Completed state**: leave the file in place. Track status via the
   `**Status**:` line under the title and the `DESIGN_NOTES.md` index.
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
