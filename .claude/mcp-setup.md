# MCP server setup notes (Claude Code only)

`/.mcp.json` at the repository root declares three project-scoped MCP
servers. These notes cover setup. CI is unaffected.

## 1. GitHub MCP (`@modelcontextprotocol/server-github`)

Structured access to PRs / issues / releases / tags / actions.

- **Requires**: `npx` (Node.js 18+) and a Personal Access Token (PAT).
- **Scopes**: `repo` (read/write) plus `read:org` is usually enough.
  Add `workflow` if you need release / tag operations.
- **Setup**: export `GITHUB_PERSONAL_ACCESS_TOKEN=<token>` in your shell
  (e.g., `.zshrc` / `.bashrc`, or via `direnv`). `.mcp.json` references
  it via `${...}` expansion.
- **When to use**:
  - "Review PR #15", "triage issues", "draft release notes", etc.
  - The `gh` CLI works too, but the MCP server returns structured
    responses that Claude can parse and act on directly.

## 2. Context7 (`@upstash/context7-mcp`)

Structured **up-to-date API documentation** for dependencies
(Julia / AtomsBase / StaticArrays / Spglib / Documenter / ...).

- **Requires**: `npx` only. No Upstash account needed (free public API).
- **When to use**:
  - "What is the latest `StaticArrays.SVector` API signature?"
  - "What does Documenter's `checkdocs = :exports` miss?"
  - Breaking-change scan when moving to Julia 1.13.
- WebFetch can do the same, but Context7 supports version-pinned queries
  and structured excerpts, which reduces back-and-forth.

## 3. arXiv MCP (`arxiv-mcp-server`)

Fetches paper abstracts / sections for references.

- **Requires**: `uvx` (i.e., `pipx`-like; `pip install uv` installs it).
  `pip install arxiv-mcp-server` also works.
- **When to use**:
  - Confirm an equation number in Drautz & Fähnle 2004 (PRB 69, 104404).
  - Pull latest revisions of Tanaka & Gohda 2025 (arXiv:2512.04458).
  - Keep `docs/src/technical_notes.md` aligned with the papers.
- Medium priority — used less often than GitHub or Context7.

## Skipping servers

Different contributors need different servers. You can leave `.mcp.json`
as-is and reject individual servers at the Claude Code approval prompt.
To disable persistently, add the server name to
`.claude/settings.local.json` under `disableMcpServers`.

## Removal policy

If a PAT is leaked, revoke it on GitHub first and `unset` the shell
environment variable. Never write tokens into `.mcp.json` — always use
`${...}` expansion.
