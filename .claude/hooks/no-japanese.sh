#!/bin/bash
# PostToolUse hook: block Edit/Write that introduces Japanese characters
# into any file committed to the repository.
#
# Policy: CLAUDE.md "Language and terminology" — everything checked in is
# English-only. Conversation with the user stays Japanese-capable.
#
# Exemptions:
#   - tools/personal/        (personal scripts, not quality-assured)
#   - .claude/bench_log.md   (historical record kept as-is)
#   - docs/specs/<slug>/...  (per-slug spec working files; Japanese allowed.
#                             docs/specs/_template/ and docs/specs/README.md
#                             stay English.)

set -e

input=$(cat)
tool_name=$(printf '%s' "$input" | python3 -c 'import sys,json; print(json.load(sys.stdin).get("tool_name",""))' 2>/dev/null || echo "")

case "$tool_name" in
  Edit|Write|MultiEdit) ;;
  *) exit 0 ;;
esac

file_path=$(printf '%s' "$input" | python3 -c 'import sys,json; print(json.load(sys.stdin).get("tool_input",{}).get("file_path",""))' 2>/dev/null || echo "")
[ -z "$file_path" ] && exit 0
[ -f "$file_path" ] || exit 0

cwd=$(pwd)
rel="${file_path#$cwd/}"

# Only enforce on files inside this repository
case "$rel" in
  /*) exit 0 ;;  # absolute path outside cwd
esac

# Skip gitignored files (Manifest.toml, .claude/settings.local.json, etc.)
# and untracked artifacts. We rely on `git check-ignore`; if git isn't
# available we still proceed to the textual check.
if command -v git >/dev/null 2>&1; then
  if git check-ignore -q "$rel" 2>/dev/null; then
    exit 0
  fi
fi

# Skip exempt subtrees and historical records.
# Note: the template and the index keep an empty body (no exit), so they fall
# through to the English-only textual check below.
case "$rel" in
  tools/personal/*) exit 0 ;;
  .claude/bench_log.md) exit 0 ;;
  docs/specs/_template/*) ;;
  docs/specs/README.md) ;;
  docs/specs/*/*) exit 0 ;;
esac

hits=$(perl -CSD -ne 'print "  line $.: $_" if /[\x{3040}-\x{30FF}\x{4E00}-\x{9FFF}]/' "$file_path" || true)
if [ -n "$hits" ]; then
  {
    echo "BLOCKED: Japanese characters detected in $rel"
    echo "(policy: CLAUDE.md Language and terminology — committed files must be English-only)"
    echo "$hits"
    echo "Re-edit the file to replace the Japanese text with English."
  } >&2
  exit 2
fi

exit 0
