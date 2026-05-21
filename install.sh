#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
target="${CODEX_HOME:-$HOME/.codex}/skills"

mkdir -p "$target"

for skill_md in "$repo_dir"/stereo-seq-*/SKILL.md; do
  skill_dir="$(dirname "$skill_md")"
  skill_name="$(basename "$skill_dir")"
  rsync -a --delete "$skill_dir"/ "$target/$skill_name"/
  echo "Installed $skill_name -> $target/$skill_name"
done

echo
echo "Done. Restart Codex to load the installed skills."
