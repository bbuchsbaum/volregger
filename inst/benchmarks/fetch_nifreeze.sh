#!/usr/bin/env bash
# Fetch the pinned nifreeze checkout used by inst/benchmarks scripts.
# Lives under vendor/nifreeze/ and is excluded from the R package build
# (.Rbuildignore) and from git (.gitignore).
#
# Usage: inst/benchmarks/fetch_nifreeze.sh [SHA]
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
target="$repo_root/vendor/nifreeze"
remote="https://github.com/nipreps/nifreeze"
sha="${1:-a3985fe5a4763f0a2eafb24e78d88fcfad60992a}"

if [ -d "$target/.git" ]; then
  current="$(git -C "$target" rev-parse HEAD)"
  if [ "$current" = "$sha" ]; then
    echo "nifreeze already at $sha"
    exit 0
  fi
  echo "Updating $target from $current to $sha"
  git -C "$target" fetch --depth 50 origin
  git -C "$target" checkout --detach "$sha"
  exit 0
fi

mkdir -p "$(dirname "$target")"
echo "Cloning $remote into $target"
git clone --filter=blob:none "$remote" "$target"
git -C "$target" checkout --detach "$sha"
echo "nifreeze pinned to $sha"
