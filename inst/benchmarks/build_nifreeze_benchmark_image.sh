#!/usr/bin/env bash
set -euo pipefail

image="${1:-volregger-nifreeze-benchmark:latest}"
platform="${VOLREGGER_NIFREEZE_PLATFORM:-linux/amd64}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

exec docker build \
  --platform "$platform" \
  -t "$image" \
  -f "$repo_root/docker/nifreeze-benchmark/Dockerfile" \
  "$repo_root/docker/nifreeze-benchmark"
