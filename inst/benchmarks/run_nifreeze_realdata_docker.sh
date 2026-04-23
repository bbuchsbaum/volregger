#!/usr/bin/env bash
set -euo pipefail

image="${VOLREGGER_NIFREEZE_IMAGE:-volregger-nifreeze-benchmark:latest}"
platform="${VOLREGGER_NIFREEZE_PLATFORM:-linux/amd64}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
forwarded_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --image)
      image="$2"
      shift 2
      ;;
    --help|-h)
      cat <<EOF
Usage:
  $(basename "$0") [--image IMAGE] [nifreeze-args...]

Runs inst/benchmarks/run_nifreeze_realdata.py inside the committed Docker image.
Build the image first with:
  inst/benchmarks/build_nifreeze_benchmark_image.sh [IMAGE]
EOF
      exit 0
      ;;
    *)
      forwarded_args+=("$1")
      shift
      ;;
  esac
done

if ! docker image inspect "$image" >/dev/null 2>&1; then
  echo "Docker image not found: $image" >&2
  echo "Build it first with: inst/benchmarks/build_nifreeze_benchmark_image.sh $image" >&2
  exit 2
fi

exec docker run --rm \
  --platform "$platform" \
  --user "$(id -u):$(id -g)" \
  -e PYTHONUNBUFFERED=1 \
  -v "$repo_root:/work" \
  -v /tmp:/tmp \
  -w /work \
  "$image" \
  "${forwarded_args[@]}"
