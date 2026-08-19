#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
fixtures="$root/crates/bellhop/tests/fixtures/golden"

"$root/tools/reference/compare-ray-linux.sh" "$fixtures/ParaBotCritical.env"
"$root/tools/reference/compare-ray-linux.sh" "$fixtures/DickinsCritical.env"
