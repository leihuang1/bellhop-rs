#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 {kraken|krakenc} CASE.env" >&2
  exit 2
fi
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
engine=$1
case_path=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
stem=$(basename "$case_path" .env)
reference="$root/target/reference/$stem-$engine/$stem"
"$root/tools/reference/run-kraken-case.sh" "$engine" "$case_path"

compare() {
  KRAKEN_FREQUENCY_SOLVER="$engine" KRAKEN_DIFFERENTIAL_ENV="$case_path" KRAKEN_DIFFERENTIAL_ROOT="$reference" \
    cargo test --release --manifest-path "$root/Cargo.toml" -p kraken --test differential_reference \
      multifrequency_fluid_matches_fresh_reference -- --ignored --exact --nocapture
}
compare
"$root/target/release/kraken" run "$case_path" --solver "$engine" --output "$reference.h5"
KRAKEN_HDF5_RESULT="$reference.h5" compare
