#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
image=${BELLHOP_REFERENCE_IMAGE:-bellhop-rs-reference:v2023.5-amd64}
commit=475108519289c6fb488b58980c644ea14eccc604
expected=f8a7a2c1e80a73431cd230a10bef5fcfc996c88889a0e1540771c3922ee2a21f
archive="$root/tools/reference/Acoustics-Toolbox-$commit.tar.gz"
trap 'rm -f "$archive" "$archive.tmp"' EXIT

actual=''
if [[ -f $archive ]]; then
  actual=$(shasum -a 256 "$archive" | awk '{print $1}')
fi
if [[ $actual != "$expected" ]]; then
  rm -f "$archive.tmp"
  curl --fail --silent --show-error --location \
    --retry 10 --retry-all-errors --connect-timeout 30 --max-time 900 \
    --output "$archive.tmp" \
    "https://codeload.github.com/oalib-acoustics/Acoustics-Toolbox/tar.gz/$commit"
  actual=$(shasum -a 256 "$archive.tmp" | awk '{print $1}')
  if [[ $actual != "$expected" ]]; then
    echo "reference archive SHA-256 mismatch: expected $expected, got $actual" >&2
    rm -f "$archive.tmp"
    exit 1
  fi
  mv "$archive.tmp" "$archive"
fi

docker build \
  --platform linux/amd64 \
  --tag "$image" \
  --file "$root/tools/reference/Dockerfile" \
  "$root/tools/reference"

rm -f "$archive"

docker run --rm --platform linux/amd64 --entrypoint /bin/sh "$image" \
  -c 'cat /opt/reference-build.txt && printf "architecture %s\n" "$(uname -m)"'
