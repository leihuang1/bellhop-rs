# Development tools

## Pinned Fortran differential reference

`reference/Dockerfile` builds the official Acoustics Toolbox `v2023.5`
BELLHOP at commit `475108519289c6fb488b58980c644ea14eccc604` for Linux x86-64.
The Debian base image, source-archive SHA-256, GNU Fortran version, and compiler
flags are fixed in the image definition.

Build the image:

```sh
tools/reference/build-image.sh
```

The build script downloads the commit archive with retries, verifies
SHA-256 `f8a7a2c1e80a73431cd230a10bef5fcfc996c88889a0e1540771c3922ee2a21f`,
and removes the temporary archive after the image is built.

Run the reference model for one case:

```sh
tools/reference/run-case.sh path/to/case.env
```

Compare an `R` run against Rust on the host, or on the authoritative pinned
Linux x86-64 Rust 1.88 environment:

```sh
tools/reference/compare-ray.sh path/to/case.env
tools/reference/compare-ray-linux.sh path/to/case.env
```

The semantic comparator checks launch angles, bounce counts, and trajectory
coordinates. It aligns isolated `1e-4 × base step` vertices because a value
within a few ulps of an SSP or boundary interface can make one compiler take
one minimum step while another reflects or changes segment immediately. Strict
aligned coordinates still use `1e-5 m`; the branch alignment window defaults
to `4.1` minimum steps and is reported separately. Override with:

```sh
BELLHOP_DIFFERENTIAL_POSITION_TOLERANCE_M=1e-6 \
BELLHOP_DIFFERENTIAL_MINIMUM_STEP_FACTOR=1 \
  tools/reference/compare-ray-linux.sh path/to/case.env
```

Run the committed critical boundary/interface cases with:

```sh
tools/reference/check-critical-rays.sh
```

Reference outputs are written below `target/reference/` and are not committed.
Ordinary parser and numerical tests use the curated fixtures under
`crates/bellhop/tests/fixtures`.
