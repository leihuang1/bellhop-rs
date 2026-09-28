# Development tools

## Pinned Fortran differential reference

`reference/Dockerfile` builds the official Acoustics Toolbox `v2023.5`
BELLHOP, KRAKEN, KRAKENC, and 2D FIELD programs at commit
`475108519289c6fb488b58980c644ea14eccc604` for Linux x86-64. The Debian base
image, source-archive SHA-256, GNU Fortran version, and compiler flags are
fixed in the image definition.

Build the image:

```sh
tools/reference/build-image.sh
```

The build script downloads the commit archive with retries, verifies
SHA-256 `f8a7a2c1e80a73431cd230a10bef5fcfc996c88889a0e1540771c3922ee2a21f`,
and removes the temporary archive after the image is built.

Run the BELLHOP reference for one case:

```sh
tools/reference/run-case.sh path/to/case.env
```

Run KRAKEN or KRAKENC on a legacy environment. If a same-stem `.flp` exists,
the helper also runs 2D FIELD; it saves the reference `.mod` and `.shd` files
under `target/reference/`. The planned Rust support matrix is in
[the KRAKEN compatibility document](../docs/kraken-compatibility.md):

```sh
tools/reference/run-kraken-case.sh kraken path/to/MunkK.env
tools/reference/run-kraken-case.sh krakenc path/to/MunkKleaky.env
```

Numerically compare a supported Rust single-fluid `.env`/`.flp` pair
against fresh KRAKEN/FIELD output (not merely a smoke test):

```sh
for case in Pekeris PekerisFiltered PekerisDense PekerisDenseLoss PekerisRefined PekerisSpline3 PekerisRigid PekerisRigidLoss PekerisHard PekerisHardBoth PekerisRigidPlane MunkLossless MunkBottomLoss MunkAnalytic SductTrapped SductPchip SductSpline; do
  tools/reference/compare-kraken.sh "crates/kraken/tests/fixtures/$case.env"
done
```

For the modes-only, derived KRAKENC Pekeris cases (faster or slower fluid
bottom, default or explicit restart option, and filtered phase-speed interval),
run the same pinned Fortran calculation and comparator used by CI:

```sh
for case in PekerisComplex PekerisComplexBlank PekerisComplexSlow PekerisComplexCLow; do
  tools/reference/run-kraken-case.sh krakenc "crates/kraken/tests/fixtures/$case.env"
  KRAKEN_COMPLEX_CASE="$case" \
  KRAKEN_COMPLEX_REFERENCE_ROOT="$PWD/target/reference/$case-krakenc/$case" \
    cargo test --release -p kraken --test differential_reference \
      complex_fluid_matches_fresh_reference -- --ignored --exact --nocapture
done
```

This compares every mode and pressure sample, including modal print precision,
mode-shape phase alignment, dimensions, and coordinate vectors. CI runs all
seventeen cases with fixed tolerances, plus the full unmodified upstream MunkK
pair and the unmodified upstream MunkAnalytic `.env` with a derived coherent
FIELD `.flp`. The original MunkAnalytic three-line `.flp` fails in v2023.5
FIELD itself and is not an accepted upstream pair. Raw small reference goldens
and their [hashes/provenance](../crates/kraken/tests/fixtures/golden/README.md)
also run in ordinary tests without Docker. MunkBottomLoss, MunkLossless and
SductTrapped are explicitly derived, reduced-grid fixtures. Unmodified MunkK
is now compared end-to-end; unmodified sduct still needs leaky-mode support.

Compare an `R` run against Rust on the host, or on the authoritative pinned
Linux x86-64 Rust 1.88 environment:

```sh
tools/reference/compare-ray.sh path/to/case.env
tools/reference/compare-ray-linux.sh path/to/case.env
```

Compare an arrivals (`A`/`a`) or pressure-field (`C`/`S`/`I`) run the same
way. The arrival comparator checks receiver counts and all eight arrival
fields at the committed single-precision storage tolerances; the pressure
comparator parses the fixed-record `.shd` layout, checks `LRecl` and the
receiver vectors, and compares every complex sample with `5e-8` absolute
pressure (override with `BELLHOP_DIFFERENTIAL_PRESSURE_TOLERANCE` and
relative tolerance `BELLHOP_DIFFERENTIAL_PRESSURE_RELATIVE_TOLERANCE`):

```sh
tools/reference/compare-arrival.sh path/to/case.env
tools/reference/compare-field.sh path/to/case.env
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

Fresh reference outputs are written below `target/reference/` and are not
committed automatically. Ordinary parser and numerical tests use curated
fixtures under `crates/bellhop/tests/fixtures` and `crates/kraken/tests/fixtures`.
