# KRAKEN Rust rewrite: compatibility target

The numerical reference is Acoustics Toolbox `v2023.5`, commit
`475108519289c6fb488b58980c644ea14eccc604`, using the pinned Linux x86-64 GNU
Fortran 12.2 environment described in [`reference.md`](reference.md).

**Status:** the first Pekeris slice is implemented in `crates/kraken`. It is
not the full matrix below; unsupported inputs in this slice are rejected.

## Products

The target is the two-dimensional normal-mode workflow:

1. `KRAKEN` computes real-eigenvalue normal modes.
2. `KRAKENC` computes complex modes.
3. `FIELD` synthesizes complex frequency-domain pressure from a mode set.

The Rust result model exposes modes and pressure fields, not BELLHOP-style ray
arrivals. Rust output will use a versioned HDF5 schema rather than the Fortran
`.mod` and `.shd` binary formats. A wideband time-domain response would be a
separate product and acceptance contract.

## Current slice

`crates/kraken` currently accepts a narrow legacy `.env`/`.flp` subset:

- one frequency and one homogeneous, lossless fluid water layer using `N`
  interpolation;
- a smooth pressure-release surface (`V`) and smooth acoustic fluid
  half-space (`A`) with zero attenuation;
- a coherent, omnidirectional line source (`X`, `O`, `C`) and one
  range-independent FIELD profile at 0 km;
- finite, in-water source/receiver depths, with FIELD depths covered by the
  mode-sampling depths in `.env`.

The mode solver uses the closed-form Pekeris dispersion relation, normalizes
mode shapes including the decaying bottom-half-space tail, and synthesizes the
line-source coherent field. The legacy parser rejects other solver, boundary,
attenuation, source, and FIELD options instead of silently changing their
meaning. The fixture is constructed for this repository; upstream `tests/PekerisRD`
is a BELLHOP case, not a KRAKEN reference. The pinned Fortran golden comparison
is `crates/kraken/tests/pekeris_reference.rs` against
`crates/kraken/tests/fixtures/Pekeris.{env,flp}`.

This is an implementation slice, not a reduction of the final support target.
The measured difference from the pinned finite-difference reference is recorded
in [`deviations.md`](deviations.md).

## Planned environment support

The planned legacy adapters accept KRAKEN `.env` and FIELD `.flp` files,
including same-stem auxiliary resources used by selected cases. A strict,
self-contained JSON adapter is also planned, following BELLHOP's
single-document input convention.

The v2023.5 KRAKEN environment reader supports these SSP interpolation options:

- `N`: N²-linear
- `C`: C-linear
- `P`: PCHIP
- `S`: cubic spline
- `A`: analytic profile

KRAKEN profiles are range-independent; range variation is represented as a
sequence of profiles for FIELD propagation. The final environment matrix
includes fluid and elastic layers, attenuation units (`N`, `F`, `M`, `m`, `W`,
`Q`, `L`), volume attenuation (`T`, `F`, `B`), multiple frequencies, and the
reference's top/bottom conditions (`V`, `R`, `A`, `F`, `P`). These cover vacuum,
rigid, half-space, tabulated-reflection, and precomputed-impedance paths.
Existing reflection tables are inputs; table generation by the separate
`BOUNCE` program is excluded.

## Planned FIELD support

The final 2D FIELD matrix includes:

- point, line, and scaled-cylindrical source geometry;
- omnidirectional and tabulated source patterns;
- coherent and incoherent mode addition;
- range-independent fields and multi-profile adiabatic or coupled-mode
  propagation;
- multiple frequencies, source depths, receiver depths/ranges, and receiver
  range offsets.

The reference rejects incoherent addition with coupled modes; Rust validation
will preserve that rule. FIELD's separate three-dimensional `FIELD3D` program
is out of scope, as are BOUNCE table generation and ray-arrival products.

## Incremental implementation and acceptance

The first slice compares modes and a range-independent coherent field against
pinned Fortran. Further representative v2023.5 cases provide feature coverage:

| Feature | Reference cases |
|---|---|
| Classic modes and FIELD | `tests/Munk/MunkK.env` + `.flp`, `tests/sduct/sductK.env` + `.flp` |
| Complex/leaky modes | `tests/MunkLeaky/MunkKwb`, `MunkKbb`, `MunkKleaky` |
| Multi-frequency | `tests/BroadBand/MunkK` |
| Adiabatic and coupled FIELD | `tests/Gulf/gulf_ad.flp`, `gulf_cm.flp` |
| Reflection inputs | `tests/TabRefCoef/neggradK_*` |

Acceptance requires differential coverage for modal wavenumbers and
attenuation, normalized/aligned mode shapes, and complex pressure-field samples.
Mode-shape comparisons account for the arbitrary sign/phase convention of
eigenvectors. Small committed goldens keep ordinary tests independent of Docker;
the pinned reference workflow covers official end-to-end cases.

## Repository shape

The implementation follows the existing BELLHOP boundaries: `crates/kraken`
contains validated cases, legacy adapters, mode solving, and FIELD. CLI and HDF5
adapters can follow once the numerical path is stable. No shared acoustics
abstraction is planned until real duplication justifies one.
