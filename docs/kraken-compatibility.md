# KRAKEN Rust rewrite: compatibility target

This document records the planned acceptance scope for a future `kraken` crate.
It is not a claim that KRAKEN support exists in this repository yet. The
numerical reference is Acoustics Toolbox `v2023.5`, commit
`475108519289c6fb488b58980c644ea14eccc604`, using the same pinned Linux x86-64
GNU Fortran 12.2 environment as the BELLHOP reference.

## Products

The Rust implementation will cover the two-dimensional normal-mode workflow:

1. `KRAKEN` computes real-eigenvalue normal modes.
2. `KRAKENC` computes complex modes.
3. `FIELD` synthesizes complex frequency-domain pressure from a mode set.

The result model exposes both mode data and pressure fields. Modes contain the
per-frequency modal wavenumber, attenuation, phase/group speeds, and sampled
mode shapes. Pressure fields contain receiver coordinates and complex pressure.
Rust output will use a versioned HDF5 schema rather than the Fortran `.mod` and
`.shd` binary formats. Those Fortran files remain reference inputs to the
pinned differential tests.

This is not a ray-arrival model: BELLHOP-style discrete arrivals are out of
scope. A wideband time-domain response would be a separate product and
acceptance contract.

## Input and environment support

The planned legacy adapters accept KRAKEN `.env` files and FIELD `.flp` files,
including same-stem auxiliary resources used by the selected cases. The modern
adapter will provide a strict, self-contained JSON case, following BELLHOP's
single-document input convention.

The v2023.5 KRAKEN environment reader supports these SSP interpolation options:

- `N`: N²-linear
- `C`: C-linear
- `P`: PCHIP
- `S`: cubic spline
- `A`: analytic profile

Unlike BELLHOP, KRAKEN's normal-mode profile is range-independent; range
variation is represented as a sequence of profiles for FIELD propagation.
The planned environment support includes fluid and elastic layers, the
reference's attenuation units (`N`, `F`, `M`, `m`, `W`, `Q`, `L`), volume
attenuation (`T`, `F`, `B`), multiple frequencies, and the top/bottom boundary
conditions implemented by the reference (`V`, `R`, `A`, `F`, `P`). The boundary
codes cover vacuum, rigid, half-space, tabulated reflection, and precomputed
impedance paths. Existing reflection tables are inputs; table generation by
the separate `BOUNCE` program is not part of this rewrite.

The input validator will reject unsupported or inconsistent combinations with
structured diagnostics; it will not silently substitute a different solver
path.

## FIELD support

The planned 2D FIELD behavior includes:

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

The first numerical slice is a single-frequency, non-attenuating fluid Pekeris
waveguide. It must parse and validate its legacy input, compute modes, synthesize
a range-independent coherent field, and compare both products with pinned
Fortran. This is an implementation slice, not a reduction of the final support
matrix.

The Pekeris input will be constructed for this project: upstream
`tests/PekerisRD` exercises BELLHOP rather than KRAKEN. Additional official
v2023.5 cases provide feature coverage:

| Feature | Reference cases |
|---|---|
| Classic modes and FIELD | `tests/Munk/MunkK.env` + `.flp`, `tests/sduct/sductK.env` + `.flp` |
| Complex/leaky modes | `tests/MunkLeaky/MunkKwb`, `MunkKbb`, `MunkKleaky` |
| Multi-frequency | `tests/BroadBand/MunkK` |
| Adiabatic and coupled FIELD | `tests/Gulf/gulf_ad.flp`, `gulf_cm.flp` |
| Reflection inputs | `tests/TabRefCoef/neggradK_*` |

Acceptance requires pinned differential coverage for modal wavenumbers and
attenuation, normalized/aligned mode shapes, and complex pressure-field samples.
Mode-shape comparisons must account for the arbitrary sign/phase convention of
eigenvectors. Small committed goldens will keep ordinary tests independent of
Docker; the pinned reference workflow will cover official end-to-end cases.

## Repository shape

The implementation will follow the existing BELLHOP boundaries: one
`crates/kraken` library for model, validated case construction, input adapters,
mode solvers, and FIELD; separate CLI and HDF5 adapters when those outputs are
ready. No shared acoustics abstraction is planned until real duplication
justifies one.
