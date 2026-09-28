# KRAKEN Rust rewrite: compatibility target

The numerical reference is Acoustics Toolbox `v2023.5`, commit
`475108519289c6fb488b58980c644ea14eccc604`, using the pinned Linux x86-64 GNU
Fortran 12.2 environment described in [`reference.md`](reference.md).

**Status:** `crates/kraken` supports range-independent, single-fluid-layer
trapped or rigid-bottom confined modes with optional bottom-half-space `W`
attenuation, and coherent
line- or point-source FIELD. It is not the full matrix below; unsupported
inputs are rejected.

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

- one frequency and one lossless, constant-density fluid water layer using
  `N` (N²-linear), `C` (sound-speed-linear), `P` (monotone PCHIP), or `S`
  (not-a-knot cubic spline) interpolation; depth-varying sound speed is
  supported, but density gradients and additional layers are not;
- a smooth vacuum (`V`) or rigid (`R`) surface and either a smooth acoustic
  fluid bottom half-space (`A`) with zero loss (`N`) or non-negative `W` bottom
  loss (dB/wavelength), or a smooth rigid (`R`) bottom with no half-space record
  or material properties. Water-column loss, other loss units, elastic layers,
  and rough boundaries are not supported;
- a coherent, omnidirectional line (`X`) or point (`R`) source (`O`, `C`),
  with one range-independent FIELD profile at 0 km;
- finite, in-water source/receiver depths, with FIELD depths covered by the
  mode-sampling depths in `.env`.

The solver uses a real acoustic finite-difference mesh with Sturm mode counts,
inverse iteration for mode shapes and group speeds on the first mesh, and up to
five Richardson-extrapolated eigenvalue meshes controlled by `RMax`. A rigid
surface includes the half-weight surface node in the Sturm and inverse-iteration
matrices; a vacuum surface fixes its pressure to zero. A rigid bottom uses a
half-weight terminal node with no half-space contribution; a fluid bottom
includes the decaying half-space contribution to normalization and the pinned
reference's first-order perturbation for bottom attenuation. Point-source FIELD
uses the reference's `sqrt(k)` modal factor and cylindrical spreading. The
former closed-form Pekeris solver remains **only in tests** as an independent
analytical cross-check. The parser rejects other boundary types, water loss,
other loss units, analytic `A` interpolation, source types, and FIELD options
instead of silently changing meaning.

Fifteen constructed fixture pairs cover Pekeris (including a three-point
spline with forced extrapolation, an alternate lossy bottom, derived rigid
surfaces with/without bottom loss, and two derived `S` waveguides with a rigid
bottom), lossless
`N` Munk, lossy-bottom `N` Munk with point-source FIELD, and derived trapped
`C/P/S` sduct profiles. The PCHIP sduct derivative also has bottom `W` loss.
**Additionally, the unmodified upstream `tests/Munk/MunkK.env` and `.flp` are
compared in CI:** 102 modes and 501,501 complex pressures. The committed
`MunkBottomLoss` golden is a *reduced-grid derivative*, not the unmodified
input. The sduct derivatives still remove the original leaky phase-speed
interval; unmodified `sductK` is not supported.
Original KRAKENC examples remain reference-only smoke tests. Upstream
`tests/PekerisRD` is a BELLHOP, not a KRAKEN, case. The same comparator in
[`tests/differential_reference.rs`](../crates/kraken/tests/differential_reference.rs)
checks committed `.mod/.shd/.prt` goldens and newly generated pinned Fortran
output in CI; provenance is recorded
[with the goldens](../crates/kraken/tests/fixtures/golden/README.md).

### Legacy syntax and limits in this slice

- Quotes, comments, comma separators, `D` exponents, explicit vectors spanning
  lines, and endpoint-subtabulated vectors (`count ≥ 3`, one or two endpoint
  values followed by `/`) are accepted. Depth vectors use upstream single
  precision; ranges and offsets use double precision. Like `ReadVector`, each
  legacy vector is sorted independently, including receiver offsets.
- Each SSP point occupies one record. Omitted trailing material values retain
  the previous point's values (Fortran defaults on the first). `/` ends that
  record, **not** the SSP; the interface depth ends the SSP. The acoustic bottom
  half-space record may similarly inherit omitted trailing values; a rigid
  bottom has **no** half-space record. In a direct `CaseDefinition`, rigid bottom
  sound speed, density, and loss must all be zero (absent material).
- Fortran null slots and repetition syntax remain unsupported and are rejected;
  this is not a general Fortran list-directed reader. All parsed numbers must
  be finite, SSP depths strictly increase from zero to the interface, and
  semantic diagnostics retain input-file records.
- Each input file is capped at 1 MiB; vectors at 100,000 entries; mesh at
  1,000,000 grid intervals; roots at 20,000 modes; mode shapes at 5,000,000
  values; all mesh searches at 300,000,000 counted operations; pressure grids
  at 1,000,000 samples and 60,000,000 modal contributions (covering full MunkK).
- `mesh_points` = 0 selects the reference's automatic base mesh (at least ten,
  ~20 points per wavelength); otherwise 10–1,000,000 is allowed if not too
  coarse. `max_range_m` = 0 uses the base mesh only; larger values control
  the reference-style extrapolation convergence criterion. If limits prevent
  convergence, the solver returns a diagnostic rather than an unverified mode.

Modes are computed in double precision and shapes stored at reference `.mod`
sampling precision. FIELD uses single-precision
wavenumbers, modal products, and accumulation, with separate range/offset phases,
following the reference's `.mod`/FIELD rounding points. Returned pressure values
are promoted to `Complex64`; that does not imply double-precision FIELD arithmetic.

This is an implementation slice, not a reduction of the final support target.
Measured differences and fixture provenance are recorded
[with the goldens](../crates/kraken/tests/fixtures/golden/README.md).

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
the pinned reference workflow numerically compares all fifteen supported
fixtures and unmodified official MunkK. Unmodified sduct remains a reference-only
smoke test; its leaky modes are not yet supported.

## Repository shape

The implementation follows the existing BELLHOP boundaries: `crates/kraken`
contains validated cases, legacy adapters, mode solving, and FIELD. CLI and HDF5
adapters can follow once the numerical path is stable. No shared acoustics
abstraction is planned until real duplication justifies one.
