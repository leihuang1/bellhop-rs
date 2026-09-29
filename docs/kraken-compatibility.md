# KRAKEN Rust rewrite: compatibility target

The numerical reference is Acoustics Toolbox `v2023.5`, commit
`475108519289c6fb488b58980c644ea14eccc604`, using the pinned Linux x86-64 GNU
Fortran 12.2 environment described in [`reference.md`](reference.md).

**Status:** `crates/kraken` supports range-independent, single-fluid-layer
trapped or rigid-bottom confined modes with optional bottom-half-space `W`
attenuation, and coherent line- or point-source FIELD. A separate
KRAKENC slice computes trapped and leaky modes and coherent line-/point-source
FIELD for lossless `N/C` water, vacuum surface and fluid bottom with optional
`W` loss, on the base mesh without extrapolation. Unmodified `MunkKleaky`,
`MunkKwb`, `MunkKbb` and `sductK` `.env/.flp` pairs pass end-to-end differential.
It is not the full matrix below; unsupported inputs are rejected.

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

The KRAKEN `load_case`/`solve` path currently accepts a narrow legacy
`.env`/`.flp` subset:

- one frequency and one lossless, constant-density fluid water layer using
  `N` (N²-linear), `C` (sound-speed-linear), `P` (monotone PCHIP), `S`
  (not-a-knot cubic spline), or `A` (the fixed 5000 m analytic Munk profile
  from upstream `misc/munk.f90`); depth-varying sound speed is supported,
  but density gradients and additional layers are not;
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
half-weight terminal node with no half-space contribution. For a constant
profile between two rigid boundaries, the zero-order plane mode lies at the
upper spectral endpoint; the search includes the rounded finite-difference
endpoint. A fluid bottom includes the decaying half-space contribution to
normalization and the pinned
reference's first-order perturbation for bottom attenuation. Analytic Munk
sampling preserves upstream's f32 intermediate `eps` and grid step even though
its sound speed is stored as f64. Point-source FIELD
uses the reference's `sqrt(k)` modal factor and cylindrical spreading. The
former closed-form Pekeris solver remains **only in tests** as an independent
analytical cross-check. Analytic `A` has no SSP point records, uses density 1,
and requires the upstream formula's 5000 m water column; it is not a general
configurable analytic profile. The parser rejects other boundary types, water
loss, other loss units, source types, and FIELD options instead of silently
changing meaning.

Seventeen constructed fixture pairs cover Pekeris (including a three-point
spline with forced extrapolation, an alternate lossy bottom, derived rigid
surfaces with/without bottom loss, two derived `S` waveguides with a rigid
bottom, and a constant `S` rigid-rigid plane mode at the spectral endpoint), lossless
`N` Munk, lossy-bottom `N` Munk with point-source FIELD, a derived analytic
`A` Munk with fixed-depth formula and `W` bottom loss, and derived trapped
`C/P/S` sduct profiles. The PCHIP sduct derivative also has bottom `W` loss.
**Additionally, the unmodified upstream `tests/Munk/MunkK.env` and `.flp` are
compared in CI:** 102 modes and 501,501 complex pressures. The committed
`MunkBottomLoss` golden is a *reduced-grid derivative*, not the unmodified
input. The analytic `A` golden is also derived. CI additionally compares the **unmodified
upstream `MunkAnalytic.env`** with its automatic mesh and extrapolation (102 modes)
and the separately derived coherent `.flp` (25 pressures). The upstream
three-line `MunkAnalytic.flp` cannot be parsed by v2023.5 FIELD itself, so the
original `.env/.flp` *pair* is not claimed as accepted. The KRAKEN sduct
derivatives still remove the original leaky phase-speed interval; **unmodified**
`sductK` is instead compared through KRAKENC below.
The separate KRAKENC `load_complex_case`/`solve_complex_modes` path currently
supports lossless `N/C`-profile water, vacuum surface, fluid half-space with
zero or `W` attenuation (including a bottom slower than the water), trapped
or leaky phase-speed intervals, and `RMax=0`. The KRAKENC legacy reader accepts both blank
and dotted restart options. The `.env` supplies modal sample depths; `.flp`
supplies FIELD geometry, used by `solve` for coherent line-/point-source
pressure or ignored by `solve_complex_modes`. The derived
`PekerisComplex.env` extends Pekeris to `cHigh=2000 m/s` and enables Fortran
restarts; `PekerisComplexBlank` uses the default restart setting and
`PekerisComplexSlow` lowers the bottom speed to 1400 m/s.
`PekerisComplexCLow` uses a 1550 m/s lower phase-speed bound and blank
restarts, returning only the three in-band modes (the excluded first root is
still used for deflation). A separate **derived** `PekerisComplexBlank.flp`
verifies nine coherent line-source pressures against committed and fresh
v2023.5 `.shd`; all modes in these cases are compared against committed and
fresh `.mod/.prt`. Scaled depth/frequency regression
checks that duplicate-root detection retains distinct modes. Fluid-bottom
sound speed must be positive for both KRAKEN and KRAKENC. The derived
`PekerisComplexGradient` and `PekerisComplexReverseGradient` change the
water SSP endpoints to 1500→1520 and 1520→1500 m/s respectively, exercising
nonconstant N²-linear water with four modes each, including one leaky mode,
against the same committed and fresh modal comparator. The **derived**
`MunkLeakyPartial` retains the first seven water SSP knots from the original
`MunkKleaky`, truncated at 1000 m: with a lossless bottom, narrower leaky
interval, 800-point mesh and six modal sampling depths, all 30 modes (nine
leaky) match committed and fresh `.mod/.prt`. `MunkLeakyPartialLoss` adds the
original 0.8 dB/wavelength `W` bottom loss; its 30 modes also pass committed
and fresh differential. Complex secant carries the shooting/deflation scaling
exponent into adjacent function comparisons. Counted root work and complex
inverse iteration then cover the **unmodified** `MunkKleaky.env/.flp`:
all 329 modes (including leaky modes) and all 150,801 complex FIELD pressure
samples pass fresh pinned `.mod/.prt/.shd` comparisons. A separate coherent
point-source `.flp` for the **derived** `MunkLeakyPartialLoss` exercises source/
receiver interpolation, positive range offsets and zero range (36 pressures)
against a committed `.shd` golden and fresh CI output. Another **derived**
`MunkLeakyPartialC.env/.flp` changes only that input's water interpolation
from `N` to `C`: all 30 modes and 36 complex pressures match committed and
fresh reference. For wide spectral intervals, complex secant seeds subsequent
roots from the previous root, as in pinned Fortran; the 300-million counted
root-work limit is unchanged. Three other **unmodified upstream** pairs also
pass fresh `.mod/.prt/.shd` comparison in CI:

| Original input | KRAKENC modes | Complex FIELD pressures |
|---|---:|---:|
| `MunkKwb.env/.flp` (cHigh=1551.91 < bottom speed 1600 m/s) | 63 | 150,801 |
| `MunkKbb.env/.flp` (cHigh=bottom speed 1600 m/s) | 102 | 150,801 |
| `sductK.env/.flp` (`C`, cHigh=100000 m/s) | 1,077 | 201,201 |

Repeated pinned runs produce identical `.mod/.shd` binaries for each original
pair; the original artifacts stay out of Git. The sduct FIELD has 216,693,477
modal contributions, below the measured 250-million-work ceiling. Water loss,
`P/S/A` complex interpolation and refined KRAKENC meshes remain unsupported. Upstream `tests/PekerisRD` is a BELLHOP, not a KRAKEN, case. The same comparator in
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
- Each `N/C/P/S` SSP point occupies one record; `A` has no SSP point records
  and requires an explicit bottom sound speed and density for a fluid bottom.
  Omitted trailing material values retain the previous point's values (Fortran
  defaults on the first). `/` ends that
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
  at 1,000,000 samples and 250,000,000 modal contributions (covering original
  sductK).
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

The supported single-fluid slices compare modes and range-independent coherent
fields against pinned Fortran. Further representative v2023.5 cases provide
feature coverage:

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
the pinned reference workflow numerically compares all seventeen constructed
KRAKEN fixtures, nine derived KRAKENC modes, three derived KRAKENC FIELD
fixtures, and unmodified upstream MunkK, MunkKleaky, MunkKwb, MunkKbb and
sductK modes and FIELD. Original `.mod/.shd` output is not committed.

## Repository shape

The implementation follows the existing BELLHOP boundaries: `crates/kraken`
contains validated cases, legacy adapters, mode solving, and FIELD. CLI and HDF5
adapters can follow once the numerical path is stable. No shared acoustics
abstraction is planned until real duplication justifies one.
