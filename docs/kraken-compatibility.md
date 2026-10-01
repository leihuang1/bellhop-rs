# KRAKEN Rust rewrite: compatibility target

The numerical reference is Acoustics Toolbox `v2023.5`, commit
`475108519289c6fb488b58980c644ea14eccc604`, using the pinned Linux x86-64 GNU
Fortran 12.2 environment described in [`reference.md`](reference.md).

**Status:** `crates/kraken` supports range-independent, single-fluid-layer
trapped or rigid-bottom confined modes with water and bottom-half-space
attenuation, and coherent line- or point-source FIELD. A separate
KRAKENC slice computes trapped and leaky modes and coherent line-/point-source
FIELD for complex `N/C/P/S` or lossless fixed-Munk-A water, vacuum surface and
fluid bottom, with up to five Richardson-extrapolated meshes. Legacy material
units N/M/m/F/W/Q/L and added T/F/B volume attenuation are supported as detailed below. Unmodified
`MunkKleaky`, `MunkKwb`, `MunkKbb`, `sductK` and `calibK` `.env/.flp` pairs
pass end-to-end differential. Single-frequency, RMax=0 KRAKENC bottom F/BRC
and P/IRC also pass the original TabRefCoef geo/brc/irc workflows, including
FIELD endpoint extension and CLI/HDF5 auxiliary-input provenance.
Discrete multi-frequency KRAKEN/KRAKENC fluid-bottom cases
are also supported; unmodified `tests/BroadBand/MunkK.env/.flp` passes at 50
and 500 Hz through **KRAKEN**. KRAKENC 500 Hz still exceeds its 300M root-work
limit and is not accepted. The `kraken` CLI now runs supported legacy pairs and writes
[KRAKEN HDF5 schema v1](kraken-output-format.md), with sequential frequency
output and atomic publication. JSON/HTTP remain unimplemented. This is not
the full matrix below; unsupported inputs are rejected.

## Products

The target is the two-dimensional normal-mode workflow:

1. `KRAKEN` computes real-eigenvalue normal modes.
2. `KRAKENC` computes complex modes.
3. `FIELD` synthesizes complex frequency-domain pressure from a mode set.

The Rust result model exposes modes and pressure fields, not BELLHOP-style ray
arrivals. Rust CLI output uses an independent versioned KRAKEN HDF5 schema,
not the Fortran `.mod` and `.shd` binary formats or BELLHOP schema v3. A wideband time-domain response would be a
separate product and acceptance contract.

## Current slice

The KRAKEN `load_case`/`solve` path currently accepts a narrow legacy
`.env`/`.flp` subset:

- one frequency and one constant-density fluid water layer using
  `N` (N²-linear), `C` (sound-speed-linear), `P` (monotone PCHIP), `S`
  (not-a-knot cubic spline), or `A` (the fixed 5000 m analytic Munk profile
  from upstream `misc/munk.f90`); depth-varying sound speed is supported,
  but density gradients and additional layers are not;
- a smooth vacuum (`V`) or rigid (`R`) surface and either a smooth acoustic
  fluid bottom half-space (`A`) with material/volume attenuation, or a smooth
  rigid (`R`) bottom with no half-space record or material properties.
  N/C/P/S water supports nonnegative absorption; analytic A remains lossless.
  Elastic layers and rough boundaries are not supported;
- a coherent, omnidirectional line (`X`) or point (`R`) source (`O`, `C`),
  with one range-independent FIELD profile at 0 km;
- finite, in-water source/receiver depths. FIELD normally interpolates the
  `.env` modal samples; it may extend either sampled endpoint by at most
  `1500 / frequency_hz` metres using the nearest two complex32 samples.
  It does not insert a zero-depth sample or clamp to an endpoint. Single-sample
  mode grids cannot extend, and depths outside the water remain rejected.

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
reference's first-mesh perturbation for water and bottom attenuation. Analytic Munk
sampling preserves upstream's f32 intermediate `eps` and grid step even though
its sound speed is stored as f64. Point-source FIELD
uses the reference's `sqrt(k)` modal factor and cylindrical spreading. The
former closed-form Pekeris solver remains **only in tests** as an independent
analytical cross-check. Analytic `A` has no SSP point records, uses density 1,
and requires the upstream formula's 5000 m water column; it is not a general
configurable analytic profile. The parser rejects unaccepted boundary types,
source types, and FIELD options instead of silently changing meaning.

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
supports complex `N/C/P/S` or lossless fixed-Munk-A water, vacuum surface, fluid
half-space with material/volume attenuation (including a bottom slower than the water), trapped
or leaky phase-speed intervals, and non-negative `RMax`. The KRAKENC legacy
reader accepts both blank and dotted restart options. The `.env` supplies modal sample depths; `.flp`
supplies FIELD geometry, used by `solve` for coherent line-/point-source
pressure or ignored by `solve_complex_modes`. The derived
`PekerisComplex.env` extends Pekeris to `cHigh=2000 m/s` and enables Fortran
restarts; `PekerisComplexBlank` uses the default restart setting. A derived
`PekerisComplexRefined.env/.flp` lowers NG to 100 and sets `RMax=1000 km`;
the pinned solver uses mesh multipliers 1/2/4, and all four modes and nine
pressures match committed and fresh reference. As in Fortran, shapes and group
speeds come from the first mesh while complex squared wavenumbers are Richardson
extrapolated. `PekerisComplexSlow` lowers the bottom speed to 1400 m/s.
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
root-work limit is unchanged. Four other **unmodified upstream** pairs also
pass fresh `.mod/.prt/.shd` comparison in CI:

| Original input | KRAKENC modes | Complex FIELD pressures |
|---|---:|---:|
| `MunkKwb.env/.flp` (cHigh=1551.91 < bottom speed 1600 m/s) | 63 | 150,801 |
| `MunkKbb.env/.flp` (cHigh=bottom speed 1600 m/s) | 102 | 150,801 |
| `sductK.env/.flp` (`C`, cHigh=100000 m/s) | 1,077 | 201,201 |
| `calibK.env/.flp` (`C`, RMax=1000 km; meshes 1/2/4) | 33 | 101,101 |

Repeated pinned runs produce identical `.mod/.shd` binaries for each original
pair; the original artifacts stay out of Git. The sduct FIELD has 216,693,477
modal contributions, below the 550-million-per-frequency work ceiling. Water loss
remains unsupported. KRAKENC `P/S/A` profile evidence is detailed below. Upstream
`tests/PekerisRD` is a BELLHOP, not a KRAKEN, case. The same comparator in
[`tests/differential_reference.rs`](../crates/kraken/tests/differential_reference.rs)
checks committed `.mod/.shd/.prt` goldens and newly generated pinned Fortran
output in CI; provenance is recorded
[with the goldens](../crates/kraken/tests/fixtures/golden/README.md).

### KRAKENC cubic and analytic water

The earlier interpolation checkpoint reused the real/lossless `Profile` for
PCHIP (`P`), not-a-knot cubic spline (`S`) and the fixed 5000 m Munk formula (`A`).
Complex N/C/P/S water and additional loss units are added by the material
checkpoint below; analytic A remains lossless. The half-space slice retains
vacuum top, constant density and A bottom, blank/dotted restarts, bounded Richardson refinement and
ordered discrete frequencies. F/P table boundaries remain **N/C only**;
complex rigid boundaries, density gradients and extra layers remain excluded.
Spline positivity and coefficient-overflow checks still run at `CaseDefinition`
validation. Both mode backends share analytic mesh sampling: the f32 step from
`5000.0 / N` is promoted to f64 separately from the physical finite-difference
step; endpoint clamping must not erase this reference rounding.

All modes, shapes, attenuation, speeds and FIELD pressures are compared at
unchanged tolerances, including actual CLI-HDF5 readback:

| Derived input | Modes | Pressures | Reference mesh multipliers |
|---|---:|---:|---|
| `MunkLeakyPartialP` / `MunkLeakyPartialS` | 30 each (9 leaky) | 36 each | 1 |
| `PekerisComplexSpline3` | 4 (1 leaky) | 9 | 1/2/4 |
| `MunkAnalyticComplex` | 102 | 25 | 1/2/4 |
| `MunkLeakyPchipBroadband` at 75/50/62.5/50 Hz | 45/30/37/30 | 36 per block | 1/2/4 |

The broadband derivative has NG=803; scaling is performed before truncation,
with the repeated 50 Hz result stored in its own frequency-index group.
The analytic derivative has NG=2003 (noninteger analytic grid step), RMax=1000 km
and eleven modal sample depths. These are explicitly derived, reduced-geometry
inputs; small raw goldens are committed and all seven single/multi-frequency
reference workflows produce identical `.mod/.shd` in three pinned runs.

The **unmodified upstream environments** `tests/Munk/MunkS.env` (S, originally
a SCOOTER fixture) and `MunkAnalytic.env` (A) are additionally run with KRAKENC:
102 modes and 25 pressures each, using a separately derived coherent `.flp`.
Neither original three-line `.flp` parses in pinned FIELD, so this is original
**environment** acceptance, not acceptance of the original input pairs or a
SCOOTER implementation. Hashes and provenance are
[with the goldens](../crates/kraken/tests/fixtures/golden/README.md).

Unmodified BroadBand/MunkK at 50 Hz also passes a local KRAKENC differential;
500 Hz hits the unchanged 300M root-work ceiling. CI explicitly checks CLI
exit 3 after the first frequency, preserves an existing output and cleans
scratch; no partial/truncated result is accepted. Full 50/500 Hz upstream
BroadBand/MunkK remains accepted through KRAKEN, not KRAKENC.

### Single-fluid material attenuation

N/C/P/S SSP nodes accept nonnegative compressional absorption. The legacy
loader converts N (neper/m), M (dB/m), F (dB/m-kHz), W (dB/wavelength), Q
(quality factor; Q=0 means no loss), L (loss parameter), and m (dB/m at freq0
with a power law and transition frequency). For m, the water header adds beta/fT
and the bottom-options record adds its own beta/fT; require beta >= 0 and fT > 0.
The fourth option character adds Thorp T, Francois–Garrison F (T, salinity, pH,
mean depth), or biological B (1..200 finite layers with positive resonance/Q).
Legacy trailing-value inheritance applies to node absorption as to sound speed.

`CaseDefinition::water_attenuation_db_per_wavelength` is empty for lossless
water or holds one effective value per SSP node **at the case frequency**.
The existing bottom attenuation field uses the same canonical unit. Both
backends share complex SSP interpolation: convert losses at the supplied nodes
before interpolation, use complex N², and interpolate P/S real and imaginary
parts independently. Reject nonfinite coefficients and cubic negative/gain
undershoot, or Im(c) > Re(c). KRAKEN uses Re(omega²/c²) for real roots and the
first-mesh volume/boundary loss perturbation; KRAKENC uses complex coefficients
through shooting, inverse iteration and normalization. Refinement, frequency
order/repetition, mode limits, work ceilings and output schema are unchanged.

Volume loss is recomputed per frequency and sampled **at SSP nodes**, not at
each finite-difference node. Thorp/FG also apply to a fluid half-space; biological
layers do not (`UpdateHSLoss` passes HUGE(depth)). Analytic Munk rejects volume
addition and has no lossy point records. F/P bottoms remain lossless N/C only;
this checkpoint does not expand table combinations or complex rigid boundaries.
Density gradients remain rejected: pinned acoustic shooting/vector/normalization
homogenize density using the top of each medium, not a physical density gradient.

Fourteen explicitly **derived** input pairs produce 27 accepted reference
workflows (13 through both backends, leaky broadband through KRAKENC only),
36 frequency blocks, 404 modes and 432 pressures. All `.mod/.shd` outputs are
byte-identical in three pinned runs. Offline goldens, API and actual CLI-HDF5
comparisons use existing tolerances; local maximum pressure error is 3.34e-8.
`WaterLossN/C/P/S` cover depth-varying complex SSPs; `WaterLossUnitN/M/F/Q/L`
cover additional material units; `WaterLossPower` covers distinct water/bottom
power laws, a transition inside 75/50/62.5/50 Hz and fractional NG scaling;
`WaterLossThorp/Fg/Bio` cover volume terms at 500 Hz, including overlapping
biological layers; `WaterLossLeaky` combines PCHIP water loss, bottom loss,
leaky modes, refinement and repeated frequencies. These are **not** acceptance
of original upstream VolAtt input pairs. Fresh CI is configured to check the
same 27 workflows through the API and actual HDF5 output.

### Tabulated KRAKENC bottoms

Bottom `F` consumes the same-stem `.brc`; `P` consumes `.irc`. This slice
requires **one frequency, RMax=0, lossless constant-density N/C single water
layer, vacuum top, smooth bottom and blank restart option**. The `B` frequency
option (even with one value), table mesh refinement, random restarts, real
KRAKEN F/P, top TRC and elastic/multilayer propagation are rejected or remain
outside this slice. There is no half-space record/material for F/P; direct
`CaseDefinition` half-space speed, density and attenuation must be zero.

Both formats require 2–100,000 ordered entries, bounded by the 1 MiB UTF-8
input limit. BRC has a count followed by one `angle magnitude phase_degrees`
record per point: finite increasing angles in 0..90 degrees, finite
nonnegative magnitudes, and finite **already unwrapped** phases. Interpolation
is linear in magnitude/phase, not complex reflection. The grazing angle uses
the real parts of complex kx/kz and reference f32 bracket selection; reflection
is zero outside the angle interval. BRC contains no frequency field; the caller
must supply a table appropriate to the single solve frequency.

IRC has a quoted title/frequency header, count and ASCII fixed-width
`(5G15.7,I5)` rows: real k², complex f, complex g, decimal scaling power.
E/D and letterless three-digit exponents are supported, not arbitrary
whitespace-delimited rewrites. The header frequency must exactly equal the
solve frequency. Real k² is finite, nonnegative and strictly increasing;
f/g are finite, powers lie in -1000..1000, and every adjacent two/three-point
window spans at most 100 powers. Brackets use real k²; up to three rescaled
points form a polynomial evaluated at complex k². Outside the real domain,
the corresponding endpoint f/g/power is retained.

Shooting, inverse iteration and normalization share the actual `(f,g,power)`
boundary; no fictitious half-space is substituted. F/P retains the f/g
normalization derivative but has no A-half-space group-speed contribution.
Critical coefficient/division and real-component grouping retain the pinned
GNU Fortran `-ffast-math` order: single-ulp differences can select a different
IRC secant root. Root work still uses the existing 20,000 roots / 300M operations
ceilings; the water-only root estimate is not a bound for table boundaries.
Singular/nonfinite boundary evaluations or normalization return diagnostics.
Reaching a lower spectral endpoint is not proof of every mathematical root.

`tools/reference/prepare-tabref.sh` copies unmodified upstream inputs and runs
**source-rebuilt pinned BOUNCE** three times to generate tables, then KRAKENC
and FIELD. The upstream directory contains no table files: these are generated
reference resources, not upstream original `.brc/.irc`. Tables and `.mod/.shd`
are identical across runs; Rust and actual CLI-HDF5 readback compare all modes,
shapes, speeds, attenuation and pressures at unchanged tolerances:

| Original TabRefCoef pair | Modes | Complex FIELD pressures |
|---|---:|---:|
| `neggradC_geo.env/.flp` | 56 | 50,601 |
| `neggradC_brc.env/.flp` + generated BRC | 54 | 50,601 |
| `neggradC_irc.env/.flp` + generated IRC | 42 | 50,601 |

Their `.env` samples 1..100 m while original `.flp` requests 0..100 m; neither
input is modified. BRC/IRC/geo results differ and are checked against their
own reference, not forced to agree. Four small **derived** 50 Hz N/C cases
with **constructed** three-point tables additionally cover line-/point-source
FIELD and nonzero IRC powers in committed goldens and fresh CI. Hashes and
source/compiler records are [with the goldens](../crates/kraken/tests/fixtures/golden/README.md).
Rust BOUNCE generation remains excluded. CLI schema-v1 additive metadata records
the exact consumed table snapshot alongside `.env/.flp`.

### Multiple frequencies

`legacy::load_frequency_cases(env, flp, ModeSolver)` returns a `Vec<Case>` in
input frequency order. Each case uses the existing `solve` (or KRAKENC
`solve_complex_modes`) and produces a separately labelled frequency-domain
result. `load_case` and `load_complex_case` reject multiple frequencies rather
than silently taking the first. The CLI/HDF5 adapter processes these cases
sequentially with cumulative output quotas; it adds no batched numerical solver,
FFT, time-domain response, JSON or HTTP endpoint.

The sixth top-option character `B` enables a frequency count/vector after
`.env` source/receiver depths. Frequencies are finite and positive; duplicates
and descending order are preserved, unlike sorted depth/range vectors. Slash-
terminated endpoint subtabulation uses double precision. The nominal `freq0`
remains `CaseDefinition::mesh_reference_frequency_hz = Some(freq0)`; `None`
uses the single case frequency. Automatic NG is resolved at freq0, then each
mesh uses `INT(NG * multiplier * frequency / freq0)`, not an already-rounded
base mesh multiplied afterward. For example NG=101 at freq0=50 Hz gives
151/303/606 intervals at 75 Hz, and 126/252/505 at 62.5 Hz. Existing `N/W`
loss semantics remain unchanged; modes and material/volume loss are recomputed
at every frequency, not obtained by scaling a previous mode set.

Two **derived** Pekeris broadband pairs have 75/50/62.5 Hz in that order,
NG=101 and RMax=1000 km: KRAKEN has 5/3/4 modes, KRAKENC has 7/4/6, with
nine line-source pressures per frequency. They pass committed and fresh
`.mod/.prt/.shd` comparisons. The **unmodified upstream**
`tests/BroadBand/MunkK.env/.flp` uses automatic `S` water sampling and `W`
bottom loss at 50/500 Hz: 102/1,023 modes, 501,501 pressures per frequency,
1,003,002 pressures total. Three pinned runs produced identical `.mod/.shd`
binaries; fresh CI compares both complete frequency blocks. Original artifacts
are not committed and do not replace the different single-frequency
`tests/Munk/MunkK` pair.

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
  half-space record may similarly inherit omitted trailing values; a rigid or
  tabulated bottom has **no** half-space record. In a direct `CaseDefinition`,
  its bottom sound speed, density, and loss must all be zero (absent material).
- Fortran null slots and repetition syntax remain unsupported and are rejected;
  this is not a general Fortran list-directed reader. All parsed numbers must
  be finite, SSP depths strictly increase from zero to the interface, and
  semantic diagnostics retain input-file records.
- Each input file is capped at 1 MiB; vectors at 100,000 entries; frequency
  count at 1,000; cloned frequency-case input vectors at 5,000,000 values.
  Per frequency: mesh at 1,000,000 grid intervals; roots at 20,000 modes;
  mode shapes at 5,000,000 values; all KRAKEN mesh searches at 2,500,000,000
  conservative operations and KRAKENC at 300,000,000 counted operations;
  pressure grids at 1,000,000 samples and 550,000,000 modal contributions.
  The larger KRAKEN/FIELD work bounds cover original BroadBand/MunkK at
  500 Hz (2,273,677,857 conservative root operations, 513,035,523 FIELD
  contributions); numerical tolerances are unchanged. The loader bounds
  cumulative input copies. The CLI additionally bounds cumulative output
  payload/file size (default 256 MiB), solving and writing one frequency at
  a time. Numerical work limits apply to each `solve`, not cumulatively
  across separate calls; the output quota is not a CPU timeout.
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
| Complex mesh extrapolation | `tests/calib/calibK.env` + `.flp` |
| Multi-frequency | `tests/BroadBand/MunkK` |
| Adiabatic and coupled FIELD | `tests/Gulf/gulf_ad.flp`, `gulf_cm.flp` |
| Reflection inputs (KRAKENC) | `tests/TabRefCoef/neggradC_*` + BOUNCE-generated tables |

Acceptance requires differential coverage for modal wavenumbers and
attenuation, normalized/aligned mode shapes, and complex pressure-field samples.
Mode-shape comparisons account for the arbitrary sign/phase convention of
eigenvectors. Small committed goldens keep ordinary tests independent of Docker;
the pinned reference workflow numerically compares all seventeen constructed
single-frequency KRAKEN fixtures, ten derived single-frequency KRAKENC mode
fixtures, four earlier derived single-frequency KRAKENC FIELD fixtures, four
cubic/analytic mode/FIELD derivatives, four table modes/FIELD fixtures, two
broadband Pekeris derivatives and a broadband PCHIP Munk derivative, and unmodified
upstream MunkK, MunkKleaky, MunkKwb, MunkKbb, sductK, calibK, BroadBand/MunkK
and all three neggradC geo/brc/irc modes and FIELD. Unmodified MunkS and
MunkAnalytic environments are additionally checked through KRAKENC with
separately derived FIELD geometry. The 27 derived water-material workflows
above also have committed goldens and configured API/CLI-HDF5 fresh comparisons.
Original `.mod/.shd` output is not committed.

## Repository shape

The implementation follows the existing BELLHOP boundaries: `crates/kraken`
contains validated cases, legacy adapters, mode solving, and FIELD. The
`kraken-cli` binary uses `bellhop-hdf5::kraken` for the supported legacy-to-HDF5
workflow; BELLHOP v3 and KRAKEN v1 schemas and result types remain independent.
JSON and HTTP adapters are still planned. No shared acoustics abstraction is
introduced.
