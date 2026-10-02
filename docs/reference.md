# Numerical reference

The compatibility oracle is the official Acoustics Toolbox release `v2023.5`:

- Repository: <https://github.com/oalib-acoustics/Acoustics-Toolbox>
- Tag: `v2023.5`
- Commit: `475108519289c6fb488b58980c644ea14eccc604`
- License: GPL-3.0

Small BELLHOP numerical goldens are committed under
`crates/bellhop/tests/fixtures/golden`. Their per-file provenance, hashes,
compiler, target, and flags are recorded alongside them. KRAKEN has seventeen
constructed single-frequency fluid cases (twelve lossless, five with bottom loss) and raw
`.mod/.shd/.prt` goldens with
[provenance and hashes](../crates/kraken/tests/fixtures/golden/README.md).
`tools/reference/compare-kraken.sh CASE.env` regenerates reference outputs and
compares mode counts, coordinates, wavenumbers, attenuation, phase/group speeds,
phase-aligned shapes, and every complex pressure sample. The same comparator
runs without Docker on committed goldens; CI also compares fresh outputs for
all seventeen cases and the unmodified upstream `MunkK.env`/`.flp` (102 modes,
501,501 complex pressures). CI also compares the unmodified upstream
`MunkAnalytic.env` with a separately derived coherent FIELD `.flp` (102 modes,
25 pressures); its original three-line `.flp` fails in pinned v2023.5 FIELD.
The unmodified official sduct has leaky modes and is **not** covered by the
derived KRAKEN trapped-mode comparison; its separate KRAKENC original-pair
comparison is described below. Measured errors and provenance accompany the goldens.
Separate KRAKENC regressions compare the **derived**
`PekerisComplex.env`, its blank-restart and slower-fluid-bottom variants
(four modes each), a `cLow=1550 m/s` blank-restart variant (three in-band
modes), a three-mesh Richardson derivative (four modes and nine pressures),
two-point graded `N`-profiles in both directions (four modes each),
and a **derived**, truncated seven-knot MunkLeaky water column with and without
0.8 dB/wavelength bottom loss (30 modes, nine leaky each) against pinned
`.mod/.prt` goldens and fresh CI output. A third seven-knot derivative changes
only `N` to `C` interpolation (30 modes). Separately derived coherent `.flp`
inputs compare nine complex line-source pressures for `PekerisComplexBlank`
and 36 point-source pressures for each of the two lossy seven-knot cases with
committed `.shd` goldens and fresh CI output. Fresh CI compares unmodified
upstream `MunkKleaky` (329 modes, 150,801 pressures), `MunkKwb` (63, 150,801),
`MunkKbb` (102, 150,801), `sductK` (1,077, 201,201) and Richardson-refined
`calibK` (33, 101,101) `.env/.flp` pairs end to end against `.mod/.prt/.shd`.
The `.env` determines modal samples; the
`.flp` determines FIELD geometry. A metamorphic
frequency/depth-scaling test checks that small distinct squared-wavenumbers
remain separate; it does not claim a fresh Fortran acceptance for scaled
inputs.

KRAKENC cubic/analytic profile acceptance reuses the lossless `Profile`, adding
**derived** `MunkLeakyPartialP/S` (30 modes including 9 leaky, 36 pressures each),
`PekerisComplexSpline3` (4 modes, 9 pressures, meshes 1/2/4), and
`MunkAnalyticComplex` (102 modes, 25 pressures, NG=2003, eleven modal depths and
meshes 1/2/4). A **derived** broadband PCHIP Munk has NG=803, 75/50/62.5/50 Hz,
45/30/37/30 modes and 36 pressures per block. Three pinned runs of all five
produce identical `.mod/.shd`; raw goldens and input/output hashes are recorded.
Analytic f32 grid-step rounding is shared by both Rust mode backends.
Unmodified `tests/Munk/MunkS.env` (a SCOOTER environment) and `MunkAnalytic.env`
also pass through KRAKENC, 102 modes/25 pressures each, with a **separately derived**
coherent `.flp`. Original three-line `.flp` fails in pinned FIELD for both;
this is not original-pair or SCOOTER acceptance. The material checkpoint below
adds complex N/C/P/S water. The smooth-boundary checkpoint adds V/R/A at either
end; P/S/A table-boundary combinations remain excluded.

Water-material acceptance uses fourteen **derived** pairs: N/C/P/S complex
water, N/M/m/F/W/Q/L material units, T/F/B volume attenuation, overlapping
biological layers, and repeated/unsorted refined frequencies. Thirteen pairs
run through both KRAKEN and KRAKENC; the leaky PCHIP partial-Munk derivative
runs through KRAKENC only: 27 workflows, 36 frequency blocks, 404 modes and
432 pressures. Three pinned runs yield identical `.mod/.shd`; all samples pass
local API and actual CLI-HDF5 comparisons at unchanged tolerances. Raw goldens,
hashes and derivation details are committed; fresh CI is configured for all 27.
Biological attenuation is sampled at SSP knots, not mesh nodes, and excluded
from half-spaces as in `UpdateHSLoss` (HUGE depth). This is not original upstream
VolAtt acceptance and does not add physical density gradients or layers.

Discrete multi-frequency regressions use two **derived** Pekeris `.env/.flp`
pairs (KRAKEN and KRAKENC): unsorted 75/50/62.5 Hz, fractional mesh scaling,
three refinement levels, `W` bottom loss, and nine pressures per frequency.
Committed `.mod/.prt/.shd` and fresh CI compare every frequency block. Fresh
CI also compares the **unmodified** `tests/BroadBand/MunkK.env/.flp` through KRAKEN:
102 modes at 50 Hz, 1,023 at 500 Hz, and 501,501 pressures at each frequency
(1,003,002 total). Three pinned runs yield identical `.mod/.shd`. The reference
reader checks frequency order, per-frequency mode counts and offsets, shapes,
printed/binary wavenumbers, group speeds and every pressure; a later-frequency
corruption regression ensures no frequency is skipped. These are frequency-
domain results, not a wideband time-domain acceptance. KRAKENC S water at
500 Hz exceeds its unchanged 300M root-work limit for this original broadband
pair. A fresh-CI CLI regression checks exit 3, no partial publication, preservation
of an existing output and scratch cleanup after 50 Hz succeeds. A local 50 Hz
KRAKENC comparison passes, but the full pair is not accepted through KRAKENC.

Smooth-boundary acceptance adds eleven **derived** pairs: all nine V/R/A top/
bottom combinations with spline water loss, a lossy rigid-rigid plane mode with
NG=68, and a radiating air top through KRAKENC only (21 workflows). Three more
**constructed** top-TRC derivatives with N/C water, cLow=1550 and A/R bottoms
add three workflows; top tables require a propagating spectral window, no
refinement/B and no simultaneous bottom table. All 24 `.mod/.shd` workflows are
byte-identical in three pinned runs and have committed goldens, hashes and
API/actual-CLI-HDF5 differential. These are not original TRC input-pair acceptance.
Automatic NG=17 rigid-plane refinement changes its mode count and remains an
explicit numerical failure, not a truncated result. An evanescent top-TRC
experiment produced a growing reference root and is explicitly outside support.

TabRefCoef acceptance adds four **derived** 50 Hz N/C single-water cases with
**constructed** BRC/IRC tables and committed `.mod/.prt/.shd`. Original
`neggradC_geo`, `neggradC_brc`, `neggradC_irc` `.env/.flp` remain unchanged:
56/54/42 modes and 50,601 pressures each, including FIELD at 0 m outside the
1..100 m modal sample interval. The upstream directory contains no tables.
`tools/reference/prepare-tabref.sh` uses source-rebuilt pinned BOUNCE (the bundled
executable/object are removed) to generate `.brc/.irc` from original `neggradB.env`,
then runs KRAKENC/FIELD three times. Tables and `.mod/.shd` must be byte-identical;
input/resource/output hashes are retained. BRC/IRC/geo are compared to their own
reference, not assumed equivalent. Pinned GNU real-component grouping prevents
IRC secant from selecting a different subset of roots; this is not a proof of
mathematical root completeness. Table multi-frequency/refinement and Rust BOUNCE
remain excluded.

The [KRAKEN CLI/HDF5 adapter](kraken-output-format.md) is independently checked
through its actual serialized output: `KRAKEN_HDF5_RESULT` makes the same test-only
comparator read `.h5` modes and FIELD instead of recomputing the in-memory
result. Fresh CI exercises both Pekeris broadband derivatives, original MunkK,
all five original fluid KRAKENC pairs, original BroadBand/MunkK, four table
derivatives, the three original TabRefCoef workflows, four cubic/analytic
derivatives, the broadband PCHIP derivative and the two unmodified MunkS/analytic
environments with derived FIELD geometry through this path,
with every frequency/sample and unchanged tolerances. Schema/metadata/units,
exact Rust-value round-trips and file/quota failure protections also run without
Docker. Existing BELLHOP schema-v3 tests remain unchanged; Rust `.h5` output is
not a Fortran golden.

BELLHOP ray trajectories
cover all `N/C/P/S/Q/A` sound-speed models. Eigenray and arrival goldens cover
Cartesian and ray-centered geometric-hat beams, Cartesian geometric-Gaussian
beams, caustics, arrival combination, and multi-depth/multi-range receiver
grids. Pressure-field goldens cover coherent, semi-coherent, and incoherent
scaling, simple-Gaussian beams, and Cartesian/ray-centered Cerveny beams.
Dedicated reflection goldens cover acousto-elastic, grain-size, and `.irc`
impedance-table amplitude and phase.

Layered-fluid acceptance adds twenty-one **derived** pairs, twenty through both
engines and a leaky broadband pair through KRAKENC: 41 workflows, 50 frequency
blocks, 530 modes and 4,902 pressures. Each `.mod/.shd` is byte-identical in
three pinned runs. The comparator traverses every fluid medium record and every
frequency/mode/FIELD sample, through both API and actual CLI-HDF5; local maximum
pressure error is 2.64e-9 at unchanged tolerances. Offline goldens and layer
metadata/budget/failure regressions run without Docker.

`OriginalLayeredDouble.env` is byte-identical to upstream `tests/TLslices/double.env`;
its `.flp` is the official shared `fieldbat.flp` selected by upstream `runtests.m`.
Pinned KRAKEN changes from 43 to 42 modes on meshes 1/2; both Rust backends
reject mode-count changes instead of clipping. `LayeredDoubleRefined` doubles
NG to 200/400/400 and passes 42 modes/501 pressures per engine; it is derived,
not acceptance of original `double`. `LayeredNormalization` removes shear only
from the original `normal.env` bottom and is also explicitly derived. Original
`normal/flused` now pass the half-space checkpoint and `elsed/ice` the finite-cap
checkpoint below; Gulf's sequence now uses multi-profile
FIELD. None is claimed as accepted by the fluid-stack stage. The wider-spectrum
three-layer experiment described in the compatibility/golden provenance remains
unaccepted; a matched subset is not advertised as full original acceptance.

Elastic-half-space acceptance adds 15 derived and three original input pairs:
27 workflows, 33 frequency blocks, 415 modes and 4,707 pressures, through API and
actual CLI/HDF5. The original scholte/normal/flused environments and official
shared fieldbat.flp are byte copies of the pinned TLslices sources; no shear is
removed and NG is unchanged. All .mod/.shd agree in three runs; local maximum
pressure error is 6.67e-8 at unchanged tolerances. The 117-record
`golden/elastic-halfspace.sha256` manifest and fresh CI source-byte comparisons
lock provenance. KRAKEN elastic top-half-space and failed slow-interface probes remain outside
acceptance. Original elsed/ice now pass the distinct finite-cap checkpoint.
See the [scope and evidence](kraken-elastic-halfspaces.md); no new parity waiver.

Homogeneous finite elasticity adds 25 pairs (23 constructed/derived, two
unmodified originals): 33 workflows, 36 blocks, 827 modes, 9,276 pressures,
triplicate identical MOD/SHD, full API and actual CLI/HDF5 comparisons at unchanged
tolerances (local max |dp| 4.22e-8). Original elsed/ice and official fieldbat.flp
retain every byte and NG; the 149-record `golden/finite-elastic.sha256` locks new
provenance. KRAKENC supports contiguous layered fluids; KRAKEN finite solids
currently require one fluid. Multi-fluid KRAKEN failures are explicitly rejected,
not added to the existing parity exception. See [finite-cap scope](kraken-finite-elastic-layers.md).

The single-profile FIELD checkpoint adds three explicitly derived representative
pairs, without expanding them into an option Cartesian product: `FieldScaled`
through KRAKEN (3 modes/9 pressures), and `FieldPattern` / `FieldIncoherent`
through KRAKENC (4/9 each). The pattern case consumes a committed seven-point
`.sbp`. Three pinned runs produced byte-identical MOD/SHD for every pair; all
modes, shapes and pressures pass direct and actual CLI-HDF5 comparison at the
unchanged tolerances, with local maximum pressure error 3.34e-8. Inputs and
reference outputs are locked by `golden/single-profile-field.sha256`; Rust HDF5
is not a golden. Multi-profile FIELD is covered by the separate checkpoint below.

The multi-profile FIELD checkpoint covers two small four-profile derivatives
and both byte-original Gulf AD/CM workflows: 24 profile blocks, 1,020 modes and
1,003,050 pressures. Full API and actual CLI/HDF5 comparisons pass unchanged
tolerances, with local maximum |dp| `4.1159031748919954e-10`. All four MOD/SHD
pairs are byte-identical across three pinned runs. Seven inputs (including the
three original Gulf files) and six small reference files are SHA-locked in
`golden/multi-profile.sha256`; original Gulf binaries are fresh-run artifacts,
not committed large goldens. CI checks original input bytes against the pinned
source tree. See [profile propagation](kraken-multi-profile-field.md).

## Pinned Linux x86-64 oracle

`tools/reference` provides the reproducible differential environment:

- Debian `bookworm-slim` image index digest:
  `sha256:abd67ffcfa541b485a3dff59865ab629aa048a6c613e639d36e7456b0b229241`
- GNU Fortran `12.2.0`
- flags: `-O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu`
- Acoustics Toolbox source archive SHA-256:
  `f8a7a2c1e80a73431cd230a10bef5fcfc996c88889a0e1540771c3922ee2a21f`
- Rust `1.88.0` bookworm image index digest:
  `sha256:af306cfa71d987911a781c37b59d7d67d934f49684058f96cf72079c3626bfe0`

The semantic trajectory comparator aligns isolated minimum-step vertices.
When a value lies within a few ulps of an SSP or boundary interface, one
compiler can take the reference's `1e-4 × base step` while another changes
segment or reflects immediately. Aligned coordinates retain the `1e-5 m`
tolerance; this discrete branch allowance is measured and reported separately.

Pinned Linux comparisons currently show:

- ParaBot: 201 rays, maximum coordinate error `1.4e-11 m`, no branch allowance
- Ellipse: 72 rays, maximum coordinate error `2.2e-11 m`, no branch allowance
- block: 50 rays, exact coordinate agreement
- DickinsBray: 501 rays, aligned error below `9e-10 m`; 16 minimum-step alignment
  operations with maximum branch displacement `0.12 m`

The committed `DickinsCritical` and `ParaBotCritical` cases preserve focused
coverage of the affected SSP-interface and curved-boundary decisions.

## Arrival and pressure differential

`compare-arrival.sh` and `compare-field.sh` extend the same pinned-Linux
procedure to arrival and pressure outputs, parsing the reference `.arr` text
and the direct-access `.shd` records directly.

- `tests/SBCX/sbcx_Arr_asc.env` (official, unmodified): 200 beams over 100
  depths × 500 ranges = 50,000 receivers ≈ 502,000 arrivals. Double-precision
  receiver-grid interpolation changes four receiver counts by one at
  influence-window edges. Matching arrivals differ by at most `3.5e-7` in
  amplitude, `2.0e-5` degrees, and `2.0e-6` seconds. Phase agrees within one
  storage quantum: the reference stores phase as single-precision degrees
  (≈ `2.2e-6` rad grid at large phases) while Rust stores single-precision
  radians, so boundary flips up to `4.2e-6` rad are admitted and reported.
- `tests/Munk/MunkB_Coh.env` (official, unmodified): 501 depths × 501 ranges =
  501,501 coherent geometric-hat pressure samples, maximum absolute error
  `2.0e-10`.
- `MunkB_Semi.env`: `9.3e-10`; `MunkB_Inc.env`: `7.3e-12`;
  `MunkB_Coh_gb.env` (geometric-Gaussian): `7.0e-10`;
  `MunkB_Coh_CervenyC.env`: 100,701 samples, `1.5e-11`;
  `MunkB_Coh_CervenyR.env`: 100,701 samples, `1.8e-12`.
- The committed `ShadedField`, `FreeFGB_grid`, `Field_CervenyF`, and
  `Field_CervenyW` goldens were generated by the same container and match
exactly (source beam pattern, Francois-Garrison attenuation, and space-filling
and WKB Cerveny widths receive numerical coverage through `cargo test`).

All seven official two-dimensional `A/a` environments parsed by the Rust
loader execute with the implemented `G/B` influence models. Three committed
arrival comparisons match receiver-by-receiver arrival counts exactly and
agree at the expected single-precision storage points.

Small pressure grids match the reference single-precision values within
`5e-8` absolute pressure. Full-grid Munk coherent, semi-coherent, incoherent,
geometric-Gaussian, and Cerveny results match the pinned container within
`1e-9` absolute pressure.
