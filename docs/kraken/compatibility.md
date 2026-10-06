# KRAKEN/KRAKENC compatibility target

The numerical reference is Acoustics Toolbox `v2023.5`, commit
`475108519289c6fb488b58980c644ea14eccc604`, using the pinned Linux x86-64 GNU
Fortran 12.2 environment described in [`reference.md`](../development/reference.md).

**Status:** `crates/kraken` supports range-independent, layered-fluid
trapped or confined modes with material/volume and fluid-half-space attenuation.
Range-independent FIELD supports line, point and scaled-cylindrical sources,
omnidirectional or tabulated patterns, and coherent or incoherent mode addition.
Multi-profile adiabatic and smooth-fluid coupled FIELD now cover both original
Gulf paths; see [profile propagation](field.md).
KRAKENC computes trapped and leaky modes and the same FIELD options for complex
`N/C/P/S` or lossless fixed-Munk-A water.
Both backends accept smooth V/R/A top and bottom combinations, with bounded
Richardson refinement; KRAKEN half-spaces must remain trapped. Legacy material
units N/M/m/F/W/Q/L and added T/F/B volume attenuation are supported as detailed below. Unmodified
`MunkKleaky`, `MunkKwb`, `MunkKbb`, `sductK` and `calibK` `.env/.flp` pairs
pass end-to-end differential. Single-frequency, RMax=0 KRAKENC bottom F/BRC
and P/IRC also pass the original TabRefCoef geo/brc/irc workflows, including
FIELD endpoint extension and CLI/HDF5 auxiliary-input provenance.
Discrete multi-frequency KRAKEN/KRAKENC fluid-bottom cases
are also supported; unmodified `tests/BroadBand/MunkK.env/.flp` passes at 50
and 500 Hz through **both engines**, including legacy/JSON CLI-HDF5, within
unchanged work ceilings. `pelagic kraken` runs supported legacy pairs and writes
[KRAKEN HDF5 schema v1](../guide/kraken-hdf5.md), with sequential frequency
output and atomic publication. [Strict self-contained JSON](json.md)
now imports/exports the same canonical cases and uses that CLI/HDF5 path; HTTP
is not provided. This is not
the full matrix below; unsupported physics is rejected. Numerical parity is
verified for the complete workflows below, not promised for every arbitrary
branch-sensitive secant spectrum.
[Elastic A half-spaces](elastic-halfspaces.md) additionally support
top/bottom boundaries through both engines over N/C/P/S fluid stacks.
Original TLslices scholte/normal/flused pass both engines, FIELD and CLI/HDF5.
Real elastic tops retain shared isolation/Brent/Solve2. Homogeneous finite solid
caps now pass the separate [finite-layer checkpoint](finite-elastic-layers.md):
both backends with contiguous layered fluids, including original elsed/ice.
KRAKEN retains pinned real-stiffness/loss rules, including the top-A compressional
perturbation but omitted elastic shear/bottom absorption; use KRAKENC for full elastic loss.

For code navigation and ownership, see the [implementation map](../development/kraken-internals.md).

## Products

The target is the two-dimensional normal-mode workflow:

1. `KRAKEN` computes real-eigenvalue normal modes.
2. `KRAKENC` computes complex modes.
3. `FIELD` synthesizes complex frequency-domain pressure from a mode set.

The Rust result model exposes modes and pressure fields, not BELLHOP-style ray
arrivals. Rust CLI output uses an independent versioned KRAKEN HDF5 schema,
not the Fortran `.mod` and `.shd` binary formats or BELLHOP schema v3. A wideband time-domain response would be a
separate product and acceptance contract.

## Current capabilities and limits

The supported input/physics slice is not the entire planned matrix below.
Each checkpoint retains its representative acceptance paths, not an inferred
Cartesian product of every option:

- Contiguous constant-density fluid stacks with independent N/C/P/S meshes,
  material/volume attenuation and smooth V/R/A boundaries. Fixed analytic Munk A
  remains single-layer and lossless. KRAKEN fluid A half-spaces stay trapped;
  KRAKENC also supports leaky intervals.
- Limited single-frequency KRAKENC N/C bottom F/BRC, P/IRC and top F/TRC
  paths; restrictions and exact resource formats remain in the
  [table/boundary checkpoints](../development/kraken-checkpoints.md#tabulated-krakenc-bottoms).
- [Elastic A half-spaces](elastic-halfspaces.md) and
  [homogeneous/depth-varying finite elastic caps](finite-elastic-layers.md),
  with the pinned real-KRAKEN stiffness/loss limitations preserved.
- [Single/multi-profile FIELD](field.md), with line/point/scaled sources,
  source patterns and coherent/incoherent addition within its published limits.
  Coupled/incoherent input and unsupported coupling grids are rejected.
- Ordered discrete frequencies, repetitions, strict [self-contained JSON](json.md)
  and actual [single-file HDF5 output](../guide/kraken-hdf5.md). Failure stops the
  sequence without publishing a successful prefix.

See [legacy syntax, precision and budgets](input.md) for the input contract.
Fluid density gradients, interleaved/pure solids, rough interfaces, unvalidated
real elastic-top/finite-solid combinations, FIELD3D, BOUNCE generation, ray
arrivals, HTTP and time-domain synthesis remain outside the accepted slice.
Unsupported combinations are rejected, not silently approximated. Successful
execution or staying under a budget is not an arbitrary-input parity certificate.

[Checkpoint history and complete evidence](../development/kraken-checkpoints.md)
retain the earlier slice boundaries, all numerical counts and historical probes.
Later checkpoints extend earlier slices without changing the original target
at `2a7a658:docs/kraken-compatibility.md`.

## Planned environment support

The planned legacy adapters accept KRAKEN `.env` and FIELD `.flp` files,
including same-stem auxiliary resources used by selected cases. The strict,
[self-contained JSON adapter](json.md) follows BELLHOP's
single-document convention, with an independent KRAKEN schema and canonical
solve-frequency materials.

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

The supported fluid slices compare modes and range-independent coherent
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
above, plus 24 derived smooth-boundary/TRC workflows and 43 derived layered-fluid
workflows, also have committed goldens and configured API/CLI-HDF5 fresh comparisons.
Original TLslices `double` additionally passes both engines with small committed
goldens; large official `.mod/.shd` outputs remain fresh-reference artifacts.

## Repository shape

The implementation follows the existing BELLHOP boundaries: `crates/kraken`
contains validated cases, legacy adapters, mode solving, and FIELD. The
`pelagic kraken` command in `cli` uses `output::kraken` for the supported legacy-to-HDF5
workflow; BELLHOP v3 and KRAKEN v1 schemas and result types remain independent.
The JSON adapter reuses those validated cases and output paths; HTTP is not
provided. No shared acoustics abstraction is introduced.
