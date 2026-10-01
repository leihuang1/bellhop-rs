# KRAKEN/KRAKENC CLI and HDF5 schema v1

This is a product adapter for the **currently accepted numerical subset** in
[kraken-compatibility.md](kraken-compatibility.md), not completion of its planned
matrix. It supports legacy `.env/.flp` inputs, KRAKEN or KRAKENC, and discrete
single/multiple-frequency modes plus coherent FIELD. JSON, HTTP, modes-only CLI
runs, reflection-table generation, ray arrivals and time-domain synthesis are
not provided. Original BroadBand/MunkK is accepted through KRAKEN; KRAKENC
now supports its S profile but still rejects 500 Hz at the unchanged work ceiling.

## CLI

```console
cargo run --release -p kraken-cli -- run path/to/case.env --output result.h5
cargo run --release -p kraken-cli -- run path/to/case.env --solver krakenc --flp geometry.flp --output complex.h5
```

The installed binary is `kraken`; `cargo install --path crates/kraken-cli`
builds it using the existing static HDF5 dependency.

- `--solver kraken|krakenc` defaults to `kraken`; the engine is not guessed.
- Both backends accept smooth V/R/A boundaries; real KRAKEN stays trapped
  at each A half-space. KRAKENC supports complex N/C/P/S or lossless fixed-Munk-A,
  with N/M/m/F/W/Q/L material units, optional T/F/B volume loss, and the existing
  bounded refinement/frequency-order contract. Analytic A remains lossless and
  rejects volume addition. Loss conversion is per solve frequency; the HDF5
  schema and exact-input provenance contract are unchanged. Original BroadBand/MunkK
  at 500 Hz still exceeds KRAKENC's 300M root-work ceiling: exit 3 preserves the
  previous output, even after 50 Hz was written to scratch. KRAKEN acceptance
  of the same pair is unchanged.
- KRAKENC `F`/`P` bottoms additionally consume same-stem `.brc`/`.irc`; this
  slice is single-frequency, RMax=0, lossless N/C water and vacuum top. Missing,
  malformed or unaccepted table combinations fail explicitly. Rust does not
  generate tables or infer an IRC solve frequency from its header.
- KRAKENC top F consumes `.trc`: single frequency, RMax=0, lossless N/C,
  blank restart/no B, smooth bottom and cLow >= last SSP-node speed.
  Evanescent top-table roots, top P and simultaneous top/bottom tables are rejected.
- `.flp` defaults to the environment's same-stem file. It remains required:
  `.env` controls modal samples; `.flp` controls FIELD geometry and mode limit.
- Output defaults to `<case-stem>.h5` in the current directory. Its parent
  directory must already exist.
- Existing destinations require `--overwrite`. Input paths and symlink aliases
  of inputs are rejected as destinations even with that flag.
- Exit codes: `0` success, `2` input/argument failure, `3` numerical failure,
  `4` output/quota/I/O failure. Numerical errors identify the zero-based
  frequency index and Hz value. Unsupported physics is rejected, not approximated.

Rust adapters may call `bellhop_hdf5::kraken::run_legacy` with input/output
paths, engine, overwrite policy and byte quota. The numerical `kraken` crate
has no production HDF5 dependency.

## Resource limits and publication

Input limits are unchanged: 1 MiB per UTF-8 file, 1,000 frequencies, and
5,000,000 total copied case-vector entries. Inputs are read once into bounded
snapshots, parsed through `legacy::load_frequency_cases_with_boundary_tables`, and
hashed from those exact bytes. Metadata is never obtained by rereading an
input after solving.

The default cumulative output limit is **268,435,456 bytes (256 MiB)**;
`--max-output-bytes N` accepts a positive override. All FIELD float32 payloads
are admitted against it before numerical work. Numeric dataset payloads are
charged cumulatively before writing; actual HDF5 file size, including metadata,
is checked after the header, each frequency flush, final summary attributes,
and HDF5 close before publication.
A rejected run publishes nothing and attempts to remove its owned scratch file;
a cleanup I/O failure is reported as a warning.

Only one frequency's `SimulationResult` and one large component-write buffer
are retained at a time. Frequency cases themselves retain the loader's bounded
input copies. Existing mesh/mode/shape/FIELD work limits still apply **per
frequency**; this output quota is not a cumulative CPU budget, timeout or
cancellation mechanism. A scratch file can temporarily exceed the physical
quota during the header or one bounded frequency's writes, before its flush check; there is
no hard byte-limited HDF5 filesystem driver.

`<output>.tmp` is exclusively created by HDF5 and is never a truncated
pre-existing file. After all frequencies succeed, HDF5 flush/close and file
sync must succeed. Without overwrite, an atomic hard link refuses a destination
that appeared during execution; with overwrite, an atomic rename replaces the
old output. Scratch/output paths share a directory/filesystem. Filesystems
without the required hard-link/rename support return an explicit output error.
An already existing scratch file is left untouched. Atomic publication is not
an all-platform power-loss durability guarantee for directory metadata.

## Root attributes and axes

This is independent of [BELLHOP schema v3](output-format.md); the common
`bellhop-hdf5` crate reuses only I/O helpers, not solver result models.

| Attribute | Type / meaning |
|---|---|
| `schema_name` | UTF-8 `kraken` |
| `schema_version` | uint32 `1` |
| `implementation` | Rust implementation/package version |
| `compatibility_reference` | pinned Acoustics Toolbox v2023.5 commit |
| `title` | `.env` title |
| `solver` | UTF-8 `kraken` or `krakenc` |
| `coordinate_convention` | range origin at source; depth positive downward |
| `frequency_count` | uint64, number of input frequencies |
| `mode_count`, `pressure_count` | uint64 totals across all frequencies |
| `max_output_bytes` | uint64 configured quota |

`/frequency_hz` is float64 `[F]`, unit `Hz`, in input order. Descending and
repeated frequencies are retained. `/inputs/env` and `/inputs/flp` have UTF-8
`filename` (the supplied path), uint64 `size_bytes`, and UTF-8 `sha256`
attributes. Bottom F/P adds `/inputs/brc` or `/inputs/irc`; top F adds `/inputs/trc`
with the same attributes, hashed from the **actual parsed snapshot**. Unused
same-stem resources are neither read nor recorded. Input contents are not embedded.
These groups and the `surface_boundary`/`bottom_boundary` attributes below are additive schema-v1
metadata; existing dataset layouts, types, units and schema identity are unchanged.
Older v1 files may omit this additive metadata.

Results live under `/frequencies/0`, `/frequencies/1`, ... using **indices,
not Hz names**. Every frequency has its own mode count and arrays; no padded
common-mode matrix or cross-frequency mode correspondence is implied.

## Per-frequency metadata

Each `/frequencies/i` group has attributes:

- `frequency_hz`: float64 actual solve frequency;
- `mesh_reference_frequency_hz`: float64 nominal frequency; for a non-broadband
  input this equals the actual frequency;
- `requested_mesh_points`: uint64 nominal NG, including `0` for automatic;
- `max_range_m`: float64 extrapolation control;
- `field_mode_limit`: uint64 requested FIELD cap (not the total stored mode count);
- `source_geometry`: UTF-8 `line` or `point`;
- `surface_boundary`: UTF-8 `V`, `R`, `A` or `F` (additive v1 metadata);
- `bottom_boundary`: UTF-8 `V`, `A`, `R`, `F` or `P` (additive v1 metadata).

## Modes

Under `/frequencies/i/modes`, let `M` be that frequency's mode count and `D`
its modal sample-depth count. All datasets are float64 and have a `unit`
attribute:

| Dataset | Shape | Unit |
|---|---|---|
| `sample_depth_m` | `[D]` | `m` |
| `horizontal_wavenumber_real`, `horizontal_wavenumber_imaginary` | `[M]` | `rad/m` |
| `phase_speed_mps`, `group_speed_mps` | `[M]` | `m/s` |
| `attenuation_nepers_per_m` | `[M]` | `neper/m` |
| `eigenfunction_real`, `eigenfunction_imaginary` | `[M,D]` | `reference_normalized` |

`eigenfunction_axis_order` is `mode,sample_depth`. `normalization` identifies
the pinned AT density-weighted normalization, including a fluid half-space
where present; eigenvectors retain their arbitrary unit phase. The
`reference_normalized` unit is not an absolute pressure unit. Stored shape
values preserve the Rust result, including its reference complex32 sampling
rounding; float64 storage does not restore precision already rounded away.
Shapes and group speeds remain first-mesh values, even when wavenumbers are
Richardson-extrapolated. Attenuation is `-Im(k)`; neither imaginary component
nor phase convention is changed during serialization.

## FIELD

Under `/frequencies/i/field`, geometry is float64:

| Dataset | Shape | Unit |
|---|---|---|
| `source_depth_m` | `[S]` | `m` |
| `receiver_depth_m` | `[Z]` | `m` |
| `receiver_range_m` | `[R]` | `m` |
| `receiver_offset_m` | `[Z]` | `m` |
| `pressure_real`, `pressure_imaginary` | `[S,Z,R]`, float32 | `1` |

`pressure_axis_order` is `source_depth,receiver_depth,receiver_range`, with
range varying fastest in row-major storage. It is **not** the flattened,
range-major BELLHOP field layout. Offsets are separate per-depth values; the
solver's effective range is range plus offset. Pressure is relative complex
pressure, not calibrated pascals. FIELD arithmetic and storage remain single
precision; writing adds no new numerical rounding.

## Acceptance

Ordinary tests run the CLI, verify schema/types/units/metadata, compare all
serialized values exactly with the Rust result, retain duplicates and varying
mode counts, and exercise overwrite/input-alias/scratch protection, input bounds,
physical/payload quotas and failure after a successful first frequency.
BELLHOP v3 tests continue unchanged.

Fresh pinned CI additionally passes CLI-produced `.h5` files back through the
same strict `.mod/.prt/.shd` comparator (`KRAKEN_HDF5_RESULT`): both derived
Pekeris broadband pairs, original single-frequency MunkK, all five original
KRAKENC pairs, and original BroadBand/MunkK (both frequencies, 1,003,002
pressures), plus four small table derivatives and all three original TabRefCoef
geo/brc/irc workflows with BOUNCE-generated resources, four cubic/analytic
KRAKENC derivatives, broadband PCHIP Munk and original MunkS/MunkAnalytic
environments with **derived** FIELD geometry, 27 water-material workflows and
24 derived smooth-boundary/TRC workflows. Every mode, shape and
pressure is checked at the existing tolerances. The table derivatives have new
Fortran goldens; no numerical tolerance is changed. Rust HDF5 artifacts
are not committed as Fortran goldens.
