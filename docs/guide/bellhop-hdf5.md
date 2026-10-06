# BELLHOP HDF5 output schema

KRAKEN/KRAKENC use a [separate schema v1](kraken-hdf5.md); this
BELLHOP v3 contract is unchanged.

Schema version 3 is written by both `pelagic bellhop run` and `POST /v1/run` through the
shared `output` crate. Local BELLHOP and KRAKEN runs share one publication
implementation, not a numerical model or schema. BELLHOP exclusively creates
`<output>.tmp`, writes schema v3, closes and syncs it, then installs the complete
file: an atomic hard link without `--overwrite`, an atomic rename with it.
Destinations appearing during the run are not clobbered without that flag.
Primary inputs, consumed auxiliary inputs and their symlink aliases cannot be
destinations. Pre-existing scratch is never truncated or removed; failure removes
only owned scratch and preserves the previous result. The output parent must
exist, with same-filesystem hard-link/rename support; this is not a directory
power-loss durability guarantee. No KRAKEN byte quota is imposed on BELLHOP.

Rust callers can use `write_hdf5_atomic` with all consumed input paths for this
policy. The existing `write_hdf5` remains a low-level path writer for trusted
private files; the HTTP implementation still uses it in a private temporary
directory and returns its bytes directly as `application/x-hdf5`.

## Root attributes

- `schema_version` (`u32`, currently `3`)
- `implementation`
- `compatibility_reference`
- `input_filename`
- `input_size_bytes` (`u64`)
- `input_sha256`
- `title`
- `frequency_hz`
- `legacy_run_options`
- `coordinate_convention`
- `warnings` (one-dimensional variable-length UTF-8 string array)

For CLI runs, input metadata describes the exact parsed primary `.env` or `.json`
snapshot, retained during loading rather than reread after validation. For HTTP runs it describes the exact request-body bytes and uses
`request.json` as the filename. Input contents and legacy auxiliary files are
not embedded. `legacy_run_options` is empty for modern JSON cases.

Each `warnings` entry is the rendered structured input diagnostic, including
source, severity, diagnostic code, message, and field. The attribute is present
as an empty array when no warnings were emitted. All source and receiver
coordinate datasets use double precision.

## Ray data

Ray trajectories use a flattened ragged-array representation under `/rays`:

- `source_depth_m`: source depths
- `ray_source_index`: source index for each ray
- `launch_angle_degrees`: launch angle for each ray, positive downward
- `point_offset`: `ray_count + 1` offsets into all point datasets
- `top_bounces`, `bottom_bounces`: final boundary-reflection counts
- `termination`: numeric termination code described by the group attribute
  `termination_codes`
- `range_m`, `depth_m`: point coordinates; depth is positive downward
- `travel_time_s`, `attenuation_time_s`: real and imaginary components of the
  complex accumulated travel time
- `amplitude`, `phase_radians`: cumulative reflection amplitude and phase

Every dataset has a `unit` attribute.

## Eigenray data

`/eigenrays` records receiver coordinates and flattened eigenray trajectories.
`receiver_offset` groups receivers by source, `eigenray_offset` groups
eigenrays by receiver, and `point_offset` groups points by eigenray. The
trajectory quantities use the same names and units as `/rays`.

## Arrival data

`/arrivals` uses `receiver_offset` to group receiver coordinates by source and
`arrival_offset` to group arrivals by receiver. Per-arrival datasets are:

- `amplitude`, `phase_radians`
- `travel_time_s`, `attenuation_time_s`
- `source_angle_degrees`, `receiver_angle_degrees`
- `top_bounces`, `bottom_bounces`

Arrival quantities intentionally use single-precision storage at the same
rounding points as Acoustics Toolbox `v2023.5`.

## Pressure-field data

`/field` stores flattened receiver samples for coherent, semi-coherent, and
incoherent calculations:

- `source_depth_m`: source depths
- `receiver_offset`: `source_count + 1` offsets grouping samples by source
- `receiver_range_m`, `receiver_depth_m`: sample coordinates
- `pressure_real`, `pressure_imaginary`: single-precision relative complex
  pressure

Samples are range-major, then depth-major for rectilinear grids. Irregular
grids contain one depth per range. Semi-coherent and incoherent results follow
the reference pressure-scaling convention and therefore have zero imaginary
components.

## Version history

- **v3:** stores source and receiver coordinates in double precision.
- **v2:** adds the root `warnings` string-array attribute and supports exact JSON
  request metadata from the HTTP service.
- **v1:** initial flattened rays, eigenrays, arrivals, and pressure fields.
