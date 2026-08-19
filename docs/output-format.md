# HDF5 output schema

Schema version 1 is used by `bellhop run`. Files are written only after a
complete simulation succeeds, using `<output>.tmp` followed by an atomic
rename.

## Root attributes

- `schema_version` (`u32`)
- `implementation`
- `compatibility_reference`
- `input_filename`
- `input_size_bytes` (`u64`)
- `input_sha256`
- `title`
- `frequency_hz`
- `legacy_run_options`
- `coordinate_convention`

Input contents are not embedded.

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
rounding points as Acoustics Toolbox `v2023.5`. Later field milestones may add
groups without changing the version-1 layouts.
