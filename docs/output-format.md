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

Every dataset has a `unit` attribute. Later field and arrival milestones will
add groups without changing the version-1 ray layout.
