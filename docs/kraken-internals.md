# KRAKEN implementation map

The public test surface remains `Case::from_definition`, the legacy/JSON
loaders, `solve`, `solve_complex_modes`, `solve_field` and `solve_frequencies`. Numerical acceptance
is still the [fixed compatibility target](kraken-compatibility.md); this map
adds no physics or new acceptance requirements.

## Input to validated Case

- `crates/kraken/src/input.rs` owns legacy acquisition and exact consumed-input
  snapshots; ENV and FLP can own different resource stems.
- `legacy.rs` owns records, inheritance, profile/frequency ordering, resource
  selection and source-location diagnostics.
- `legacy/material.rs` owns **Legacy materials**: raw absorption and power laws
  travel with each fluid, solid or half-space, never in canonical loss fields.
  Its `case_definition` selects one frequency and converts the whole material
  stack using `attenuation::db_per_wavelength`. Storage accounting happens
  before frequency copies. Biological loss uses node depths, but half-space
  conversion retains the reference's HUGE-depth exclusion.
- `json.rs` imports canonical definitions directly, without legacy conversion.
- `layers::validate` owns fluid SSP/loss rules, interpolation validation and
  minimum-speed calculation through the existing layer iterator. It retains
  per-layer diagnostic names and aggregate ordering, including existing repeated
  interpolation diagnostics. Valid additional profiles are constructed once
  during Case validation, not once for checking and again for their minima.
- `Case::from_definition` retains topology, spectral, geometry and remaining
  boundary/material checks. Validated Cases stay immutable.

## Validated Case to modes and FIELD

`solver::solve_modes` chooses `modes.rs` (KRAKEN) or `complex_modes.rs`
(KRAKENC). Both use `Profile`, finite-layer meshes and elastic impedances.
Their root searches, precision, work limits and deflation arithmetic remain
separate Implementations.

`refinement.rs` owns **Modal refinement** behind `seed` and `accept`:

- mesh multipliers remain 1, 2, 4, 8, 16;
- raw root history is separate from the Richardson table, and includes
  cLow-excluded roots used for KRAKENC deflation;
- mesh two still uses the original scan; Neville seeds start on mesh three;
- first-mesh shapes/group speeds are retained; ordinary fluid and KRAKENC
  mode-count changes are rejected, while real finite solids retain surviving
  first-mesh data when Solve2 reduces its search bound;
- KRAKEN extrapolates real k² and retains first-mesh loss; KRAKENC extrapolates
  complex k². Standard arithmetic bounds share bookkeeping, not a new numerical
  trait or plugin Interface.

`solve_frequencies` keeps real Solve2's run-local bound across ordered blocks,
without retaining all results; the first error ends the iterator. Independent
`solve`/`solve_field` calls reset the bound. HDF5 legacy/JSON execution and API
differential checks share this Interface; no global solver state is used.

FIELD synthesis and multi-profile propagation still own their intentionally
different precision/operation grouping. The
[known three-layer refinement gap](kraken-layered-refinement-gap.md) is unchanged.

`tests/refactoring.rs` characterizes the Case/load/solve Interfaces; existing
behavior, JSON round-trip, material, FIELD, CLI/HDF5 and fresh-reference tests
remain the acceptance surface. No helper-only test replaces those checks.
