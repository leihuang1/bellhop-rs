# Intentional compatibility deviations

Only intentional differences that affect numerical or discrete behavior belong here. Ordinary porting defects are fixed directly and covered by tests.

## Initial KRAKEN Pekeris slice

`crates/kraken` uses the analytical Pekeris dispersion relation and mode
normalization for homogeneous-fluid cases. Pinned KRAKEN uses a finite-difference
mesh; its `.mod` eigenvectors come from the first mesh, while the reported
eigenvalues use mesh extrapolation. The following are measured differences
against the pinned reference, not a bound for arbitrary environments or meshes:

| Fixture | Mesh | Modes / pressures | Max wavenumber error (m^-1) | Max shape error | Max pressure error |
|---|---:|---:|---:|---:|---:|
| Pekeris | 1000 | 3 / 9 | 2.3e-11 | 5.7e-7 | 1.3e-6 |
| PekerisFiltered | 1000 | 2 / 18 | 2.0e-11 | 5.7e-7 | 1.3e-7 |
| PekerisDense | 4000 | 7 / 30 | 4.3e-11 | 2.2e-7 | 3.9e-7 |

Group speed comparison is limited by `.prt` print precision of `0.01 m/s`.
Mode-shape comparisons align arbitrary unit phase but do not rescale amplitudes.
The tolerances remain `5e-10 m^-1` for printed wavenumbers, `1e-6` for shapes,
and `2e-6` for complex pressure.

The initial all-double FIELD implementation exceeded the existing pressure
limit on PekerisDense (`3.7e-6` at the first failing sample). Restoring the
reference's single-precision wavenumber/products/accumulation and separate
range/offset exponentials reduced the maximum to `3.9e-7`; the tolerance was
not widened. Geometry, interpolation, and analytical mode results remain double
precision. FIELD pressures are promoted from single precision in the public
`Complex64` result.

Committed raw Fortran goldens and fresh CI reference runs use the same
comparator. General profiles still require a different solver path and their
own differential results; the analytical solver is also an independent
baseline for that future work.

## Double-precision source and receiver geometry

bellhop-rs stores and interpolates source and receiver depths in `f64` instead
of reproducing Acoustics Toolbox's single-precision position grids. This keeps
all geometric coordinates at the solver's native precision, preserves modern
JSON inputs, and avoids rejecting decimal coordinates that lie exactly on an
`f64` water-column boundary.

The pinned `sbcx_Arr_asc` comparison changes four of 50,000 receiver arrival
counts by one at influence-window edges. Receivers with matching counts remain
within `3.5e-7` amplitude, `2.0e-5` degrees, and `2.0e-6` seconds of the
single-precision reference. Critical rays and the 501,501-sample Munk coherent
field retain their previous tolerances. Arrival values and pressure remain
single precision at their documented output rounding points.

## Precalculated `.irc` bottom reflection

Acoustics Toolbox BELLHOP `v2023.5` reads a bottom `P` option and loads the
same-stem `.irc` file through `misc/RefCoef.f90`, but
`bellhop.f90::Reflect2D` has no `P` branch and terminates with an unknown
boundary condition when the first such reflection occurs.

bellhop-rs implements the evidently intended behavior: it uses
`RefCoef.f90::InterpolateIRC`'s power-scaled quadratic interpolation and the
complex reflection-coefficient formula used by `Kraken/bounce.f90`. A reduced
BOUNCE-generated `.irc`/`.brc` pair verifies amplitude and phase. This makes a
traditional input useful where the pinned BELLHOP oracle is internally
incomplete.

The top-boundary `P` option remains rejected because the `.irc` generation and
loading path defines a bottom impedance only.
