# Intentional compatibility deviations

Only intentional differences that affect numerical or discrete behavior belong here. Ordinary porting defects are fixed directly and covered by tests.

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
