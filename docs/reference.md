# Numerical reference

The compatibility oracle is the official Acoustics Toolbox release `v2023.5`:

- Repository: <https://github.com/oalib-acoustics/Acoustics-Toolbox>
- Tag: `v2023.5`
- Commit: `475108519289c6fb488b58980c644ea14eccc604`
- License: GPL-3.0

Small numerical goldens are committed under
`crates/bellhop/tests/fixtures/golden`. Their per-file provenance, hashes,
compiler, target, and flags are recorded alongside them. Ray trajectories
cover all `N/C/P/S/Q/A` sound-speed models. Eigenray and arrival goldens cover
Cartesian and ray-centered geometric-hat beams, Cartesian geometric-Gaussian beams, caustics,
arrival combination, and multi-depth/multi-range receiver grids. Pressure-field
goldens cover coherent, semi-coherent, and incoherent scaling, simple-Gaussian
beams, and Cartesian/ray-centered Cerveny beams.

A local full-fan comparison has also exercised the seven two-dimensional ray
cases used by the official MATLAB test scripts. Range-independent PCHIP and
range-dependent quadrilateral cases agree point-for-point to approximately
`3.1e-10 m`; curved and piecewise boundary cases expose compiler-sensitive
near-zero boundary-crossing branches and remain under investigation.

All seven official two-dimensional `A/a` environments parsed by the Rust
loader execute with the implemented `G/B` influence models. Three committed
arrival comparisons match receiver-by-receiver arrival counts exactly and
agree at the expected single-precision storage points. Larger exploratory
cases differ by one or two arrivals out of approximately 35,000 to 502,000;
these differences follow the same near-boundary floating-point branches and
remain under investigation.

Small pressure grids match the reference single-precision values within
`5e-8` absolute pressure. Exploratory full-grid Munk coherent results and
Cerveny `F/M/W` width calculations also match at sampled single-precision
storage points on the local reference build.

The authoritative differential-test container will additionally pin a Linux
Fortran compiler, compiler flags, x86-64 target, and container digest. Those
container values have not yet been selected.
