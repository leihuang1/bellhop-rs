# Numerical reference

The compatibility oracle is the official Acoustics Toolbox release `v2023.5`:

- Repository: <https://github.com/oalib-acoustics/Acoustics-Toolbox>
- Tag: `v2023.5`
- Commit: `475108519289c6fb488b58980c644ea14eccc604`
- License: GPL-3.0

Small `R`-mode golden trajectories are committed under
`crates/bellhop/tests/fixtures/golden`. Their per-file provenance, hashes,
compiler, target, and flags are recorded alongside them. They cover all
`N/C/P/S/Q/A` sound-speed models.

A local full-fan comparison has also exercised the seven two-dimensional ray
cases used by the official MATLAB test scripts. Range-independent PCHIP and
range-dependent quadrilateral cases agree point-for-point to approximately
`3.1e-10 m`; curved and piecewise boundary cases expose compiler-sensitive
near-zero boundary-crossing branches and remain under investigation.

The authoritative differential-test container will additionally pin a Linux
Fortran compiler, compiler flags, x86-64 target, and container digest. Those
container values have not yet been selected.
