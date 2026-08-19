# Legacy input compatibility

The input reader targets the two-dimensional environment syntax accepted by Acoustics Toolbox `v2023.5`.

Implemented:

- Fortran list-directed records used by legacy inputs
- quoted and unquoted fields, commas, `!` comments, `/` termination, repetition, null values, and `D` exponents
- legacy vector subtabulation and sorting
- `N/C/P/S/Q/A` SSP selection
- supported boundary, attenuation, run, receiver-grid, source, and beam options
- development-mode single-beam selection
- range-dependent `.ssp` matrices
- short and long, piecewise-linear and curvilinear `.ati`/`.bty` boundary data
- `.brc`/`.trc` magnitude and phase tables
- `.sbp` source beam patterns
- same-directory, same-stem auxiliary-file resolution through `load_case`
- semantic validation without invoking the numerical solver

Not yet implemented:

- precalculated internal `.irc` reflection tables
- numerical execution

Three-dimensional options are rejected explicitly.

The complete case loader has also been exercised against 69 two-dimensional environments invoked by the official `v2023.5` MATLAB test scripts, including all auxiliary files they request. `tests/MunkRot/MunkRot.env` is intentionally rejected because its run type explicitly selects three-dimensional behavior.
