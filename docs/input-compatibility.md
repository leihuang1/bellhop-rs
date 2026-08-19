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

Numerically implemented:

- deterministic single-threaded `R` ray traces
- `N/C/P/S/Q/A` sound-speed evaluation and derivatives
- reference modified polygon/midpoint integration with interface and boundary step reduction
- flat, piecewise-linear, curvilinear, and range-dependent boundaries
- vacuum, rigid, acousto-elastic, grain-size, and tabulated reflection conditions
- source beam-pattern weighting and reflection amplitude/phase accumulation
- `E` eigenray and `A/a` arrival calculations for Cartesian/ray-centered geometric-hat (`G/g`) and Cartesian geometric-Gaussian (`B`) beams
- simple-Gaussian (`S`) eigenray detection
- point and line source scaling, rectilinear and irregular receiver grids, caustic phases, arrival combination, and strongest-arrival retention

Not yet implemented:

- precalculated internal `.irc` reflection tables and `W` table generation
- Cerveny (`C/R`) influence models and simple-Gaussian (`S`) arrival influence
- `C`, `S`, and `I` field calculations

Three-dimensional options are rejected explicitly.

The complete case loader has also been exercised against 69 two-dimensional environments invoked by the official `v2023.5` MATLAB test scripts, including all auxiliary files they request. `tests/MunkRot/MunkRot.env` is intentionally rejected because its run type explicitly selects three-dimensional behavior.
