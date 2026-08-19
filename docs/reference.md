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
Cartesian and ray-centered geometric-hat beams, Cartesian geometric-Gaussian
beams, caustics, arrival combination, and multi-depth/multi-range receiver
grids. Pressure-field goldens cover coherent, semi-coherent, and incoherent
scaling, simple-Gaussian beams, and Cartesian/ray-centered Cerveny beams.
Dedicated reflection goldens cover acousto-elastic, grain-size, and `.irc`
impedance-table amplitude and phase.

## Pinned Linux x86-64 oracle

`tools/reference` provides the reproducible differential environment:

- Debian `bookworm-slim` image index digest:
  `sha256:abd67ffcfa541b485a3dff59865ab629aa048a6c613e639d36e7456b0b229241`
- GNU Fortran `12.2.0`
- flags: `-O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu`
- Acoustics Toolbox source archive SHA-256:
  `f8a7a2c1e80a73431cd230a10bef5fcfc996c88889a0e1540771c3922ee2a21f`
- Rust `1.88.0` bookworm image index digest:
  `sha256:af306cfa71d987911a781c37b59d7d67d934f49684058f96cf72079c3626bfe0`

The semantic trajectory comparator aligns isolated minimum-step vertices.
When a value lies within a few ulps of an SSP or boundary interface, one
compiler can take the reference's `1e-4 × base step` while another changes
segment or reflects immediately. Aligned coordinates retain the `1e-5 m`
tolerance; this discrete branch allowance is measured and reported separately.

Pinned Linux comparisons currently show:

- ParaBot: 201 rays, maximum coordinate error `1.4e-11 m`, no branch allowance
- Ellipse: 72 rays, maximum coordinate error `2.2e-11 m`, no branch allowance
- block: 50 rays, exact coordinate agreement
- DickinsBray: 501 rays, aligned error below `9e-10 m`; 16 minimum-step alignment
  operations with maximum branch displacement `0.12 m`

The committed `DickinsCritical` and `ParaBotCritical` cases preserve focused
coverage of the affected SSP-interface and curved-boundary decisions.

All seven official two-dimensional `A/a` environments parsed by the Rust
loader execute with the implemented `G/B` influence models. Three committed
arrival comparisons match receiver-by-receiver arrival counts exactly and
agree at the expected single-precision storage points. Larger exploratory
cases can differ by one or two arrivals out of approximately 35,000 to
502,000; these differences follow the same discrete near-boundary branches.

Small pressure grids match the reference single-precision values within
`5e-8` absolute pressure. Exploratory full-grid Munk coherent results and
Cerveny `F/M/W` width calculations also match at sampled single-precision
storage points.
