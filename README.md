# bellhop-rs

A modern Rust implementation of the two-dimensional BELLHOP underwater-acoustics model.

The compatibility baseline is Acoustics Toolbox `v2023.5`. The project loads complete legacy 2D cases, including their `.ssp`, `.ati`, `.bty`, `.brc`, `.trc`, `.irc`, and `.sbp` inputs. It also provides a strict, [self-contained JSON format](docs/json-input.md) for the CLI and HTTP service. The deterministic compatibility solver implements `R` ray traces; `E` eigenrays and `A/a` arrivals with geometric-hat and Cartesian geometric-Gaussian beams; and `C/S/I` pressure fields with geometric-hat, Cartesian geometric-Gaussian, simple-Gaussian, and Cartesian/ray-centered Cerveny beams. All `N/C/P/S/Q/A` sound-speed models and supported 2D boundaries are available.

## CLI

```console
cargo run -p bellhop-cli -- validate path/to/case.env
cargo run -p bellhop-cli -- export path/to/case.env > case.json
cargo run -p bellhop-cli -- run case.json --output result.h5
```

`validate` and `run` accept legacy `.env` or modern `.json` cases. `export`
resolves legacy auxiliary files and emits one canonical JSON document. `run`
writes a [versioned HDF5 result](docs/output-format.md) through the shared
`bellhop-hdf5` crate, using the same exclusive-scratch, atomic publication policy
as KRAKEN. Existing outputs require `--overwrite`; consumed inputs and their
symlink aliases cannot be destinations.

## HTTP service

```console
RUST_LOG=info cargo run -p bellhop-server --release
curl -H 'Content-Type: application/json' \
  --data-binary @examples/field-g.json \
  http://localhost:8080/v1/run \
  --output result.h5
```

The multi-client service accepts modern JSON only. `/v1/run` returns the full
HDF5 result, while `/v1/arrivals` and `/v1/field` return focused versioned JSON
arrival and pressure-field results for their corresponding run kinds. It
includes bounded blocking execution, timeouts, server-controlled simulation
limits, optional Bearer
authentication, tracing, health checks, generated OpenAPI at `/openapi.json`,
and Swagger UI at `/docs`. See the [HTTP service guide](docs/api.md).

The reference-unsupported ray-centered geometric-Gaussian path remains
unavailable. BELLHOP v2023.5's `W` boundary option is rejected because the
reference advertises table generation but contains neither a generator nor the
layer inputs needed to produce a table.

## KRAKEN/KRAKENC

The workspace also contains a separate KRAKEN rewrite: `kraken` currently covers
validated range-independent `N/C/P/S` fluid stacks and single-layer fixed analytic
Munk `A` profiles, smooth V/R/A outer boundaries, material/volume attenuation,
and FIELD with line, point or scaled-cylindrical sources,
omnidirectional/tabulated patterns, and coherent/incoherent mode addition.
[Multi-profile adiabatic and smooth-fluid coupled FIELD](docs/kraken-multi-profile-field.md)
now cover both original Gulf propagation paths. A
separate, narrow KRAKENC path supports complex modes and the same FIELD options for
complex `N/C/P/S` or lossless fixed-Munk-A fluid stacks with Richardson
mesh extrapolation, covering trapped and leaky spectral intervals. Legacy
N/M/m/F/W/Q/L loss units and T/F/B volume loss are supported; density gradients
within a fluid layer remain excluded. Fluid
interfaces may have density, sound-speed and loss jumps; each layer has its own
mesh, with all budgets shared across the stack. Both backends accept smooth V/R/A top and
bottom combinations; KRAKENC also accepts a propagating, single-frequency N/C
top F/TRC subset (cLow >= last-node speed), with no refinement or simultaneous bottom table. Unmodified
MunkKleaky, MunkKwb, MunkKbb, sductK and refined calibK `.env/.flp` pairs are
accepted end to end (1,077 modes and 201,201 complex pressures for sductK).
Cubic/analytic KRAKENC derivatives and unmodified MunkS/MunkAnalytic environments
with separately derived FIELD geometry also pass modes, FIELD and CLI/HDF5.
Single-frequency, RMax=0 KRAKENC F/BRC and P/IRC bottoms also accept original
TabRefCoef geo/brc/irc inputs with generated reference tables, including all
151,803 pressures; auxiliary-file snapshots are recorded in CLI/HDF5 provenance.
Twenty-one derived layered-fluid pairs pass 41 API/CLI-HDF5 workflows (530 modes,
4,902 pressures). Byte-original TLslices `double` and official `fieldbat.flp`
add two accepted workflows (84 modes, 1,002 pressures). Both engines retain the
surviving first-mesh data when refinement reduces the spectrum from 43 to 42;
the denser-mesh derivative remains separately labelled.
A [known wide three-layer KRAKENC refinement gap](docs/kraken-layered-refinement-gap.md)
remains outside acceptance as a non-blocking release exception: Rust can return
five modes where the pinned reference returns four. CLI success is not a guarantee
of arbitrary-input reference parity; the extra mode may affect coherent FIELD.
[Elastic half-spaces](docs/kraken-elastic-halfspaces.md) cover top/bottom A boundaries
through both engines: 37 full API/CLI-HDF5 workflows, including original TLslices
`scholte/normal/flused` (473 modes, 5,715 pressures). Real tops retain pinned
shared isolation/Brent/Solve2 and normalization: top compressional loss acts,
while elastic shear/bottom absorption is omitted. Use KRAKENC for full elastic loss.
[Finite elastic caps](docs/kraken-finite-elastic-layers.md) now pass
50 full API/CLI-HDF5 workflows, including original `elsed/ice` (933 modes,
10,536 pressures). Both engines support contiguous fluid stacks; KRAKEN retains
its reference real-stiffness/loss limitation. `solve_frequencies` preserves
ordered finite-elastic search state across a run.
Depth-varying finite cp/cs/density/P/S loss also passes 23 complete workflows,
158 modes and 1,827 pressures, including ordered frequency conversion,
self-contained JSON and additive HDF5 material-profile metadata.
Analytic Munk and F/P/TRC boundaries remain single-layer.
`legacy::load_frequency_cases` supports discrete multi-frequency KRAKEN/KRAKENC
runs; original `BroadBand/MunkK` passes through both engines at 50 and 500 Hz
(102/1,023 modes and 1,003,002 complex pressures per engine), including
legacy/JSON CLI-HDF5. KRAKENC uses about 153M root work at 500 Hz, below its
unchanged 300M ceiling; late failures still publish no partial result.
This is not a time-domain response. Its
full 2D acceptance target and explicit exclusions are in
[the KRAKEN compatibility matrix](docs/kraken-compatibility.md).

```console
cargo run --release -p kraken-cli -- run crates/kraken/tests/fixtures/PekerisBroadband.env --output broadband.h5
cargo run --release -p kraken-cli -- run crates/kraken/tests/fixtures/PekerisComplexBroadband.env --solver krakenc --output complex.h5
cargo run --release -p kraken-cli -- export crates/kraken/tests/fixtures/PekerisBroadband.env > case.json
cargo run --release -p kraken-cli -- run case.json --output modern.h5
```

The `kraken` binary reads `.env` plus the same-stem `.flp` (or `--flp PATH`),
the selected `.brc/.irc` bottom or `.trc` top resource for KRAKENC, and a
same-stem `.sbp` when the FIELD input requests a source pattern, then writes
[KRAKEN HDF5 schema v1](docs/kraken-output-format.md). Multi-profile ENV sequences
retain every profile's modes and synthesize one FIELD grid per frequency. Frequencies
are solved and written sequentially in input order, including duplicates;
existing outputs require `--overwrite`. The default cumulative output quota is
256 MiB (`--max-output-bytes`). This adapter does **not** expand the numerical
matrix. It also accepts [strict self-contained JSON v1](docs/kraken-json-input.md),
with all resources inline and exact parsed-byte provenance. HTTP and time-domain
products are not provided; BELLHOP schema v3 is unchanged.

## Workspace

- `bellhop`: legacy/JSON input models and deterministic solver
- `bellhop-hdf5`: separate BELLHOP v3 and KRAKEN v1 HDF5 adapters
- `bellhop-cli`: local validation, conversion, and simulation
- `bellhop-server`: synchronous JSON/HDF5 HTTP service
- `kraken`: layered-fluid 2D KRAKEN/KRAKENC modes and FIELD
- `kraken-cli`: legacy/JSON KRAKEN/KRAKENC inputs to atomic, bounded HDF5 output

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
