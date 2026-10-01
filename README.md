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
`bellhop-hdf5` crate, using a temporary file and atomic rename. Existing outputs
require `--overwrite`.

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
and coherent line- or point-source FIELD. A separate, narrow
KRAKENC path supports complex modes and coherent line-/point-source FIELD for
complex `N/C/P/S` or lossless fixed-Munk-A water over a fluid bottom with Richardson
mesh extrapolation, covering trapped and leaky spectral intervals. Legacy
N/M/m/F/W/Q/L loss units and T/F/B volume loss are supported; density gradients within a layer and elastic media remain excluded. Fluid
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
4,902 pressures). The original TLslices `double` is still rejected for changing
mesh mode counts; its denser-mesh derivative is accepted, not relabelled original.
Analytic Munk and F/P/TRC boundaries remain single-layer.
`legacy::load_frequency_cases` supports discrete multi-frequency KRAKEN/KRAKENC
runs; original `BroadBand/MunkK` passes through KRAKEN at 50 and 500 Hz (102/1,023 modes,
1,003,002 total complex pressures). KRAKENC 500 Hz still hits its unchanged
300M root-work ceiling; no partial result is published. This is not a time-domain response. Its
full 2D acceptance target and explicit exclusions are in
[the KRAKEN compatibility matrix](docs/kraken-compatibility.md).

```console
cargo run --release -p kraken-cli -- run crates/kraken/tests/fixtures/PekerisBroadband.env --output broadband.h5
cargo run --release -p kraken-cli -- run crates/kraken/tests/fixtures/PekerisComplexBroadband.env --solver krakenc --output complex.h5
```

The `kraken` binary reads `.env` plus the same-stem `.flp` (or `--flp PATH`),
and the selected `.brc/.irc` bottom or `.trc` top resource for KRAKENC, then writes [KRAKEN HDF5 schema v1](docs/kraken-output-format.md). Frequencies
are solved and written sequentially in input order, including duplicates;
existing outputs require `--overwrite`. The default cumulative output quota is
256 MiB (`--max-output-bytes`). This adapter does **not** expand the numerical
matrix or provide JSON/HTTP or time-domain products. BELLHOP schema v3 is unchanged.

## Workspace

- `bellhop`: legacy/JSON input models and deterministic solver
- `bellhop-hdf5`: separate BELLHOP v3 and KRAKEN v1 HDF5 adapters
- `bellhop-cli`: local validation, conversion, and simulation
- `bellhop-server`: synchronous JSON/HDF5 HTTP service
- `kraken`: layered-fluid 2D KRAKEN/KRAKENC modes and coherent FIELD
- `kraken-cli`: supported legacy KRAKEN/KRAKENC pairs to atomic, bounded HDF5 output

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
