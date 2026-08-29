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

The multi-client service accepts modern JSON only and returns HDF5 directly.
It includes bounded blocking execution, timeouts, server-controlled simulation
limits, optional Bearer authentication, tracing, health checks, generated
OpenAPI at `/openapi.json`, and Swagger UI at `/docs`. See the
[HTTP service guide](docs/api.md).

The reference-unsupported ray-centered geometric-Gaussian path remains
unavailable. BELLHOP v2023.5's `W` boundary option is rejected because the
reference advertises table generation but contains neither a generator nor the
layer inputs needed to produce a table.

## Workspace

- `bellhop`: legacy/JSON input models and deterministic solver
- `bellhop-hdf5`: shared HDF5 schema writer
- `bellhop-cli`: local validation, conversion, and simulation
- `bellhop-server`: synchronous JSON/HDF5 HTTP service

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
