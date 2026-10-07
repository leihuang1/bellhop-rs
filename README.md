# Pelagic

Two-dimensional underwater acoustics in Rust: independent BELLHOP ray/beam
and KRAKEN/KRAKENC normal-mode solvers, with Acoustics Toolbox `v2023.5` as the
fixed compatibility reference. Supported physics and limitations are documented
[for BELLHOP](docs/bellhop/compatibility.md) and
[for KRAKEN/KRAKENC](docs/kraken/compatibility.md); successful execution is not
an arbitrary-input numerical-parity certificate.

## Install

Windows x64 CLI builds are attached to new [GitHub Releases](https://github.com/leihuang1/pelagic/releases)
as `pelagic-vX.Y.Z-windows-x64.zip`, with a SHA-256 checksum file. Extract the ZIP
and run `pelagic.exe`; HDF5 and the MSVC runtime are statically linked. The ZIP
includes license notices, the exact source commit and dependency source links.
`pelagic-server` is not included.

To build from source, use Rust 1.88 or later and a C/C++ build toolchain with
CMake for static HDF5. The HTTP build also downloads Swagger UI assets.

```sh
cargo install --path crates/cli
cargo install --path crates/server
```

The executables are `pelagic` and `pelagic-server`. The five workspace crates
are internal (`publish = false`); there are no old CLI aliases or compatibility
packages.

## Compute

```sh
pelagic bellhop validate examples/field-g.json
pelagic bellhop run examples/field-g.json --output results/field
pelagic kraken run crates/kraken/tests/fixtures/Pekeris.env --output results/pekeris --format both
pelagic kraken run crates/kraken/tests/fixtures/PekerisComplexBlank.env --solver krakenc --output results/complex
pelagic kraken export crates/kraken/tests/fixtures/Pekeris.env > case.json
```

Both solvers accept legacy inputs and their own strict, self-contained JSON
schema. Computation defaults to **native output**; `--format hdf5|both` selects
HDF5 or both from one solve. All formats write result directories and safely
update only verified owned artifacts, preserving unrelated files and consumed
inputs. See the [CLI guide](docs/guide/cli.md) and
[native layouts/publication contract](docs/guide/native-output.md).

## HTTP

```sh
RUST_LOG=info pelagic-server
curl --fail-with-body -H 'Content-Type: application/json' \
  --data-binary @examples/field-g.json http://localhost:8080/v1/run --output result.h5
```

The independent HTTP adapter supports BELLHOP JSON only. Its routes and
HDF5/JSON responses are unchanged. See [HTTP configuration and limits](docs/guide/http.md).

## Documentation

[Documentation index](docs/README.md): input contracts, separate HDF5 schemas,
capabilities, architecture, pinned-reference reproduction and historical evidence.
[CONTEXT.md](CONTEXT.md) defines the domain language.

## License

GPL-3.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
