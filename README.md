# Pelagic

Two-dimensional underwater acoustics in Rust: independent BELLHOP ray/beam
and KRAKEN/KRAKENC normal-mode solvers, with Acoustics Toolbox `v2023.5` as the
fixed compatibility reference. Supported physics and limitations are documented
[for BELLHOP](docs/bellhop/compatibility.md) and
[for KRAKEN/KRAKENC](docs/kraken/compatibility.md); successful execution is not
an arbitrary-input numerical-parity certificate.

## Install

Requires Rust 1.88 or later and a C/C++ build toolchain with CMake for static
HDF5. The HTTP build also downloads Swagger UI assets.

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
pelagic bellhop run examples/field-g.json --output field.h5
pelagic kraken run crates/kraken/tests/fixtures/Pekeris.env --output modes.h5
pelagic kraken run crates/kraken/tests/fixtures/PekerisComplexBlank.env --solver krakenc --output complex.h5
pelagic kraken export crates/kraken/tests/fixtures/Pekeris.env > case.json
```

Both solvers accept legacy inputs and their own strict, self-contained JSON
schema. This release still writes **HDF5 only**, to a single file; existing
outputs require `--overwrite`. Consumed inputs and their symlink aliases cannot
be destinations. See the [CLI guide](docs/guide/cli.md).

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

GPL-3.0-or-later. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
