# bellhop-rs

A modern Rust implementation of the two-dimensional BELLHOP underwater-acoustics model.

The compatibility baseline is Acoustics Toolbox `v2023.5`. The project loads complete legacy 2D cases, including their `.ssp`, `.ati`, `.bty`, `.brc`, `.trc`, and `.sbp` inputs. The deterministic compatibility solver currently implements `R` ray-trace runs with all `N/C/P/S/Q/A` sound-speed models and supported 2D boundaries.

## Commands

```console
cargo run -p bellhop-cli -- validate path/to/case.env
cargo run -p bellhop-cli -- run path/to/ray-case.env --output result.h5
```

`run` writes a [versioned HDF5 result](docs/output-format.md) using a temporary file and atomic rename. Existing outputs require `--overwrite`. Field, eigenray, and arrival calculations are not implemented yet.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
