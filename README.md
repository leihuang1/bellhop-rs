# bellhop-rs

A modern Rust implementation of the two-dimensional BELLHOP underwater-acoustics model.

The compatibility baseline is Acoustics Toolbox `v2023.5`. The project loads complete legacy 2D cases, including their `.ssp`, `.ati`, `.bty`, `.brc`, `.trc`, and `.sbp` inputs. The deterministic compatibility solver implements `R` ray traces; `E` eigenrays and `A/a` arrivals with geometric-hat and Cartesian geometric-Gaussian beams; and `C/S/I` pressure fields with geometric-hat, Cartesian geometric-Gaussian, simple-Gaussian, and Cartesian/ray-centered Cerveny beams. All `N/C/P/S/Q/A` sound-speed models and supported 2D boundaries are available.

## Commands

```console
cargo run -p bellhop-cli -- validate path/to/case.env
cargo run -p bellhop-cli -- run path/to/ray-case.env --output result.h5
```

`run` writes a [versioned HDF5 result](docs/output-format.md) using a temporary file and atomic rename. Existing outputs require `--overwrite`. Precalculated `.irc` tables, `W` table generation, and reference-unsupported ray-centered geometric-Gaussian beams remain unavailable.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
