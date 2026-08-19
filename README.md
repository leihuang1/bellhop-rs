# bellhop-rs

A modern Rust implementation of the two-dimensional BELLHOP underwater-acoustics model.

The compatibility baseline is Acoustics Toolbox `v2023.5`. The project currently loads and validates complete legacy 2D cases, including their `.ssp`, `.ati`, `.bty`, `.brc`, `.trc`, and `.sbp` inputs. It does not yet contain a numerical solver.

## Current command

```console
cargo run -p bellhop-cli -- validate path/to/case.env
```

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
