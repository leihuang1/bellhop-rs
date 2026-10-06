# Computing with Pelagic

Install from the repository with `cargo install --path crates/cli`, or replace
`pelagic` below with `cargo run --release -p cli --`. Use `pelagic --help` and
each subcommand's `--help` for the complete argument list.

## BELLHOP

```sh
pelagic bellhop validate path/to/case.env
pelagic bellhop export path/to/case.env > case.json
pelagic bellhop validate case.json
pelagic bellhop run case.json --output result.h5
```

`validate` checks inputs without solving. `export` resolves consumed legacy
resources and emits canonical self-contained JSON on stdout; warnings go to
stderr. JSON input never discovers neighboring auxiliary files.

See [input compatibility](../bellhop/compatibility.md),
[JSON schema](../bellhop/json.md) and [HDF5 v3](bellhop-hdf5.md).

## KRAKEN/KRAKENC

```sh
pelagic kraken run path/to/case.env --output result.h5
pelagic kraken run path/to/case.env --solver krakenc --flp geometry.flp --output complex.h5
pelagic kraken export path/to/case.env --solver krakenc > case.json
pelagic kraken run case.json --output modern.h5
```

Legacy input needs ENV plus FIELD geometry: `--flp` defaults to the same-stem
FLP. ENV owns boundary-table resources; FLP owns the source-pattern stem.
Legacy `--solver` defaults to `kraken`. JSON declares its engine; an explicit
`--solver` must match it, and `--flp` is invalid. Frequencies retain input order
and repetitions; solving/writing stops on the first failure.

See [fixed capabilities and exclusions](../kraken/compatibility.md),
[JSON schema](../kraken/json.md) and [HDF5 v1, quota and publication](kraken-hdf5.md).

## Output and failures

This release preserves the existing single-file HDF5 behavior. `--output`
names a file, defaulting to `<case-stem>.h5` in the current directory. Its
parent must exist. Existing files require explicit `--overwrite`; consumed
inputs and their symlink aliases are protected even with that flag. Owned
scratch is removed on failure; pre-existing scratch is never overwritten.
KRAKEN's default cumulative output quota remains 256 MiB, controlled by
`--max-output-bytes`; BELLHOP has no KRAKEN byte quota.

Exit codes are `0` success, `2` input/argument failure, `3` numerical failure,
and `4` output failure. Completion does not certify arbitrary-input oracle
parity. Native output formats, result directories and default replacement are
not implemented here.

The HTTP adapter is a separate `pelagic-server` executable; it supports only
BELLHOP and keeps its existing [HTTP protocol](http.md).
