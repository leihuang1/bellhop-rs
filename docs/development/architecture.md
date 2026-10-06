# Pelagic architecture

Five internal, unpublished crates:

```text
cli (pelagic) ──> bellhop, kraken, output
server (pelagic-server) ──> bellhop, output
output ──> bellhop, kraken
```

`bellhop` and `kraken` are independent solver Modules. Each exposes validated
Cases and its own results; no common Case, result model, JSON schema, solver
trait or real/complex root engine is introduced.

```text
crates/{bellhop,kraken}/src/
  lib.rs          Interface exports
  model.rs        unvalidated input types
  case.rs         validation and immutable Case invariants
  diagnostic.rs   diagnostics
  result.rs       numerical products
  input/          legacy/JSON parsing and acquisition
  solver/         solver-local numerical Implementation
```

Root exports remain convenient Rust Interfaces; their implementations live in
the responsibility directories. Files follow actual concepts rather than line
counts. Shared layout does not imply symmetric arithmetic: KRAKEN/KRAKENC
root searches and FIELD precision stages retain their numerical Locality.
Legacy material conversion lives with KRAKEN input; its solver consumes
canonical solve-frequency material. See the [implementation map](kraken-internals.md).

`cli/src/{bellhop,kraken}.rs` retain solver-specific argument handling behind
one dispatch entry. `server` is an independent HTTP Adapter, still BELLHOP-only.
Neither introduces solver-independent transport models.

`output/src/{bellhop,kraken}.rs` own HDF5 v3/v1 separately. Only existing HDF5
metadata I/O and single-file publication policy are shared. Publication, input
protection, sequential failure handling, quotas and rollback behavior are
unchanged. The `implementation` metadata identifies Pelagic; dataset schemas,
precision and numerical contents are unchanged.

This first structural PR does not implement native output, `--format`, result
directories, artifact manifests or default replacement. Those are a separate
Adapter capability and acceptance boundary, not a side effect of moving files.
The original numerical/input target, references, tolerances, budgets and legacy
fixture bytes remain unchanged. Naming is not evidence of increased Depth;
the Leverage here is one CLI and predictable responsibility navigation, not a
new universal core.
