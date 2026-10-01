# Wide three-layer KRAKENC refinement: open parity gap

**Not fixed; not an accepted workflow.** PR #25 remains draft. No solver,
reference, tolerance or accepted golden was changed by this diagnosis.

## Reproducer and expected contract

Use `LayeredFluidThreeWide.env/.flp`, changing only RMax from 0 to 1000 km.
The three finite layers have NG=101/83/127, with cLow/cHigh=1400/1800 m/s.
Pinned Acoustics Toolbox v2023.5 (`475108519289c6fb488b58980c644ea14eccc604`),
GNU Fortran 12.2.0 with the documented fast-math flags, finds five modes on
mesh 1 but only four on mesh 2. Rust finds five on both meshes. This is **not**
a fifth root introduced by Richardson extrapolation or HDF5 serialization.

If Rust reproduced that search, its existing mode-count-change guard would
reject the run; silently truncating five modes to four is not the contract.
The isolated, intentionally ignored reproducer currently **fails**:

```sh
cargo test --release -p kraken --test layered \
  wide_three_layer_refinement_rejects_the_pinned_count_change \
  -- --ignored --exact --nocapture
```

The ignored test records unfinished work; it does not count as passing validation.
The existing full differential comparator independently fails `reference=4,
Rust=5` through both API and actual CLI-produced HDF5.

## Localization

Three hypotheses were checked: mesh/search-state differences, floating-point
coefficient/interface arithmetic, and deflation scaling/termination differences.

- Mesh-2 initial spectral seed is identical:
  `x = 0.047895828732633719 + 0i` (x is k²).
- All 625 real mesh-2 coefficients are bit-identical. Imaginary coefficients
  differ by at most `2.97e-21`; this does not establish bitwise equivalence of
  the complete operators, but rules out a gross mesh/material discrepancy here.
- First four roots agree to small floating-point errors. The fourth raw root's
  real part differs by about `1.67e-15`. `Solve2` multiplies that root by the
  f32 literal `1.00001`, then `RootFinderSecant` evaluates two nearby points.
- At the fifth search, Fortran's first real secant step goes from about
  `0.035248` to `0.029282` and eventually converges outside the interval:
  `x ≈ 0.0273422640556015 - 0.000981650507728i`.
  Rust instead first steps to about `0.034250` and converges inside it:
  `x ≈ 0.032081981205861 - 0.0000628866847424i`.
- Holding the first four deflated roots and the two evaluation points identical
  makes the real interface transfers bit-identical at those probe points and
  the complex deflated dispersion evaluations agree to relative < `1e-12`.
  The natural trajectories do not hold these points identical. Their initial
  secant denominator subtracts values of order `8.8e10` whose real parts differ
  by only hundreds/thousands; small evaluation errors change the next step.
- Scaling powers remain zero throughout this fifth-root probe; a scaling-loop
  difference is not the cause of this particular divergence.

This localizes the count divergence to a roundoff-sensitive, deflated secant
trajectory. It does **not** prove there are no other porting defects or that
either solver finds every mathematical root.

## Controlled reference-only experiment

`tools/reference/probe-layered-refinement.py` copies the source from the pinned
image and rebuilds it in a disposable container. It changes only the real part
of the mesh-2, mode-5 initial guess, **after** the tolerance has been computed.
No material, mesh, deflated root, stopping tolerance or other initial guess is
changed. The installed oracle/image and all accepted goldens remain untouched.

| Reference executable | Seed change | Final modes |
|---|---:|---:|
| Unmodified installed oracle | none | 4 |
| Instrumented control | 0 ULP | 4 |
| Diagnostic rebuild | +256 ULP | 5 |
| Diagnostic rebuild | -256 ULP | 5 |

Control `.mod/.shd` are byte-identical to the unmodified oracle. The perturbation
is `1.7763568394002505e-15` in x, relative `5.04e-14`, or 0.806% of the unchanged
`2.2030040619372846e-13` secant stopping tolerance. Both altered-reference fifth
roots have phase speed about 1753.9575 m/s, inside the 1800 m/s bound. Both
five-mode/63-pressure products match the unmodified Rust API and actual CLI-HDF5
at the existing strict tolerances (max pressure error `1.4725502861e-9`).
**These altered runs are diagnostic evidence, not fixed-oracle acceptance.**

Reproduce from the repository root with the pinned image already built:

```sh
# Optionally set BELLHOP_REFERENCE_IMAGE to a locally available pinned image ID.
python3 tools/reference/probe-layered-refinement.py
# The script refuses to overwrite target/reference/layered-refinement-gap.
root="$PWD/target/reference/layered-refinement-gap"
cargo build --release -p kraken-cli
target/release/kraken run "$root/pinned.env" --solver krakenc --output "$root/rust.h5"
KRAKEN_FREQUENCY_SOLVER=krakenc KRAKEN_DIFFERENTIAL_ENV="$root/pinned.env" \
KRAKEN_DIFFERENTIAL_ROOT="$root/pinned" KRAKEN_HDF5_RESULT="$root/rust.h5" \
  cargo test --release -p kraken --test differential_reference \
    multifrequency_fluid_matches_fresh_reference -- --ignored --exact --nocapture
# Expected failure: mode count, reference 4 / Rust 5.
# For diagnostic comparisons only: change reference root to plus256 or minus256.
# Unset KRAKEN_HDF5_RESULT to compare the API instead.
```

## Attempted fix and remaining decision

Changing the bottom k² arithmetic from `(omega/c)^2` to the reference expression
`omega²/c²` reduced a boundary-evaluation discrepancy but still returned five
Rust modes. It was reverted. Temporary Rust instrumentation and frozen probe
constants were also removed; no perturbation was added to production code.

There is no justified root-deletion or case-specific rejection rule from this
experiment. The Rust fifth root is also found by the reference when a tiny seed
perturbation changes its search path. Reproducing the unmodified four-root path
would require further control of floating-point trajectories, not a missing
interface condition; selecting arithmetic solely because it happens to skip
this root would be a fragile, platform-sensitive workaround.

Under the current fixed-oracle parity requirement, the blocker remains open.
Accepting this known search limitation, or adopting an independently validated
root-completeness criterion, requires an explicit scope/acceptance decision;
neither has been done here. The 41 already accepted workflows remain separate.
