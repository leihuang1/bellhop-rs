# Wide three-layer KRAKENC refinement: open parity gap

**Not fixed; not an accepted workflow.** Following the user's scope review,
this documented exception no longer blocks review/release of PR #25's validated
fluid-stack subset. It is not a numerical fix or a general waiver for other
parity failures. No solver, reference, tolerance or accepted golden was changed.

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

## Attempted fix

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

## Scope decision and impact assessment

The user authorized a scope adjustment if its practical impact on this stage
was limited. Based on the existing evidence, accept a **release-gate exception
for this documented workflow**, not acceptance of its numerical result. The
41 validated workflows remain the strict acceptance set. Arbitrary complex
fluid stacks are not certified as fixed-oracle-equivalent merely because their
inputs parse and the solver finishes.

### Physical impact is not established as small

Pinned `KrakenField/EvaluateMod.f90::Evaluate` forms the coherent pressure as a
sum of modal contributions proportional to
`C_m * phi_m(z) * exp(-i*k_m*r) / sqrt(k_m)` for a point source, followed by
cylindrical spreading. Adding/removing one root can change interference and
transmission-loss minima; mode-count fraction is not a pressure-error estimate.

The diagnostic fifth root has `-Im(k) ≈ 1.7555e-4 /m`, so its attenuation-only
amplitude factor is `exp(-1.7555e-4*r)` and its e-folding distance is about
5.7 km. `LayeredFluidThreeWide.flp` samples 0.5, 1 and 2 km, where that factor
is roughly 0.92, 0.84 and 0.70. These are analytical attenuation factors, **not**
measured fractions of the total pressure: source/receiver coupling, modal
normalization, phase and spreading also matter. RMax=1000 km controls mesh
refinement; it is not the FIELD receiver range and cannot justify neglecting
this mode. No new pressure-error calculation was performed for this review.

The tiny seed perturbation explains why strict search-path parity is fragile;
it does not bound the physical effect, prove Rust is universally more accurate,
or prove this behavior occurs only in the one known fixture.

### Impact on this delivery is bounded

- No previously accepted input or capability is removed. Independent refined
  three-layer, wide base-mesh, layered loss/broadband/leaky, V/R/A and interface
  workflows retain their existing complete comparisons. This is not a waiver
  for all three-layer refinement, all wide spectra or all KRAKENC failures.
- API, CLI/HDF5 schema, atomic publication, quotas, work ceilings, physics and
  the original `double` count-change rejection are unchanged. Layer-density
  gradients, elasticity and multi-profile FIELD remain out of scope.
- The known case still may return five modes and CLI exit 0, without a warning
  or a parity-certification flag. This exception is documentary, **not** runtime
  rejection/protection. Users requiring exact legacy reproduction must exclude
  this workflow and independently validate unverified configurations.
- Changing cHigh to 1700 or RMax to 0 creates a different calculation. Their
  accepted fixtures are not automatic substitutes for a user's wide refined
  problem. HDF5 provenance identifies inputs; it does not certify parity.

### Acceptance rules retained

1. All 41 accepted workflows must continue to match the **unmodified** pinned
   oracle, for every reference mode and pressure, at unchanged tolerances.
2. This wide refined workflow stays outside acceptance. Keep its diagnosis,
   reference-only perturbation tool and explicitly ignored failing reproducer;
   none may be relabelled a passing test or included in acceptance totals.
3. Do not trim reference/Rust roots, promote perturbed outputs to goldens, or
   replace the fixed oracle with an unvalidated root-completeness criterion.
4. Any new parity failure, regression in an accepted workflow, or evidence of
   broader practical impact requires a separate scope review; no blanket waiver.
5. Promotion of this workflow requires fresh unmodified-oracle API/CLI-HDF5
   evidence and resolution of the count-change contract, or a separately
   approved numerical acceptance contract. The current exception does neither.

This permits PR #25 to leave draft for normal review; it does not authorize
merging or claim a complete KRAKENC rewrite. This assessment used the existing
source, fixtures and diagnostic evidence only; no new code or tests were run.
