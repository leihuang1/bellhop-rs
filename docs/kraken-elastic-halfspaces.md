# Elastic half-space checkpoint

This is the first elastic capability block, **not finite elastic-layer support**.
The oracle remains unmodified Acoustics Toolbox v2023.5, commit
`475108519289c6fb488b58980c644ea14eccc604`, Linux x86-64 / GNU Fortran 12.2.0,
with the existing pinned flags and unchanged numerical tolerances.

## Supported boundary subset

- N/C/P/S constant-density finite fluid stacks over an elastic A bottom:
  KRAKEN and KRAKENC modes, coherent FIELD and actual CLI/HDF5.
- Elastic A top, including simultaneous elastic A bottom: **KRAKENC only**.
- Positive shear speed, positive density and bulk modulus (`cp² > 4/3 cs²`),
  finite nonnegative compressional/shear losses with `Im(c) <= Re(c)`.
- Existing frequency ordering, independent fluid meshes, refinement and shared
  budgets remain in force. No new root-work allowance or partial publication.
- Source/receiver/modal depths remain inside the finite **fluid** stack. These
  are pressure eigenfunctions, not elastic displacement/strain samples.
- Analytic Munk, F/P/TRC combinations, finite solids, roughness and multi-profile
  FIELD are not enabled by this block. Unsupported combinations fail explicitly.

`Boundary::ElasticHalfSpace { shear_sound_speed_mps,
shear_attenuation_db_per_wavelength }` supplies shear material; the existing
surface/bottom sound-speed, density and loss fields supply compressional
material. Both losses are canonical solve-frequency dB/wavelength. Legacy A
records with cs>0 select this variant; conversion is repeated per frequency,
including shear power-law/volume loss. A fluid cs=0 record still forbids shear
loss. No elastic boundary is replaced with a fictitious fluid half-space.

The implementation ports `BCImpedance{c}Mod`'s pressure/normal-velocity impedance
and its normalization derivative. Real KRAKEN adds the elastic boundary's
mode-count contribution and caps effective cHigh at the shear speed; the
reference's 0.85*cMin adjustment permits interface-wave roots. KRAKENC retains
complex P/S radiation roots and bounded deflated secant search. Root completeness
and arbitrary-input parity are not inferred from successful fixture comparisons.

**Reference limitation retained:** KRAKEN's real elastic boundary ignores
elastic material attenuation, even when the input records nonzero P/S loss.
Fluid volume loss still acts. Use KRAKENC for elastic attenuation. HDF5 records
requested losses and `*_elastic_attenuation_model=reference_real|complex` so this
is not mistaken for implemented elastic-loss perturbation. Group speeds follow
the pinned compressional half-space slowness formula; matching them is not an
independent validation of full elastic energy/group dispersion.

## Complete acceptance evidence

18 input pairs: **15 derived** and **3 unmodified original environments with
their official shared FIELD geometry**. They produce 27 workflows, 33 frequency
blocks, 415 modes and 4,707 pressures. Every mode/sample/pressure is compared
through the API and actual CLI/HDF5. Each `.mod/.shd` is byte-identical in three
independent pinned runs. Local maximum pressure error is 6.67e-8. No mode is
trimmed and no numeric tolerance increased; shear-only loss uses the existing
lossy imaginary-wavenumber tolerance, not the lossless exact-zero comparison.

| Inputs | Engines | Modes per frequency | Pressures per engine |
|---|---|---|---:|
| ElasticHalfBottomN/C/P/S | both, each input | 5 | 63 |
| ElasticHalfTopN/C/P/S | KRAKENC, each input | 6 | 63 |
| ElasticHalfBothN/C/P/S | KRAKENC, each input | 5 | 63 |
| ElasticHalfLeaky | KRAKENC | 5 | 63 |
| ElasticHalfPower, 75/50/62.5/50 Hz | both | 7/5/6/5 | 252 |
| ElasticHalfShearOnly, coherent line source | both | 5 | 63 |
| OriginalElasticScholte | both | 45 | 501 |
| OriginalElasticNormal | both | 44 | 501 |
| OriginalElasticFlused | both | 46 | 501 |

Original environments are byte copies of `tests/TLslices/{scholte,normal,flused}.env`.
All three FLPs are byte copies of `tests/TLslices/fieldbat.flp`, selected by the
unmodified upstream `runtests.m`. In particular original normal retains bottom
cs=2000; it is independent of the previous shear-removed `LayeredNormalization`.
Fresh CI checks the input bytes against the installed pinned source before
running the complete comparisons. See [golden provenance](../crates/kraken/tests/fixtures/golden/README.md)
and the 117-record `golden/elastic-halfspace.sha256` manifest. Rust HDF5 is
readback evidence, never a replacement golden.

Ordinary tests check material validation/source locations, per-frequency shear
conversion, repeated frequencies, reference real-loss omission/cutoff, HDF5
metadata and preservation/cleanup on input, numerical and quota failures.
The merged layered-review diagnostic issue is also locked down: an additional
layer's interpolation error no longer falsely blames the first water layer.

## Work deliberately left open

**KRAKEN elastic top is rejected at Case validation / CLI exit 2**, not returned
as an unverified successful result. On the constructed ElasticHalfTopN input,
the initial real port's mode count jumps 3→1 at adjacent floating-point x around
0.04252483222633925; its first putative root fails inverse iteration. Pinned
Fortran itself counts five roots on its first meshes but publishes four after
refinement, whereas KRAKENC publishes six. The current per-root bisection cannot
be assumed equivalent to Fortran's coupled interval isolation/Brent/Solve2
trajectory. This motivates the explicit restriction; it does not establish
which search found every mathematical root. Temporary Rust tracing was removed.
No additional release-parity exception is introduced.

A separate constructed slow-shear Scholte probe (cp=3000, cs=1000, 100 m water)
reached a group speed near 17344.62 m/s in both ports, while the pinned G14.6
print rounded to 17344.6. It did not pass the existing 0.005 comparison. A cp=1800
probe also failed the strict checks (KRAKEN inverse iteration; KRAKENC print
precision). Neither probe is a committed/accepted fixture; the comparator was
not relaxed and neither proves low-speed interface-wave acceptance.

The subsequent [homogeneous finite-cap block](kraken-finite-elastic-layers.md)
now validates compound-matrix transfer and original `elsed/ice` through both
engines, FIELD and actual CLI/HDF5. Graded/interleaved solids and KRAKEN
elastic-top-half-space isolation remain later work. Multi-profile FIELD, JSON/HTTP,
BOUNCE generation and time-domain synthesis remain later work. The existing
[wide layered refinement exception](kraken-layered-refinement-gap.md) is unchanged
and is not a waiver for new elastic failures.
