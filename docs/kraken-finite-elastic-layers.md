# Homogeneous finite elastic-layer checkpoint

Finite solids are now propagated with the pinned five-component compound matrix,
not replaced by acoustic layers or half-spaces. This bounded block covers smooth,
homogeneous solid caps above/below a contiguous N/C/P/S fluid stack:

- **KRAKENC:** multiple fluid layers and multiple finite solid caps, with complex
  P/S attenuation. Outer solid boundary must be V, R or elastic A.
- **KRAKEN:** one finite fluid layer with homogeneous solid caps; multiple-fluid
  finite-elastic combinations are explicitly rejected at Case validation / CLI
  exit 2. The previous multi-fluid elastic-half-space capability is unchanged.
- Finite cp/cs/density/loss are constant within each solid; jumps between solids
  are allowed, with full compound state carried across each interface. Roughness,
  graded solids, solids interleaved with fluids, pure-solid cases, analytic Munk,
  F/P/TRC combinations and depth-local biological loss in solids remain excluded.
- Modal/source/receiver samples remain **fluid pressure**, at absolute depths
  inside the fluid interval. No displacement/strain sampling in solids is exposed.
- An acoustic outer A behind finite solids is rejected: the pinned BC routine
  does not initialize its five-component compound state in that branch.

`ElasticLayer` carries absolute bottom depth, cp, cs, density, canonical
solve-frequency P/S dB/wavelength losses and nominal NG. `top_elastic_layers`
starts at 0; `bottom_elastic_layers` starts at `fluid_bottom_depth_m()`. The
existing water fields still define the first fluid, now possibly below a cap.
`fluid_top_depth_m()` and `fluid_bottom_depth_m()` bound acoustic samples;
`total_depth_m()` is the physical last finite interface, excluding the half-space.
Legacy SSP inheritance is retained; a fluid after a solid must explicitly reset
cs to zero where inherited records would otherwise retain shear material.

Each solid uses its own shear-based automatic mesh or explicit NG. Frequency
scaling precedes truncation, including every refinement level. The reference
Euler first step, modified midpoint, terminal filter and decimal scaling exponent
are retained. The first physical medium's h supplies Neville seeds. KRAKEN
finite-solid cases use deflated real Solve2 secants instead of acoustic inertia;
KRAKENC reuses its complex secant, first-mesh shapes/group speeds and Richardson
roots. The first two raw mesh scans and later Neville seeds stay distinct from
extrapolated eigenvalues. Neither root cropping nor normalization rescaling is
used to disguise differential failures.

## Loss semantics and budgets

KRAKEN finite-solid stiffness uses **Re(c²)**, unlike its elastic-half-space
formula's **Re(c)²**. Changing finite-solid loss can therefore change the real
spectrum and pressure, but no elastic absorption perturbation is added. Original
ice has nonzero P/S losses yet zero real-KRAKEN modal attenuation; zeroing those
losses changes its spectrum. KRAKENC includes complex stiffness and absorption.
The original half-space limitation and cs cutoff are unchanged. Group speeds
follow the reference's fluid/compressional-half-space slowness expression, not an
independent full-elastic energy/group-dispersion validation.

Limits are unchanged and shared with fluid media: at most 500 total finite media,
1M total mesh intervals per frequency/level, 100K total conceptual profile/loss
entries, 20K roots, 5M sampled shapes, 1M FIELD samples, 5M copied-input entries,
1000 frequencies. Parsing counts all raw SSP records before compressing uniform
solids; copied solid material costs seven entries per frequency. Root work counts
five compound-component steps per solid node against the existing 2.5B KRAKEN /
300M KRAKENC limits. These are logical work/storage bounds, not wall-time or native
HDF5 descriptor guarantees. Atomic publication/quota/input-protection contracts
remain unchanged.

## Complete reference evidence

Acoustics Toolbox v2023.5 commit `475108519289c6fb488b58980c644ea14eccc604`, Linux
x86-64 / GNU Fortran 12.2.0, the same flags and unmodified numerical tolerances.
25 input pairs (23 derived/constructed plus two original pairs) give **33
workflows, 36 frequency blocks, 827 modes and 9,276 pressures**. Every mode, shape,
FIELD sample and material attribute is compared through API and actual CLI/HDF5.
All `.mod/.shd` are byte-identical in three independent pinned runs. Local max
pressure error: 4.22e-8. The 149-record `golden/finite-elastic.sha256` locks inputs
and generated artifacts; Rust HDF5 is never a golden.

- OriginalElasticSediment is byte-identical `tests/TLslices/elsed.env` (5000 m
  fluid, finite 100 m cp1400/cs700/rho1.5 sediment, elastic A bottom). Both engines:
  46 modes / 501 pressures, with every interface mode retained.
- OriginalElasticIce is byte-identical `tests/TLslices/ice.env` (30 m cp3000/cs1400
  lossy solid cap, fluid 30..5000 m, acoustic A bottom). Both: 44 / 501.
- Both official FLPs are byte-identical `fieldbat.flp`, selected by upstream
  `runtests.m`. The six FiniteSingleIce/Sediment C/P/S derivatives change only title
  and interpolation. They do not replace the independently accepted originals.
- FiniteElasticBottom/Top/Both N/C/P/S, shear-only, power-law and V/R bottom cases
  exercise independently meshed/dense layered fluids and homogeneous solids in
  KRAKENC. The four-frequency power case preserves 75/50/62.5/50 Hz order/repeats.
- FiniteElasticStack adds two different solid media on each side of two fluids;
  no density averaging or interpolation crosses solid/material interfaces.

## Open numerical scope

The initial multi-fluid finite-solid **KRAKEN** port exposed new secant trajectory
failures: TopN changed 5→6 roots between meshes; ShearOnly published six versus
pinned three; Power at 62.5 Hz retained seven versus pinned six. These are not
accepted workflows and are not covered by the old layered release exception.
Those inputs now explicitly fail validation for KRAKEN rather than publish a
known unverified result; KRAKENC versions independently pass full references.
This restriction reflects unvalidated search parity, not a physical prohibition
in the original Fortran. It does not claim which search found all mathematical
roots. No temporary trace, reference seed patch, tolerance increase or work-budget
increase is promoted into goldens.

Next: validate KRAKEN coupled multi-fluid/solid search before relaxing its guard,
then graded finite solids / remaining KRAKEN elastic-top-half-space isolation.
Multi-profile FIELD, JSON/HTTP, BOUNCE generation, elastic displacement output and
time-domain synthesis remain later blocks. A successful arbitrary CLI run is
not a parity certificate.
