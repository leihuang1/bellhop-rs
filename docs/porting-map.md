# Fortran to Rust porting map

| Reference | Rust | Status |
|---|---|---|
| `Bellhop/ReadEnvironmentBell.f90::ReadEnvironment` | `bellhop::legacy::env` | Parser verified |
| `Bellhop/ReadEnvironmentBell.f90::ReadTopOpt` | `bellhop::legacy::env` | Parser verified |
| `misc/SourceReceiverPositions.f90` | `bellhop::legacy::env` | Parser verified |
| `Bellhop/angleMod.f90::ReadRayElevationAngles` | `bellhop::legacy::env` | Parser verified |
| `misc/subtabulate.f90` | `bellhop::legacy::env` | Parser verified |
| `Bellhop/sspMod.f90::Quad` input loading | `bellhop::legacy::auxiliary` | Parser verified |
| `Bellhop/bdryMod.f90::ReadATI` | `bellhop::legacy::auxiliary` | Parser verified |
| `Bellhop/bdryMod.f90::ReadBTY` | `bellhop::legacy::auxiliary` | Parser verified |
| `misc/RefCoef.f90::ReadReflectionCoefficient` (`.brc`/`.trc`) | `bellhop::legacy::auxiliary` | Parser verified |
| `misc/beampattern.f90::ReadPat` | `bellhop::legacy::auxiliary` | Parser verified |
| `Bellhop/Step.f90::Step2D` | `bellhop::solver::integrator::step_2d` | Golden verified |
| `Bellhop/Step.f90::ReduceStep2D` | `bellhop::solver::integrator::reduce_step` | Golden verified |
| `Bellhop/sspMod.f90` `N/C/P/S/Q/A` models | `bellhop::solver::ssp` | Golden verified |
| `Bellhop/bdryMod.f90::ComputeBdryTangentNormal` | `bellhop::solver::boundary` | Differentially exercised |
| `bellhop.f90::TraceRay2D` | `bellhop::solver::trace_ray` | Golden verified |
| `bellhop.f90::Reflect2D` | `bellhop::solver::reflection::reflect_2d` | Differentially exercised |
| `Bellhop/WriteRay.f90` | `bellhop::solver::{RayTrajectory, RayPoint}` | In-memory equivalent; legacy writer intentionally omitted |
| `Bellhop/influence.f90` | — | Not started |
