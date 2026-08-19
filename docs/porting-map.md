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
| `Bellhop/Step.f90` | — | Not started |
| `Bellhop/sspMod.f90` numerical models | — | Not started |
| `Bellhop/influence.f90` | — | Not started |
