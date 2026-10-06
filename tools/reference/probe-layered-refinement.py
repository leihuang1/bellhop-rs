#!/usr/bin/env python3
"""Diagnostic only: perturb one pinned KRAKENC secant seed, not the oracle/goldens."""
import io
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tarfile

repo = Path(__file__).resolve().parents[2]
out = repo / "target/reference/layered-refinement-gap"
out.mkdir(parents=True, exist_ok=False)  # Keep old evidence rather than overwrite it.
image = os.environ.get("BELLHOP_REFERENCE_IMAGE", "pelagic-reference:v2023.5-amd64")
docker = ["docker", "run", "--rm", "--platform", "linux/amd64"]
archive = subprocess.run(
    docker + ["--entrypoint", "tar", image, "cf", "-", "-C",
              "/opt/acoustics-toolbox/Kraken", "krakenc.f90"],
    check=True, capture_output=True,
).stdout
with tarfile.open(fileobj=io.BytesIO(archive)) as files:
    source = files.extractfile("krakenc.f90").read().decode()
fixtures = repo / "crates/kraken/tests/fixtures"
env = (fixtures / "LayeredFluidThreeWide.env").read_text()
needle = "1400.0 1800.0\n0.0\n"
assert env.count(needle) == 1
env = env.replace(needle, "1400.0 1800.0\n1000.0\n")
call = "       CALL RootFinderSecant( x, Tolerance, Iteration, MaxIteration, ErrorMessage, Funct )"
assert source.count(call) == 1
for stem, ulps in [("pinned", None), ("control", 0), ("plus256", 256), ("minus256", -256)]:
    (out / f"{stem}.env").write_text(env)
    shutil.copyfile(fixtures / "LayeredFluidThreeWide.flp", out / f"{stem}.flp")
    if ulps is not None:
        # Do not change tolerance, coefficients, deflated roots or any other seed.
        patch = f"""       IF (iSet == 2 .AND. mode == 5) THEN
          WRITE(*,'(A,3ES26.17)') '[DEBUG-gap] original seed/tolerance', x, Tolerance
          x = x + ({ulps}.0D0) * SPACING(DBLE(x))
          WRITE(*,'(A,2ES26.17)') '[DEBUG-gap] perturbed seed', x
       END IF
{call}
       WRITE(*,'(A,2I6,2ES26.17)') '[DEBUG-gap] root', iSet, mode, x"""
        (out / f"{stem}.f90").write_text(source.replace(call, patch))

# The installed oracle and Docker image stay unchanged. Rebuild only in this
# disposable container; altered outputs are diagnostic artifacts, never goldens.
subprocess.run(docker + ["--volume", f"{out}:/work", "--workdir", "/work",
    "--entrypoint", "/bin/sh", image, "-c", """
set -eu
/usr/local/bin/krakenc-fortran pinned
/usr/local/bin/field-fortran pinned
for variant in control plus256 minus256; do
  cp /work/$variant.f90 /opt/acoustics-toolbox/Kraken/krakenc.f90
  rm -f /opt/acoustics-toolbox/Kraken/krakenc.o
  make -C /opt/acoustics-toolbox/Kraken krakenc.exe FC=gfortran \
    FFLAGS='-O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu -I../misc' \
    > /work/$variant.build.log 2>&1
  /opt/acoustics-toolbox/Kraken/krakenc.exe $variant
  /usr/local/bin/field-fortran $variant
 done
"""], check=True)
for ext in ["mod", "shd"]:
    assert (out / f"pinned.{ext}").read_bytes() == (out / f"control.{ext}").read_bytes(), ext
for stem, expected in [("pinned", 4), ("control", 4), ("plus256", 5), ("minus256", 5)]:
    data = (out / f"{stem}.mod").read_bytes()
    record = 4 * struct.unpack_from("<I", data)[0]
    # This fixture has one frequency, three media, eight mode-sample depths.
    assert struct.unpack_from("<4i", data, 84) == (1, 3, 8, 8)
    modes = struct.unpack_from("<i", data, 5 * record)[0]
    assert modes == expected, (stem, modes)
    print(f"{stem}: {modes} modes")
print(f"Diagnostic evidence: {out}; this does NOT establish pinned parity.")
