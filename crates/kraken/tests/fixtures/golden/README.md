# Pekeris reference artifacts

Generated only by Acoustics Toolbox v2023.5, commit
`475108519289c6fb488b58980c644ea14eccc604`, not by the Rust solver.
The pinned `tools/reference/Dockerfile` uses Linux x86-64, GNU Fortran
12.2.0 (Debian 12.2.0-14+deb12u1), and
`-O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu`.

| Case | Coverage | Modes / pressure samples |
|---|---|---|
| Pekeris | original 50 Hz / 100 m baseline, NG=1000 | 3 / 9 |
| PekerisFiltered | equivalent SSP slash spellings, phase-speed selection, FIELD limit=1, two interpolated source depths, tilted receiver array | 2 / 18 |
| PekerisDense | 75 Hz / 120 m, different speeds/densities, NG=4000, intermediate SSP point, surface/interface receivers, zero range | 7 / 30 |

Regenerate and compare from the repository root:

```sh
for case in Pekeris PekerisFiltered PekerisDense; do
  tools/reference/compare-kraken.sh "crates/kraken/tests/fixtures/$case.env"
done
```

Fresh outputs are in `target/reference/$case-kraken/`. To update goldens,
review the numerical changes first, copy that case's `.mod`, `.shd`, and
`.prt` here, trim trailing whitespace from `.prt` lines, and update these hashes.
The binary files are unmodified; only text whitespace is normalized. `.prt`
contains nondeterministic CPU timings; comparisons read its numerical modal
table, not the entire text.
The binary layouts are interpreted only by the test-only reader in
`../../pekeris_reference.rs`.

The same comparator runs against committed goldens in ordinary tests and
fresh Fortran files in CI. It checks mode/grid counts, coordinates, frequency,
printed and binary wavenumbers, attenuation, phase/group speed, unit-phase
aligned mode shapes (without amplitude rescaling), and every complex pressure.
Tolerances remain `5e-10` for printed k, `1e-6` for shapes, and `2e-6` for
pressure. Binary k uses `2e-8` because `.mod` stores complex32; group speed
uses `0.005` because `.prt` prints only two decimal places.

## SHA-256

Paths below are relative to the parent `fixtures/` directory.

```text
affd23e8f6e0096498a18ae82ccab87f66d430faa9b5ea652238d88a2eaf23fa  Pekeris.env
38527d0ab6d5360f3ab424b75c8a0c81c375cebd446886a640f639540a03c7d7  Pekeris.flp
a5b9f14bd2d15b9667bedd22a7c9ecf1af2101000a6236013968556b0de32a9c  PekerisDense.env
699f29da292e40a739bab73cf0675c7ecd8aa7f3cf4c808e4b54883e2894845d  PekerisDense.flp
e77a71442d8c519363726e672310a9ab8fdd629ab2672552aba0dd97b504160b  PekerisFiltered.env
2f6e92633dece71a5ff413d40829cf95b4cb243609bc0ee22aa56a410f27e67e  PekerisFiltered.flp
6e31524a1cd0a7feb0adb25646b8aba6dd5fce7d76d80e3d79e3cb3009a16d28  golden/Pekeris.mod
063fc54d0442a79c8512a803828a078fb647a3af024a6ab45d874dba5b6db2af  golden/Pekeris.prt
cf0c549a05dba48085038a202daf422267f89001ef119bfbfa177661b0767628  golden/Pekeris.shd
01471b22453d649c7cf6a2b022ce912891eaa672864f0c5cd0f4a34f92bf7722  golden/PekerisDense.mod
6c2d059a0d63122d9517e9725cf5f14bde0a132b7e6c5cf3c097a72895825328  golden/PekerisDense.prt
d389df249da65b4eac0073640fc73b152c4239e4ea9f1c1829316bd3fe0099b6  golden/PekerisDense.shd
3130ee57251111c43ac0e75ae0eb785057b4377db8ed8b6666ed2e82425d82ff  golden/PekerisFiltered.mod
ccff31b052f385edd1a68c0a7bbe3214259f9c5d4d06993858d6ad5e7712091d  golden/PekerisFiltered.prt
869ed1f696aa63648b20cc1812bf791c416386d573396e9b4f7508352ecda586  golden/PekerisFiltered.shd
```
