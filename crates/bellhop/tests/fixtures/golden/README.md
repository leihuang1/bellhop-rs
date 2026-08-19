# Numerical golden fixtures

`MunkP_one_ray.env` is derived from Acoustics Toolbox `v2023.5`
`tests/Munk/MunkB_ray.env` at commit
`475108519289c6fb488b58980c644ea14eccc604`. Modifications are limited to the
title, one source depth, one launch angle, and a 10.1 km trace box.

The other `.env` files are small, project-authored or reduced official cases
covering every `N/C/P/S/Q/A` SSP algorithm. Their `.ray` files and
`MunkP_one_ray.ray` were generated from the
unmodified `v2023.5` Fortran source using GNU Fortran 14.2.0 on Apple arm64
with:

```text
-O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu
```

SHA-256:

```text
3ecccefc3f47428de53cf07672cb549cadbae425a7dc66e00909b4b209124b8c  Analytic_one_ray.env
be89ebaa6a6d237a1451f6cf4bcd0fc26f6bebd370b182205322256877830120  Analytic_one_ray.ray
afc2bd68348802ed2d7fd3ff02c8eaafdc408221af5059cf59a122c6b3d274c3  CLinear_one_ray.env
0205066cdc9a4b771c6dd1271c5fa4cd623a01339f9cdf0876a6c44f19aff272  CLinear_one_ray.ray
0bb8fde32081d577b0a07718cf94488760e7c580c6b34149f4a05911b1030d59  MunkP_one_ray.env
a0f1c76e4b11bd9075a091e4a31986b8bcc6e8fffd1dd9c60a12afeb89680b32  MunkP_one_ray.ray
ebf1f6d0daa5a93e46cf39be4643b0fa04fef2a67c72a003264200bc9ab6d04e  N2_one_ray.env
1f67bdca4175f5fa79d16b42dd2fcaa7f9c2564ab493d3587f72001ad2dac8cb  N2_one_ray.ray
b4aee66522de2009c869c6a4decc41cd8e5c4ab9631753c0adc31804835fac77  Quadrilateral_one_ray.env
ea8ef08e631e0ec27e1f4d4eeeba6f5708a54d6d9cf6480c30509231b5b31c0a  Quadrilateral_one_ray.ssp
0df1f9168e791b00ff18412be355e792938d591e9bd3ea21f6d46bee82a3100e  Quadrilateral_one_ray.bty
4b0bbfcf50c93d6a90e7e4c282be20f41709646bb8bc7e00de775d9350a1d069  Quadrilateral_one_ray.ray
991081418e0c06cb16fcb47e28913756dcb721b2dc4b1462db8eedd49a015dc9  Spline_one_ray.env
51ca6fcb102dbf688ea9f769c8f88a5ca40c83e98e30fa299e5160227ed07524  Spline_one_ray.ray
```

This small golden catches computation-order regressions. The eventual pinned
Linux differential container remains the authoritative cross-platform oracle.
