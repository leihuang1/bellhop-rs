# Numerical golden fixtures

`MunkP_one_ray.env` is derived from Acoustics Toolbox `v2023.5`
`tests/Munk/MunkB_ray.env` at commit
`475108519289c6fb488b58980c644ea14eccc604`. Modifications are limited to the
title, one source depth, one launch angle, and a 10.1 km trace box.

The other `.env` files are small, project-authored, reduced, or unmodified
official cases. They cover every `N/C/P/S/Q/A` SSP algorithm, Cartesian
geometric-hat eigenrays and arrivals, and Cartesian geometric-Gaussian
arrivals, including a ray-centered geometric-hat case. The `.ray` and `.arr` files were generated from the unmodified
`v2023.5` Fortran source using GNU Fortran 14.2.0 on Apple arm64 with:

```text
-O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu
```

SHA-256:

```text
3ecccefc3f47428de53cf07672cb549cadbae425a7dc66e00909b4b209124b8c  Analytic_one_ray.env
be89ebaa6a6d237a1451f6cf4bcd0fc26f6bebd370b182205322256877830120  Analytic_one_ray.ray
afc2bd68348802ed2d7fd3ff02c8eaafdc408221af5059cf59a122c6b3d274c3  CLinear_one_ray.env
0205066cdc9a4b771c6dd1271c5fa4cd623a01339f9cdf0876a6c44f19aff272  CLinear_one_ray.ray
14553411d9ecb2c21c07457cf6fb4077782c0c06eae18949a90f05281fc60f1f  GeoHat_arrival.env
ea53e614c83bf394eed92bd88274e3e1be3bf9e23bd798b0f7c33c99d13b29e9  GeoHat_arrival.arr
f5f6fbd0be0ccfee58771d3de66ccf6989cb404a2751150165ad9613765ba1cd  GeoHat_eigen.env
6caeb6b8dd16b636cd331583796506a3f35c468a20afe00ffc7903199f213d7f  GeoHat_eigen.ray
44ec956e7619ef9290827c86f529dc6cfbbb15b9abf5253f1c3a45cbf838762e  GeoHatRay_arrival.env
2d20874e6961c36e2cb3817ea12e8f8a89e73b5fe191c63570703e85eee96c2c  GeoHatRay_arrival.arr
0bb8fde32081d577b0a07718cf94488760e7c580c6b34149f4a05911b1030d59  MunkP_one_ray.env
a0f1c76e4b11bd9075a091e4a31986b8bcc6e8fffd1dd9c60a12afeb89680b32  MunkP_one_ray.ray
0b4ea0c627430c32e10e972e829b606efb12ecd75a5c400316fec4f227b2d31a  MunkB_Arr.env
fb5630ab36dab8ebc433dfe137ac213d2a730d6dd3fa0e72528186a9607d976a  MunkB_Arr.arr
ebf1f6d0daa5a93e46cf39be4643b0fa04fef2a67c72a003264200bc9ab6d04e  N2_one_ray.env
1f67bdca4175f5fa79d16b42dd2fcaa7f9c2564ab493d3587f72001ad2dac8cb  N2_one_ray.ray
b4aee66522de2009c869c6a4decc41cd8e5c4ab9631753c0adc31804835fac77  Quadrilateral_one_ray.env
ea8ef08e631e0ec27e1f4d4eeeba6f5708a54d6d9cf6480c30509231b5b31c0a  Quadrilateral_one_ray.ssp
0df1f9168e791b00ff18412be355e792938d591e9bd3ea21f6d46bee82a3100e  Quadrilateral_one_ray.bty
4b0bbfcf50c93d6a90e7e4c282be20f41709646bb8bc7e00de775d9350a1d069  Quadrilateral_one_ray.ray
991081418e0c06cb16fcb47e28913756dcb721b2dc4b1462db8eedd49a015dc9  Spline_one_ray.env
51ca6fcb102dbf688ea9f769c8f88a5ca40c83e98e30fa299e5160227ed07524  Spline_one_ray.ray
aa2ed61ea679acd8e54fbb5122a1f612d100282caf340c4bdd39ef32de059343  calibBarr.env
fca078e8b082443cac28aefc53da9cc795379a8010ce7316f6f269c26b968d5c  calibBarr.arr
```

This small golden catches computation-order regressions. The eventual pinned
Linux differential container remains the authoritative cross-platform oracle.
