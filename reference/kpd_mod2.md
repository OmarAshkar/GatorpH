# KPD Model with Linear Parameterization

## Usage

``` r
kpd_mod2(
  edk50 = 10,
  kde = 0.1,
  kd = 0.5,
  ks = 3.5,
  gamma = 1,
  eta.edk50 = 0.1,
  eta.kde = 0.1,
  eta.kd = 0.1,
  eta.ks = 0.1,
  sigma_add = 0.01
)
```

## Arguments

- edk50:

  Dose producing 50

  kdeElimination rate constant (1/h) from the virtual compartment KDE
  (default is 0.1).

  kdElimination rate constant (1/h) from the response compartment
  (default is 0.5).

  ksZero-order production rate constant (default is 3.5).

  gammaHill coefficient (default is 1).

  eta.edk50Inter-individual variability (IIV) on edk50 (default is 0.1).

  eta.kdeIIV on kde (default is 0.1).

  eta.kdIIV on kd (default is 0.1).

  eta.ksIIV on ks (default is 0.1).

  sigma_addAdditive residual error (default is 0.01).

An rxode2 model object representing the KPD model with linear
parameterization. KPD Model with Linear Parameterization Omar I.
Elashkar
