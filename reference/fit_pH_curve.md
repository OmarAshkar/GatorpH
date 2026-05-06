# Fit pH Curve using NLME

Fits a pharmacodynamic model to pH data using nonlinear mixed-effects
modeling. If \`stratify\` is TRUE, fits separate models for each group.

## Usage

``` r
fit_pH_curve(
  data,
  model,
  stratify = FALSE,
  estmethod = "focei",
  cov_params = c("kd", "kde", "edk50", "gamma"),
  cov_fixedeffects = c("t.kd", "t.kde", "t.edk50", "t.gamma"),
  include_gamma = TRUE
)
```

## Arguments

- data:

  Data frame containing pH data. Must include columns: pH, time, id,
  group.

- model:

  rxode2/nlmixr2 model to fit.

- stratify:

  Logical indicating whether to fit separate models for each group
  (default is FALSE).

- estmethod:

  Estimation method to use (default is "focei").

- cov_params:

  Character vector of covariate parameters to include in the model
  (default is c("kd", "kde", "edk50", "gamma")).

- cov_fixedeffects:

  Character vector of fixed effect names corresponding to the covariate
  parameters (default is c("t.kd", "t.kde", "t.edk50", "t.gamma")).

- include_gamma:

  logical indicating whether to include the gamma parameter in the model
  (default is TRUE).

## Value

nlmixr2 fit object.

## Author

Omar I. Elashkar
