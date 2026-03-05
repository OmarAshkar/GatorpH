# Fit pH Curve using NLME

Fits a pharmacodynamic model to pH data using nonlinear mixed-effects
modeling. If \`stratify\` is TRUE, fits separate models for each group.

## Usage

``` r
fit_pH_curve(
  data,
  model,
  amt,
  stratify = FALSE,
  estmethod = "focei",
  dose_time = 5
)
```

## Arguments

- data:

  Data frame containing pH data. Must include columns: pH, time, id,
  group.

- model:

  rxode2/nlmixr2 model to fit.

- amt:

  Dose amount to administer at time 0.

- stratify:

  Logical indicating whether to fit separate models for each group
  (default is FALSE).

- estmethod:

  Estimation method to use (default is "focei").

- dose_time:

  Time (positive) of baseline to add before time 0 (default is 5).

## Value

nlmixr2 fit object.

## Author

Omar I. Elashkar
