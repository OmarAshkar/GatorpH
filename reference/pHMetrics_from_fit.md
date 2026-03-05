# Get pH Metrics from NLME Fit

Get pH Metrics from NLME Fit

## Usage

``` r
pHMetrics_from_fit(
  x,
  ph_threshold = 5.4,
  time_start = 0,
  time_end = 50,
  step = 0.1,
  dose = 100,
  onlymean = TRUE,
  plot = FALSE,
  stratify_by = "None"
)
```

## Arguments

- x:

  nlmixr2 fit object.

- ph_threshold:

  pH threshold for calculations (default is 5.4).

- time_start:

  Start time for simulation.

- time_end:

  End time for simulation.

- step:

  Time step for simulation.

- dose:

  Dose amount to administer at time 0.

- onlymean:

  Logical indicating whether to simulate only the mean profile (default
  is TRUE).

- plot:

  Logical indicating whether to plot the simulated pH profiles (default
  is FALSE).

- stratify_by:

  Stratification option for plotting ("None", "Subject", "Group").

## Value

A data frame containing pH metrics calculated from the simulated pH
profiles.

## Author

Omar I. Elashkar
