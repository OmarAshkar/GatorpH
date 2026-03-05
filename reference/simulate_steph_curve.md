# Simulate Oral pH-Time Curve

Simulates oral pH-time curves based on a given PK-PD model.

## Usage

``` r
simulate_steph_curve(
  model,
  time = seq(0, 50, by = 0.1),
  dose = 100,
  baseline_time = -5,
  nsub = 1,
  baseline = 7,
  step = 0.1,
  group = "A",
  ignoreBSV = TRUE,
  ignoreRUV = TRUE
)
```

## Arguments

- model:

  rxode2 model object representing the PK-PD model.

- time:

  Vector of time points for simulation (default is seq(0, 50, by =
  0.1)).

- dose:

  Dose amount (default is 100).

- baseline_time:

  Time of baseline to adjust time points (default is -5).

- nsub:

  Number of subjects to simulate (default is 1).

- baseline:

  Baseline pH value. Default is 7. Must be either 1 or same number as
  number of subjects.

- step:

  Time step for simulation (default is 0.1).

- group:

  Group identifier for the subjects (default is "A").

- ignoreBSV:

  Logical indicating whether to ignore between-subject variability
  (default is TRUE).

- ignoreRUV:

  Logical indicating whether to ignore residual unexplained variability
  (default is TRUE).

## Value

A data frame containing the simulated pH-time data.

## Author

Omar I. Elashkar
