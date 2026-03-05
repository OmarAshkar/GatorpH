# Run direct estimation for pH data

Runs direct estimation calculations for pH data, including minimum pH,
time to minimum pH, time under pH threshold, and AUC under pH threshold.

## Usage

``` r
run_direct_estimation(
  x,
  ph_threshold = 5.4,
  time_start = 0,
  time_end = 50,
  method = "linear"
)
```

## Arguments

- x:

  Data frame containing pH data.

- ph_threshold:

  pH threshold for calculations.

- time_start:

  Start time for area and time under pH calculation.

- time_end:

  End time for area and time under pH calculation.

- method:

  Method for AUC calculation ("linear", "log", "linear_down_log_up").

## Value

A data frame with subject and group pH parameters
