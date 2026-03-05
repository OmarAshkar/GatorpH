# Integrate Area Under pH Threshold

Calculates the area under the pH curve below a specified pH threshold
using the trapezoidal rule. It assumes the flipped shape

## Usage

``` r
integratepHArea(
  x,
  y,
  ph_threshold = 5.4,
  time_start = 0,
  time_end = 50,
  method = "linear",
  interpolate = TRUE
)
```

## Arguments

- x:

  Vector of time points.

- y:

  Vector of pH values corresponding to the time points.

- ph_threshold:

  pH threshold (default is 5.4).

- time_start:

  Start time for integration (default is 0).

- time_end:

  End time for integration (default is 50).

## Value

Numeric value representing the area under the pH curve below the
threshold.

## Author

Omar I. Elashkar
