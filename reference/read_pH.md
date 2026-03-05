# Read pH Data

Reads pH data from a CSV or Excel file and performs basic validation.

## Usage

``` r
read_pH(file_path, baseline = 7, baseline_time = -5)
```

## Arguments

- file_path:

  Path to the data file.

- baseline:

  Baseline pH value to assign for each subject if not present in the
  data (default is 7).

- baseline_time:

  Time of baseline to adjust time points (default is -5). If NULL, no
  adjustment is made.

## Value

A data frame containing the pH data.

## Author

Omar I. Elashkar
