# Plot pH vs Time with Threshold

Generates a plot of pH values over time for different subjects/groups,
highlighting a specified pH threshold.

## Usage

``` r
plot_pH_time(
  res,
  ph_threshold = 5.4,
  show_id = TRUE,
  stratify_by = "None",
  showAvg = FALSE,
  showDosing = FALSE
)
```

## Arguments

- res:

  Data frame containing pH data.

- ph_threshold:

  pH threshold to highlight on the plot (default is 5.4).

- show_id:

  Logical indicating whether to show subject IDs in the legend (default
  is TRUE).

- stratify_by:

  Variable to stratify the plot by (default is "None"). Options are
  "None", "Subject", or "Group".

- showAvg:

  Logical indicating whether to show the average pH curve across
  subjects/groups (default is FALSE).

- showDosing:

  Logical indicating whether to show dosing times on the plot (default
  is FALSE).

## Value

A ggplot2 object representing the pH vs time plot.

## Author

Omar I. Elashkar
