#' KPD Model
#' @param edk50 Dose producing 50% of Emax (default is 10).
#' @param kde Elimination rate constant (1/h) from the virtual compartment KDE (default is 0.1).
#' @param kd Elimination rate constant (1/h) from the response compartment (default is 0.5).
#' @param ks Zero-order production rate constant (default is 3.5).
#' @param gamma Hill coefficient (default is 1).
#' @param eta.edk50 Inter-individual variability (IIV) on edk50 (default is 0.1).
#' @param eta.kde IIV on kde (default is 0.1).
#' @param eta.kd IIV on kd (default is 0.1).
#' @param eta.ks IIV on ks (default is 0.1).
#' @param sigma_add Additive residual error (default is 0.01).
#' @return An rxode2 model object representing the KPD model.
#' @author Omar I. Elashkar
#' @export
kpd_mod <- function(
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
) {
  rxode2::ini({
    t.edk50 = log(edk50) # dose producing 50% of Emax
    t.kde = log(kde) # elimination rate constant (1/h) from the virtual compartment KDE
    t.kd = log(kd)
    t.ks = log(ks)
    t.gamma = log(gamma) # Hill coefficient

    eta.edk50 ~ eta.edk50
    eta.kde ~ eta.kde
    eta.kd ~ eta.kd
    eta.ks ~ eta.ks
    eta.gamma ~ 0.0 # No IIV on Gamma

    sigma_add <- sigma_add
  })

  rxode2::model({
    edk50 <- exp(t.edk50 + eta.edk50)
    kde <- exp(t.kde + eta.kde)
    kd <- exp(t.kd + eta.kd)
    ks <- exp(t.ks + eta.ks)
    gamma <- exp(t.gamma + eta.gamma)

    # depot(0) = amt
    # ks <- baseline * kd
    resp(0) = ks / kd # Initial condition for response

    d / dt(depot) = -kde * depot

    IR = kde * depot # drug elimiation rate
    d / dt(resp) = ks *
      (1 - (IR)**gamma / (edk50**gamma + (IR)**gamma)) -
      kd * resp

    resp ~ add(sigma_add)

    auc(0) = 0
    # area under the curve below pH
    if (resp < 5.4) {
      d / dt(auc) = resp # accumulate AUC when below threshold
    }

    # time below pH
    # time_under_ph(0) = 0
    # if (resp < 5.4) {
    #     d / dt(time_under_ph) = 0.1
    # }
  })
}


#' KPD Model with Linear Parameterization
#' @param edk50 Dose producing 50% of Emax (default is 10).
#' @param kde Elimination rate constant (1/h) from the virtual compartment KDE (default is 0.1).
#' @param kd Elimination rate constant (1/h) from the response compartment (default is 0.5).
#' @param ks Zero-order production rate constant (default is 3.5).
#' @param gamma Hill coefficient (default is 1).
#' @param eta.edk50 Inter-individual variability (IIV) on edk50 (default is 0.1).
#' @param eta.kde IIV on kde (default is 0.1).
#' @param eta.kd IIV on kd (default is 0.1).
#' @param eta.ks IIV on ks (default is 0.1).
#' @param sigma_add Additive residual error (default is 0.01).
#' @return An rxode2 model object representing the KPD model with linear parameterization.
#' @author Omar I. Elashkar
#' @export
kpd_mod2 <- function(
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
) {
  rxode2::ini({
    t.edk50 = log(edk50) # dose producing 50% of Emax
    t.kde = log(kde) # elimination rate constant (1/h) from the virtual compartment KDE
    t.kd = log(kd)
    t.ks = log(ks)
    t.gamma = log(gamma) # Hill coefficient

    eta.edk50 ~ eta.edk50
    eta.kde ~ eta.kde
    eta.kd ~ eta.kd
    eta.ks ~ eta.ks
    eta.gamma ~ 0.0 # No IIV on Gamma

    sigma_add <- sigma_add
  })

  rxode2::model({
    edk50 <- exp(t.edk50 + eta.edk50)
    kde <- exp(t.kde + eta.kde)
    kd <- exp(t.kd + eta.kd)
    ks <- exp(t.ks + eta.ks)
    gamma <- exp(t.gamma + eta.gamma)

    # depot(0) = amt
    # ks <- baseline * kd
    resp(0) = ks / kd # Initial condition for response

    d / dt(depot) = -kde * depot
    # IR = kde * depot
    ## A / A50 + A
    d / dt(resp) = ks *
      (1 - depot**gamma / (edk50**gamma + depot**gamma)) -
      kd * resp

    resp ~ add(sigma_add)
  })
}


getnoVarIds <- function(simRes) {
  if (is.null(simRes$id)) {
    simRes <- dplyr::mutate(simRes, id = 1)
  }
  novar <- simRes |>
    dplyr::group_by(.data$id) |>
    dplyr::summarize(sd = sd(resp)) |>
    dplyr::filter(sd == 0) |>
    dplyr::pull("id")
  nanids <- simRes |>
    dplyr::group_by(.data$id) |>
    dplyr::summarize(anynan = any(is.na(resp))) |>
    dplyr::filter(anynan) |>
    dplyr::pull("id")
  unique(c(novar, nanids))
}

#' Calculate AUC
#' Linear to tmin, log after tmin
#' @noRd
get_auc_linear_down_log_up <- function(time, conc) {
  # Check input validity
  if (length(time) != length(conc)) {
    stop("time and conc must have same length")
  }
  if (any(diff(time) <= 0)) {
    stop("time must be strictly increasing")
  }

  auc <- 0
  t_min <- calc_t_min(data.frame(time = time, pH = conc))$t_min

  for (i in 2:length(time)) {
    dt <- time[i] - time[i - 1]
    c1 <- conc[i - 1]
    c2 <- conc[i]

    # log-up case (after t-min, exponential-like)
    if (t_min < time[i] && c1 > 0 && c2 > 0 && c1 != c2) {
      auc <- auc + (c2 - c1) / log(c2 / c1) * dt
    } else {
      # linear-down or equal case
      auc <- auc + 0.5 * (c1 + c2) * dt
    }
  }
  auc
}

get_auc_linear <- function(time, conc) {
  # Check input validity
  if (length(time) != length(conc)) {
    stop("time and conc must have same length")
  }
  # if (any(diff(time) <= 0)) {
  #   stop("time must be strictly increasing")
  # }

  auc <- 0
  for (i in 2:length(time)) {
    dt <- time[i] - time[i - 1]
    c1 <- conc[i - 1]
    c2 <- conc[i]

    # linear case
    auc <- auc + 0.5 * (c1 + c2) * dt
  }
  auc
}

get_auc_log <- function(time, conc) {
  # Check input validity
  if (length(time) != length(conc)) {
    stop("time and conc must have same length")
  }
  if (any(diff(time) <= 0)) {
    stop("time must be strictly increasing")
  }

  auc <- 0
  for (i in 2:length(time)) {
    dt <- time[i] - time[i - 1]
    c1 <- conc[i - 1]
    c2 <- conc[i]

    # log case
    if (c1 > 0 && c2 > 0 && c1 != c2) {
      auc <- auc + (c2 - c1) / log(c2 / c1) * dt
    } else {
      # fallback to linear if non-positive concentrations
      # https://onlinehelp.certara.com/phoenix/8.3/topics/Partial_area_calculation.htm
      auc <- auc + 0.5 * (c1 + c2) * dt
    }
  }
  auc
}

check_time_varying <- function(x, group = "id", column) {
  baseline_vary <- x |>
    dplyr::group_by(.data[[group]]) |>
    dplyr::summarize(
      n_baselines = dplyr::n_distinct(.data[[column]]),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$n_baselines > 1)

  nrow(baseline_vary) > 0
}

#' Read pH Data
#' @description Reads pH data from a CSV or Excel file and performs basic validation.
#' @param file_path Path to the data file.
#' @param dose_time Time of dose administration (default is 1).
#' @param amt Dose amount (default is 100).
#' @return A data frame containing the pH data.
#' @author Omar I. Elashkar
#' @export
read_pH <- function(file_path, dose_time = 1, amt = 100) {
  checkmate::assertFileExists(file_path)

  file_ext <- tools::file_ext(file_path)
  checkmate::assertChoice(file_ext, choices = c("csv", "xls", "xlsx"))
  # checkmate::assertNumeric(baseline, lower = 0, upper = 14)
  checkmate::assertNumeric(dose_time, finite = TRUE, min.len = 1)
  checkmate::assertNumeric(amt, finite = TRUE, min.len = 1)

  if (file_ext %in% c("xls", "xlsx")) {
    dat <- readxl::read_excel(file_path)
    dat <- as.data.frame(dat)
  } else if (file_ext == "csv") {
    dat <- read.csv(file_path)
  } else {
    stop("Unsupported file type. Please provide a CSV or Excel file.")
  }

  # remove empty rows (trailing empty)
  dat <- dat |> janitor::remove_empty(which = c("rows"))

  if (is.null(dat$group)) {
    dat$group <- "default"
  }

  if (is.null(dat$flowrate)) {
    dat$flowrate <- NA
  }

  if (is.null(dat$buffering)) {
    dat$buffering <- NA
  }

  if (check_time_varying(dat, group = "id", column = "group")) {
    stop(
      "Group varies within subject(s). Please ensure unique group for each subject ID in the data."
    )
  }
  if (check_time_varying(dat, group = "id", column = "flowrate")) {
    stop(
      "Flowrate varies within subject(s). Please ensure unique flowrate for each subject ID in the data."
    )
  }
  if (check_time_varying(dat, group = "id", column = "buffering")) {
    stop(
      "Buffering varies within subject(s). Please ensure unique buffering for each subject ID in the data."
    )
  }

  adm_df <- data.frame(time = dose_time, adm = 1, amt = amt)
  stopifnot(nrow(adm_df) == length(dose_time))

  dat <- dat |>
    dplyr::group_by(.data$id, .data$group) |>
    dplyr::mutate(adm = 0)

  dat <- split(dat, dat$id) |>
    lapply(function(df) {
      df <- dplyr::bind_rows(df, adm_df) |>
        mutate(id = df$id[1], group = df$group[1]) |>
        dplyr::arrange(time) |>
        dplyr::mutate(amt = ifelse(.data$adm == 1, amt, NA))
    })
  dat <- do.call(rbind, dat)

  # if (is.null(dat$baseline)) {
  #   dat$baseline <- baseline
  # } else{
  #   # make sure not NA
  #   if(any(is.na(dat$baseline))){
  #     stop("Baseline cannot be NA. Please fill in missing baseline values or set baseline parameter.")
  #   }

  #   # make sure not time varying
  #   baseline_vary <- check_time_varying(dat, group = "id", column = "baseline")
  #   if (baseline_vary) {
  #     stop("Baseline varies within subject(s). Please ensure baseline is constant for each subject or set baseline parameter.")
  #   }

  # }

  # if(!is.null(baseline_time)){
  #   dat <- split(dat, dat$id) |> lapply(function(df) {
  #       predose = data.frame(
  #         time = baseline_time,
  #         pH = df$baseline[2],
  #         id = df$id[1],
  #         group = df$group[1],
  #         baseline = df$baseline[1],
  #         flowrate = df$flowrate[1],
  #         buffering = df$buffering[1]
  #       )
  #       df <- dplyr::bind_rows(predose, df)
  #       df
  #     })
  #   dat <- do.call(rbind, dat) |>
  #     dplyr::mutate(time = time + abs(baseline_time)) |>
  #     dplyr::arrange(id, time)
  # }

  dat$group_code <- factor_to_numeric(dat$group)
  dat <- dat |>
    # dplyr::select(-"baseline") |>
    dplyr::mutate(id = as.numeric(id)) |>
    janitor::remove_empty(which = c("rows", "cols")) # recheck remove empty

  check_data(dat)
  dat
}


split_pH_data <- function(x) {
  split(x, x$adm)
}

#' Check pH Data
#' @description Validates the structure and content of pH data.
#' @param x Data frame containing pH data.
#' @return TRUE if the data is valid; otherwise, an error is raised.
#' @author Omar I. Elashkar
#' @noRd
check_data <- function(x, sim = FALSE) {
  checkmate::assertDataFrame(x)
  checkmate::assertNames(
    names(x),
    must.include = c("id", "pH", "time", "group", "group_code", "amt", "adm")
  )
  if (!sim) {
    checkmate::assertNumeric(x$pH, lower = 0, upper = 14)
  }
  checkmate::assertNumeric(x$time, finite = TRUE)
  # check id and group are nonempty
  checkmate::assertIntegerish(x$id, any.missing = FALSE, min.len = 1)
  checkmate::assertNumeric(x$id, any.missing = FALSE, min.len = 1)

  checkmate::assertVector(x$group, any.missing = FALSE, min.len = 1)
  checkmate::assertNumeric(x$group_code, any.missing = FALSE, min.len = 1)

  checkmate::assertNumeric(x$time, finite = TRUE)
  checkmate::assertNumeric(x$amt, any.missing = TRUE)

  if (!isTRUE(all.equal(sort(unique(x$adm)), c(0, 1)))) {
    stop(
      "Administration indicator (adm) must be binary (0 or 1). Found values: ",
      paste(unique(x$adm), collapse = ", ")
    )
  }

  # ensure group, group_code, and baseline do not vary within each subject
  if ("group" %in% names(x)) {
    if (check_time_varying(x, group = "id", column = "group")) {
      stop(
        "Group varies within subject(s). Please ensure group is constant for each subject."
      )
    }
  }

  if ("group_code" %in% names(x)) {
    if (check_time_varying(x, group = "id", column = "group_code")) {
      stop(
        "Group code varies within subject(s). Please ensure group code is constant for each subject."
      )
    }
  }

  # if ("baseline" %in% names(x)) {
  #   baseline_vary <- x |>
  #     dplyr::group_by(.data$id) |>
  #     dplyr::summarize(n_baselines = dplyr::n_distinct(.data$baseline), .groups = "drop") |>
  #     dplyr::filter(.data$n_baselines > 1)

  #   if (nrow(baseline_vary) > 0) {
  #     stop("Baseline varies within subject(s): ", paste(baseline_vary$id, collapse = ", "))
  #   }
  # }

  # check if more than 10 groups
  if ("group" %in% names(x)) {
    n_groups <- length(unique(x$group))
    if (n_groups > 5) {
      stop(
        "GatorpH cannot handle more than 10 groups. Found ",
        n_groups,
        " groups."
      )
    }
  }

  TRUE
}

#' Check if pH Crosses Threshold on both upper and lower bound
#' @description Checks if the pH data crosses a specified threshold on both sides of the minimum pH value.
#' @param xtime Vector of time points.
#' @param xph Vector of pH values corresponding to the time points.
#' @param ph_threshold pH threshold (default is 5.4).
#' @return A logical vector of length 2 indicating whether the pH crosses the threshold before and after the minimum pH point.
#' @author Omar I. Elashkar
#' @noRd
check_crossing <- function(xtime, xph, ph_threshold = 5.4) {
  minpH <- min(xph, na.rm = TRUE)
  timeMin <- xtime[which.min(xph)]
  # no cross at all
  if (minpH >= ph_threshold) {
    return(FALSE)
  }

  # first point below threshold
  timeFirst <- xtime[which(xph < ph_threshold)[1]]

  # last point below threshold
  timeLast <- xtime[which(xph < ph_threshold)[length(which(
    xph < ph_threshold
  ))]]

  # cross before min pH
  if (any(xph[xtime <= timeFirst] >= ph_threshold, na.rm = TRUE)) {
    before <- TRUE
  } else {
    before <- FALSE
  }
  if (any(xph[xtime >= timeLast] >= ph_threshold, na.rm = TRUE)) {
    after <- TRUE
  } else {
    after <- FALSE
  }

  c(before, after)
}

#' Check if first and last pH points passed as crossing the threshold
#' Must have interpolated values for support
#' @param xtime Vector of time points.
#' @param xph Vector of pH values corresponding to the time points.
#' @param ph_threshold pH threshold (default is 5.4).
#' @return A logical vector of length 2 indicating whether the first and last points cross the threshold.
#' @author Omar I. Elashkar
#' @noRd
check_below_threshold <- function(
  xtime,
  xph,
  ph_threshold = 5.4,
  start_time = 0,
  end_time = 50
) {
  firstTimePoint <- find_closest(xtime, start_time, threshold = 0.2)
  lastTimePoint <- find_closest(xtime, end_time, threshold = 0.2)

  firstpHPoint <- xph[xtime == firstTimePoint][1]
  lastpHPoint <- xph[xtime == lastTimePoint][1]

  before <- (ph_threshold - firstpHPoint) >= 0
  after <- (ph_threshold - lastpHPoint) >= 0

  c(isTRUE(before), isTRUE(after))
}

find_closest <- function(vec, target, threshold = 0.2) {
  # Calculate absolute differences
  diffs <- abs(vec - target)

  # Find the index of the smallest difference
  best_index <- which.min(diffs)

  # Check if the smallest difference is within the threshold
  if (diffs[best_index] <= threshold) {
    return(vec[best_index])
  } else {
    return(NA) # Returns NA if no value is close enough
  }
}

interpolate_pH <- function(xtime, xph) {
  time_points <- seq(min(xtime), max(xtime), by = 0.1)
  time_points <- time_points[
    time_points >= min(xtime) & time_points <= max(xtime)
  ]

  newpH <- approx(xtime, xph, xout = time_points)$y
  new <- data.frame(time = time_points, pH = newpH)
  new
}

#' Integrate Area Under pH Threshold
#' Calculates the area under the pH curve below a specified pH threshold using the trapezoidal rule.
#' It assumes the flipped shape
#' @param x Vector of time points.
#' @param y Vector of pH values corresponding to the time points.
#' @param ph_threshold pH threshold (default is 5.4).
#' @param time_start Start time for integration (default is 0).
#' @param time_end End time for integration (default is 50).
#' @return Numeric value representing the area under the pH curve below the threshold.
#' @author Omar I. Elashkar
#' @export
integratepHArea <- function(
  x,
  y,
  ph_threshold = 5.4,
  time_start = 0,
  time_end = 50,
  method = "linear",
  interpolate = TRUE,
  add_support_points = FALSE,
  plot = FALSE
) {
  checkmate::assertNumeric(x)
  checkmate::assertNumeric(y)
  checkmate::assertNumber(ph_threshold, lower = 0, upper = 14)
  checkmate::assertNumber(time_start, finite = TRUE)
  checkmate::assertNumber(time_end, finite = TRUE)
  checkmate::assertChoice(
    method,
    choices = c("linear", "log", "linear_down_log_up")
  )

  tmpdata <- data.frame(time = x, pH = y)

  # Interpolation is needed as there is usually no support points on the threshold
  if (interpolate) {
    tmpdata <- interpolate_pH(tmpdata$time, tmpdata$pH)
    stopifnot(nrow(tmpdata) >= length(x))
  }

  tmpdata <- tmpdata |>
    dplyr::filter(time >= time_start & time <= time_end)

  message(
    "Crossing check: ",
    paste(
      check_crossing(tmpdata$time, tmpdata$pH, ph_threshold = ph_threshold),
      collapse = ", "
    )
  )

  # The method will add points exactly at the start and end times if there are any existing data points below the threshold plus no crossing with threshold at this time point
  if (add_support_points) {
    blwThreshold <- check_below_threshold(
      tmpdata$time,
      tmpdata$pH,
      ph_threshold = ph_threshold,
      start_time = time_start,
      end_time = time_end
    )
    checkCross <- check_crossing(
      tmpdata$time,
      tmpdata$pH,
      ph_threshold = ph_threshold
    )
    nsupport <- 0
    len_before <- nrow(tmpdata)
    for (i in 1:2) {
      if (blwThreshold[i] & !checkCross[i]) {
        nsupport <- nsupport + 1
        tmpdata <- rbind(
          tmpdata,
          data.frame(
            time = ifelse(i == 1, time_start, time_end),
            pH = ph_threshold + 0.04
          )
        )
        message(paste(
          "Added support point at time",
          ifelse(i == 1, time_start, time_end),
          "with pH =",
          ph_threshold + 0.04
        ))
      }
    }
    tmpdata <- tmpdata[order(tmpdata$time), ]

    stopifnot(nsupport >= 0 & nsupport <= 2)
    stopifnot(nrow(tmpdata) == len_before + nsupport)

    message(paste(
      "Added",
      nsupport,
      "support points at start and end times for integration."
    ))

    message(
      "Crossing check: ",
      paste(
        check_crossing(tmpdata$time, tmpdata$pH, ph_threshold = ph_threshold),
        collapse = ", "
      )
    )
  }

  if (plot) {
    print(
      ggplot(tmpdata, aes(x = time, y = pH)) +
        geom_line() +
        geom_hline(
          yintercept = ph_threshold,
          linetype = "dashed",
          color = "red"
        ) +
        labs(
          title = "pH over Time with Threshold",
          x = "Time",
          y = "pH",
          subtitle = paste(
            "interpolate = ",
            interpolate,
            ", add_support_points = ",
            add_support_points
          )
        )
    )
  }

  if (
    all(check_crossing(tmpdata$time, tmpdata$pH, ph_threshold = ph_threshold))
  ) {
    tmpdata <- tmpdata |>
      dplyr::filter(
        abs(.data$pH - ph_threshold) <= 0.02 | .data$pH < ph_threshold
      ) # keep only points below threshold or close to it

    x <- tmpdata$time
    y <- tmpdata$pH

    y <- ph_threshold - y # flip
    if (method == "linear") {
      res <- get_auc_linear(x, y)
    } else if (method == "log") {
      res <- get_auc_log(x, y)
    } else if (method == "linear_down_log_up") {
      res <- get_auc_linear_down_log_up(x, y)
    }

    if (length(res) == 0) {
      res <- NA
    }
  } else {
    res <- NA
  }
  res
}


#' Calculate Time Under pH Threshold
#' @description Calculates the total time each subject/group spends below a specified pH threshold.
#' @param x Data frame containing pH data.
#' @param ph_threshold pH threshold (default is 5.4).
#' @return A data frame with the total time under the pH threshold for each subject/group.
#' @author Omar I. Elashkar
#' @export
calc_time_under_pH <- function(x, ph_threshold = 5.4) {
  check_data(x)

  x <- split_pH_data(x)[["0"]]
  time_under_ph_func <- function(xtime, xph, ph_threshold) {
    newdata <- interpolate_pH(xtime, xph)
    # filter only points below pH threshold
    newdata <- newdata[newdata$pH < ph_threshold, ]

    if (min(newdata$time) == Inf | max(newdata$time) == -Inf) {
      return(list(start_time = NA, end_time = NA))
    }

    list(
      start_time = min(newdata$time, na.rm = TRUE),
      end_time = max(newdata$time, na.rm = TRUE)
    )
  }

  x |>
    dplyr::group_by(.data$id, .data$group) |>
    dplyr::summarise(
      time_under_pH = time_under_ph_func(.data$time, .data$pH, ph_threshold)
    ) |>
    dplyr::mutate(start_time = .data$time_under_pH$start_time) |>
    dplyr::mutate(end_time = .data$time_under_pH$end_time) |>
    dplyr::mutate(
      time_under_pH = .data$end_time - .data$start_time
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct()
}

calc_area_under_pH <- function(
  x,
  ph_threshold = 5.4,
  time_start = 0,
  time_end = 50,
  method = "linear",
  add_support_points = FALSE,
  plot = FALSE # print plot for diagnostics
) {
  check_data(x)

  x <- split_pH_data(x)[["0"]]

  checkMinMax <- check_min_max(x)

  min_time <- checkMinMax[1]
  max_time <- checkMinMax[2]

  if (min_time > time_start | max_time < time_end) {
    stop(paste(
      "Data does not cover the full integration interval. Consider adjusting time_start and time_end parameters to",
      round(min_time, 2),
      "and",
      round(max_time, 2),
      "respectively."
    ))
  }

  # plot_pH_time(x)
  x |>
    dplyr::select("id", "group", "time", "pH") |>
    dplyr::group_by(.data$id, .data$group) |>
    dplyr::summarise(
      area_under_pH_no_interpolation = integratepHArea(
        x = .data$time,
        y = .data$pH,
        ph_threshold = ph_threshold,
        time_start = time_start,
        time_end = time_end,
        method = method,
        interpolate = FALSE,
        add_support_points = add_support_points,
        plot = plot
      ),
      area_under_pH = integratepHArea(
        x = time,
        y = pH,
        ph_threshold = ph_threshold,
        time_start = time_start,
        time_end = time_end,
        method = method,
        interpolate = TRUE,
        add_support_points = add_support_points,
        plot = plot
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(across(
      .cols = c("area_under_pH_no_interpolation", "area_under_pH"),
      \(x) ifelse(x < 0, 0, x)
    )) |>
    dplyr::mutate(
      area_label = paste0(
        "Area < pH_{",
        ph_threshold,
        " \\{",
        round(time_start, 2),
        ", ",
        round(time_end, 2),
        "\\} }"
      )
    ) |>
    dplyr::mutate(dur_label = paste0("T < pH_{", ph_threshold, "}"))
}

calc_pH_min <- function(x) {
  check_data(x)
  x <- split_pH_data(x)[["0"]]
  x |>
    dplyr::group_by(.data$id, .data$group) |>
    dplyr::summarise(pH_min = min(.data$pH, na.rm = TRUE), .groups = "drop")
}

calc_t_min <- function(x) {
  check_data(x)
  x <- split_pH_data(x)[["0"]]

  x |>
    dplyr::group_by(.data$id, .data$group) |>
    dplyr::summarise(t_min = .data$time[which.min(.data$pH)], .groups = "drop")
}

rxensure <- function(mod) {
  if (inherits(mod, "function")) {
    mod <- mod()
  }
  mod
}

#' Simulate Oral pH-Time Curve
#' @description Simulates oral pH-time curves based on a given PK-PD model.
#' @param model rxode2 model object representing the PK-PD model.
#' @param time Vector of time points for simulation (default is seq(0, 50, by = 0.1)).
#' @param dose Dose amount (default is 100).
#' @param dose_time Time of dose administration (default is 5).
#' @param nsub Number of subjects to simulate (default is 1).
#' @param baseline Baseline pH value. Default is 7. Must be either 1 or same number as number of subjects.
#' @param step Time step for simulation (default is 0.1).
#' @param group Group identifier for the subjects (default is "A").
#' @param ignoreBSV Logical indicating whether to ignore between-subject variability (default is TRUE).
#' @param ignoreRUV Logical indicating whether to ignore residual unexplained variability (default is TRUE).
#' @param include_gamma logical indicating whether to include the gamma parameter in the simulation (default is TRUE).
#' @return A data frame containing the simulated pH-time data.
#' @author Omar I. Elashkar
#' @noRd
simulate_steph_curve <- function(
  model,
  time = seq(0, 50, by = 0.1),
  dose = 100,
  dose_time = 5,
  nsub = 1,
  baseline = 7,
  step = 0.1,
  group = "A",
  ignoreBSV = TRUE,
  ignoreRUV = TRUE,
  include_gamma = TRUE
) {
  checkmate::assertNumber(dose, finite = TRUE)
  checkmate::assertNumber(nsub, lower = 1)
  model <- rxensure(model)
  checkmate::assertNumeric(dose_time, lower = 0, upper = Inf)
  checkmate::assertNumeric(baseline, lower = 0, upper = 14, null.ok = TRUE)

  if (!include_gamma) {
    model <- remove_gamma(model)
  }

  if (!is.null(baseline)) {
    inidf <- model$iniDf
    tkd <- inidf[inidf$name == "t.kd", "est"] |> exp()
    kd.sd <- inidf[inidf$name == "eta.kd", "est"] |> sqrt()
    # tks <- inidf[inidf$name == "t.ks", "est"] |> exp()
    tks <- mean(baseline) * tkd
    ikd_vals <- tkd * exp(rnorm(nsub, mean = 0, sd = kd.sd))
    checkmate::assertNumber(tkd, finite = TRUE)
    checkmate::assertNumber(kd.sd, finite = TRUE)
    iks_vals <- baseline * ikd_vals

    eta_ks_vals <- log(iks_vals / tks)
    eta_kd_vals <- log(ikd_vals / tkd)

    inidf <- inidf[!(inidf$name %in% c("eta.ks", "eta.kd")), ]
    inidf[inidf$name == "t.ks", "est"] <- log(tks)
    model$iniDf <- inidf
  }

  ev <- rxode2::et(amt = dose, cmt = "depot", time = abs(dose_time))
  ev <- ev |>
    rxode2::et(time = unique(c(0, time))) |>
    rxode2::et(id = seq(1, nsub))
  if (ignoreBSV) {
    # resp ~ add(sigma_add)
    model <- rxode2::zeroRe(model, which = "omega")
  }
  if (ignoreRUV) {
    model <- rxode2::zeroRe(model, which = "sigma")
  }

  basecovariates <- c("ks", "kd", "kde", "edk50")
  if (include_gamma) {
    basecovariates <- c(basecovariates, "gamma")
  }
  model <- parse_covariate(
    1,
    model,
    parameters = basecovariates,
    fixed_effects = paste0("t.", basecovariates)
  )

  group <- group
  group_code <- factor_to_numeric(group)

  if (!is.null(baseline)) {
    covdf <- data.frame(
      id = seq(1, nsub),
      eta.ks = eta_ks_vals,
      eta.kd = eta_kd_vals,
      group_code = group_code
    )
  } else {
    covdf <- data.frame(
      id = seq(1, nsub),
      group_code = group_code
    )
  }
  sim <- rxode2::rxSolve(model, events = ev, iCov = covdf, addDosing = TRUE)

  if (nsub == 1) {
    message("only one subject simulated, setting id to 1")
    sim <- sim |> dplyr::mutate(id = 1)
  }
  idx <- 1
  repeat {
    problematic_ids <- getnoVarIds(sim)

    if (length(problematic_ids) == 0 | idx > 10) {
      if (idx == 10) {
        browser()
      }
      break
    }
    message(
      "Regenerating subjects with no variability in pH response: ",
      paste(problematic_ids, collapse = ", ")
    )
    sim2 <- rxode2::rxSolve(
      model,
      events = ev |> dplyr::filter(.data$id %in% problematic_ids),
      iCov = covdf |> dplyr::filter(.data$id %in% problematic_ids),
      nSub = length(problematic_ids)
    )
    # replace only the problematic ids with regenerated subjects
    sim2$id <- rep(
      problematic_ids,
      each = nrow(sim2) / length(problematic_ids)
    )

    sim <- sim |>
      dplyr::filter(!(.data$id %in% problematic_ids)) |>
      dplyr::bind_rows(sim2)
    idx <- idx + 1
  }

  if (length(getnoVarIds(sim)) != 0) {
    stop(
      "Some subjects still have no variability in pH response after regeneration attempts."
    )
  }

  sim <- as.data.frame(sim)

  sim$id <- as.numeric(sim$id)
  stopifnot(sum(is.na(sim$id)) == 0)

  sim$group <- group
  sim$group_code <- group_code
  sim$pH <- sim$resp

  sim$group_code <- factor_to_numeric(sim$group)

  sim
}

sim_to_pH_data <- function(x) {
  x <- x |>
    dplyr::rename("adm" = "evid") |>
    dplyr::mutate(group_code = factor_to_numeric(.data$group)) |>
    dplyr::select("id", "group", "group_code", "time", "pH", "adm", "amt")
  check_data(x)
  x
}


plot_pkpd_curve <- function(res) {
  pkplot <- ggplot2::ggplot(
    res,
    ggplot2::aes(x = time, y = central, group = id)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Time (min)", y = "Drug Concentration (mg/mL)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  pdplot <- plot_pH_time(res)

  patchwork::wrap_plots(list(pkplot, pdplot)) + patchwork::plot_layout(ncol = 1)
}

#' Plot pH vs Time with Threshold
#' @description Generates a plot of pH values over time for different subjects/groups, highlighting a specified pH threshold.
#' @param res Data frame containing pH data.
#' @param ph_threshold pH threshold to highlight on the plot (default is 5.4).
#' @param show_id Logical indicating whether to show subject IDs in the legend (default is TRUE).
#' @param stratify_by Variable to stratify the plot by (default is "None"). Options are "None", "Subject", or "Group".
#' @param showAvg Logical indicating whether to show the average pH curve across subjects/groups (default is FALSE).
#' @param showDosing Logical indicating whether to show dosing times on the plot (default is FALSE).
#' @return A ggplot2 object representing the pH vs time plot.
#' @author Omar I. Elashkar
#' @export
plot_pH_time <- function(
  res,
  ph_threshold = 5.4,
  show_id = TRUE,
  stratify_by = "None",
  showAvg = FALSE,
  showDosing = FALSE
) {
  check_data(res)
  checkmate::assertChoice(
    stratify_by,
    choices = c("None", "Subject", "Group")
  )

  obs_data <- split_pH_data(res)[["0"]]
  obs_data$id <- as.factor(obs_data$id)
  plt <- ggplot2::ggplot(
    obs_data,
    ggplot2::aes(x = time, y = pH, color = id)
  ) +
    ggplot2::geom_line(aes(group = id)) +
    ggplot2::labs(x = "Time", y = "pH") +
    ggplot2::geom_hline(
      aes(yintercept = ph_threshold, color = "threshold"),
      linetype = "dashed"
    ) +
    ggplot2::geom_point()
  if (!show_id) {
    plt <- plt + ggplot2::theme(legend.position = "none")
  }
  if (showAvg & stratify_by != "Subject") {
    grouped <- stratify_by == "Group"

    restmp <- res |> dplyr::mutate(id = as.numeric(as.character(id)))
    avgdata <- curve_averaging(restmp, interpolate = TRUE, grouped = grouped)

    if (stratify_by != "Group") {
      avgdata$group <- "Average"
    }
    plt <- plt +
      ggplot2::geom_line(
        data = avgdata,
        aes(x = time, y = mean, color = group),
        linewidth = 1.5
      )
  }
  if (showDosing) {
    dosing_data <- split_pH_data(res)[["1"]] |>
      dplyr::select(time, id, group) |>
      dplyr::distinct()

    minpH <- min(obs_data$pH, na.rm = TRUE)
    dosing_data <- dosing_data |> dplyr::mutate(minpH = minpH)

    plt <- plt +
      ggplot2::geom_segment(
        data = dosing_data,
        aes(
          x = time,
          xend = time,
          y = minpH,
          yend = -Inf,
          color = "Administration"
        ),
        linetype = "solid",
        lineend = "round",
        linejoin = "round",
        arrow = ggplot2::arrow(
          length = ggplot2::unit(0.2, "cm"),
          ends = "last"
        ),
        linewidth = 1
      )
  }

  if (stratify_by == "Group") {
    plt <- plt + ggplot2::facet_wrap(~group)
  }
  if (stratify_by == "Subject") {
    plt <- plt + ggplot2::facet_wrap(~id)
  }

  unique_ids <- unique(obs_data$id)
  id_colors <- scales::hue_pal()(length(unique_ids))
  id_color_map <- setNames(id_colors, unique_ids)
  clr_scales <- c(
    "Administration" = "red",
    "Average" = "blue",
    "threshold" = "black",
    id_color_map
  )
  plt +
    ggplot2::scale_color_manual(values = clr_scales, na.value = "grey") +
    ggplot2::theme_minimal()
}


#' Read Digitized pH Data from WebPlotDigitizer
#' @noRd
digitizeread <- function(x) {
  df_raw <- read_csv(x, col_names = FALSE)

  headers <- df_raw[1:2, ] |>
    t() |>
    as.data.frame() |>
    fill(V1, .direction = "down") |>
    dplyr::mutate(name = paste(V1, V2, sep = "_")) |>
    pull(name)

  df <- df_raw[-c(1, 2), ]
  colnames(df) <- headers

  df <- df |> dplyr::mutate(across(everything(), as.numeric))

  res <- df |>
    pivot_longer(
      cols = everything(),
      names_to = c("group", ".value"),
      names_pattern = "(.*)_(X|Y)"
    ) |>
    group_by(.data$group) |>
    dplyr::mutate(id = cur_group_id()) |>
    dplyr::ungroup() |>
    dplyr::rename(time = X, pH = Y)
  # ensure no value more than 16 and not less than -1
  stopifnot(all(res$pH <= 15, na.rm = TRUE))
  stopifnot(all(res$pH >= -1, na.rm = TRUE))

  # round up to 0 if less than 0
  res <- res |>
    dplyr::mutate(pH = ifelse(pH < 0, 0, pH)) |>
    dplyr::mutate(pH = ifelse(pH > 14, 14, pH))
  res
}


#' Fit pH Curve using NLME
#' @description Fits a pharmacodynamic model to pH data using nonlinear mixed-effects modeling. If `stratify` is TRUE, fits separate models for each group.
#' @param data Data frame containing pH data. Must include columns: pH, time, id, group.
#' @param model rxode2/nlmixr2 model to fit.
#' @param stratify Logical indicating whether to fit separate models for each group (default is FALSE).
#' @param estmethod Estimation method to use (default is "focei").
#' @param cov_params Character vector of covariate parameters to include in the model (default is c("kd", "kde", "edk50", "gamma")).
#' @param cov_fixedeffects Character vector of fixed effect names corresponding to the covariate parameters (default is c("t.kd", "t.kde", "t.edk50", "t.gamma")).
#' @param include_gamma logical indicating whether to include the gamma parameter in the model (default is TRUE).
#' @return nlmixr2 fit object.
#' @author Omar I. Elashkar
#' @export
fit_pH_curve <- function(
  data,
  model,
  # amt,
  stratify = FALSE,
  estmethod = "focei",
  # dose_time = 5,
  cov_params = c("kd", "kde", "edk50", "gamma"),
  cov_fixedeffects = c("t.kd", "t.kde", "t.edk50", "t.gamma"),
  include_gamma = TRUE,
   nprint = 50
) {
  if ("ks" %in% cov_params && "kd" %in% cov_params) {
    stop("cov_params cannot contain both 'ks' and 'kd'")
  }

  check_data(data, sim = TRUE)
  checkmate::assertLogical(stratify, len = 1)
  checkmate::assertChoice(
    estmethod,
    choices = c("focei", "saem", "bobyqa", "uobyqa", "fo", "foce")
  )
  # checkmate::assertNumber(amt, finite = TRUE)
  # checkmate::assertNumber(dose_time, finite = TRUE, upper = Inf, lower = 0)

  model <- rxensure(model)
  if (!include_gamma) {
    model <- remove_gamma(model)
  }

  # preserve group and group_code information
  group_info <- data |>
    dplyr::select(id, group, group_code) |>
    dplyr::distinct()

  data <- data |>
    dplyr::mutate(cmt = NA_character_) |>
    dplyr::rename(evid = "adm")

  # data <- split(data, data$id) |>
  #   lapply(function(df) {
  #     df |>
  #       dplyr::filter(time != dose_time) |>
  #       add_row(
  #         pH = NA,
  #         id = df$id[1],
  #         group = df$group[1],
  #         group_code = df$group_code[1],
  #         cmt = "depot",
  #         evid = 1,
  #         amt = amt,
  #         time = dose_time,
  #         .before = 1
  #       )
  #   })
  # data <- do.call(rbind, data)

  # assert correct dose time for predose for each subject exist
  # tmpdat <- data |> dplyr::group_by(.data$id) |>
  #   dplyr::filter(any(.data$time == dose_time)) |>
  #   dplyr::summarise(n_predose = sum(.data$time == dose_time, na.rm = TRUE), .groups = "drop") |>
  #   dplyr::ungroup() |>
  #   dplyr::filter(n_predose == 0)

  # if(nrow(tmpdat) > 0){
  #   stop("Missing predose time point for subject(s): ", paste(tmpdat$id, collapse = ", "))
  # }

  data <- data |>
    dplyr::rename(DV = "pH") |>
    select("id", "group", "group_code", "time", "DV", "evid", "amt", "cmt") |>
    dplyr::arrange(id, time)

  uniqueids <- unique(data$id)

  if (length(uniqueids) > 1 && length(unique(data$group_code)) > 1) {
    model <- parse_covariate(
      unique(data$group_code),
      model,
      parameters = cov_params,
      fixed_effects = cov_fixedeffects
    )
  } else {
    warning(
      "Only one group or one subject in the data. Not including group covariate in the model."
    )
  }

  if (estmethod %in% c("bobyqa", "uobyqa", "fo", "foce", "focei")) {
    ctrl <- eval(parse(text=paste0("nlmixr2est::", estmethod, "Control", sep = "")))
    finalFit <- nlmixr2est::nlmixr2(
      zeroRe(model, which = "omega"),
      data,
      est = estmethod, 
      control = ctrl(print = nprint)
    )
  } else {
    if (length(uniqueids) == 1) {
      stop(
        "Only one subject in the data. Cannot fit NLME model. Switch to naive method"
      )
    }

    # finalFit <- nlmixr2est::nlmixr2(
    #   rxode2::zeroRe(model),
    #   data,
    #   est = "uobyqa"
    # )

    # model <- finalFit |>
    #   ini(eta.ks = unfix) |>
    #   ini(eta.kd = unfix) |>
    #   ini(eta.kde=unfix) |>
    #   ini(eta.edk50 = unfix) |>
    #   ini(sigma_add = unfix) |>

    #   ini(eta.ks = 0.1) |>
    #   ini(eta.kd = 0.1) |>
    #   ini(eta.kde=0.1) |>
    #   ini(eta.edk50 = 0.1) |>
    #   ini(sigma_add = 0.1)

    if (estmethod == "focei") {
      ctrl <- nlmixr2est::foceiControl(
        maxOuterIterations = 100,
        maxInnerIterations = 100,
        print = nprint
      )
    } else if (estmethod == "saem") {
      ctrl <- nlmixr2est::saemControl(print = nprint)
    }

    finalFit <- nlmixr2est::nlmixr2(
      model,
      data,
      est = estmethod,
      control = ctrl
    )
  }

  finalFit |>
    dplyr::mutate(ID = as.numeric(as.character(.data$ID))) |>
    # read groups
    dplyr::left_join(group_info, by = c("ID" = "id")) |>
    dplyr::rename(id = "ID")
}

fit_param_table <- function(fit) {
  fit$parFixed
}

#' Format pH Metrics Table
#' @description Formats a data frame containing pH metrics for presentation, including formatting numbers, adding labels, and organizing columns into spanners.
#' @param x Data frame containing pH metrics to format.
#' @return A gt table object representing the formatted pH metrics table.
#' @author Omar I. Elashkar
#' @export
format_metrics_tab <- function(x) {
  x <- x |>
    dplyr::select(
      -c("area_under_pH_no_interpolation", "start_time", "end_time")
    ) |>
    dplyr::mutate(area_label = paste0("$$", area_label, "$$")) |>
    dplyr::mutate(dur_label = paste0("$$", dur_label, "$$")) |>

    gt::gt() |>
    gt::fmt_number(columns = everything(), decimals = 2) |>
    gt::fmt_markdown(columns = c("area_label", "dur_label")) |>
    gt::cols_label(
      id = "ID",
      group = "Group",
      pH_min = "pH_{min}",
      t_min = "T_{min}",
      time_under_pH = "T < pH_{c}",
      area_under_pH = "Area < pH_{c}",
      area_label = "Area Label",
      dur_label = "Duration Label",
      dplyr::matches("ks") ~ "KS",
      dplyr::matches("kd") ~ "KD",
      dplyr::matches("kde") ~ "KDE",
      dplyr::matches("edk50") ~ "EDK50",
      dplyr::matches("gamma") ~ "Gamma",
      .fn = \(x) gt::md(paste0("$$", x, "$$"))
    )

  x <- x |>
    gt::tab_spanner(
      label = "Model Parameters",
      columns = c(
        matches("ks"),
        matches("kd"),
        matches("kde"),
        matches("edk50"),
        matches("gamma")
      )
    ) |>
    gt::tab_spanner(
      label = "pH Metrics",
      columns = c("area_under_pH", "time_under_pH", "pH_min", "t_min")
    )

  x
}


#' Plot Observed vs Predicted pH
#' @description Generates a plot of observed pH values versus predicted pH values from a fitted model, including a reference line for perfect predictions.
#' @param fit nlmixr2 fit object containing the observed and predicted values.
#' @return A ggplot2 object representing the observed vs predicted pH plot.
#' @author Omar I. Elashkar
#' @export
fit_obs_vs_ipred_plot <- function(fit) {
  as.data.frame(fit) |>
    ggplot2::ggplot(aes(x = DV, y = IPRED)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
    ggplot2::labs(x = "Observed pH", y = "Predicted pH") +
    ggplot2::theme_minimal()
}

#' Plot Individual pH Profiles
#' @description Generates a plot of individual pH profiles over time for each subject/group, including observed points and predicted lines from a fitted model.
#' @param fit nlmixr2 fit object containing the observed and predicted values.
#' @return A ggplot2 object representing the individual pH profiles plot.
#' @author Omar I. Elashkar
#' @export
fit_individual_plot <- function(fit) {
  as.data.frame(fit) |>
    dplyr::mutate(id = as.factor(id)) |>
    ggplot2::ggplot(aes(x = TIME, y = DV, color = id)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(aes(y = IPRED)) +
    ggplot2::facet_wrap(~id) +
    ggplot2::labs(x = "Time", y = "pH") +
    ggplot2::theme_minimal()
}


#' Run direct estimation for pH data
#' @description Runs direct estimation calculations for pH data, including minimum pH, time to minimum pH, time under pH threshold, and AUC under pH threshold.
#' @param x Data frame containing pH data.
#' @param ph_threshold pH threshold for calculations.
#' @param time_start Start time for area and time under pH calculation.
#' @param time_end End time for area and time under pH calculation.
#' @param method Method for AUC calculation ("linear", "log", "linear_down_log_up").
#' @return A data frame with subject and group pH parameters
#' @export
run_direct_estimation <- function(
  x,
  ph_threshold = 5.4,
  time_start = 0,
  time_end = 50,
  method = "linear",
  add_support_points = FALSE,
  use_baseline = TRUE
) {
  check_data(x)

  phmin <- calc_pH_min(x)
  tmin <- calc_t_min(x)
  time_under_ph <- calc_time_under_pH(
    x,
    ph_threshold = ph_threshold
  )
  auc <- calc_area_under_pH(
    x,
    ph_threshold = ph_threshold,
    time_start = time_start,
    time_end = time_end,
    method = method,
    add_support_points = add_support_points
  )
  res <- phmin |>
    dplyr::left_join(tmin, by = c("id", "group")) |>
    dplyr::left_join(time_under_ph, by = c("id", "group")) |>
    dplyr::left_join(auc, by = c("id", "group"))

  class(res) <- c("direct_estimation", "QuantPH", class(res))
  res
  res
}

#' Get pH Metrics from NLME Fit
#' @param x nlmixr2 fit object.
#' @param time_start Start time for simulation.
#' @param time_end End time for simulation.
#' @param step Time step for simulation.
#' @param dose Dose amount to administer at time 0.
#' @param plot Logical indicating whether to plot the simulated pH profiles (default is FALSE).
#' @param ph_threshold pH threshold for calculations (default is 5.4).
#' @param stratify_by Stratification option for plotting ("None", "Subject", "Group").
#' @return A data frame containing pH metrics calculated from the simulated pH profiles.
#' @author Omar I. Elashkar
#' @export
pHMetrics_from_fit <- function(
  x,
  ph_threshold = 5.4,
  time_start = 0,
  time_end = 50,
  step = 0.1,
  stratify_by = "None",
  include_gamma = TRUE,
  add_support_points = FALSE
) {
  checkmate::assertClass(x, "nlmixr2FitCore")
  model <- x$finalUi
  uniqueids <- unique(x$data$id)
  nsub <- length(uniqueids)
  time <- seq(time_start, time_end, by = step)
  estMethod <- x$est
  onlymean <- estMethod == "uobyqa" | estMethod == "bobyqa"

  # fit_individual_plot(x)

  dose_time <- x$origData |>
    dplyr::filter(.data$evid == 1) |>
    dplyr::pull("time") |>
    unique()

  dose <- x$origData |>
    dplyr::filter(.data$evid == 1) |>
    dplyr::pull("amt") |>
    unique()

  fixedparamdf <- x$parFixedDf |>
    tibble::rownames_to_column("name") |>
    dplyr::select(name, Estimate) |>
    tidyr::pivot_wider(names_from = "name", values_from = "Estimate") |>
    dplyr::select(dplyr::starts_with("t."))

  has_covariates <- any(grepl("cov_", rownames(x$parFixedDf)))
  if (has_covariates) {
    if (length(unique(x$origData$group_code)) > 1) {
      icovDf <- as.data.frame(x) |>
        dplyr::select(
          "id",
          dplyr::starts_with("group"),
          dplyr::starts_with("eta."),
          dplyr::starts_with("cov_")
        ) |>
        dplyr::mutate(group_code = factor_to_numeric(.data$group)) |>
        tidyr::pivot_longer(
          cols = dplyr::starts_with("cov_"),
          names_to = "covariate",
          values_to = "value"
        ) |>
        dplyr::mutate(
          covariate = paste0(.data$covariate, "_", .data$group_code)
        ) |>
        dplyr::distinct() |>
        tidyr::pivot_wider(
          names_from = "covariate",
          values_from = "value",
          values_fill = 0
        ) |>
        dplyr::distinct() |>
        dplyr::mutate(id = as.numeric(as.character(.data$id)))
    } else {
      icovDf <- as.data.frame(x) |>
        dplyr::select(
          "id",
          dplyr::starts_with("group"),
          dplyr::starts_with("eta."),
          dplyr::starts_with("cov_")
        ) |>
        dplyr::mutate(group_code = factor_to_numeric(.data$group)) |>
        # add covariate manually
        mutate(
          cov_ks = 0,
          cov_edk50 = 0,
          cov_kd = 0,
          cov_kde = 0,
          cov_gamma = 0
        ) |>
        tidyr::pivot_longer(
          cols = dplyr::starts_with("cov_"),
          names_to = "covariate",
          values_to = "value"
        ) |>
        dplyr::mutate(
          covariate = paste0(.data$covariate, "_", .data$group_code)
        ) |>
        dplyr::distinct() |>
        tidyr::pivot_wider(
          names_from = "covariate",
          values_from = "value",
          values_fill = 0
        ) |>
        dplyr::distinct() |>
        dplyr::mutate(id = as.numeric(as.character(.data$id)))
    }
  } else {
    # if no covariates, just get ID and group information
    icovDf <- as.data.frame(x) |>
      dplyr::select(
        "id",
        dplyr::starts_with("group"),
        dplyr::starts_with("eta.")
      ) |>
      dplyr::mutate(group_code = factor_to_numeric(.data$group)) |>
      dplyr::distinct() |>
      dplyr::mutate(id = as.numeric(as.character(.data$id)))
  }

  new_mod <- model |> rxode2::zeroRe(which = "sigma")

  ev <- rxode2::et(amt = dose, cmt = "depot", time = dose_time) |> # use dose_time here
    rxode2::et(time = time)

  newini <- new_mod$iniDf |>
    dplyr::filter(!grepl("eta\\.", .data$name), !grepl("cov_", .data$name))

  # assert fixed effects parameter are final
  newini <- newini |>
    dplyr::mutate(
      est = dplyr::case_when(
        name == "t.edk50" ~ fixedparamdf$t.edk50,
        name == "t.kde" ~ fixedparamdf$t.kde,
        name == "t.kd" ~ fixedparamdf$t.kd,
        name == "t.ks" ~ fixedparamdf$t.ks,
        name == "t.gamma" ~ fixedparamdf$t.gamma,
        TRUE ~ est
      )
    )

  rxode2::ini(new_mod) <- newini

  if (onlymean) {
    icovDf <- icovDf |>
      # dplyr::mutate(across(starts_with("eta."), ~ 0)) |>
      dplyr::mutate(
        eta.edk50 = 0,
        eta.kde = 0,
        eta.kd = 0,
        eta.ks = 0,
        eta.gamma = 0
      ) |> # ensure regressors
      dplyr::distinct()
  }

  if (!include_gamma) {
    icovDf <- icovDf |>
      dplyr::mutate(eta.gamma = 0)
  }

  ids <- as.numeric(as.character(unique(icovDf$id)))
  ev <- rxode2::et(amt = dose, cmt = "depot", time = dose_time) |> # use dose_time here
    rxode2::et(time = time) |>
    rxode2::et(id = ids)

  simRes <- rxSolve(
    new_mod,
    iCov = icovDf, # TODO add flowrate, buffering, substance
    events = ev,
    addDosing = TRUE
  )

  # fix any subjects with no variability in pH response
  idx <- 1
  repeat {
    problematic_ids <- getnoVarIds(simRes)

    if (length(problematic_ids) == 0 | idx > 10) {
      if (idx == 10) {
        browser()
      }
      break
    }
    message(
      "Regenerating subjects with no variability in pH response: ",
      paste(problematic_ids, collapse = ", ")
    )
    sim2 <- rxSolve(
      new_mod,
      iCov = icovDf |> dplyr::filter(.data$id %in% problematic_ids),
      events = ev |> dplyr::filter(.data$id %in% problematic_ids)
    )
    # replace only the problematic ids with regenerated subjects
    simRes <- simRes |>
      dplyr::filter(!(.data$id %in% problematic_ids)) |>
      dplyr::bind_rows(sim2)
    idx <- idx + 1
  }

  stopifnot(length(getnoVarIds(simRes)) == 0)

  if (is.null(simRes$id)) {
    simRes <- dplyr::mutate(simRes, id = x$origData$id[1])
  }
  simRes <- as.data.frame(simRes) |>
    dplyr::group_by(.data$id) |>
    dplyr::filter(sum(.data$resp, na.rm = TRUE) > 0) |>
    dplyr::ungroup() |>
    dplyr::rename(pH = "resp")

  # simRes$id <- as.factor(simRes$sim.id)

  # Get original groups from x$origData
  orig_groups <- icovDf |>
    dplyr::select("id", "group", "group_code") |>
    dplyr::distinct() |>
    dplyr::mutate(id = as.numeric(as.character(.data$id)))

  simRes <- simRes |>
    dplyr::select(
      -dplyr::starts_with("group"),
      -dplyr::starts_with("group_code")
    ) |>
    dplyr::left_join(orig_groups, by = c("id" = "id")) |>
    dplyr::rename(adm = "evid")

  # no matter if pooled or not, all original groups must have calculations
  stopifnot(length(unique(simRes$group)) == length(unique(x$origData$group)))
  stopifnot(unique(simRes$group) == unique(x$origData$group))

  plt <- plot_pH_time(
    simRes,
    show_id = FALSE,
    stratify_by = stratify_by,
    ph_threshold = ph_threshold,
    showDosing = TRUE
  ) +
    labs(
      subtitle = paste0(
        ifelse(onlymean, "Mean Profile", "Individual Profiles"),
        " from method ",
        estMethod
      )
    )

  if (!onlymean) {
    originalData <- nlme::getData(x) |> dplyr::filter(.data$evid == 0)
    plt <- plt +
      ggplot2::geom_point(
        data = originalData,
        aes(x = .data$time, y = .data$DV),
        color = "red",
        size = 1,
        shape = 4
      )
  }

  derivedDf <- run_direct_estimation(
    simRes,
    ph_threshold = ph_threshold,
    time_start = time_start,
    time_end = time_end,
    add_support_points = add_support_points
  ) |>
    as.data.frame()

  if (onlymean) {
    derivedDf <- derivedDf |>
      # replace actual id with placeholder
      dplyr::mutate(
        id = ifelse(onlymean, ".", as.numeric(as.character(.data$id)))
      )
  }

  stopifnot(nrow(derivedDf) == nrow(icovDf))

  if (onlymean) {
    paramsdf <- as.data.frame(x) |>
      dplyr::mutate(id = ".")
  } else {
    paramsdf <- as.data.frame(x) |>
      dplyr::mutate(id = as.numeric(as.character(.data$id)))
  }

  param_cols <- c("edk50", "kde", "kd", "ks")
  if ("gamma" %in% names(paramsdf)) {
    param_cols <- c(param_cols, "gamma")
  }

  derivedDf <- dplyr::left_join(
    derivedDf,
    paramsdf |>
      dplyr::select(dplyr::all_of(c("id", param_cols, "group"))) |>
      dplyr::group_by(.data$id, .data$group) |>
      dplyr::summarize(
        across(dplyr::all_of(param_cols), mean),
        .groups = "keep"
      ) |>
      dplyr::ungroup() |>
      dplyr::distinct(),
    by = c("id" = "id", "group" = "group")
  )

  if (onlymean) {
    derived_summarize_cols <- c(
      "edk50",
      "kde",
      "kd",
      "ks",
      "pH_min",
      "t_min",
      "start_time",
      "end_time",
      "time_under_pH",
      "area_under_pH",
      "area_under_pH_no_interpolation"
    )

    if ("gamma" %in% names(derivedDf)) {
      derived_summarize_cols <- c(
        "edk50",
        "kde",
        "kd",
        "ks",
        "gamma",
        "pH_min",
        "t_min",
        "start_time",
        "end_time",
        "time_under_pH",
        "area_under_pH",
        "area_under_pH_no_interpolation"
      )
    }

    derivedDf <- derivedDf |>
      dplyr::group_by(.data$group, .data$area_label, .data$dur_label) |>
      dplyr::summarize(
        across(dplyr::all_of(derived_summarize_cols), mean),
        .groups = "keep"
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(id = ".")
    stopifnot(nrow(derivedDf) == length(unique(x$origData$group_code)))
  } else {
    stopifnot(nrow(derivedDf) == length(unique(x$origData$id)))
  }
  stopifnot(sort(unique(derivedDf$group)) == sort(unique(x$origData$group)))

  list(
    derivedDf = derivedDf,
    simRes = simRes,
    plot = plt
  )
}

curve_averaging <- function(
  x,
  na.rm = TRUE,
  interpolate = TRUE,
  xout = seq(0, 50, by = 0.1),
  grouped = FALSE
) {
  check_data(x)

  df <- x |>
    dplyr::select(
      time,
      pH,
      id,
      dplyr::all_of(if (grouped) "group" else character())
    )

  if (!grouped) {
    df <- df |>
      dplyr::mutate(group = "All")
  }

  if (interpolate) {
    df <- df |>
      group_by(.data$id, .data$group) |>
      dplyr::reframe(
        x_new = xout,
        y_new = approx(time, pH, xout = xout)$y
      ) |>
      dplyr::filter(!is.na(.data$y_new)) |> # Remove rows with NA values
      dplyr::ungroup() |>
      dplyr::rename(pH = "y_new", time = "x_new")
  }

  avg_df <- df |>
    dplyr::group_by(.data$time, .data$group) |>
    dplyr::summarise(
      mean = mean(.data$pH, na.rm = na.rm),
      sd = sd(.data$pH, na.rm = na.rm),
      n = sum(!is.na(.data$pH)),
      .groups = "drop"
    )

  avg_df
}

plot_curve_averaging <- function(x) {
  group <- length(unique(x$group)) > 1
  plt <- ggplot2::ggplot(x, ggplot2::aes(x = .data$time, y = .data$mean)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = mean - sd,
        ymax = mean + sd
      ),
      alpha = 0.2,
      fill = "lightblue"
    )

  if (group) {
    plt <- plt + ggplot2::facet_wrap(~group)
  }

  plt
}


avg_to_pHdata <- function(x) {
  checkmate::assertDataFrame(x)
  checkmate::assertNames(
    names(x),
    must.include = c("time", "mean")
  )

  data.frame(
    time = x$time,
    pH = x$mean,
    id = factor(1),
    group = factor("A")
  )
}


get_nsub <- function(x) {
  length(unique(x$id))
}

#' convert factor, char or numeric to numeric
#' @description Converts a factor, character, or numeric vector to numeric.
#' Returns the same numeric vector regardless of the input type (factor or character).
#' @param x A factor, character, or numeric vector.
#' @return A numeric vector.
#' @details
#' If input is numeric, returns as-is.
#' If input is factor, converts to character then to numeric (preserving original numeric values).
#' If input is character, converts to factor then to numeric (assigning numeric codes to unique values).
#' @author Omar I. Elashkar
#' @noRd
factor_to_numeric <- function(x) {
  if (is.numeric(x)) {
    return(x)
  } else if (is.factor(x)) {
    return(as.numeric(x))
  } else if (is.character(x)) {
    return(as.numeric(as.factor(x)))
  } else {
    stop("Input must be a factor, character, or numeric vector.")
  }
}

#' Multi-level Categorical Covariate Effects into NLME Model
#' @description Adds multi-level categorical covariate effects into an NLME model by modifying the model's initial parameter estimates and model code to include covariate effects for specified parameters.
#' @param groups A vector of group identifiers (numeric, factor, or character) for each subject. The first unique value is treated as the reference group.
#' @param model An rxode2 model object to be modified.
#' @param cov_name The base name for the covariate parameters to be added (default is "group_code"). The actual parameter names will be constructed as "cov_parameter_group".
#' @param parameters A vector of parameter names (e.g., c("edk50", "kde", "kd", "ks", "gamma")) for which covariate effects should be added.
#' @param fixed_effects A vector of fixed effect names corresponding to the parameters (e.g., c("t.edk50", "t.kde", "t.kd", "t.ks", "t.gamma")). These should be mu-referenced fixed effects that represent the average effect across groups. The function will add covariate effects on top of these fixed effects for each non-reference group.
#' @return A modified rxode2 model object with added covariate effects for the specified parameters based on the provided groups.
#' @details
#' The fixed effects must be mu-referenced.
#' @noRd
#' @author Omar I. Elashkar
parse_covariate <- function(
  groups,
  model = kpd_mod,
  cov_name = "group_code",
  parameters = c("edk50", "kde", "kd", "ks", "gamma"),
  fixed_effects = c("t.edk50", "t.kde", "t.kd", "t.ks", "t.gamma")
) {
  model <- rxensure(model)

  if (is.null(parameters) && is.null(fixed_effects)) {
    return(model)
  }
  if (xor(is.null(parameters), is.null(fixed_effects))) {
    stop(
      "`parameters` and `fixed_effects` must both be NULL or both be non-NULL"
    )
  }

  checkmate::assertIntegerish(groups, lower = 0)
  stopifnot(length(parameters) == length(fixed_effects))
  groups <- factor_to_numeric(groups)

  n_groups <- unique(length(groups))
  if (n_groups > 5) {
    stop("Currently, GatorpH cannot handle 10 groups")
  }

  n_parameters <- n_groups - 1 # reference group has no covariate effect

  # cov parameters formula: cov..parameter..group
  cov_parameters <- expand.grid(
    parameter = parameters,
    group = unique(groups)[-1] # exclude reference group
  ) |>
    dplyr::mutate(cov_name = paste0("cov_", .data$parameter, "_", .data$group))

  oini <- model$iniDf
  ntheta <- max(as.numeric(oini$ntheta), na.rm = TRUE)
  ini <- dplyr::bind_rows(
    oini,
    cov_parameters |>
      dplyr::select(cov_name) |>
      dplyr::mutate(ntheta = ntheta + dplyr::row_number()) |>
      dplyr::mutate(
        name = .data$cov_name,
        est = 0.5,
        upper = Inf,
        lower = -Inf,
        fix = FALSE
      ) |>
      dplyr::select(name, ntheta, est, lower, upper, fix)
  ) |>
    dplyr::arrange(ntheta, neta1, neta2, name)

  # new ini
  new_mod <- model

  # cov ifelse block by group
  ref_group <- unique(groups)[1]
  other_groups <- unique(groups)[-1]

  cov_ifelse_lines <- c()

  add_covariate_to_param_line <- function(lines, param) {
    assign_pattern <- paste0("^\\s*", param, "\\s*<-\\s*(.*)$")
    idx <- grep(assign_pattern, lines, perl = TRUE)
    if (length(idx) == 0) {
      stop("Could not find assignment starting with '", param, " <-'")
    }

    has_symbol <- function(expr, symbol_name) {
      if (is.symbol(expr)) {
        return(identical(as.character(expr), symbol_name))
      }
      if (!is.call(expr)) {
        return(FALSE)
      }
      any(vapply(
        as.list(expr)[-1],
        has_symbol,
        logical(1),
        symbol_name = symbol_name
      ))
    }

    inject_covariate_effect <- function(expr, param) {
      cov_name <- paste0("cov_", param)
      eta_name <- paste0("eta.", param)

      if (has_symbol(expr, cov_name)) {
        return(expr)
      }

      cov_expr <- substitute(1 - x, list(x = as.name(cov_name)))

      flatten_plus_terms <- function(node) {
        if (
          is.call(node) &&
            identical(as.character(node[[1]]), "+") &&
            length(node) == 3
        ) {
          c(flatten_plus_terms(node[[2]]), flatten_plus_terms(node[[3]]))
        } else {
          list(node)
        }
      }

      rebuild_plus_terms <- function(terms) {
        if (length(terms) == 1) {
          return(terms[[1]])
        }
        Reduce(function(lhs, rhs) call("+", lhs, rhs), terms)
      }

      if (has_symbol(expr, eta_name)) {
        terms <- flatten_plus_terms(expr)
        eta_idx <- which(vapply(
          terms,
          function(term) {
            is.symbol(term) && identical(as.character(term), eta_name)
          },
          logical(1)
        ))

        if (length(eta_idx) == 0) {
          return(call("+", expr, cov_expr))
        }

        for (k in rev(eta_idx)) {
          terms <- append(terms, list(cov_expr), after = k - 1)
        }

        return(rebuild_plus_terms(terms))
      }

      call("+", expr, cov_expr)
    }

    for (i in idx) {
      rhs_txt <- sub(assign_pattern, "\\1", lines[i], perl = TRUE)
      rhs_expr <- parse(text = rhs_txt)[[1]]

      if (
        is.call(rhs_expr) &&
          identical(as.character(rhs_expr[[1]]), "exp") &&
          length(rhs_expr) >= 2
      ) {
        rhs_expr[[2]] <- inject_covariate_effect(rhs_expr[[2]], param)
      } else {
        rhs_expr <- inject_covariate_effect(rhs_expr, param)
      }

      rhs_new <- paste(deparse(rhs_expr, width.cutoff = 500L), collapse = " ")
      lines[i] <- paste0(param, " <- ", rhs_new)
    }

    lines
  }

  # Reference group: all covariate effects = 1
  for (param in parameters) {
    cov_ifelse_lines <- c(
      cov_ifelse_lines,
      paste0("if(", cov_name, " == ", ref_group, ") {"),
      paste0("  cov_", param, " <- 1"),
      "}"
    )
  }

  # Other groups: covariate effects by parameter
  for (grp in other_groups) {
    for (param in parameters) {
      cov_ifelse_lines <- c(
        cov_ifelse_lines,
        paste0("if(", cov_name, " == ", grp, ") {"),
        paste0("  cov_", param, " <- cov_", param, "_", grp),
        "}"
      )
    }
  }

  omodel_lines <- rxode2::modelExtract(new_mod, endpoint = NA)

  for (param in parameters) {
    omodel_lines <- add_covariate_to_param_line(omodel_lines, param)
  }

  # Combine cov_ifelse_lines and omodel_lines, then evaluate as model code
  model_code <- c(cov_ifelse_lines, omodel_lines)

  # Build model code string and evaluate with base R
  new_mod <- rxode2::rxUiDecompress(new_mod)
  new_mod$lstExpr <- as.list(str2lang(paste0(
    "{",
    paste(model_code, collapse = "\n"),
    "}"
  ))[-1])
  new_mod$ini <- ini
  new_mod <- new_mod$fun()
  rxode2::ini(new_mod) <- ini
  new_mod
}

remove_gamma <- function(model) {
  model <- rxensure(model)
  model <- model |>
    rxode2::ini(t.gamma = log(1), eta.gamma = 0) |>
    rxode2::ini(t.gamma = fix)

  model
}


theme_academic <- function(base_size = 18, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Center and bold the title for clarity
      plot.title = element_text(
        size = rel(1.2),
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      # Enhance axis appearance
      axis.title = element_text(size = rel(1.1), face = "bold"),
      axis.text = element_text(size = rel(1.0), color = "black"),
      axis.line = element_line(linewidth = 0.8, color = "black"),

      # Legend formatting
      legend.title = element_text(size = rel(1.0), face = "bold"),
      legend.text = element_text(size = rel(0.9)),
      legend.position = "right",
      # Remove any remaining clutter
      panel.background = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' Calculate Observed vs Predicted Metrics
#' @description Calculates metrics comparing observed values (DV) to predicted values (IPRED) from a fitted model, including Mean Absolute Percentage Error (MAPE), Root Mean Squared Error (RMSE), and Mean Absolute Error (MAE).
#' @param fit A data frame containing at least the columns DV (observed values) and IPRED (predicted values).
#' @return A data frame with calculated metrics: MAPE, RMSE, and MAE.
#' @author Omar I. Elashkar
#' @export
fit_obs_metrics <- function(fit){
  mape <- mean(abs(fit$DV - fit$IPRED) / abs(fit$DV), na.rm = TRUE)
  rmse <- sqrt(mean((fit$DV - fit$IPRED)^2, na.rm = TRUE))
  mae <- mean(abs(fit$DV - fit$IPRED), na.rm = TRUE)
  r_squared <- 1 - sum((fit$DV - fit$IPRED)^2, na.rm = TRUE) / sum((fit$DV - mean(fit$DV, na.rm = TRUE))^2, na.rm = TRUE)
  
  data.frame(
    MAPE = mape,
    RMSE = rmse,
    MAE = mae, 
    R_squared = r_squared
  )
}


check_min_max <- function(x){
  check_min_max <- x |> 
    dplyr::group_by(.data$id, .data$group) |>
    dplyr::summarise(
      min_time = min(.data$time, na.rm = TRUE),
      max_time = max(.data$time, na.rm = TRUE),
      .groups = "drop"
    ) 
    # find lastest min and earliest max across groups
  min_time <- max(check_min_max$min_time, na.rm = TRUE)
  max_time <- min(check_min_max$max_time, na.rm = TRUE)

  c(min_time, max_time)
}