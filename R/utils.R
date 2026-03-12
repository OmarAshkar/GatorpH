kpd_mod <- function(edk50 = 0.5, kde = 0.5, kd = 0.5, ks = 0.5, 
    gamma = 1, 
    eta.edk50 = 0.1, eta.kde = 0.1, eta.kd = 0.1, eta.ks = 0.1, sigma_add = 0.01){
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
        eta.gamma ~ 0.1

        sigma_add <- sigma_add
    })

    rxode2::model({

        edk50 <- exp(t.edk50  + eta.edk50)
        kde <- exp(t.kde +  eta.kde)
        kd <- exp(t.kd +  eta.kd)
        ks <- exp(t.ks +  eta.ks)
        gamma <- exp(t.gamma + eta.gamma)


        # depot(0) = amt
        # ks <- baseline * kd
        resp(0) = ks/kd # Initial condition for response

        d/dt(depot) = -kde * depot

        IR = kde * depot # drug elimiation rate 
        d/dt(resp) = ks  * (1 - (IR)**gamma / (edk50**gamma + (IR)**gamma)) - kd*resp 

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


# mlx/Ooi parameterization
kpd_mod2 <- function(edk50 = 0.5, kde = 0.5, kd = 0.5, ks = 0.5, gamma = 1, 
    eta.edk50 = 0.1, eta.kde = 0.1, eta.kd = 0.1, eta.ks = 0.1, sigma_add = 0.01){
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
        eta.gamma ~ 0.1

        sigma_add <- sigma_add
    })

    rxode2::model({

        edk50 <- exp(t.edk50 +  eta.edk50)
        kde <- exp(t.kde + eta.kde)
        kd <- exp(t.kd + eta.kd)
        ks <- exp(t.ks + eta.ks)
        gamma <- exp(t.gamma + eta.gamma)

        # depot(0) = amt
        # ks <- baseline * kd
        resp(0) = ks/kd # Initial condition for response

        d/dt(depot) = -kde * depot
        IR = kde * depot
        d/dt(resp) = ks * (1 - depot**gamma/(edk50**gamma + depot**gamma) ) - kd*resp

        resp ~ add(sigma_add)

    })

}

turnover_stim_breakdown_mod <- function(){
    rxode2::ini({
        # ---- PK PARAMETERS ----
        ka = 1 # absorption rate constant (1/h)
        kel = 0.5 # elimination rate constant (1/h)
        V = fix(1.1) # volume of saliva (mL)

        # ---- PD PARAMETERS ----
        R0 = fix(7) # Baseline pH

        t.kout = 0.1 # first-order loss rate of response

        Emax = 0.6 # maximum stimulatory effect on breakdown

        EC50 = 0.3 # concentration for 50% of Emax

        sigma_add <- 0.01

    })

    rxode2::model({
        kout <- t.kout
        # ---- PK MODEL ----
        d / dt(depot) = -ka * depot # Absorption from depot
        d / dt(central) = ka * depot - kel * central # Elimination from central compartment
        cp = central / V # Plasma concentration

        # ---- PD MODEL ----
        kin = R0 * kout # Zero-order production rate of response
        resp(0) = R0 # Initial condition for response
        d / dt(resp) = kin - kout * resp * (1 + Emax * cp / (EC50 + cp))

        resp ~ add(sigma_add)

        auc(0) = 0
        # area under the curve below pH
        if (resp < 5.4) {
            d / dt(auc) = resp # accumulate AUC when below threshold
        }

        # time below pH
        time_under_ph(0) = 0
        if (resp < 5.4) {
            d / dt(time_under_ph) = 0.1
        }
    })
}


getnoVarIds <- function(simRes){
  if(is.null(simRes$id)){
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

check_time_varying <- function(x, group = "id", column){

    baseline_vary <- x |>
      dplyr::group_by(.data[[group]]) |>
      dplyr::summarize(n_baselines = dplyr::n_distinct(.data[[column]]), .groups = "drop") |>
      dplyr::filter(.data$n_baselines > 1)

    nrow(baseline_vary) > 0
}

#' Read pH Data
#' @description Reads pH data from a CSV or Excel file and performs basic validation.
#' @param file_path Path to the data file.
#' @param baseline Baseline pH value to assign for each subject if not present in the data (default is 7).
#' @param baseline_time Time of baseline to adjust time points (default is -5). If NULL, no adjustment is made.
#' @return A data frame containing the pH data.
#' @author Omar I. Elashkar
#' @export
read_pH <- function(file_path, baseline = 7, baseline_time = -5) {
  checkmate::assertFileExists(file_path)

  file_ext <- tools::file_ext(file_path)
  checkmate::assertChoice(file_ext, choices = c("csv", "xls", "xlsx"))
  checkmate::assertNumeric(baseline, lower = 0, upper = 14)
  checkmate::assertNumber(baseline_time, finite = TRUE, null.ok = TRUE, upper = 0, lower = -Inf) 

  if (file_ext %in% c("xls", "xlsx")) {
    dat <- readxl::read_excel(file_path)
    dat <- as.data.frame(dat)
  } else if (file_ext == "csv") {
    dat <- read.csv(file_path)
  } else {
    stop("Unsupported file type. Please provide a CSV or Excel file.")
  }

  if (is.null(dat$baseline)) {
    dat$baseline <- baseline
  } else{
    # make sure not NA
    if(any(is.na(dat$baseline))){
      stop("Baseline cannot be NA. Please fill in missing baseline values or set baseline parameter.")
    }

    # make sure not time varying
    baseline_vary <- check_time_varying(dat, group = "id", column = "baseline")
    if (baseline_vary) {
      stop("Baseline varies within subject(s). Please ensure baseline is constant for each subject or set baseline parameter.")
    }

  }
  if(is.null(dat$flowrate)){
    dat$flowrate <- NA
  }
  if(is.null(dat$buffering)){
    dat$buffering <- NA
  }
  
  if(is.null(dat$group)){
    dat$group <- "default"
  } 
  

  if(!is.null(baseline_time)){
    dat <- split(dat, dat$id) |> lapply(function(df) {
        predose = data.frame(
          time = baseline_time,
          pH = df$baseline[2],
          id = df$id[1],
          group = df$group[1],
          baseline = df$baseline[1],
          flowrate = df$flowrate[1],
          buffering = df$buffering[1]
        )
        df <- dplyr::bind_rows(predose, df)
        df
      })
    dat <- do.call(rbind, dat) |>
      dplyr::mutate(time = time + abs(baseline_time)) |> 
      dplyr::arrange(id, time)
  }
  
  dat$group_code <- factor_to_numeric(dat$group)
  dat <- dat |> 
    dplyr::select(-"baseline") |>
    dplyr::mutate(id = as.numeric(id)) |> 
    janitor::remove_empty(which = c("rows", "cols"))


  check_data(dat)
  dat
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
    must.include = c("id", "pH", "time", "group", "group_code")
  )
  if (!sim) {
    checkmate::assertNumeric(x$pH, lower = 0, upper = 14)
  }
  checkmate::assertNumeric(x$time, finite = TRUE)
  # check id and group are nonempty
  checkmate::assertIntegerish(x$id, any.missing = FALSE, min.len = 1)
  checkmate::assertNumeric(x$id, any.missing = FALSE, min.len = 1)

  checkmate::assertVector(x$group, any.missing = FALSE, min.len = 1)
  checkmate::assertNumeric(x$time, finite = TRUE)

  # ensure group, group_code, and baseline do not vary within each subject
  if ("group" %in% names(x)) {
    if(check_time_varying(x, group = "id", column = "group")){
      stop("Group varies within subject(s). Please ensure group is constant for each subject.")
    }
    
  }

  if ("group_code" %in% names(x)) {
    if(check_time_varying(x, group = "id", column = "group_code")){
      stop("Group code varies within subject(s). Please ensure group code is constant for each subject.")
    }

  }

  if ("baseline" %in% names(x)) {
    baseline_vary <- x |>
      dplyr::group_by(.data$id) |>
      dplyr::summarize(n_baselines = dplyr::n_distinct(.data$baseline), .groups = "drop") |>
      dplyr::filter(.data$n_baselines > 1)
    
    if (nrow(baseline_vary) > 0) {
      stop("Baseline varies within subject(s): ", paste(baseline_vary$id, collapse = ", "))
    }
  }

  # check if more than 10 groups
  if ("group" %in% names(x)) {
    n_groups <- length(unique(x$group))
    if (n_groups > 5) {
      stop("GatorpH cannot handle more than 10 groups. Found ", n_groups, " groups.")
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
  # cross before min pH
  if (any(xph[xtime < timeMin] >= ph_threshold, na.rm = TRUE)) {
    before <- TRUE
  } else {
    before <- FALSE
  }
  if (any(xph[xtime > timeMin] >= ph_threshold, na.rm = TRUE)) {
    after <- TRUE
  } else {
    after <- FALSE
  }

  c(before, after)
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
#' @description Calculates the area under the pH curve below a specified pH threshold using the trapezoidal rule.
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
  interpolate = TRUE
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

  if (interpolate) {
    newdata <- interpolate_pH(x, y)
    stopifnot(nrow(newdata) >= length(x))
    x <- newdata$time
    y <- newdata$pH
  }
  tmpdata <- data.frame(x = x, y = y) |>
    dplyr::filter(x >= time_start & x <= time_end)

  if (all(check_crossing(tmpdata$x, tmpdata$y, ph_threshold = ph_threshold))) {
    tmpdata <- tmpdata |>
      dplyr::filter(abs(y - ph_threshold) <= 0.1 | y < ph_threshold) # keep only points below threshold or close to it

    x <- tmpdata$x
    y <- tmpdata$y

    y <- ph_threshold - y # flip
    if (method == "linear") {
      res <- get_auc_linear(x, y)
    } else if (method == "log") {
      res <- get_auc_log(x, y)
    } else if (method == "linear_down_log_up") {
      res <- get_auc_linear_down_log_up(x, y)
    }

    if(length(res) == 0){
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

  time_under_ph_func <- function(xtime, xph, ph_threshold) {
    newdata <- interpolate_pH(xtime, xph)
    # filter only points below pH threshold
    newdata <- newdata[newdata$pH < ph_threshold, ]

    if(min(newdata$time) == Inf | max(newdata$time) == -Inf){
      return(list(start_time = NA, end_time = NA))
    }

    list(start_time = min(newdata$time, na.rm = TRUE), end_time = max(newdata$time, na.rm = TRUE))
  }

  x |>
    dplyr::group_by(.data$id, .data$group) |>
    dplyr::summarise(
      time_under_ph = time_under_ph_func(.data$time, .data$pH, ph_threshold)
    ) |>
    dplyr::mutate(start_time = .data$time_under_ph$start_time) |>
    dplyr::mutate(end_time = .data$time_under_ph$end_time) |>
    dplyr::mutate(
      time_under_ph = .data$end_time - .data$start_time
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct()


}

calc_area_under_pH <- function(
  x,
  ph_threshold = 5.4,
  time_start = 0,
  time_end = 50,
  method = "linear"
) {
  check_data(x)

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
        interpolate = FALSE
      ),
      area_under_pH = integratepHArea(
        x = time,
        y = pH,
        ph_threshold = ph_threshold,
        time_start = time_start,
        time_end = time_end,
        method = method,
        interpolate = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(across(.cols = c("area_under_pH_no_interpolation", "area_under_pH"),
                  \(x) ifelse(x < 0, 0, x))) |>
    dplyr::mutate(auc = paste0("AUC_", time_start, ",", time_end))
}

calc_pH_min <- function(x) {
  check_data(x)
  x |>
    dplyr::group_by(.data$id, .data$group) |>
    dplyr::summarise(pH_min = min(.data$pH, na.rm = TRUE), .groups = "drop")
}

calc_t_min <- function(x) {
  check_data(x)

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
#' @param baseline_time Time of baseline to adjust time points (default is -5).
#' @param nsub Number of subjects to simulate (default is 1).
#' @param baseline Baseline pH value. Default is 7. Must be either 1 or same number as number of subjects.
#' @param step Time step for simulation (default is 0.1).
#' @param group Group identifier for the subjects (default is "A").
#' @param ignoreBSV Logical indicating whether to ignore between-subject variability (default is TRUE).
#' @param ignoreRUV Logical indicating whether to ignore residual unexplained variability (default is TRUE).
#' @param include_gamma logical indicating whether to include the gamma parameter in the simulation (default is TRUE).
#' @return A data frame containing the simulated pH-time data.
#' @author Omar I. Elashkar
simulate_steph_curve <- function(
  model,
  time = seq(0, 50, by = 0.1),
  dose = 100,
  baseline_time = -5,
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
  checkmate::assertNumeric(baseline_time, lower = -Inf, upper = 0)
  checkmate::assertNumeric(baseline, lower = 0, upper = 14, null.ok = TRUE)
  
  if(!include_gamma){
    model <- remove_gamma(model)
  }

  if(!is.null(baseline)){
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

  ev <- rxode2::et(amt = dose, cmt = "depot", time = abs(baseline_time))
  ev <- ev |>
    rxode2::et(time = unique(c(0, time + abs(baseline_time)))) |>
    rxode2::et(id = seq(1, nsub))
  if (ignoreBSV) {
        # resp ~ add(sigma_add)
    model <- rxode2::zeroRe(model, which = "omega")
  }
  if (ignoreRUV) {
    model <- rxode2::zeroRe(model, which = "sigma")
  }

  basecovariates <- c("ks", "kd", "kde", "edk50")
  if(include_gamma){
    basecovariates <- c(basecovariates, "gamma")
  }
  model <- parse_covariate(1, model, parameters = basecovariates, 
    fixed_effects = paste0("t.", basecovariates))

  group <- group
  group_code <- factor_to_numeric(group)
  
  if(!is.null(baseline)){
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
  sim <- rxode2::rxSolve(model, events = ev, iCov = covdf)

  if (nsub == 1) {
    message("only one subject simulated, setting id to 1")
    sim <- sim |> dplyr::mutate(id = 1)
  }
  idx <- 1
  repeat {
    problematic_ids <- getnoVarIds(sim)

    if(length(problematic_ids) == 0 | idx > 10) {
      if(idx == 10){browser()}
      break
    }
    message(
      "Regenerating subjects with no variability in pH response: ",
      paste(problematic_ids, collapse = ", ")
    )
    sim2 <- rxode2::rxSolve(model, 
      events = ev |> dplyr::filter(.data$id %in% problematic_ids),
      iCov = covdf |> dplyr::filter(.data$id %in% problematic_ids),
      nSub = length(problematic_ids))
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

  if(length(getnoVarIds(sim)) != 0){
    stop("Some subjects still have no variability in pH response after regeneration attempts.")
  }
  
  sim <- as.data.frame(sim)

  sim$id <- as.numeric(sim$id)
  stopifnot(sum(is.na(sim$id)) == 0 )

  sim$group <- group
  sim$group_code <- group_code
  sim$pH <- sim$resp

  sim$group_code <- factor_to_numeric(sim$group)

  sim
}

summarize_sim <- function(res) {
  res |>
    group_by(.data$id, .data$group) |>
    summarize(aucUnderpH = max(auc), timeUnderpH = max(time_under_ph))
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
#' @return A ggplot2 object representing the pH vs time plot.
#' @author Omar I. Elashkar
#' @export
plot_pH_time <- function(
  res,
  ph_threshold = 5.4,
  show_id = TRUE,
  stratify_by = "None", 
  showAvg = FALSE
) {
  checkmate::assertChoice(
    stratify_by,
    choices = c("None", "Subject", "Group")
  )

  res$id <- as.factor(res$id)
  plt <- ggplot2::ggplot(
    res,
    ggplot2::aes(x = time, y = pH, color = id)
  ) +
    ggplot2::geom_line(aes(group = id)) +
    ggplot2::labs(x = "Time", y = "pH") +
    ggplot2::geom_hline(
      yintercept = ph_threshold,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::geom_point() + 
    scale_color_viridis_d()
  if (!show_id) {
    plt <- plt + ggplot2::theme(legend.position = "none")
  }
  if(showAvg & stratify_by != "Subject"){
    grouped <- stratify_by == "Group"
    
    restmp <- res |> dplyr::mutate(id = as.numeric(as.character(id)))
    avgdata <- curve_averaging(restmp, interpolate = TRUE, grouped = grouped)

    if(stratify_by != "Group"){
      avgdata$group <- "Average"
    }
    plt <- plt + ggplot2::geom_line(data = avgdata, 
      aes(x = time, y = mean, color = group), linewidth = 1.5)
  }
  
  if (stratify_by == "Group") {
    plt <- plt + ggplot2::facet_wrap(~group)
  }
  if(stratify_by == "Subject"){
    plt <- plt + ggplot2::facet_wrap(~id)
  }

  plt
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
#' @param amt Dose amount to administer at time 0.
#' @param stratify Logical indicating whether to fit separate models for each group (default is FALSE).
#' @param estmethod Estimation method to use (default is "focei").
#' @param dose_time Time (positive) of baseline to add before time 0 (default is 5). 
#' @param include_gamma logical indicating whether to include the gamma parameter in the model (default is TRUE).
#' @return nlmixr2 fit object.
#' @author Omar I. Elashkar
#' @export
fit_pH_curve <- function(data, model, amt, stratify = FALSE, estmethod = "focei", dose_time = 5,
  cov_params = c("kd", "kde", "edk50", "gamma"), cov_fixedeffects = c("t.kd", "t.kde", "t.edk50", "t.gamma"), include_gamma = TRUE
) {

  if ("ks" %in% cov_params && "kd" %in% cov_params) {
    stop("cov_params cannot contain both 'ks' and 'kd'")
  }
  
  check_data(data, sim = TRUE)
  checkmate::assertNumber(amt, finite = TRUE)
  checkmate::assertLogical(stratify, len = 1)
  checkmate::assertChoice(estmethod, choices = c("focei", "saem", "bobyqa", "uobyqa"))
  checkmate::assertNumber(dose_time, finite = TRUE, upper = Inf, lower = 0)
  
  model <- rxensure(model)
  if(!include_gamma){
    model <- remove_gamma(model)
  }
  
  # preserve group and group_code information
  group_info <- data |>
    dplyr::select(id, group, group_code) |>
    dplyr::distinct()

  data <- data |>
    dplyr::mutate(evid = 0, amt = NA_integer_, cmt = NA_character_)

  data <- split(data, data$id) |>
    lapply(function(df) {
      df |> 
        dplyr::filter(time != dose_time) |> 
        add_row(
          pH = NA,
          id = df$id[1],
          group = df$group[1],
          group_code = df$group_code[1],
          cmt = "depot",
          evid = 1,
          amt = amt,
          time = dose_time,
          .before = 1
        )
    })

  data <- do.call(rbind, data) 

  # check data
  check_data(data)
  # assert correct dose time for predose for each subject exist 
  tmpdat <- data |> dplyr::group_by(.data$id) |>
    dplyr::filter(any(.data$time == dose_time)) |>
    dplyr::summarise(n_predose = sum(.data$time == dose_time, na.rm = TRUE), .groups = "drop") |>
    dplyr::ungroup() |>
    dplyr::filter(n_predose == 0)
  if(nrow(tmpdat) > 0){
    stop("Missing predose time point for subject(s): ", paste(tmpdat$id, collapse = ", "))
  }
  

  data <-   data |> 
    dplyr::rename( DV = "pH") |> 
    select("id", "group", "group_code", "time", "DV", "evid", "amt", "cmt") |>
    dplyr::arrange(id, time)

  uniqueids <- unique(data$id)



  if(length(uniqueids) > 1 && length(unique(data$group_code)) > 1){
    model <- parse_covariate(unique(data$group_code), model, parameters= cov_params, 
      fixed_effects = cov_fixedeffects)
  } else{
    warning("Only one group or one subject in the data. Not including group covariate in the model.")
  }

  if (estmethod == "bobyqa" || estmethod == "uobyqa") {
    finalFit <- nlmixr2est::nlmixr2(
      zeroRe(model, which = "omega"),
      data,
      est = estmethod
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

    if(estmethod == "focei"){
      ctrl <- nlmixr2est::foceiControl(
        maxOuterIterations = 100,
        maxInnerIterations = 100
      )
    } else if(estmethod == "saem"){
      ctrl <- nlmixr2est::saemControl()
    } 

    finalFit <- nlmixr2est::nlmixr2(
      model, 
      data,
      est = estmethod, 
      control = ctrl
    )
  }

  finalFit |> 
    dplyr::mutate(ID = as.numeric(as.character(.data$ID)))|>
    # readd groups
    dplyr::left_join(group_info, by = c("ID" = "id"))
}

fit_param_table <- function(fit) {
  fit$parFixed
}

fit_format_param_table <- function(fit) {
  fit$parFixed
}

fit_obs_vs_ipred_plot <- function(fit) {
  as.data.frame(fit) |>
    ggplot2::ggplot(aes(x = DV, y = IPRED)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
    ggplot2::labs(x = "Observed pH", y = "Predicted pH") +
    ggplot2::theme_minimal()
}

fit_individual_plot <- function(fit) {
  as.data.frame(fit) |>
    dplyr::mutate(ID = as.factor(ID)) |>
    ggplot2::ggplot(aes(x = TIME, y = DV, color = ID)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(aes(y = IPRED)) +
    ggplot2::facet_wrap(~ID) +
    ggplot2::labs(x = "Time", y = "pH") +
    ggplot2::theme_minimal()
}



#' Summarize Direct Estimation Results
#' @description Provides summary statistics for direct estimation results.
#' @param object Data frame containing direct estimation results with class 'QuantPH'.
#' @return A data frame with summary statistics for AUC, time under pH threshold, minimum pH, and time to minimum pH.
#' @author Omar I. Elashkar
#' @method summary QuantPH
#' @export
summary.QuantPH <- function(object, ...) {
  object |>
    dplyr::group_by(.data$group) |>
    dplyr::summarize(
      mean_auc = mean(.data$area_under_pH, na.rm = TRUE),
      median_auc = median(.data$area_under_pH, na.rm = TRUE),
      min_auc = min(.data$area_under_pH, na.rm = TRUE),
      max_auc = max(.data$area_under_pH, na.rm = TRUE),
      sd_auc = sd(.data$area_under_pH, na.rm = TRUE),

      mean_time_under_pH = mean(.data$time_under_pH, na.rm = TRUE),
      median_time_under_pH = median(.data$time_under_pH, na.rm = TRUE),
      min_time_under_pH = min(.data$time_under_pH, na.rm = TRUE),
      max_time_under_pH = max(.data$time_under_pH, na.rm = TRUE),
      sd_time_under_pH = sd(.data$time_under_pH, na.rm = TRUE),

      mean_pH_min = mean(.data$pH_min, na.rm = TRUE),
      median_pH_min = median(.data$pH_min, na.rm = TRUE),
      min_pH_min = min(.data$pH_min, na.rm = TRUE),
      max_pH_min = max(.data$pH_min, na.rm = TRUE),

      mean_t_min = mean(.data$t_min, na.rm = TRUE),
      median_t_min = median(.data$t_min, na.rm = TRUE),
      min_t_min = min(.data$t_min, na.rm = TRUE),
      max_t_min = max(.data$t_min, na.rm = TRUE),
    )
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
  method = "linear"
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
    time_end = time_end
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
  dose = 100,
  plot = FALSE,
  stratify_by = "None",
  include_gamma = TRUE
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
    dplyr::pull("time") |> unique()
  baseline_time <- -dose_time


  fixedparamdf <- x$parFixedDf |>
    tibble::rownames_to_column("name") |>
    dplyr::select(name, Estimate) |>
    tidyr::pivot_wider(names_from = "name", values_from = "Estimate") |> 
    dplyr::select(dplyr::starts_with("t."))

  has_covariates <- any(grepl("cov_", rownames(x$parFixedDf)))
  if(has_covariates){
    if(length(unique(x$origData$group_code)) > 1){
      icovDf <- as.data.frame(x) |> 
        dplyr::select("ID", dplyr::starts_with("group"), dplyr::starts_with("eta."), dplyr::starts_with("cov_")) |> 
        dplyr::mutate(group_code = factor_to_numeric(.data$group)) |>
        tidyr::pivot_longer(cols = dplyr::starts_with("cov_"), 
          names_to = "covariate", values_to = "value") |>
        dplyr::mutate(covariate = paste0(.data$covariate, "_", .data$group_code)) |> 
        dplyr::distinct() |>
        tidyr::pivot_wider(names_from = "covariate", values_from = "value", values_fill = 0) |>
        dplyr::distinct() |>
        dplyr::mutate(ID = as.numeric(as.character(.data$ID)))
    } else {
      icovDf <- as.data.frame(x) |> 
        dplyr::select("ID", dplyr::starts_with("group"), dplyr::starts_with("eta."), dplyr::starts_with("cov_")) |> 
        dplyr::mutate(group_code = factor_to_numeric(.data$group)) |>
        # add covariate manually 
        mutate(cov_ks = 0, cov_edk50 = 0, cov_kd = 0, cov_kde = 0, cov_gamma = 0) |>
        tidyr::pivot_longer(cols = dplyr::starts_with("cov_"), 
          names_to = "covariate", values_to = "value") |>
        dplyr::mutate(covariate = paste0(.data$covariate, "_", .data$group_code)) |> 
        dplyr::distinct() |>
        tidyr::pivot_wider(names_from = "covariate", values_from = "value", values_fill = 0) |>
        dplyr::distinct() |>
        dplyr::mutate(ID = as.numeric(as.character(.data$ID)))
    }
  } else { # if no covariates, just get ID and group information
    icovDf <- as.data.frame(x) |> 
      dplyr::select("ID", dplyr::starts_with("group"), dplyr::starts_with("eta.")) |> 
      dplyr::mutate(group_code = factor_to_numeric(.data$group)) |>
      dplyr::distinct() |>
      dplyr::mutate(ID = as.numeric(as.character(.data$ID))) 
  }

  new_mod <- model |> rxode2::zeroRe(which = "sigma")
  
  ev <- rxode2::et(amt = dose, cmt = "depot", time = dose_time) |> # use dose_time here
      rxode2::et(time = time) 

  newini <- new_mod$iniDf |> 
    dplyr::filter(!grepl("eta\\.", .data$name), !grepl("cov_", .data$name)) 

  # assert fixed effects parameter are final 
  newini <- newini |> 
    dplyr::mutate(est = dplyr::case_when(
      name == "t.edk50" ~ fixedparamdf$t.edk50,
      name == "t.kde" ~ fixedparamdf$t.kde,
      name == "t.kd" ~ fixedparamdf$t.kd,
      name == "t.ks" ~ fixedparamdf$t.ks,
      name == "t.gamma" ~ fixedparamdf$t.gamma,
      TRUE ~ est
    ))

  rxode2::ini(new_mod) <- newini

  if (onlymean) {
    icovDf <- icovDf |> 
      # dplyr::mutate(across(starts_with("eta."), ~ 0)) |>
      dplyr::mutate(eta.edk50 = 0, eta.kde = 0, eta.kd = 0, eta.ks = 0, eta.gamma = 0) |> # ensure regressors
      dplyr::distinct() 
  } 

  if(!include_gamma){
    icovDf <- icovDf |> 
      dplyr::mutate(eta.gamma = 0) 
  }

  
  ids <- as.numeric(as.character(unique(icovDf$ID)))
  ev <- rxode2::et(amt = dose, cmt = "depot", time = dose_time) |> # use dose_time here
    rxode2::et(time = time) |>
    rxode2::et(id = ids)
  
  simRes <- rxSolve(
    new_mod,
    iCov = icovDf, # TODO add flowrate, buffering, substance
    events = ev
  )

  # fix any subjects with no variability in pH response
  idx <- 1
  repeat{
    problematic_ids <- getnoVarIds(simRes)

    if(length(problematic_ids) == 0 | idx > 10) {
      if(idx == 10){browser()}
      break
    }
    message(
      "Regenerating subjects with no variability in pH response: ",
      paste(problematic_ids, collapse = ", ")
    )
    sim2 <- rxSolve(
      new_mod,
      iCov = icovDf |> dplyr::filter(.data$ID %in% problematic_ids),
      events = ev |> dplyr::filter(.data$id %in% problematic_ids)
    )
    # replace only the problematic ids with regenerated subjects
    simRes <- simRes |>
      dplyr::filter(!(.data$id %in% problematic_ids)) |>
      dplyr::bind_rows(sim2)
    idx <- idx + 1
  }

  stopifnot(length(getnoVarIds(simRes)) == 0)
  
  if(is.null(simRes$id)){
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
    dplyr::select("ID", "group", "group_code") |>
    dplyr::distinct() |> 
    dplyr::mutate(id = as.numeric(as.character(.data$ID))) |>
    dplyr::select(-"ID")
  
  simRes <- simRes |>
    dplyr::select(-dplyr::starts_with("group"), -dplyr::starts_with("group_code")) |>
    dplyr::left_join(orig_groups, by = c("id" = "id")) 

  # no matter if pooled or not, all original groups must have calculations
  stopifnot(length(unique(simRes$group)) == length(unique(x$origData$group)))  
  stopifnot(unique(simRes$group) == unique(x$origData$group))
  
  if (plot) {
    plt <- plot_pH_time(
      simRes,
      show_id = FALSE,
      stratify_by = stratify_by,
      ph_threshold = ph_threshold
    ) +
      labs(
        subtitle = paste0(ifelse(onlymean, "Mean Profile", "Individual Profiles"),
          " from method ", estMethod)
      )
      if(!onlymean){
        originalData <- nlme::getData(x) |> dplyr::filter(.data$evid == 0)
        plt <- plt + ggplot2::geom_point(
          data = originalData, 
          aes(x = .data$time, y = .data$DV),
          color = "red",
          size = 1, 
          shape = 4
        )
      }

    print(plt)
  } else {
    plt <- NA
  }

  derivedDf <- run_direct_estimation(simRes, ph_threshold = ph_threshold) |> 
    as.data.frame() |> 
    dplyr::rename(ID = "id")

  
  if(onlymean){ 
      derivedDf <- derivedDf |> 
        # replace actual id with placeholder
        dplyr::mutate(ID = ifelse(onlymean, 
        ".", 
        as.numeric(as.character(.data$ID)))) 
    }

  stopifnot(nrow(derivedDf) == nrow(icovDf))

  if(onlymean){
    paramsdf <- as.data.frame(x) |> 
      dplyr::mutate(ID = ".")
  } else {
    paramsdf <- as.data.frame(x) |> 
      dplyr::mutate(ID = as.numeric(as.character(.data$ID)))
  }

  param_cols <- c("edk50", "kde", "kd", "ks")
  if ("gamma" %in% names(paramsdf)) {
    param_cols <- c(param_cols, "gamma")
  }
  
  derivedDf <-  dplyr::left_join(
      derivedDf,
      paramsdf |>
        dplyr::select(dplyr::all_of(c("ID", param_cols, "group"))) |>
        dplyr::group_by(.data$ID, .data$group) |>
        dplyr::summarize(across(dplyr::all_of(param_cols), mean), .groups = "keep") |>
        dplyr::ungroup() |>
        dplyr::distinct(),
      by = c("ID" = "ID", "group" = "group")
    )
  
      

  if(onlymean){
    derived_summarize_cols <- c("edk50", "kde", "kd", "ks", "pH_min", "t_min", "start_time", "end_time", "time_under_ph", "area_under_pH", "area_under_pH_no_interpolation")

    if ("gamma" %in% names(derivedDf)) {
      derived_summarize_cols <- c("edk50", "kde", "kd", "ks", "gamma", "pH_min", "t_min", "start_time", "end_time", "time_under_ph", "area_under_pH", "area_under_pH_no_interpolation")
    }
    
    derivedDf <- derivedDf |> 
      dplyr::group_by(.data$group, .data$auc) |>
      dplyr::summarize(across(dplyr::all_of(derived_summarize_cols), mean), .groups = "keep") |>
      dplyr::ungroup() |>
      dplyr::mutate(ID = ".")
    stopifnot(nrow(derivedDf) == length(unique(x$origData$group_code)))
  } else{
    stopifnot(nrow(derivedDf) == length(unique(x$origData$id)))
  }
  stopifnot(sort(unique(derivedDf$group)) == sort(unique(x$origData$group)))

  derivedDf
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
      dplyr::select(time, pH, id, dplyr::all_of(if (grouped) "group" else character()))

    if (!grouped) {
      df <- df |>
        dplyr::mutate(group = "All")
    }

    if(interpolate){
      df <- df |>
        group_by(.data$id, .data$group) |>
        dplyr::reframe(
          x_new = xout,
          y_new = approx(time, pH, xout = xout)$y
        ) |>
        dplyr::filter(!is.na(.data$y_new)) |> # Remove rows with NA values
        dplyr::ungroup() |> 
        dplyr::rename(pH = .data$y_new, time = .data$x_new)
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
#' - If input is numeric, returns as-is
#' - If input is factor, converts to character then to numeric (preserving original numeric values)
#' - If input is character, converts to factor then to numeric (assigning numeric codes to unique values)
#' @author Omar I. Elashkar
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
parse_covariate <- function(groups, 
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
    stop("`parameters` and `fixed_effects` must both be NULL or both be non-NULL")
  }
  
  checkmate::assertIntegerish(groups, lower = 0)
  stopifnot(length(parameters) == length(fixed_effects))
  groups <- factor_to_numeric(groups)

  n_groups <- unique(length(groups))
  if(n_groups > 5){
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
  ini <- dplyr::bind_rows(oini,
    cov_parameters |> 
      dplyr::select(cov_name) |> 
      dplyr::mutate(ntheta = ntheta + dplyr::row_number()) |>
      dplyr::mutate(name = .data$cov_name, est = 0.5, upper = Inf, lower = -Inf, fix = FALSE) |>
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
      any(vapply(as.list(expr)[-1], has_symbol, logical(1), symbol_name = symbol_name))
    }

    inject_covariate_effect <- function(expr, param) {
      cov_name <- paste0("cov_", param)
      eta_name <- paste0("eta.", param)

      if (has_symbol(expr, cov_name)) {
        return(expr)
      }

      cov_expr <- substitute(1 - x, list(x = as.name(cov_name)))

      flatten_plus_terms <- function(node) {
        if (is.call(node) && identical(as.character(node[[1]]), "+") && length(node) == 3) {
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
          function(term) is.symbol(term) && identical(as.character(term), eta_name),
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

      if (is.call(rhs_expr) && identical(as.character(rhs_expr[[1]]), "exp") && length(rhs_expr) >= 2) {
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
    cov_ifelse_lines <- c(cov_ifelse_lines,
      paste0("if(", cov_name, " == ", ref_group, ") {"),
      paste0("  cov_", param, " <- 1"),
      "}")
  }
  
  # Other groups: covariate effects by parameter
  for (grp in other_groups) {
    for (param in parameters) {
      cov_ifelse_lines <- c(cov_ifelse_lines,
        paste0("if(", cov_name, " == ", grp, ") {"),
        paste0("  cov_", param, " <- cov_", param, "_", grp),
        "}")
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
  new_mod$lstExpr <- as.list(str2lang(paste0("{", paste(model_code, collapse="\n"), "}"))[-1])
  new_mod$ini <- ini
  new_mod <- new_mod$fun()
  rxode2::ini(new_mod) <- ini
  new_mod
}

remove_gamma <- function(model){
  model <- rxensure(model)
  model <- model |> 
    rxode2::ini(t.gamma = log(1), eta.gamma = 0) |> 
    rxode2::ini(t.gamma = fix, eta.gamma = fix) 

  model
}