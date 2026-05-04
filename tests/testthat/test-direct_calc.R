test_that("larssen", {
    dat <- read_pH(system.file("extdata", "Larsen1997_restingpH.csv", package = "GatorpH"),  
        dose_time = 0)

    dat <- dat |> dplyr::filter(id == 1)

    plot_pH_time(dat) |> expect_no_error()

    res <- calc_area_under_pH(dat, ph_threshold = 7, method = "linear", add_support_points = TRUE)$area_under_pH
    expect_true(all.equal(res, 29.275, tol = 0.005))
})

test_that("validated with trapz", {
    res <- simulate_steph_curve(kpd_mod(0.8, 0.8, 0.1, 0.8), nsub = 1)
    plot_pH_time(res)
    
    newres <- res |> dplyr::filter(pH < 5.4)

    # incorrect way, for testing purposes only
    x <- pracma::trapz(newres$time, newres$pH)
    y <- get_auc_linear(newres$time, newres$pH)
    expect_true(all.equal(x, y))

    # correct way
    x <- pracma::trapz(newres$time, (newres$pH - 5.4))
    y <- pracma::trapz(newres$time, (5.4 - newres$pH))
    z <- integratepHArea(res$time, res$pH, ph_threshold = 5.4, interpolate = FALSE)
    z <- integratepHArea(res$time, res$pH, ph_threshold = 5.4, interpolate = TRUE)
    expect_true(all.equal(y, z, tol = 0.03))

    run_direct_estimation(res) 
})

test_that("direct_estimation multiple groups", {
    suppressWarnings({
        res <- simulate_steph_curve(kpd_mod(0.8, 0.8, 0.1, 0.8), nsub = 10, group = "A")
        res2 <- simulate_steph_curve(kpd_mod(0.8, 0.3, 0.1, 0.8), nsub = 10, group = "B")
    })
    res2$id <- as.factor(as.numeric(as.character(res2$id)) + 10)
    res <- dplyr::bind_rows(res, res2)

    plot_pH_time(res, stratify_by = "Group") |> expect_no_error()

    x <- run_direct_estimation(res) |> expect_no_error()
    x$group |> expect_equal(rep(c("A", "B"), each = 10))
})



test_that("check_crossing works", {

    # double crossing
    time <- c(0, 1, 2, 3, 4, 5)
    pH <- c(6, 5.5, 5, 4.5, 5, 6)
    
    all(check_crossing(time, pH)) |> expect_true()
    plot(time, pH)

    res <- interpolate_pH(time, pH)

    expect_equal(max(res$time), max(time))
    expect_equal(min(res$time),  min(time))

    
    # no crossing
    time <- c(0, 1, 2, 3, 4)
    pH <- c(7, 6.5, 6, 5, 5.3)
    all(check_crossing(time, pH)) |> expect_false()

    plot(time, pH)
    res <- interpolate_pH(time, pH)
    expect_equal(max(res$time), max(time))
    expect_equal(min(res$time),  min(time))

    # single crossing (start only)
    time <- c(0, 1, 2, 3, 4)
    pH <- c(5, 5.5, 6, 6.5, 7)
    check_crossing(time, pH) |> expect_equal(c(FALSE, TRUE))

    # single crossing (end only)
    time <- c(0, 1, 2, 3, 4)
    pH <- c(7, 6.5, 6, 5.5, 5)
    check_crossing(time, pH) |> expect_equal(c(TRUE, FALSE))
    plot(time, pH)
})

test_that("integratepHArea works with interpolation", {
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.3, ks = 4.9, kd = 0.7),
        nsub = 1, time = c(0, 5, 10, 30))

    plot_pH_time(dat)

    integratepHArea(
        dat$time,
        dat$pH,
        ph_threshold = 5.7,
        interpolate = TRUE
    )

    res <- calc_area_under_pH(
        dat,
        ph_threshold = 5.7,
        method = "linear"
    )
    # expect_true(is.na(res$area_under_pH_no_interpolation))
    expect_true(!is.na(res$area_under_pH))
    expect_true(res$auc == "AUC_0,50")
})



test_that("integratepHArea works with interpolation and support", {
    dat <- simulate_steph_curve(
        kpd_mod(5, 0.3, ks = 3, kd = 0.7),
        nsub = 1, time = c(0, 5, 10, 30), baseline = 5)

    plot_pH_time(dat)

    
    integratepHArea(
        dat$time,
        dat$pH,
        ph_threshold = 5.7,
        interpolate = TRUE, 
        add_support_points = FALSE
    ) |> expect_equal(NA)

    is.na(integratepHArea(
        dat$time,
        dat$pH,
        ph_threshold = 5.7,
        interpolate = TRUE, 
        add_support_points = TRUE
    )) |> expect_false()

    res <- calc_area_under_pH(
        dat,
        ph_threshold = 5.7,
        method = "linear"
    )
    # expect_true(is.na(res$area_under_pH_no_interpolation))
    expect_true(is.na(res$area_under_pH))
    expect_true(res$auc == "AUC_0,50")
    
    res <- calc_area_under_pH(
        dat,
        ph_threshold = 5.7,
        method = "linear", 
        add_support_points = TRUE
    )
    # expect_true(is.na(res$area_under_pH_no_interpolation))
    expect_true(!is.na(res$area_under_pH))
    expect_true(res$auc == "AUC_0,50")
})



test_that("integratepHArea with support Fadel", {
    dat <- system.file("extdata", "Fadel_2013.csv", package = "GatorpH") 

    # as mentioned in the paper
    dat <- read_pH(dat, dose_time = 0)

    plot_pH_time(dat, 6.5)

    res <- run_direct_estimation(dat, ph_threshold = 6.5)
    expect_true(all(!is.na(res$area_under_pH)))
    expect_true(all(!is.na(res$time_under_ph)))

    
    long_interval <- calc_area_under_pH(
        dat,
        ph_threshold = 7,
        time_start = 0, 
        time_end = 50,
        method = "linear", 
        add_support_points = TRUE,
        plot = TRUE
    )$area_under_pH 
    all(is.na(long_interval)) |> expect_true()
    
    exact_interval <- calc_area_under_pH(
        dat,
        ph_threshold = 7,
        time_start = 0, 
        time_end = 30,
        method = "linear", 
        add_support_points = TRUE,
        plot = TRUE
    )$area_under_pH 
    all(!is.na(exact_interval)) |> expect_true()
    
    short_interval <- calc_area_under_pH(
        dat,
        ph_threshold = 7,
        time_start = 0, 
        time_end = 25,
        method = "linear", 
        add_support_points = TRUE,
        plot = TRUE
    )$area_under_pH 
    all(!is.na(short_interval)) |> expect_true()

    
    long_interval <- calc_area_under_pH(
        dat,
        ph_threshold = 7,
        method = "linear", 
        # add_support_points = TRUE,
        plot = TRUE
    )$area_under_pH 



    res <- run_direct_estimation(dat, ph_threshold = 7, add_support_points = TRUE)
    expect_true(all(!is.na(res$area_under_pH)))
    expect_true(all(!is.na(res$time_under_ph)))


})
