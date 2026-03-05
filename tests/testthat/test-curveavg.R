test_that("simulate multiple subject curve avg", {
    dat <- simulate_steph_curve(
        kpd_mod(3, 0.7, ks = 0.5, kd = 0.1),
        nsub = 10, 
        baseline = 7,
        time = c(0, 5 , 10, 15, 20, 30, 50),
        dose = 100,
        baseline_time = -5,
        ignoreBSV = FALSE,
        group = "A"
        )
        
    dat2 <- simulate_steph_curve(
        kpd_mod(3, 0.7, ks = 0.5, kd = 0.1),
        nsub = 10, 
        baseline = 7,
        time = c(0, 5 , 10, 15, 20, 30, 50),
        dose = 100,
        baseline_time = -5,
        ignoreBSV = FALSE, 
        group = "B"
        )
    dat2 <- dat2 |> dplyr::mutate(id = id + 10) 
    dat <- dplyr::bind_rows(dat, dat2)  |>
        dplyr::mutate(group_code = as.numeric(as.factor(.data$group)))
    plot_pH_time(dat, showAvg = TRUE) |> expect_no_error()


    res <- curve_averaging(dat, grouped = TRUE)
    # run_direct_estimation(avg_to_pHdata(res))
    plot_curve_averaging(res) |> expect_no_error()
})


test_that("simulate multiple subject curve avg no interp", {
    dat <- simulate_steph_curve(
        kpd_mod(3, 0.7, ks = 0.5, kd = 0.1),
        nsub = 10, 
        baseline = 7,
        time = c(0, 10, 15, 20, 30, 50),
        ignoreBSV = FALSE
        )
    plot_pH_time(dat) 
    res <- curve_averaging(dat, interpolate = FALSE)
    plot_curve_averaging(res)  |> expect_no_error()
})
