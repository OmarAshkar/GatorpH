test_that("make_covariates",{
    parse_covariate(c(1, 2, 3), model = kpd_mod()) |> expect_no_error()
    parse_covariate(c(1), model = kpd_mod()) |> expect_no_error()

    factor_to_numeric(c("A", "B", "C", "A")) |> expect_equal(c(1, 2, 3, 1))
    
    # no covariates
    parse_covariate(c(1:3), model = kpd_mod(), 
        parameters = NULL, 
        fixed_effects = NULL
    ) |> expect_no_error()

})

test_that("simulate single subject", {
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 1,
        baseline = 7,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    ) |>
        expect_no_error()
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 1,
        baseline = 1:10,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    ) |>
        expect_error()
})

test_that("simulate multiple subject", {
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 1,
        baseline = rnorm(10, 7, 0.1),
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    ) |>
        expect_error()

    suppressWarnings({
        dat <- simulate_steph_curve(
            kpd_mod(0.8, 0.8, 0.1, 0.8),
            nsub = 10,
            baseline = rnorm(10, 7, 0.1),
            time = c(0, 10, 15, 20, 30),
            ignoreBSV = FALSE
        ) |>
            expect_no_error()
    })
})

test_that("Naive pool single sub", {
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 1,
        baseline = 7,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    )

    plot_pH_time(dat)

    fit <- fit_pH_curve(
        dat,
        model = kpd_mod(0.8, 0.8, 0.1, 0.8),
        amt = 100,
        estmethod = "uobyqa", 
        cov_params = "ks",
        cov_fixedeffects = "t.ks"
    )

    fit_individual_plot(fit) |> expect_no_error()

    pHdatMean <- pHMetrics_from_fit(
        fit,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        dose = 100,
        plot = TRUE
    )
    expect_true(pHdatMean$area_under_pH > 0)
    expect_true(pHdatMean$ID == ".")
    expect_true(pHdatMean$group == "A")
    expect_true(pHdatMean$edk50 > 0)
    expect_true(pHdatMean$kde > 0)
    expect_true(pHdatMean$kd > 0)
    expect_true(pHdatMean$ks > 0)

})




test_that("Naive pool multiple subs, multiple groups", {
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 30,
        baseline = rnorm(30, 7, 0.1),
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    )
    dat2 <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 30,
        baseline = rnorm(30, 7, 0.1),
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    ) |> dplyr::mutate(group = "B", group_code = 2) |> 
        dplyr::mutate(id = id + 30)

    dat <- dplyr::bind_rows(dat, dat2)

    plot_pH_time(dat)

    fit <- fit_pH_curve(
        dat,
        model = kpd_mod(0.8, 0.8, 0.1, 0.8),
        amt = 100,
        estmethod = "uobyqa"
    )

    fit_individual_plot(fit) |> expect_no_error()

    pHdatMean <- pHMetrics_from_fit(
        fit,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        dose = 100,
        plot = TRUE
    )
    expect_true(nrow(pHdatMean) == 2)
    expect_true(all(pHdatMean$area_under_pH > 0))
    expect_true(all(pHdatMean$ID == c(".", ".")))
    expect_true(all(unique(pHdatMean$group) == c("A", "B")))
    expect_true(length(unique(pHdatMean$edk50)) == 2)
})

test_that("NLME fit 1", {
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8, 1, 0.1, 0.1, 0.1, 0.1, 0.01),
        nsub = 30,
        baseline = rnorm(30, mean = 7, sd = 0.0),
        baseline_time = -5,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    )
    dat2 <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 30,
        baseline = rnorm(30, 7, 0),
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    ) |> dplyr::mutate(group = "B", group_code = 2) |> 
        dplyr::mutate(id = as.factor(as.numeric(as.character(id)) + 30))

    dat <- dplyr::bind_rows(dat, dat2)
    plot_pH_time(dat, stratify_by = "Group")

    fit <- fit_pH_curve(
        dat,
        model = kpd_mod(0.8, 0.8, 0.1, 0.8, 1, 0.1, 0.1, 0.1, 0.1, 0.01),
        amt = 100,
        estmethod = "saem",
        dose_time = 5
    )

    fit_individual_plot(fit) + ggplot2::facet_wrap(~group) 
    fit_obs_vs_ipred_plot(fit)

    expect_true(all.equal(
        unique(fit$ks / fit$kd),
        rep(7, 60),
        tolerance = 0.05
    ))

    fit$parFixedDf
    # all.equal(unname(fit$theta), c(log(0.8), log(0.8), log(0.1), 1, 0.01))
    suppressWarnings(
        pHdatMean <- pHMetrics_from_fit(
            fit,
            time_start = 0,
            time_end = 50,
            step = 0.1,
            dose = 100,
            plot = TRUE
        )
    )

    expect_true(all(pHdatMean$area_under_pH > 0))
    expect_true(all(pHdatMean$pH_min > 0))
    expect_true(all(pHdatMean$t_min > 0))
    expect_true(all(pHdatMean$time_under_ph > 0))
    
    expect_true(all(unique(pHdatMean$group) == c("A", "B")))
    expect_true(length(unique(pHdatMean$edk50)) == 2)

    suppressWarnings(
        pHdatFull <- pHMetrics_from_fit(
            fit,
            time_start = 0,
            time_end = 50,
            step = 0.1,
            dose = 100,
            plot = TRUE
        )
    )

    expect_true(any(pHdatFull$area_under_pH > 0))
    expect_true(any(pHdatFull$pH_min > 0))
    expect_true(any(pHdatFull$t_min > 0))
    expect_true(any(pHdatFull$time_under_ph > 0))
    expect_false(any(is.infinite(pHdatFull$area_under_pH)))
    expect_false(any(is.infinite(pHdatFull$pH_min)))
    expect_false(any(is.infinite(pHdatFull$t_min)))
    expect_false(any(is.infinite(pHdatFull$time_under_ph)))

    expect_true(all(!is.na(pHdatFull$edk50)))
    expect_true(all(!is.na(pHdatFull$kde)))
    expect_true(all(!is.na(pHdatFull$kd)))
    expect_true(all(!is.na(pHdatFull$ks)))

})


test_that("NLME fit mean2", {
    dat <- simulate_steph_curve(
        kpd_mod(
            edk50 = 0.835,
            kde = 2.07,
            kd = 2.1,
            ks = 1.56,
            1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.01
        ),
        nsub = 50,
        baseline = rnorm(50, mean = 7, sd = 0.0),
        baseline_time = -5,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    )
    plot_pH_time(dat, stratify_by = "Subject")

    fit <- fit_pH_curve(
        dat,
        model = kpd_mod(
            edk50 = 0.835,
            kde = 2.07,
            kd = 2.1,
            ks = 1.56,
            1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.01
        ),
        amt = 100,
        estmethod = "focei",
        dose_time = 5
    )
    fit_individual_plot(fit)

    expect_true(all.equal(
        unique(fit$ks / fit$kd),
        rep(7, 30),
        tolerance = 0.05
    ))
    fit$parFixedDf
    # all.equal(unname(fit$theta), c(log(0.8), log(0.8), log(0.1), 1, 0.01))
    suppressWarnings(
        pHdatMean <- pHMetrics_from_fit(
            fit,
            time_start = 0,
            time_end = 50,
            step = 0.1,
            dose = 100,
            plot = TRUE
        )
    )

    expect_true(pHdatMean$area_under_pH > 0)
    expect_true(pHdatMean$pH_min > 0)
    expect_true(pHdatMean$t_min > 0)
    expect_true(pHdatMean$time_under_ph > 0)
})

test_that("NLME fit full2", {
    dat <- simulate_steph_curve(
        kpd_mod(
            edk50 = 0.835,
            kde = 2.07,
            kd = 2.1,
            ks = 1.56,
            1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.01
        ),
        nsub = 50,
        baseline = rnorm(50, mean = 7, sd = 0.0),
        baseline_time = -5,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    )
    plot_pH_time(dat, stratify_by = "Subject")

    fit <- fit_pH_curve(
        dat,
        model = kpd_mod(
            edk50 = 0.835,
            kde = 2.07,
            kd = 2.1,
            ks = 1.56,
            1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.01
        ),
        amt = 100,
        estmethod = "focei",
        dose_time = 5
    )
    fit_individual_plot(fit)

    expect_true(all.equal(
        unique(fit$ks / fit$kd),
        rep(7, 30),
        tolerance = 0.05
    ))
    fit$parFixedDf
    # all.equal(unname(fit$theta), c(log(0.8), log(0.8), log(0.1), 1, 0.01))

    suppressWarnings(
        pHdatFull <- pHMetrics_from_fit(
            fit,
            time_start = 0,
            time_end = 50,
            step = 0.1,
            dose = 100,
            plot = TRUE
        )
    )

    expect_true(any(pHdatFull$area_under_pH > 0))
    expect_true(any(pHdatFull$pH_min > 0))
    expect_true(any(pHdatFull$t_min > 0))
    expect_true(any(pHdatFull$time_under_ph > 0))
    expect_false(any(is.infinite(pHdatFull$area_under_pH)))
    expect_false(any(is.infinite(pHdatFull$pH_min)))
    expect_false(any(is.infinite(pHdatFull$t_min)))
    expect_false(any(is.infinite(pHdatFull$time_under_ph)))
})


test_that("nlme birkhed", {
    d <- read_pH(system.file(
        "extdata",
        "Birkhed1_copy.csv",
        package = "GatorpH"
    ))

    # naive fit
    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        amt = 100,
        estmethod = "bobyqa",
        dose_time = 5,
        cov_params = "ks",
        cov_fixedeffects = "t.ks",
        include_gamma = TRUE
    )

    fit_individual_plot(fit) 

    suppressWarnings(
        pHdatMean <- pHMetrics_from_fit(
            fit,
            ph_threshold = 6.5,
            time_start = 0,
            time_end = 50,
            step = 0.1,
            dose = 100,
            plot = TRUE
        )
    )
    nrow(pHdatMean) |> expect_equal(2)
    expect_false(any(is.na(pHdatMean$area_under_pH)))
    expect_false(any(is.na(pHdatMean$edk50)))

    ## nlme 
    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        amt = 100,
        estmethod = "saem",
        cov_params = "ks",
        cov_fixedeffects = "t.ks",
        dose_time = 5
    )
    pHdatnlme <- pHMetrics_from_fit(
        fit,
        ph_threshold = 6.5,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        dose = 100,
        plot = TRUE
    )

    expect_false(any(is.na(pHdatnlme$area_under_pH)))
    expect_false(any(is.na(pHdatnlme$edk50)))

    

})

test_that("fejeskov", {
    d <- read_pH(system.file(
        "extdata",
        "Fejerskov Data_GatorpH Test.csv",
        package = "GatorpH"
    ), baseline_time = 0)

    expect_false(any(is.na(d$pH)))

    # test naive fit 
    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        amt = 100,
        estmethod = "uobyqa",
        dose_time = 0+0.01
    )

    fit_individual_plot(fit)
    suppressWarnings(
        pHdatMean <- pHMetrics_from_fit(
            fit,
            ph_threshold = 5.5,
            time_start = 0,
            time_end = 50,
            step = 0.1,
            dose = 100,
            plot = TRUE
        )
    )
    nrow(pHdatMean) |> expect_equal(2)
    expect_false(any(is.na(pHdatMean$area_under_pH)))
    expect_false(any(is.na(pHdatMean$edk50)))

    ## nlme
    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        amt = 100,
        estmethod = "saem",
        dose_time = 0+0.01,
        covmethod = "linFim"
    )
    fit_individual_plot(fit)
    pHdatnlme <- pHMetrics_from_fit(
        fit,
        ph_threshold = 6.5,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        dose = 100,
        plot = TRUE
    )
    expect_false(any(is.na(pHdatnlme$area_under_pH)))
    expect_false(any(is.na(pHdatnlme$edk50)))
})

test_that("remove_gamma", {
    mod <- kpd_mod()
    mod_no_gamma <- remove_gamma(mod)
    expect_false(any(grepl("gamma", names(mod_no_gamma))))

    
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 10,
        baseline = 1,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE, 
        include_gamma= FALSE
    ) 

    fit <- fit_pH_curve(
        dat,
        model = mod_no_gamma,
        amt = 100,
        estmethod = "bobyqa",
        cov_params = "ks",
        cov_fixedeffects = "t.ks", 
        include_gamma = FALSE
    )

    res <- pHMetrics_from_fit(
        fit,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        dose = 100,
        plot = TRUE, 
        include_gamma = FALSE
    )  |> expect_no_error()

    
    fit <- fit_pH_curve(
        dat,
        model = mod_no_gamma,
        amt = 100,
        estmethod = "saem",
        cov_params = "ks",
        cov_fixedeffects = "t.ks", 
        include_gamma = FALSE
    )

    res <- pHMetrics_from_fit(
        fit,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        dose = 100,
        plot = TRUE, 
        include_gamma = FALSE
    )  |> expect_no_error()


    ## 
    
    d <- read_pH(system.file(
        "extdata",
        "Birkhed1_copy.csv",
        package = "GatorpH"
    ))

    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        amt = 100,
        estmethod = "bobyqa",
        dose_time = 5,
        cov_params = "ks",
        cov_fixedeffects = "t.ks",
        include_gamma = FALSE
    )

    pHdatnlme <- pHMetrics_from_fit(
        fit,
        ph_threshold = 6.5,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        dose = 100,
        plot = TRUE, 
        include_gamma = FALSE
    )

    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        amt = 100,
        estmethod = "saem",
        dose_time = 5,
        cov_params = c(),
        cov_fixedeffects = c(),
        include_gamma = FALSE
    )
    pHdatnlme <- pHMetrics_from_fit(
        fit,
        ph_threshold = 6.5,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        dose = 100,
        plot = TRUE, 
        include_gamma = FALSE
    )

    
})
