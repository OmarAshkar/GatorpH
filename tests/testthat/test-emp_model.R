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
    ) |> sim_to_pH_data()

    plot_pH_time(dat)

    fit <- fit_pH_curve(
        dat,
        model = kpd_mod(0.8, 0.8, 0.1, 0.8),
        estmethod = "bobyqa", 
        cov_params = "ks",
        cov_fixedeffects = "t.ks"
    )

    fit_individual_plot(fit) |> expect_no_error()

    pHdatMean <- pHMetrics_from_fit(
        fit,
        time_start = 0,
        time_end = 50,
        step = 0.1
    )

    expect_no_error(pHdatMean$plot)

    expect_true(pHdatMean$derivedDf$area_under_pH > 0)
    expect_true(pHdatMean$derivedDf$id == ".")
    expect_true(pHdatMean$derivedDf$group == "A")
    expect_true(pHdatMean$derivedDf$edk50 > 0)
    expect_true(pHdatMean$derivedDf$kde > 0)
    expect_true(pHdatMean$derivedDf$kd > 0)
    expect_true(pHdatMean$derivedDf$ks > 0)

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

    dat <- dplyr::bind_rows(dat, dat2) |> sim_to_pH_data()

    plot_pH_time(dat)

    fit <- fit_pH_curve(
        dat,
        model = kpd_mod(0.8, 0.8, 0.1, 0.8),
        estmethod = "bobyqa"
    )

    fit_individual_plot(fit) |> expect_no_error()

    pHdatMean <- pHMetrics_from_fit(
        fit,
        time_start = 0,
        time_end = 50,
        step = 0.1
    )
    expect_true(nrow(pHdatMean$derivedDf) == 2)
    expect_true(all(pHdatMean$derivedDf$area_under_pH > 0))
    expect_true(all(pHdatMean$derivedDf$id == c(".", ".")))
    expect_true(all(unique(pHdatMean$derivedDf$group) == c("A", "B")))
    expect_true(length(unique(pHdatMean$derivedDf$edk50)) == 2)
})

test_that("NLME fit 1", {
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8, 1, 0.1, 0.1, 0.1, 0.1, 0.01),
        nsub = 30,
        baseline = rnorm(30, mean = 7, sd = 0.0),
        dose_time = 1,
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
        dplyr::mutate(id = id + 30)

    dat <- dplyr::bind_rows(dat, dat2) |> sim_to_pH_data()
    plot_pH_time(dat, stratify_by = "Group")

    fit <- fit_pH_curve(
        dat,
        model = kpd_mod(0.8, 0.8, 0.1, 0.8, 1, 0.1, 0.1, 0.1, 0.1, 0.01),
        estmethod = "saem"
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
            step = 0.1
        )
    )

    expect_true(all(pHdatMean$derivedDf$area_under_pH > 0))
    expect_true(all(pHdatMean$derivedDf$pH_min > 0))
    expect_true(all(pHdatMean$derivedDf$t_min > 0))
    expect_true(all(pHdatMean$derivedDf$time_under_ph > 0))

    expect_true(all(unique(pHdatMean$derivedDf$group) == c("A", "B")))
    expect_true(length(unique(pHdatMean$derivedDf$edk50)) == 2)

    suppressWarnings(
        pHdatFull <- pHMetrics_from_fit(
            fit,
            time_start = 0,
            time_end = 50,
            step = 0.1
        )
    )

    expect_true(any(pHdatFull$derivedDf$area_under_pH > 0))
    expect_true(any(pHdatFull$derivedDf$pH_min > 0))
    expect_true(any(pHdatFull$derivedDf$t_min > 0))
    expect_true(any(pHdatFull$derivedDf$time_under_ph > 0))
    expect_false(any(is.infinite(pHdatFull$derivedDf$area_under_pH)))
    expect_false(any(is.infinite(pHdatFull$derivedDf$pH_min)))
    expect_false(any(is.infinite(pHdatFull$derivedDf$t_min)))
    expect_false(any(is.infinite(pHdatFull$derivedDf$time_under_ph)))

    expect_true(all(!is.na(pHdatFull$derivedDf$edk50)))
    expect_true(all(!is.na(pHdatFull$derivedDf$kde)))
    expect_true(all(!is.na(pHdatFull$derivedDf$kd)))
    expect_true(all(!is.na(pHdatFull$derivedDf$ks)))

}


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
        dose_time = 1,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    ) |> sim_to_pH_data()
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
        estmethod = "focei"
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
            step = 0.1
        )
    )

    expect_true(pHdatMean$derivedDf$area_under_pH > 0)
    expect_true(pHdatMean$derivedDf$pH_min > 0)
    expect_true(pHdatMean$derivedDf$t_min > 0)
    expect_true(pHdatMean$derivedDf$time_under_ph > 0)
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
        dose_time = 1,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE
    ) |> sim_to_pH_data()
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
        estmethod = "focei"
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
            step = 0.1
        )
    )

    expect_true(any(pHdatFull$derivedDf$area_under_pH > 0))
    expect_true(any(pHdatFull$derivedDf$pH_min > 0))
    expect_true(any(pHdatFull$derivedDf$t_min > 0))
    expect_true(any(pHdatFull$derivedDf$time_under_ph > 0))
    expect_false(any(is.infinite(pHdatFull$derivedDf$area_under_pH)))
    expect_false(any(is.infinite(pHdatFull$derivedDf$pH_min)))
    expect_false(any(is.infinite(pHdatFull$derivedDf$t_min)))
    expect_false(any(is.infinite(pHdatFull$derivedDf$time_under_ph)))
})


test_that("nlme birkhed with gamma", {
    d <- read_pH(system.file(
        "extdata",
        "Birkhed1_copy.csv",
        package = "GatorpH"
    ), dose_time = 0.5)
    plot_pH_time(d, stratify_by = "Subject", showDosing = TRUE)

    # naive fit 
    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        estmethod = "bobyqa",
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
            include_gamma = TRUE
        )
    )
    nrow(pHdatMean$derivedDf) |> expect_equal(2)
    expect_false(any(is.na(pHdatMean$derivedDf$area_under_pH)))
    expect_false(any(is.na(pHdatMean$derivedDf$edk50)))

    ## nlme 
    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        estmethod = "saem",
        cov_params = "ks",
        cov_fixedeffects = "t.ks", 
        include_gamma = TRUE
    )
    pHdatnlme <- pHMetrics_from_fit(
        fit,
        ph_threshold = 6.5,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        include_gamma = TRUE
    )

    expect_false(any(is.na(pHdatnlme$derivedDf$area_under_pH)))
    expect_false(any(is.na(pHdatnlme$derivedDf$edk50)))

})

test_that("fejeskov", {
    d <- read_pH(system.file(
        "extdata",
        "Fejerskov Data_GatorpH Test.csv",
        package = "GatorpH"
    ), dose_time = 1)


    # test naive fit 
    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        estmethod = "uobyqa"
    )

    fit_individual_plot(fit)
    suppressWarnings(
        pHdatMean <- pHMetrics_from_fit(
            fit,
            ph_threshold = 5.5,
            time_start = 0,
            time_end = 50,
            step = 0.1
        )
    )
    nrow(pHdatMean$derivedDf) |> expect_equal(2)
    expect_false(any(is.na(pHdatMean$derivedDf$area_under_pH)))
    expect_false(any(is.na(pHdatMean$derivedDf$edk50)))

    ## nlme
    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        estmethod = "saem"
    )
    fit_individual_plot(fit)
    pHdatnlme <- pHMetrics_from_fit(
        fit,
        ph_threshold = 6.5,
        time_start = 0,
        time_end = 50,
        step = 0.1
    )
    expect_false(any(is.na(pHdatnlme$derivedDf$area_under_pH)))
    expect_false(any(is.na(pHdatnlme$derivedDf$edk50)))
})

test_that("remove_gamma", {
    mod <- kpd_mod()
    mod_no_gamma <- remove_gamma(mod)
    expect_false(any(grepl("gamma", names(mod_no_gamma))))

    
    dat <- simulate_steph_curve(
        kpd_mod(0.8, 0.8, 0.1, 0.8),
        nsub = 10,
        dose_time = 1,
        time = c(0, 10, 15, 20, 30),
        ignoreBSV = FALSE, 
        include_gamma= FALSE
    )  |> sim_to_pH_data()

    fit <- fit_pH_curve(
        dat,
        model = mod_no_gamma,
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
        include_gamma = FALSE
    )  |> expect_no_error()

    
    fit <- fit_pH_curve(
        dat,
        model = mod_no_gamma,
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
        estmethod = "bobyqa",
        cov_params = "ks",
        cov_fixedeffects = "t.ks",
        include_gamma = FALSE
    )

    fit$objDf
    pHdatnlme <- pHMetrics_from_fit(
        fit,
        ph_threshold = 6.5,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        include_gamma = FALSE
    )

    fit <- fit_pH_curve(
        d,
        model = kpd_mod(),
        estmethod = "saem",
        cov_params = c(),
        cov_fixedeffects = c(),
        include_gamma = FALSE
    )
    
    pHdatnlme <- pHMetrics_ fsrom_fit(
        fit,
        ph_threshold = 6.5,
        time_start = 0,
        time_end = 50,
        step = 0.1,
        include_gamma = FALSE
    )

    
})
