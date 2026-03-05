test_that("simulation single baseline", {
  dat <- simulate_steph_curve(
    kpd_mod(0.8, 0.8, ks = 4.9, kd = 0.7),
    nsub = 10,
    baseline = 7,
    time = c(0, 5, 10, 15, 20, 30),
    dose = 100,
    baseline_time = -5,
    ignoreBSV = FALSE
  )

  phvals <- dat |> dplyr::filter(time == 0) |> dplyr::pull(pH)
  all.equal(phvals, rep(7, 10), tolerance = 0.01) |> expect_true()
})


test_that("simulation multiple baseline", {
  skip_on_cran()
  baseline_vec <- rnorm(10, mean = 7, sd = 0.1)
  dat <- simulate_steph_curve(
    kpd_mod(0.8, 0.8, ks = 4.9, kd = 0.7),
    nsub = 10,
    baseline = baseline_vec,
    time = c(0, 5, 10, 15, 20, 30),
    dose = 100,
    baseline_time = -5,
    ignoreBSV = FALSE
  )
  plot_pH_time(dat)

  phvals <- dat |> dplyr::filter(time == 0) |> dplyr::pull(pH)
  all.equal(phvals, baseline_vec) |> expect_true()
})
