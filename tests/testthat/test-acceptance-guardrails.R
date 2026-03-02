test_that("volreg is deterministic for fixed input", {
  set.seed(101)
  x <- array(rnorm(16 * 16 * 8 * 6), dim = c(16, 16, 8, 6))
  ctl <- utils::modifyList(volreg_control(), list(max_iter = 6L, n_samples = 3000L, n_restarts = 1L))

  fit1 <- volreg(x, resample = FALSE, qc = TRUE, control = ctl)
  fit2 <- volreg(x, resample = FALSE, qc = TRUE, control = ctl)

  expect_equal(fit1$motion, fit2$motion, tolerance = 1e-12)
  expect_equal(fit1$qc$fd, fit2$qc$fd, tolerance = 1e-12)
  expect_equal(fit1$qc$dvars, fit2$qc$dvars, tolerance = 1e-12)
})

test_that("volreg runtime stays within baseline envelope on small fixture", {
  skip_on_cran()

  set.seed(102)
  x <- array(rnorm(24 * 24 * 10 * 12), dim = c(24, 24, 10, 12))
  ctl <- utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 2500L))

  elapsed <- system.time({
    fit <- volreg(x, resample = TRUE, qc = TRUE, control = ctl)
    expect_equal(dim(fit$motion), c(12, 6))
  })[["elapsed"]]

  # Generous guardrail for local CPU variability.
  expect_lt(unname(elapsed), 10)
})
