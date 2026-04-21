test_that("motion transform helpers recover exact translation displacement", {
  mask <- array(TRUE, dim = c(3, 3, 3))
  est <- matrix(c(1, 0, 0, 0, 0, 0), nrow = 1)
  truth <- array(diag(4), dim = c(4, 4, 1))

  disp <- displacements_within_mask(mask, est[1, ], reference = truth[, , 1], pixdim = c(2, 2, 2))
  summ <- motion_displacement_summary(mask, est, truth, pixdim = c(2, 2, 2))

  expect_true(all(abs(disp - 1) < 1e-8))
  expect_equal(summ$disp_median, 1, tolerance = 1e-8)
  expect_equal(summ$disp_p95, 1, tolerance = 1e-8)
  expect_equal(summ$disp_max, 1, tolerance = 1e-8)
  expect_equal(summ$fd_est, 1, tolerance = 1e-8)
  expect_equal(summ$fd_truth, 0, tolerance = 1e-8)
  expect_equal(summ$fd_error, 1, tolerance = 1e-8)
})

test_that("invert_motion returns the inverse homogeneous transform", {
  motion <- matrix(c(1.2, -0.5, 0.8, 0.1, -0.05, 0.08), nrow = 1)
  mats <- motion_matrices(motion)
  inv_mats <- invert_motion(motion)

  ident <- mats[, , 1] %*% inv_mats[, , 1]
  expect_equal(ident, diag(4), tolerance = 1e-10)
})

test_that("volreg_profile exposes named benchmark configurations", {
  prof <- volreg_profile(
    "ic_nuisance_full",
    overrides = list(n_samples = 1800L, stencil_n = 1200L)
  )

  expect_s3_class(prof, "volreg_profile")
  expect_equal(prof$name, "ic_nuisance_full")
  expect_equal(prof$engine, "rigid_robust")
  expect_true("nuisance" %in% prof$components)
  expect_true("template_basis" %in% prof$components)
  expect_true(isTRUE(prof$control$use_ic_stencil))
  expect_true(isTRUE(prof$control$whiten_modes))
  expect_true(isTRUE(prof$control$use_nuisance_marginalization))
  expect_true(isTRUE(prof$control$nuisance_include_template))
  expect_equal(prof$control$n_samples, 1800L)
  expect_equal(prof$control$stencil_n, 1200L)
})

test_that("fit summaries and tabular output expose richer QC fields", {
  set.seed(200)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 5), dim = c(8, 8, 4, 5))

  fit <- volreg(
    x,
    resample = FALSE,
    qc = TRUE,
    control = utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 1800L))
  )

  sm <- summary(fit)
  tab <- as.data.frame(fit)

  expect_s3_class(sm, "summary.volregger_fit")
  expect_equal(sm$n_frames, 5)
  expect_true("fd_mean" %in% names(sm))
  expect_true("fit_failures" %in% names(sm))
  expect_true(all(c("cost_drop", "motion_spike", "fit_failure", "censor_suggest") %in% names(tab)))
})
