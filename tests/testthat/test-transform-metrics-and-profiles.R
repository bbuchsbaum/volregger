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

test_that("fast_native profile targets the native speed frontier", {
  prof <- volreg_profile("fast_native")

  expect_s3_class(prof, "volreg_profile")
  expect_equal(prof$engine, "rigid_robust")
  expect_true(all(c("ic_stencil", "whiten", "fast_native") %in% prof$components))
  expect_true(isTRUE(prof$control$use_ic_stencil))
  expect_true(isTRUE(prof$control$whiten_modes))
  expect_false(isTRUE(prof$control$use_nuisance_marginalization))
  expect_equal(prof$control$stencil_n, 1000L)
  expect_equal(prof$control$capture_topk, 4L)
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

test_that("real-data benchmark truth metrics expose a broad scorecard", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  mask <- array(TRUE, dim = c(3, 3, 3))
  base <- array(seq_len(prod(dim(mask))), dim = dim(mask))
  clean <- array(NA_real_, dim = c(dim(mask), 4))
  for (i in seq_len(dim(clean)[4])) {
    clean[, , , i] <- base + i
  }
  attr(clean, "pixdim") <- c(2, 2, 2)

  truth <- array(rep(diag(4), dim(clean)[4]), dim = c(4, 4, dim(clean)[4]))
  metrics <- env$compute_benchmark_metrics(
    input_run = clean,
    corrected = clean,
    mask = mask,
    motion_mats = truth,
    truth = truth,
    clean = clean
  )

  expected <- c(
    "disp_p95",
    "fd_error",
    "trans_abs_p95",
    "rot_abs_deg_p95",
    "rel_frob_p95",
    "rel_rot_deg_p95",
    "mae_to_truth",
    "nrmse_to_truth",
    "residual_dvars_to_truth",
    "temporal_corr_to_truth",
    "edge_ratio_to_truth"
  )

  expect_true(all(expected %in% names(metrics)))
  expect_true(all(env$benchmark_metric_names() %in% names(env$empty_benchmark_metrics())))
  expect_equal(metrics$mae_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$nrmse_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$residual_dvars_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$trans_abs_p95, 0, tolerance = 1e-12)
  expect_equal(metrics$rot_abs_deg_p95, 0, tolerance = 1e-12)
  expect_equal(metrics$rel_frob_p95, 0, tolerance = 1e-12)
  expect_equal(metrics$temporal_corr_to_truth, 1, tolerance = 1e-12)
  expect_equal(metrics$edge_ratio_to_truth, 1, tolerance = 1e-12)
})
