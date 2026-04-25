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

test_that("16 x 6 matrix representation is parsed as six transforms", {
  mats <- array(rep(diag(4), 6), dim = c(4, 4, 6))
  mats[1, 4, ] <- seq_len(6)
  flat <- matrix(as.numeric(mats), nrow = 16, ncol = 6)

  parsed <- motion_matrices(flat)

  expect_equal(dim(parsed), c(4L, 4L, 6L))
  expect_equal(parsed, mats)
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

test_that("fast_native_pyramid profile enables bounded coarse-to-fine fitting", {
  prof <- volreg_profile("fast_native_pyramid")

  expect_s3_class(prof, "volreg_profile")
  expect_true("pyramid" %in% prof$components)
  expect_true(isTRUE(prof$control$use_pyramid))
  expect_equal(prof$control$downsample, c(4L, 2L, 1L))
  expect_equal(prof$control$max_iter, c(4L, 3L, 4L))
  expect_equal(prof$control$n_samples, c(1000L, 2000L, 8000L))
})

test_that("fast_native_parallel profile enables opt-in frame parallelism", {
  prof <- volreg_profile("fast_native_parallel")

  expect_s3_class(prof, "volreg_profile")
  expect_true("parallel_frames" %in% prof$components)
  expect_true(isTRUE(prof$control$parallel_frames))
  expect_true(isTRUE(prof$control$use_ic_stencil))
  expect_true(isTRUE(prof$control$whiten_modes))
})

test_that("bare high-level APIs use fast_fmri while dense baseline stays opt-in", {
  set.seed(201)
  x0 <- array(rnorm(6 * 6 * 4), dim = c(6, 6, 4))
  x <- array(rep(as.numeric(x0), 4), dim = c(6, 6, 4, 4))

  fit_default <- volreg(x, resample = FALSE, qc = FALSE)
  expect_equal(fit_default$control$profile, "fast_fmri")
  expect_true(isTRUE(fit_default$control$use_ic_stencil))
  expect_true(isTRUE(fit_default$control$whiten_modes))
  expect_equal(fit_default$control$edge_exclude_frac, 0.05)
  expect_true(isTRUE(fit_default$control$template_refresh_valid_only))

  ctl <- utils::modifyList(volreg_control(), list(max_iter = 1L, n_samples = 500L))
  fit_manual <- volreg(x, resample = FALSE, qc = FALSE, control = ctl)
  expect_null(fit_manual$control$profile)
  expect_true(isTRUE(fit_manual$control$use_ic_stencil))

  fit_dense <- volreg(x, resample = FALSE, qc = FALSE, profile = "dense_baseline")
  expect_equal(fit_dense$control$profile, "dense_baseline")
  expect_false(isTRUE(fit_dense$control$use_ic_stencil))
  expect_false(isTRUE(fit_dense$control$whiten_modes))

  est_default <- estimate_motion(x)
  expect_equal(est_default$warm_start$passes, 2L)
})

test_that("fast_fmri profile adds fMRI-specific edge and template safeguards", {
  prof <- volreg_profile("fast_fmri")

  expect_s3_class(prof, "volreg_profile")
  expect_equal(prof$engine, "rigid_robust")
  expect_true(all(c("fast_native", "edge_exclude", "valid_template_refresh", "parallel_frames") %in% prof$components))
  expect_true(isTRUE(prof$control$use_ic_stencil))
  expect_true(isTRUE(prof$control$whiten_modes))
  expect_equal(prof$control$stencil_n, 1000L)
  expect_equal(prof$control$edge_exclude_frac, 0.05)
  expect_true(isTRUE(prof$control$template_refresh_valid_only))
  expect_true(isTRUE(prof$control$parallel_frames))
})

test_that("IC-polish profile keeps IC setup and enables bounded analytic polish", {
  prof <- volreg_profile("ic_nuisance_polish")

  expect_s3_class(prof, "volreg_profile")
  expect_true(all(c("ic_stencil", "whiten", "nuisance", "adaptive_huber", "parallel_frames", "ic_polish") %in% prof$components))
  expect_true(isTRUE(prof$control$use_ic_stencil))
  expect_true(isTRUE(prof$control$adaptive_huber))
  expect_equal(prof$control$huber_k, 0.15)
  expect_equal(prof$control$max_iter, 12L)
  expect_true(isTRUE(prof$control$parallel_frames))
  expect_true(isTRUE(prof$control$use_ic_polish))
  expect_false(isTRUE(prof$control$use_dense_polish))
  expect_equal(prof$control$polish_iter, 3L)
  expect_equal(prof$control$polish_n_samples, 2000L)
})

test_that("adaptive-Huber profile keeps IC setup and scales robust threshold", {
  prof <- volreg_profile("ic_nuisance_adaptive")

  expect_s3_class(prof, "volreg_profile")
  expect_true(all(c("ic_stencil", "whiten", "nuisance", "adaptive_huber") %in% prof$components))
  expect_true(isTRUE(prof$control$use_ic_stencil))
  expect_true(isTRUE(prof$control$adaptive_huber))
  expect_false(isTRUE(prof$control$use_dense_polish))
})

test_that("tuned adaptive-Huber profile uses a conservative residual-scale multiplier", {
  prof <- volreg_profile("ic_nuisance_adaptive_tuned")

  expect_s3_class(prof, "volreg_profile")
  expect_true(all(c("adaptive_huber", "huber_0.15", "iter_12", "parallel_frames") %in% prof$components))
  expect_true(isTRUE(prof$control$adaptive_huber))
  expect_true(isTRUE(prof$control$template_refresh_valid_only))
  expect_equal(prof$control$huber_k, 0.15)
  expect_equal(prof$control$max_iter, 12L)
  expect_true(isTRUE(prof$control$parallel_frames))
})

test_that("guarded adaptive-Huber profile enables common-objective candidate selection", {
  prof <- volreg_profile("ic_nuisance_adaptive_guarded")

  expect_s3_class(prof, "volreg_profile")
  expect_true(all(c("adaptive_huber", "huber_0.5", "guarded") %in% prof$components))
  expect_true(isTRUE(prof$control$adaptive_huber))
  expect_true(isTRUE(prof$control$adaptive_huber_guard))
  expect_equal(prof$control$huber_k, 0.5)
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

test_that("real-data benchmark thread resolver supports auto and explicit counts", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  expect_equal(env$resolve_benchmark_nthreads("1"), 1L)
  expect_equal(env$resolve_benchmark_nthreads("2"), 2L)
  expect_true(env$resolve_benchmark_nthreads("auto") >= 1L)
  expect_true(env$resolve_benchmark_nthreads("auto") <= 4L)
  expect_error(env$resolve_benchmark_nthreads("0"), "positive integer")
})

test_that("real-data benchmark scalar metadata follows estimator last-value convention", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  expect_equal(env$benchmark_scalar_last(c(30L, 20L, 12L), NA_integer_), 12L)
  expect_equal(env$benchmark_scalar_last(numeric(), 1), 1)
  expect_equal(env$benchmark_scalar_last(NULL, "default"), "default")
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
    "residual_dvars_ratio_to_truth",
    "temporal_corr_to_truth",
    "fd_abs_corr_to_truth",
    "edge_ratio_to_truth",
    "edge_error_to_truth",
    "highfreq_residual_ratio",
    "nrmse_inner_to_truth",
    "nrmse_boundary_to_truth",
    "pose_jitter_error",
    "truth_image_fidelity_loss"
  )

  expect_true(all(expected %in% names(metrics)))
  expect_true(all(env$benchmark_metric_names() %in% names(env$empty_benchmark_metrics())))
  expect_equal(metrics$mae_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$nrmse_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$residual_dvars_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$residual_dvars_ratio_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$trans_abs_p95, 0, tolerance = 1e-12)
  expect_equal(metrics$rot_abs_deg_p95, 0, tolerance = 1e-12)
  expect_equal(metrics$rel_frob_p95, 0, tolerance = 1e-12)
  expect_equal(metrics$pose_jitter_error, 0, tolerance = 1e-12)
  expect_equal(metrics$temporal_corr_to_truth, 1, tolerance = 1e-12)
  expect_equal(metrics$edge_ratio_to_truth, 1, tolerance = 1e-12)
  expect_equal(metrics$edge_error_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$highfreq_residual_ratio, 0, tolerance = 1e-12)
  expect_equal(metrics$nrmse_inner_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$nrmse_boundary_to_truth, 0, tolerance = 1e-12)
  expect_equal(metrics$truth_image_fidelity_loss, 0, tolerance = 1e-12)
})

test_that("real-data benchmark FD correlation uses frame-to-frame displacement", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  truth_pose <- cbind(c(0, 1.0, 2.5, 2.0, 3.5), 0, 0, 0, 0, 0)
  est_pose <- cbind(c(0, 1.1, 2.3, 2.1, 3.2), 0, 0, 0, 0, 0)
  truth <- motion_matrices(truth_pose)
  est <- motion_matrices(est_pose)

  dyn <- env$motion_dynamics_error_summary(est, truth)
  expected_step_corr <- stats::cor(
    env$framewise_disp_from_transform(est),
    env$framewise_disp_from_transform(truth)
  )
  expected_abs_corr <- stats::cor(fd_from_transform(est), fd_from_transform(truth))

  expect_equal(dyn$fd_corr_to_truth, expected_step_corr, tolerance = 1e-12)
  expect_equal(dyn$fd_abs_corr_to_truth, expected_abs_corr, tolerance = 1e-12)
})

test_that("real-data benchmark scorecard ranks truth metrics by direction", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  metric_names <- env$benchmark_metric_names()
  make_row <- function(method, elapsed_sec, disp_p95, corr_to_truth, edge_ratio_to_truth) {
    metric_vals <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    metric_vals[["disp_p95"]] <- disp_p95
    metric_vals[["corr_to_truth"]] <- corr_to_truth
    metric_vals[["edge_ratio_to_truth"]] <- edge_ratio_to_truth
    cbind(
      data.frame(
        method = method,
        dataset_mode = "truth_injected",
        scenario = "low_motion",
        method_model = method,
        elapsed_sec = elapsed_sec,
        estimate_sec = elapsed_sec,
        apply_sec = NA_real_,
        report_sec = 0,
        peak_rss_mb = NA_real_,
        temp_bytes_written = NA_real_,
        status = "ok",
        error = NA_character_,
        stringsAsFactors = FALSE
      ),
      as.data.frame(as.list(metric_vals), stringsAsFactors = FALSE)
    )
  }

  raw <- rbind(
    make_row("volregger", elapsed_sec = 2, disp_p95 = 0.5, corr_to_truth = 0.99, edge_ratio_to_truth = 1.02),
    make_row("nifreeze", elapsed_sec = 4, disp_p95 = 1.5, corr_to_truth = 0.95, edge_ratio_to_truth = 0.70)
  )
  scorecard <- env$make_benchmark_scorecard(raw)
  summary <- env$summarize_benchmark_scorecard(scorecard)
  claim_summary <- env$summarize_benchmark_scorecard(scorecard, claim_only = TRUE)

  expect_true(all(c("elapsed_sec", "disp_p95", "corr_to_truth", "edge_ratio_to_truth") %in% scorecard$metric))
  expect_true(all(c("role", "claim_metric") %in% names(scorecard)))
  expect_true(scorecard$claim_metric[scorecard$metric == "disp_p95"])
  expect_false(scorecard$claim_metric[scorecard$metric == "tsnr_ratio"])
  expect_equal(scorecard$role[scorecard$metric == "tsnr_ratio"], "diagnostic_smoothing")
  expect_identical(scorecard$comparison[scorecard$metric == "disp_p95"], "win")
  expect_identical(scorecard$comparison[scorecard$metric == "corr_to_truth"], "win")
  expect_identical(scorecard$comparison[scorecard$metric == "edge_ratio_to_truth"], "win")
  expect_gt(sum(summary$wins), 0)
  expect_gt(sum(claim_summary$wins), 0)
  expect_equal(sum(summary$losses), 0)
  expect_equal(sum(claim_summary$losses), 0)
})

test_that("real-data benchmark no-regression gate includes protected diagnostics", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  metric_names <- env$benchmark_metric_names()
  make_row <- function(method, fd_corr_to_truth, tsnr_ratio) {
    metric_vals <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    metric_vals[["disp_p95"]] <- if (identical(method, "volregger")) 0.8 else 1.0
    metric_vals[["fd_error"]] <- if (identical(method, "volregger")) 0.6 else 0.9
    metric_vals[["fd_corr_to_truth"]] <- fd_corr_to_truth
    metric_vals[["pose_jitter_error"]] <- if (identical(method, "volregger")) 0.02 else 0.04
    metric_vals[["nrmse_to_truth"]] <- if (identical(method, "volregger")) 0.08 else 0.10
    metric_vals[["edge_error_to_truth"]] <- if (identical(method, "volregger")) 0.04 else 0.08
    metric_vals[["truth_image_fidelity_loss"]] <- if (identical(method, "volregger")) 0.12 else 0.20
    metric_vals[["tsnr_ratio"]] <- tsnr_ratio
    cbind(
      data.frame(
        method = method,
        dataset_mode = "truth_injected",
        scenario = "low_motion",
        method_model = method,
        elapsed_sec = NA_real_,
        estimate_sec = NA_real_,
        apply_sec = NA_real_,
        report_sec = NA_real_,
        peak_rss_mb = NA_real_,
        temp_bytes_written = NA_real_,
        status = "ok",
        error = NA_character_,
        stringsAsFactors = FALSE
      ),
      as.data.frame(as.list(metric_vals), stringsAsFactors = FALSE)
    )
  }

  raw <- rbind(
    make_row("volregger", fd_corr_to_truth = 0.97, tsnr_ratio = 3.0),
    make_row("afni_3dvolreg", fd_corr_to_truth = 0.98, tsnr_ratio = 4.0)
  )
  gate <- env$make_no_regression_gate(raw)
  summary <- env$summarize_no_regression_gate(gate)

  expect_equal(
    sort(gate$metric),
    sort(c(
      "disp_p95",
      "fd_error",
      "fd_corr_to_truth",
      "pose_jitter_error",
      "nrmse_to_truth",
      "edge_error_to_truth",
      "truth_image_fidelity_loss",
      "tsnr_ratio"
    ))
  )
  expect_identical(gate$comparison[gate$metric == "fd_corr_to_truth"], "tie")
  expect_identical(gate$comparison[gate$metric == "tsnr_ratio"], "loss")
  expect_equal(summary$regressions, 1)
})

test_that("real-data benchmark tSNR context exposes fidelity tradeoffs", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  metric_names <- env$benchmark_metric_names()
  make_row <- function(method, tsnr_ratio, fidelity_loss, edge_error, disp_p95) {
    metric_vals <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    metric_vals[["tsnr_ratio"]] <- tsnr_ratio
    metric_vals[["truth_image_fidelity_loss"]] <- fidelity_loss
    metric_vals[["edge_error_to_truth"]] <- edge_error
    metric_vals[["disp_p95"]] <- disp_p95
    metric_vals[["fd_error"]] <- if (identical(method, "volregger")) 0.8 else 2.0
    metric_vals[["fd_corr_to_truth"]] <- if (identical(method, "volregger")) 0.98 else 0.99
    metric_vals[["pose_jitter_error"]] <- if (identical(method, "volregger")) 0.02 else 0.05
    metric_vals[["nrmse_to_truth"]] <- if (identical(method, "volregger")) 0.09 else 0.20
    metric_vals[["blur_proxy_ratio"]] <- if (identical(method, "volregger")) 1.05 else 1.00
    cbind(
      data.frame(
        method = method,
        dataset_mode = "truth_injected",
        scenario = "low_motion",
        method_model = method,
        elapsed_sec = NA_real_,
        estimate_sec = NA_real_,
        apply_sec = NA_real_,
        report_sec = NA_real_,
        peak_rss_mb = NA_real_,
        temp_bytes_written = NA_real_,
        status = "ok",
        error = NA_character_,
        stringsAsFactors = FALSE
      ),
      as.data.frame(as.list(metric_vals), stringsAsFactors = FALSE)
    )
  }

  raw <- rbind(
    make_row("volregger", tsnr_ratio = 3.5, fidelity_loss = 0.15, edge_error = 0.06, disp_p95 = 1.0),
    make_row("afni_3dvolreg", tsnr_ratio = 3.7, fidelity_loss = 0.30, edge_error = 0.10, disp_p95 = 2.0)
  )
  ctx <- env$make_tsnr_context(raw)

  expect_true(ctx$tsnr_loss)
  expect_true(ctx$truth_fidelity_win)
  expect_true(ctx$edge_win)
  expect_true(ctx$displacement_win)
  expect_equal(ctx$tsnr_delta, -0.2)
  expect_equal(ctx$fidelity_loss_margin, 0.15)
})

test_that("interpolation guardrail summary reports deltas from baseline", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  raw <- data.frame(
    interp = c("lanczos", "linear", "lanczos", "linear"),
    scenario = rep(c("low_motion", "hard_motion_plus_nuisance"), each = 2),
    dataset_mode = "truth_injected",
    method = "volregger",
    status = "ok",
    elapsed_sec = c(2.0, 1.0, 3.0, 1.5),
    apply_sec = c(0.8, 0.4, 1.2, 0.6),
    disp_p95 = c(0.9, 1.1, 3.0, 3.4),
    truth_image_fidelity_loss = c(0.15, 0.25, 0.5, 0.8),
    stringsAsFactors = FALSE
  )

  summary <- env$summarize_interpolation_guardrail(raw, baseline_interp = "lanczos")

  low_linear <- summary[summary$scenario == "low_motion" & summary$interp == "linear", ]
  expect_equal(low_linear$apply_sec_vs_lanczos, -0.4)
  expect_equal(low_linear$disp_p95_vs_lanczos, 0.2)
  expect_equal(low_linear$truth_image_fidelity_loss_vs_lanczos, 0.1)
  expect_true(file.exists(testthat::test_path("..", "..", "inst", "benchmarks", "benchmark_interpolation_guardrail.R")))
})

test_that("real-data benchmark claim summary spans multiple truth scenarios", {
  env <- new.env(parent = globalenv())
  source(testthat::test_path("..", "..", "inst", "benchmarks", "realdata_benchmark_lib.R"), local = env)

  metric_names <- env$benchmark_metric_names()
  make_row <- function(method, scenario, elapsed_sec, disp_p95, mae_to_truth, tsnr_ratio) {
    metric_vals <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    metric_vals[["disp_p95"]] <- disp_p95
    metric_vals[["mae_to_truth"]] <- mae_to_truth
    metric_vals[["tsnr_ratio"]] <- tsnr_ratio
    cbind(
      data.frame(
        method = method,
        dataset_mode = "truth_injected",
        scenario = scenario,
        method_model = method,
        elapsed_sec = elapsed_sec,
        estimate_sec = elapsed_sec,
        apply_sec = NA_real_,
        report_sec = 0,
        peak_rss_mb = NA_real_,
        temp_bytes_written = NA_real_,
        status = "ok",
        error = NA_character_,
        stringsAsFactors = FALSE
      ),
      as.data.frame(as.list(metric_vals), stringsAsFactors = FALSE)
    )
  }

  raw <- rbind(
    make_row("volregger", "low_motion", 2, 0.5, 10, 2.0),
    make_row("nifreeze", "low_motion", 4, 1.5, 30, 1.2),
    make_row("volregger", "hard_motion_plus_nuisance", 3, 0.8, 12, 2.5),
    make_row("nifreeze", "hard_motion_plus_nuisance", 5, 3.0, 40, 1.4)
  )

  claim_summary <- env$summarize_benchmark_scorecard(
    env$make_benchmark_scorecard(raw),
    claim_only = TRUE
  )

  expect_equal(
    sort(unique(claim_summary$scenario)),
    sort(c("low_motion", "hard_motion_plus_nuisance"))
  )
  expect_equal(sum(claim_summary$losses), 0)
  expect_gt(sum(claim_summary$wins), 0)
})
