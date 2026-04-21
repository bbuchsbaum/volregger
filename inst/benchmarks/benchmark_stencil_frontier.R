#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("volregger", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("Need either installed volregger or devtools to load local package", call. = FALSE)
    }
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(volregger)
  }
})

make_phantom <- function(nx = 24, ny = 24, nz = 8) {
  gx <- array(rep(seq(-1, 1, length.out = nx), each = ny * nz), dim = c(nx, ny, nz))
  gy <- array(rep(rep(seq(-1, 1, length.out = ny), each = nz), times = nx), dim = c(nx, ny, nz))
  gz <- array(rep(seq(-1, 1, length.out = nz), times = nx * ny), dim = c(nx, ny, nz))

  a <- exp(-((gx - 0.25)^2 + (gy + 0.2)^2 + (gz)^2) / 0.12)
  b <- 0.7 * exp(-((gx + 0.45)^2 + (gy - 0.35)^2 + (gz + 0.1)^2) / 0.06)
  c <- 0.10 * gx - 0.05 * gy + 0.03 * gz
  120 * (a + b + c)
}

build_run <- function(base, nt = 20, noise_sd = 2, seed = 1L) {
  set.seed(seed)
  run <- array(rep(as.numeric(base), nt), dim = c(dim(base), nt))
  run <- run + array(rnorm(length(run), sd = noise_sd), dim = dim(run))

  truth_motion <- matrix(0, nrow = nt, ncol = 6)
  truth_motion[c(4, 8, 13, 17), 1] <- c(1.2, -1.8, 2.1, -1.4)
  truth_motion[c(6, 11, 15), 2] <- c(0.8, -0.9, 1.0)
  truth_motion[c(9, 18), 6] <- c(1.5, -1.2) * pi / 180

  moved <- apply_motion(run, truth_motion, interp = "lanczos")
  list(run = moved, truth_pose = truth_motion)
}

pose_stats <- function(est_par, truth_pose, iters) {
  truth_correction <- -truth_pose
  err <- est_par - truth_correction
  trans_norm <- sqrt(rowSums(err[, 1:3, drop = FALSE]^2))
  rot_norm <- sqrt(rowSums(err[, 4:6, drop = FALSE]^2))
  c(
    trans_rmse = sqrt(mean(err[, 1:3]^2)),
    rot_rmse = sqrt(mean(err[, 4:6]^2)),
    fail_rate = mean(trans_norm > 2.0 | rot_norm > (3 * pi / 180)),
    mean_iters = mean(iters)
  )
}

benchmark_frontier <- function() {
  base <- make_phantom()
  sim <- build_run(base, nt = 20, noise_sd = 2, seed = 42)
  eval_mask <- base > stats::quantile(base, probs = 0.35)

  methods <- list(
    list(name = "dense-9000", use_ic = FALSE, n_points = 9000L),
    list(name = "uniform-2500", use_ic = FALSE, n_points = 2500L),
    list(name = "stencil-500", use_ic = TRUE, n_points = 500L),
    list(name = "stencil-1000", use_ic = TRUE, n_points = 1000L),
    list(name = "stencil-2000", use_ic = TRUE, n_points = 2000L),
    list(name = "stencil-4000", use_ic = TRUE, n_points = 4000L)
  )

  rows <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    m <- methods[[i]]

    profile <- volreg_profile(
      if (isTRUE(m$use_ic)) "ic_stencil" else "dense_baseline",
      overrides = list(
        max_iter = 8L,
        n_restarts = 2L,
        n_samples = as.integer(m$n_points),
        robust_template = FALSE,
        template_refresh = FALSE,
        stencil_n = as.integer(m$n_points)
      )
    )
    ctl <- profile$control

    tm <- system.time({
      est <- estimate_motion(sim$run, ref = "index", ref_index = 1, control = ctl)
    })

    st <- pose_stats(est$par, sim$truth_pose, est$iters)
    truth_disp <- motion_displacement_summary(
      mask = eval_mask,
      motion = est$par,
      truth = invert_motion(sim$truth_pose)
    )

    rows[[i]] <- data.frame(
      method = m$name,
      profile = profile$name,
      components = paste(profile$components, collapse = "+"),
      use_ic_stencil = isTRUE(m$use_ic),
      n_points = as.integer(m$n_points),
      elapsed_sec = unname(tm[["elapsed"]]),
      trans_rmse = st[["trans_rmse"]],
      rot_rmse = st[["rot_rmse"]],
      mask_disp_median = stats::median(truth_disp$disp_median),
      mask_disp_p95 = stats::median(truth_disp$disp_p95),
      mask_disp_max = stats::median(truth_disp$disp_max),
      fd_error_median = stats::median(truth_disp$fd_error),
      fail_rate = st[["fail_rate"]],
      mean_iters = st[["mean_iters"]],
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

res <- benchmark_frontier()
out_csv <- file.path("inst", "benchmarks", "stencil_frontier_results.csv")
write.csv(res, out_csv, row.names = FALSE)

cat("Saved:", out_csv, "\n")
print(res)
