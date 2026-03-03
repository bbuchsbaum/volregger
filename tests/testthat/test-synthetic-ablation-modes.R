.deg <- function(x) x * pi / 180

.pose_err <- function(est, truth_pose, idx = seq_len(nrow(est$par))) {
  err <- est$par[idx, , drop = FALSE] + truth_pose[idx, , drop = FALSE]
  c(
    trans_rmse = sqrt(mean(err[, 1:3]^2)),
    rot_rmse = sqrt(mean(err[, 4:6]^2)),
    fail_rate = mean(est$flags[idx] != 0)
  )
}

test_that("ablation: robust template is safer on bad-reference corruption", {
  seeds <- 1:3
  trans_ratio <- numeric(length(seeds))
  rot_ratio <- numeric(length(seeds))

  for (ii in seq_along(seeds)) {
    s <- seeds[ii]
    subject <- make_continuous_subject(seed = 1200 + s)
    nt <- 11
    tt <- seq(0, 1, length.out = nt)
    poses <- cbind(
      1.2 * sin(2 * pi * tt),
      0.8 * cos(2 * pi * tt),
      0.4 * sin(3 * pi * tt),
      .deg(1.4) * sin(2 * pi * tt),
      .deg(1.0) * cos(2 * pi * tt),
      .deg(1.2) * sin(3 * pi * tt)
    )
    poses <- sweep(poses, 2, poses[1, ], FUN = "-")

    sim <- simulate_continuous_run(
      subject = subject,
      nx = 16,
      ny = 16,
      nz = 8,
      nt = nt,
      poses = poses,
      noise_sd = 2.0,
      drift_amp = 0.02,
      bias_amp = 0.05,
      supersample = 1L,
      seed = 5000 + s
    )

    run <- sim$run
    mid <- (nt + 1) %/% 2
    run[4:10, 5:12, 3:6, mid] <- 0.15 * run[4:10, 5:12, 3:6, mid]
    run[2:15, , , mid] <- run[2:15, , , mid] + 0.20 * run[1:14, , , mid]

    ctl_fixed <- utils::modifyList(
      volreg_control(),
      list(max_iter = 8L, n_samples = 2500L, n_restarts = 2L, robust_template = FALSE, template_refresh = FALSE)
    )
    ctl_robust <- utils::modifyList(
      volreg_control(),
      list(max_iter = 8L, n_samples = 2500L, n_restarts = 2L, robust_template = TRUE, template_refresh = TRUE)
    )

    e_fixed <- estimate_motion(run, ref = "middle", control = ctl_fixed)
    e_robust <- estimate_motion(run, ref = "middle", control = ctl_robust)

    idx <- setdiff(seq_len(nt), mid)
    st_fixed <- .pose_err(e_fixed, poses, idx)
    st_rob <- .pose_err(e_robust, poses, idx)

    trans_ratio[ii] <- st_rob[["trans_rmse"]] / st_fixed[["trans_rmse"]]
    rot_ratio[ii] <- st_rob[["rot_rmse"]] / st_fixed[["rot_rmse"]]
  }

  expect_true(median(rot_ratio) < 1.0, paste0("median rot ratio=", signif(median(rot_ratio), 4)))
  expect_true(max(trans_ratio) < 1.35, paste0("max trans ratio=", signif(max(trans_ratio), 4)))
})

test_that("ablation: overlap taper preserves overlap under partial-FOV motion", {
  seeds <- 1:3
  overlap_diff <- numeric(length(seeds))

  for (ii in seq_along(seeds)) {
    s <- seeds[ii]
    subject <- make_continuous_subject(seed = 1300 + s)
    nt <- 10
    tt <- seq(0, 1, length.out = nt)

    poses <- cbind(
      6.0 * sin(2 * pi * tt),
      0.3 * cos(2 * pi * tt),
      6.0 * sin(2 * pi * tt),
      .deg(6) * sin(2 * pi * tt),
      .deg(1) * cos(2 * pi * tt),
      .deg(2) * sin(3 * pi * tt)
    )
    poses <- sweep(poses, 2, poses[1, ], FUN = "-")

    sim <- simulate_continuous_run(
      subject = subject,
      nx = 14,
      ny = 14,
      nz = 7,
      nt = nt,
      poses = poses,
      noise_sd = 2.0,
      drift_amp = 0,
      bias_amp = 0.06,
      supersample = 1L,
      seed = 7400 + s
    )

    ctl_taper <- utils::modifyList(
      volreg_control(),
      list(max_iter = 6L, n_samples = 1800L, n_restarts = 2L, use_overlap_taper = TRUE)
    )
    ctl_notaper <- utils::modifyList(
      volreg_control(),
      list(max_iter = 6L, n_samples = 1800L, n_restarts = 2L, use_overlap_taper = FALSE)
    )

    e_taper <- estimate_motion(sim$run, ref = "index", ref_index = 1, control = ctl_taper)
    e_notaper <- estimate_motion(sim$run, ref = "index", ref_index = 1, control = ctl_notaper)

    overlap_diff[ii] <- mean(e_taper$overlap) - mean(e_notaper$overlap)
  }

  expect_true(median(overlap_diff) > 0.15, paste0("median overlap diff=", signif(median(overlap_diff), 4)))
})

test_that("ablation: packet-aware resampling reduces intra-TR residual mismatch", {
  ratios <- numeric(3)

  for (s in seq_len(3)) {
    subject <- make_continuous_subject(seed = 1500 + s)
    nt <- 12
    nx <- 16; ny <- 16; nz <- 8
    tt <- seq(0, 1, length.out = nt)

    poses <- cbind(
      1.2 * sin(2 * pi * tt),
      0.6 * cos(2 * pi * tt),
      0.3 * sin(3 * pi * tt),
      .deg(1.2) * sin(2 * pi * tt),
      .deg(0.9) * cos(2 * pi * tt),
      .deg(1.1) * sin(3 * pi * tt)
    )
    poses <- sweep(poses, 2, poses[1, ], FUN = "-")

    slopes <- matrix(0, nrow = nt, ncol = 6)
    for (t in seq_len(nt - 1)) {
      slopes[t, ] <- poses[t + 1, ] - poses[t, ]
    }
    slopes[nt, ] <- slopes[nt - 1, ]

    slice_times <- seq(0, 1, length.out = nz)

    sim <- simulate_continuous_run(
      subject = subject,
      nx = nx,
      ny = ny,
      nz = nz,
      nt = nt,
      poses = poses,
      pose_slopes = slopes,
      slice_times = slice_times,
      noise_sd = 0.8,
      drift_amp = 0,
      bias_amp = 0.02,
      supersample = 1L,
      seed = 9000 + s
    )

    truth_corr <- -poses
    y_vol <- apply_motion(sim$run, truth_corr, interp = "lanczos")
    y_pkt <- apply_motion(sim$run, truth_corr, interp = "lanczos", slice_times = slice_times)

    sim0 <- simulate_continuous_run(
      subject = subject,
      nx = nx,
      ny = ny,
      nz = nz,
      nt = nt,
      poses = matrix(0, nrow = nt, ncol = 6),
      noise_sd = 0,
      drift_amp = 0,
      bias_amp = 0,
      supersample = 2L,
      seed = 9500 + s
    )

    mse_vol <- mean((y_vol - sim0$run)^2)
    mse_pkt <- mean((y_pkt - sim0$run)^2)
    ratios[s] <- mse_pkt / mse_vol
  }

  expect_true(median(ratios) < 0.90, paste0("median packet/volume mse ratio=", signif(median(ratios), 4)))
})
