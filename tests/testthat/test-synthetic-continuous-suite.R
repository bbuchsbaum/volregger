.deg <- function(x) x * pi / 180

.estimate_with_ref1 <- function(run, control_overrides = list()) {
  ctl <- utils::modifyList(
    volreg_control(),
    list(
      max_iter = 8L,
      n_samples = 2500L,
      n_restarts = 2L,
      robust_template = FALSE,
      template_refresh = FALSE
    )
  )
  if (length(control_overrides) > 0) {
    ctl <- utils::modifyList(ctl, control_overrides)
  }

  estimate_motion(
    run,
    ref = "index",
    ref_index = 1,
    control = ctl
  )
}

.pose_error_stats <- function(est_par, truth_pose) {
  truth_correction <- -truth_pose
  err <- est_par - truth_correction
  list(
    trans_rmse = sqrt(mean(err[, 1:3]^2)),
    rot_rmse = sqrt(mean(err[, 4:6]^2)),
    trans_p95 = unname(stats::quantile(rowSums(abs(err[, 1:3])), 0.95)),
    rot_p95 = unname(stats::quantile(rowSums(abs(err[, 4:6])), 0.95))
  )
}

test_that("continuous null: zero-motion run does not induce large spurious pose", {
  subject <- make_continuous_subject(seed = 1001)
  nt <- 10
  poses <- matrix(0, nrow = nt, ncol = 6)

  sim <- simulate_continuous_run(
    subject = subject,
    nx = 16,
    ny = 16,
    nz = 8,
    nt = nt,
    poses = poses,
    noise_sd = 2.0,
    drift_amp = 0.04,
    bias_amp = 0.08,
    supersample = 1L,
    seed = 2001
  )

  est <- .estimate_with_ref1(sim$run)

  expect_true(stats::median(abs(est$par[, 1])) < 0.45)
  expect_true(stats::median(abs(est$par[, 2])) < 0.45)
  expect_true(stats::median(abs(est$par[, 3])) < 0.45)
  expect_true(stats::median(abs(est$par[, 4])) < .deg(2.0))
  expect_true(stats::median(abs(est$par[, 5])) < .deg(2.0))
  expect_true(stats::median(abs(est$par[, 6])) < .deg(2.0))
})

test_that("continuous rigid accuracy: smooth trajectory recovered within tolerance across noise", {
  subject <- make_continuous_subject(seed = 1002)
  nt <- 12
  tt <- seq(0, 1, length.out = nt)

  pose <- cbind(
    1.8 * sin(2 * pi * tt),
    0.9 * cos(2 * pi * tt),
    0.5 * sin(4 * pi * tt),
    .deg(1.8) * sin(2 * pi * tt),
    .deg(1.2) * cos(2 * pi * tt),
    .deg(1.5) * sin(3 * pi * tt)
  )
  pose <- sweep(pose, 2, pose[1, ], FUN = "-")

  for (nsd in c(1.0, 3.0)) {
    sim <- simulate_continuous_run(
      subject = subject,
      nx = 16,
      ny = 16,
      nz = 8,
      nt = nt,
      poses = pose,
      noise_sd = nsd,
      drift_amp = 0.02,
      bias_amp = 0.05,
      supersample = 1L,
      seed = 2100 + as.integer(10 * nsd)
    )

    est <- .estimate_with_ref1(sim$run)
    st <- .pose_error_stats(est$par, pose)

    expect_true(st$trans_rmse < 0.95, paste0("noise=", nsd, " trans_rmse=", signif(st$trans_rmse, 4)))
    expect_true(st$rot_rmse < .deg(3.0), paste0("noise=", nsd, " rot_rmse(rad)=", signif(st$rot_rmse, 4)))
    expect_true(st$trans_p95 < 3.0, paste0("noise=", nsd, " trans_p95=", signif(st$trans_p95, 4)))
  }
})

test_that("continuous capture-range: capture boost expands success region", {
  subject <- make_continuous_subject(seed = 1003)
  nx <- 16; ny <- 16; nz <- 8
  dx_grid <- c(-6, 0, 6)
  rz_grid <- .deg(c(-8, 0, 8))

  score_mode <- function(capture_boost) {
    ok <- 0L
    total <- 0L
    for (dx in dx_grid) {
      for (rz in rz_grid) {
        total <- total + 1L
        poses <- matrix(0, nrow = 2, ncol = 6)
        poses[2, 1] <- dx
        poses[2, 6] <- rz

        sim <- simulate_continuous_run(
          subject = subject,
          nx = nx,
          ny = ny,
          nz = nz,
          nt = 2,
          poses = poses,
          noise_sd = 1.0,
          drift_amp = 0,
          bias_amp = 0.04,
          supersample = 1L,
          seed = as.integer(3000 + 10 * dx + 100 * rz)
        )

        est <- .estimate_with_ref1(sim$run, control_overrides = list(n_restarts = 2L, capture_boost = capture_boost))
        ex <- est$par[2, 1] + dx
        er <- est$par[2, 6] + rz
        if (abs(ex) < 1.5 && abs(er) < .deg(4.0)) {
          ok <- ok + 1L
        }
      }
    }
    c(ok = ok, total = total)
  }

  s0 <- score_mode(FALSE)
  s2 <- score_mode(TRUE)

  expect_true(s2[["ok"]] >= s0[["ok"]], paste0("capture ok=", s2[["ok"]], " baseline ok=", s0[["ok"]]))
  expect_true(s2[["ok"]] >= 2)
})
