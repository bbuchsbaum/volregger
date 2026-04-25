.deg <- function(x) x * pi / 180

.run_ic_whiten_fixture <- function(seed, whiten_modes) {
  subject <- make_continuous_subject(seed = 3300 + seed)
  nt <- 10
  tt <- seq(0, 1, length.out = nt)

  poses <- cbind(
    2.5 * sin(2 * pi * tt),
    1.4 * cos(2 * pi * tt),
    1.0 * sin(3 * pi * tt),
    .deg(4.5) * sin(2 * pi * tt),
    .deg(4.0) * cos(2 * pi * tt),
    .deg(4.2) * sin(3 * pi * tt)
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
    bias_amp = 0.06,
    supersample = 1L,
    seed = 4300 + seed
  )

  ctl <- utils::modifyList(
    volreg_control(),
    list(
      max_iter = 8L,
      n_samples = 2000L,
      n_restarts = 2L,
      use_ic_stencil = TRUE,
      stencil_n = 800L,
      stencil_bx = 3L,
      stencil_by = 3L,
      stencil_bz = 2L,
      capture_boost = TRUE,
      capture_trans_mm = 9,
      capture_rot_deg = 9,
      capture_topk = 6L,
      whiten_modes = isTRUE(whiten_modes)
    )
  )

  est <- estimate_motion(sim$run, ref = "index", ref_index = 1, control = ctl)
  err <- est$par + poses

  c(
    trans_rmse = sqrt(mean(err[, 1:3]^2)),
    rot_rmse = sqrt(mean(err[, 4:6]^2)),
    fail_rate = mean(est$flags != 0),
    mean_iters = mean(est$iters)
  )
}

.run_nuisance_fixture <- function(seed, use_nuisance) {
  subject <- make_continuous_subject(seed = 6200 + seed)
  nt <- 12
  tt <- seq(0, 1, length.out = nt)

  poses <- cbind(
    2.0 * sin(2 * pi * tt),
    1.1 * cos(2 * pi * tt),
    0.8 * sin(3 * pi * tt),
    .deg(2.5) * sin(2 * pi * tt),
    .deg(2.0) * cos(2 * pi * tt),
    .deg(2.2) * sin(3 * pi * tt)
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
    drift_amp = 0.06,
    bias_amp = 0.12,
    supersample = 1L,
    seed = 7200 + seed
  )

  ctl <- utils::modifyList(
    volreg_control(),
    list(
      max_iter = 8L,
      n_samples = 2200L,
      n_restarts = 2L,
      use_ic_stencil = TRUE,
      stencil_n = 1800L,
      whiten_modes = TRUE,
      capture_boost = TRUE,
      use_nuisance_marginalization = isTRUE(use_nuisance),
      nuisance_include_template = FALSE,
      nuisance_orthogonalize = TRUE,
      nuisance_ridge = 1e-4
    )
  )

  est <- estimate_motion(sim$run, ref = "index", ref_index = 1, control = ctl)
  err <- est$par + poses

  c(
    trans_rmse = sqrt(mean(err[, 1:3]^2)),
    rot_rmse = sqrt(mean(err[, 4:6]^2)),
    fail_rate = mean(est$flags != 0),
    mean_iters = mean(est$iters)
  )
}

.capture_success_count <- function(capture_trans_mm, capture_rot_deg, capture_topk) {
  subject <- make_continuous_subject(seed = 2555)
  dx_grid <- c(-10, -6, -2, 2, 6, 10)
  rz_grid <- .deg(c(-12, 12))

  n_ok <- 0L

  for (dx in dx_grid) {
    for (rz in rz_grid) {
      poses <- matrix(0, nrow = 2, ncol = 6)
      poses[2, 1] <- dx
      poses[2, 6] <- rz

      sim <- simulate_continuous_run(
        subject = subject,
        nx = 18,
        ny = 18,
        nz = 9,
        nt = 2,
        poses = poses,
        noise_sd = 1.0,
        drift_amp = 0,
        bias_amp = 0.03,
        supersample = 1L,
        seed = 8200 + as.integer(40 * dx + 200 * rz)
      )

      ctl <- utils::modifyList(
        volreg_control(),
        list(
          max_iter = 8L,
          n_samples = 1800L,
          n_restarts = 2L,
          use_ic_stencil = TRUE,
          stencil_n = 800L,
          stencil_bx = 3L,
          stencil_by = 3L,
          stencil_bz = 2L,
          whiten_modes = TRUE,
          capture_boost = TRUE,
          capture_trans_mm = capture_trans_mm,
          capture_rot_deg = capture_rot_deg,
          capture_topk = capture_topk
        )
      )

      est <- estimate_motion(sim$run, ref = "index", ref_index = 1, control = ctl)
      expect_true(all(is.finite(est$par)))
      expect_true(all(is.finite(est$cost_final)))

      ex <- est$par[2, 1] + dx
      er <- est$par[2, 6] + rz
      if (abs(ex) < 2.0 && abs(er) < .deg(5.0)) {
        n_ok <- n_ok + 1L
      }
    }
  }

  n_ok
}

test_that("ic whitening is numerically stable and non-inferior on coupled motion", {
  seeds <- 1:6
  base <- t(vapply(seeds, function(s) .run_ic_whiten_fixture(s, FALSE), numeric(4)))
  white <- t(vapply(seeds, function(s) .run_ic_whiten_fixture(s, TRUE), numeric(4)))

  expect_true(all(is.finite(base)))
  expect_true(all(is.finite(white)))

  colnames(base) <- colnames(white) <- c("trans_rmse", "rot_rmse", "fail_rate", "mean_iters")

  trans_ratio <- median(white[, "trans_rmse"] / base[, "trans_rmse"])
  rot_ratio <- median(white[, "rot_rmse"] / base[, "rot_rmse"])
  fail_delta <- median(white[, "fail_rate"] - base[, "fail_rate"])
  iter_ratio <- median(white[, "mean_iters"] / base[, "mean_iters"])

  expect_true(trans_ratio <= 1.08, paste0("trans ratio=", signif(trans_ratio, 4)))
  expect_true(rot_ratio <= 1.05, paste0("rot ratio=", signif(rot_ratio, 4)))
  expect_true(fail_delta <= 0.15, paste0("fail delta=", signif(fail_delta, 4)))
  expect_true(iter_ratio <= 1.12, paste0("iter ratio=", signif(iter_ratio, 4)))
})

test_that("nuisance marginalization improves bias/drift robustness on IC stencil path", {
  seeds <- 1:6
  base <- t(vapply(seeds, function(s) .run_nuisance_fixture(s, FALSE), numeric(4)))
  nuis <- t(vapply(seeds, function(s) .run_nuisance_fixture(s, TRUE), numeric(4)))

  expect_true(all(is.finite(base)))
  expect_true(all(is.finite(nuis)))

  colnames(base) <- colnames(nuis) <- c("trans_rmse", "rot_rmse", "fail_rate", "mean_iters")

  trans_ratio <- median(nuis[, "trans_rmse"] / base[, "trans_rmse"])
  rot_ratio <- median(nuis[, "rot_rmse"] / base[, "rot_rmse"])
  fail_delta <- median(nuis[, "fail_rate"] - base[, "fail_rate"])
  iter_ratio <- median(nuis[, "mean_iters"] / base[, "mean_iters"])

  expect_true(trans_ratio <= 0.98, paste0("trans ratio=", signif(trans_ratio, 4)))
  expect_true(rot_ratio <= 0.95, paste0("rot ratio=", signif(rot_ratio, 4)))
  expect_true(fail_delta <= 0.10, paste0("fail delta=", signif(fail_delta, 4)))
  expect_true(iter_ratio <= 1.20, paste0("iter ratio=", signif(iter_ratio, 4)))
})

test_that("expanded capture settings stay stable on very-large motion cases", {
  ok_default <- .capture_success_count(
    capture_trans_mm = 8,
    capture_rot_deg = 8,
    capture_topk = 4L
  )
  ok_expanded <- .capture_success_count(
    capture_trans_mm = 14,
    capture_rot_deg = 14,
    capture_topk = 12L
  )

  expect_true(ok_expanded >= ok_default - 1L, paste0("default=", ok_default, " expanded=", ok_expanded))
  expect_true(ok_expanded >= 2L, paste0("expanded=", ok_expanded))
})
