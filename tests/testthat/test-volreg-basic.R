test_that("volreg returns expected structure on identity-like input", {
  set.seed(1)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 5), dim = c(8, 8, 4, 5))

  fit <- volreg(x, resample = TRUE, qc = TRUE)

  expect_s3_class(fit, "volregger_fit")
  expect_equal(dim(fit$corrected), dim(x))
  expect_equal(dim(fit$motion), c(5, 6))
  expect_equal(dim(fit$mats), c(16, 5))
  expect_equal(length(fit$qc$fd), 5)
  expect_true(all(abs(fit$motion[, 1:3]) < 1e-6))
  expect_true(all(abs(fit$motion[, 4:6]) < 1e-12))
})

test_that("apply_motion returns same shape", {
  set.seed(2)
  x <- array(rnorm(6 * 6 * 4 * 4), dim = c(6, 6, 4, 4))
  motion <- matrix(0, nrow = 4, ncol = 6)

  y <- apply_motion(x, motion)

  expect_equal(dim(y), dim(x))
})

test_that("estimate_motion detects known translation perturbation", {
  x0 <- array(0, dim = c(9, 9, 5))
  x0[5, 5, 3] <- 5
  x0[4, 5, 3] <- 2
  x0[6, 5, 3] <- 1

  x <- array(0, dim = c(9, 9, 5, 3))
  x[, , , 1] <- x0
  x[, , , 2] <- x0
  x[, , , 3] <- x0

  # Shift frame 2 by +1 in x direction.
  x[2:9, , , 2] <- x0[1:8, , ]

  est <- estimate_motion(x, control = utils::modifyList(volreg_control(), list(max_iter = 8L, n_samples = 2000L)))

  expect_equal(dim(est$par), c(3, 6))
  # Middle-frame reference implies the shifted frame is near zero and
  # neighboring frames carry the compensating translation.
  expect_true(abs(est$par[2, 1]) < 0.3)
  expect_true(max(abs(est$par[c(1, 3), 1])) > 0.3)
  expect_true(any(est$cost_final[c(1, 3)] < est$cost_init[c(1, 3)]))
})

test_that("robust_mean reference anchors to the consensus frame", {
  x0 <- array(0, dim = c(9, 9, 5))
  x0[5, 5, 3] <- 5
  x0[4, 5, 3] <- 2
  x0[6, 5, 3] <- 1

  x <- array(0, dim = c(9, 9, 5, 3))
  x[, , , 1] <- x0
  x[, , , 2] <- x0
  x[, , , 3] <- x0

  x[2:9, , , 2] <- x0[1:8, , ]

  ctl <- utils::modifyList(volreg_control(), list(max_iter = 8L, n_samples = 2000L))
  est_mid <- estimate_motion(x, ref = "middle", control = ctl)
  est_rob <- estimate_motion(x, ref = "robust_mean", control = ctl)

  expect_true(abs(est_mid$par[2, 1]) < 0.3)
  expect_true(abs(est_rob$par[2, 1]) > 0.3)
  expect_true(abs(est_rob$par[2, 1]) > max(abs(est_rob$par[c(1, 3), 1])))
  expect_true(mean(abs(est_rob$par[c(1, 3), 1])) < mean(abs(est_mid$par[c(1, 3), 1])))
})

test_that("apply_motion applies translation components", {
  x <- array(0, dim = c(7, 7, 3, 1))
  x[4, 4, 2, 1] <- 1
  motion <- matrix(0, nrow = 1, ncol = 6)
  motion[1, 1] <- 1

  y <- apply_motion(x, motion, interp = "linear")

  peak_lin <- which.max(y[, , 2, 1])
  peak_idx <- arrayInd(peak_lin, .dim = c(7, 7))
  expect_equal(unname(peak_idx[1, 1]), 5)
  expect_equal(unname(peak_idx[1, 2]), 4)
})

test_that("apply_motion applies rotation components", {
  x <- array(0, dim = c(7, 7, 3, 1))
  # point to the right of center
  x[6, 4, 2, 1] <- 1
  motion <- matrix(0, nrow = 1, ncol = 6)
  motion[1, 6] <- pi / 2

  y <- apply_motion(x, motion, interp = "linear")

  peak_lin <- which.max(y[, , 2, 1])
  peak_idx <- arrayInd(peak_lin, .dim = c(7, 7))
  # Expected around upper side after +90 deg around center
  expect_true(peak_idx[1, 1] <= 5)
  expect_true(peak_idx[1, 2] >= 5)
})

test_that("apply_motion interpolation modes are distinct on fractional shifts", {
  set.seed(10)
  x <- array(rnorm(8 * 8 * 4 * 1), dim = c(8, 8, 4, 1))
  motion <- matrix(0, nrow = 1, ncol = 6)
  motion[1, 1] <- 0.5
  motion[1, 2] <- -0.25

  y_lin <- apply_motion(x, motion, interp = "linear")
  y_bsp <- apply_motion(x, motion, interp = "bspline")
  y_lan <- apply_motion(x, motion, interp = "lanczos")

  expect_equal(dim(y_lin), dim(x))
  expect_equal(dim(y_bsp), dim(x))
  expect_equal(dim(y_lan), dim(x))

  d1 <- max(abs(y_lin - y_bsp))
  d2 <- max(abs(y_lin - y_lan))
  expect_true(d1 > 0)
  expect_true(d2 > 0)
})

test_that("slice-aware packet interpolation changes output when offsets supplied", {
  set.seed(11)
  x <- array(rnorm(8 * 8 * 4 * 3), dim = c(8, 8, 4, 3))
  motion <- rbind(
    c(0, 0, 0, 0, 0, 0),
    c(1, 0, 0, 0, 0, 0),
    c(2, 0, 0, 0, 0, 0)
  )

  y_vol <- apply_motion(x, motion, interp = "linear")
  y_pkt <- apply_motion(
    x, motion,
    interp = "linear",
    slice_times = c(0.0, 0.33, 0.66, 1.0)
  )

  expect_true(max(abs(y_vol - y_pkt)) > 0)
})

test_that("rigid_spline engine returns spline payload", {
  set.seed(12)
  x <- array(rnorm(8 * 8 * 4 * 4), dim = c(8, 8, 4, 4))
  fit <- volreg(
    x,
    engine = "rigid_spline",
    slice_times = c(0.0, 0.33, 0.66, 1.0),
    resample = TRUE
  )

  expect_false(is.null(fit$pose_spline))
  expect_equal(dim(fit$pose_spline$par_smooth), c(4, 6))
  expect_equal(length(fit$qc$packet_correction_mag), 4)
})

test_that("framewise_disp and dvars behave consistently", {
  x <- array(0, dim = c(5, 5, 3, 4))
  x[, , , 2] <- 1
  x[, , , 3] <- 2
  x[, , , 4] <- 3

  motion <- rbind(
    c(0, 0, 0, 0, 0, 0),
    c(1, 0, 0, 0, 0, 0),
    c(1, 1, 0, 0, 0, 0),
    c(1, 1, 1, 0, 0, 0)
  )

  fd <- framewise_disp(motion)
  dv <- dvars(x, robust = FALSE)

  expect_equal(length(fd), 4)
  expect_equal(fd[1], 0)
  expect_true(all(fd[-1] >= 0))
  expect_equal(length(dv), 4)
  expect_true(is.na(dv[1]))
  expect_true(all(is.finite(dv[-1])))
})

test_that("repair_spin_history returns structured correction outputs", {
  set.seed(20)
  x <- array(rnorm(8 * 8 * 4 * 5), dim = c(8, 8, 4, 5))
  motion <- matrix(0, nrow = 5, ncol = 6)
  motion[, 3] <- c(0, 0.5, 1.2, 0.7, 0.2)

  rep <- repair_spin_history(x, motion)

  expect_equal(dim(rep$corrected), dim(x))
  expect_equal(dim(rep$gains), c(4, 5))
  expect_equal(dim(rep$spike_mask), c(4, 5))
  expect_equal(length(rep$through_plane), 5)
})
