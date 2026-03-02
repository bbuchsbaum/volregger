make_phantom <- function(nx = 24, ny = 24, nz = 8) {
  gx <- array(rep(seq(-1, 1, length.out = nx), each = ny * nz), dim = c(nx, ny, nz))
  gy <- array(rep(rep(seq(-1, 1, length.out = ny), each = nz), times = nx), dim = c(nx, ny, nz))
  gz <- array(rep(seq(-1, 1, length.out = nz), times = nx * ny), dim = c(nx, ny, nz))

  a <- exp(-((gx - 0.25)^2 + (gy + 0.2)^2 + (gz)^2) / 0.12)
  b <- 0.7 * exp(-((gx + 0.45)^2 + (gy - 0.35)^2 + (gz + 0.1)^2) / 0.06)
  c <- 0.1 * gx - 0.05 * gy + 0.03 * gz
  a + b + c
}

shift_volume_int <- function(v, dx, dy, dz) {
  d <- dim(v)
  out <- array(0, dim = d)
  for (k in seq_len(d[3])) {
    for (j in seq_len(d[2])) {
      for (i in seq_len(d[1])) {
        ii <- i - dx
        jj <- j - dy
        kk <- k - dz
        if (ii >= 1 && ii <= d[1] && jj >= 1 && jj <= d[2] && kk >= 1 && kk <= d[3]) {
          out[i, j, k] <- v[ii, jj, kk]
        }
      }
    }
  }
  out
}

make_shifted_run <- function(base, shifts) {
  nt <- nrow(shifts)
  x <- array(0, dim = c(dim(base), nt))
  for (t in seq_len(nt)) {
    x[, , , t] <- shift_volume_int(base, shifts[t, 1], shifts[t, 2], shifts[t, 3])
  }
  x
}

test_that("synthetic contract: estimator returns finite bounded diagnostics", {
  base <- make_phantom()
  shifts <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(-1, 0, 1), c(1, -1, 0), c(0, 0, -1))
  x <- make_shifted_run(base, shifts)

  ctl <- utils::modifyList(volreg_control(), list(max_iter = 8L, n_samples = 5000L, n_restarts = 2L))
  est <- estimate_motion(x, control = ctl)

  expect_equal(dim(est$par), c(nrow(shifts), 6))
  expect_true(all(is.finite(est$par)))
  expect_true(all(is.finite(est$cost_init)))
  expect_true(all(is.finite(est$cost_final)))
  expect_true(all(est$overlap >= 0 & est$overlap <= 1))
  expect_true(all(est$iters >= 0))
})

test_that("synthetic oracle: objective improves on known shifted phantom run", {
  base <- make_phantom()
  shifts <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(-1, 0, 1), c(1, -1, 0), c(0, 0, -1))
  x <- make_shifted_run(base, shifts)

  ctl <- utils::modifyList(volreg_control(), list(max_iter = 8L, n_samples = 5000L, n_restarts = 2L))
  est <- estimate_motion(x, control = ctl)

  improved <- est$cost_final <= est$cost_init
  expect_true(mean(improved) >= 0.9)
  expect_true(mean(est$cost_final) < mean(est$cost_init))
})

test_that("metamorphic: global rigid shift of all frames keeps FD near zero", {
  set.seed(300)
  base <- make_phantom()
  nt <- 6
  x0 <- array(rep(as.numeric(base), nt), dim = c(dim(base), nt))

  x_shifted <- array(0, dim = dim(x0))
  for (t in seq_len(nt)) {
    x_shifted[, , , t] <- shift_volume_int(x0[, , , t], dx = 2, dy = -1, dz = 1)
  }

  ctl <- utils::modifyList(volreg_control(), list(max_iter = 6L, n_samples = 3500L))
  e0 <- estimate_motion(x0, control = ctl)
  e1 <- estimate_motion(x_shifted, control = ctl)

  fd0 <- framewise_disp(e0$par)
  fd1 <- framewise_disp(e1$par)
  expect_true(max(fd0) < 1e-8)
  expect_true(max(fd1) < 1e-8)
})

test_that("property fuzz: random small-shift runs remain stable and improve", {
  base <- make_phantom()
  nt <- 7
  ctl <- utils::modifyList(volreg_control(), list(max_iter = 7L, n_samples = 4000L, n_restarts = 2L))

  for (s in 1:10) {
    set.seed(400 + s)
    shifts <- cbind(sample(-1:1, nt, replace = TRUE), sample(-1:1, nt, replace = TRUE), sample(-1:1, nt, replace = TRUE))
    shifts[1, ] <- c(0, 0, 0)
    x <- make_shifted_run(base, shifts)

    est <- estimate_motion(x, control = ctl)

    expect_true(all(is.finite(est$par)))
    expect_true(all(is.finite(est$cost_final)))
    expect_true(mean(est$cost_final <= est$cost_init) >= 0.85)
  }
})

test_that("adversarial: dropout corruption yields finite fit and nonzero warning flags", {
  set.seed(500)
  x <- array(rnorm(20 * 20 * 8 * 8), dim = c(20, 20, 8, 8))
  x[, , , 5] <- 0
  x[, , 1:2, 6] <- 0

  est <- estimate_motion(x, control = utils::modifyList(volreg_control(), list(max_iter = 6L, n_samples = 3000L, n_restarts = 2L)))

  expect_true(all(is.finite(est$cost_final)))
  expect_true(any(est$flags != 0))
})

test_that("regression: non-finite input is rejected early", {
  x <- array(0, dim = c(8, 8, 4, 3))
  x[1, 1, 1, 2] <- NA_real_

  expect_error(
    estimate_motion(x),
    regexp = "non-finite",
    fixed = FALSE
  )
})
