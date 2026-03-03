make_vignette_phantom <- function(nx = 16, ny = 16, nz = 8) {
  grid <- expand.grid(x = seq_len(nx), y = seq_len(ny), z = seq_len(nz))
  d <- sqrt((grid$x - 8.5)^2 + (grid$y - 8.5)^2 + (grid$z - 4.5)^2)
  x <- array(0, dim = c(nx, ny, nz))
  x[d <= 5] <- 100

  d2 <- sqrt((grid$x - 5.5)^2 + (grid$y - 11)^2 + (grid$z - 3.5)^2)
  x[d2 <= 2] <- 140

  for (i in seq_len(nx)) {
    x[i, , ] <- x[i, , ] + 1.2 * i
  }
  x
}

replicate_run <- function(base, nt) {
  array(rep(as.numeric(base), nt), dim = c(dim(base), nt))
}

warp_single_frame <- function(vol3d, pose6) {
  arr4 <- array(vol3d, dim = c(dim(vol3d), 1L))
  out <- apply_motion(arr4, matrix(pose6, nrow = 1), interp = "lanczos")
  out[, , , 1]
}

build_sparse_motion_run <- function(base, noise_sd, seed, injected_motion) {
  set.seed(seed)

  nt <- nrow(injected_motion)
  run <- replicate_run(base, nt)
  run <- run + array(rnorm(length(run), sd = noise_sd), dim = dim(run))

  for (t in seq_len(nt)) {
    if (any(injected_motion[t, ] != 0)) {
      run[, , , t] <- warp_single_frame(run[, , , t], injected_motion[t, ])
    }
  }

  run
}

evaluate_sparse_case <- function(base,
                                 noise_sd,
                                 seed,
                                 injected_motion,
                                 improve_frames = integer(),
                                 control_overrides = list()) {
  run <- build_sparse_motion_run(base, noise_sd, seed, injected_motion)

  ctl <- utils::modifyList(
    volreg_control(),
    list(max_iter = 8L, n_samples = 2500L, n_restarts = 2L)
  )
  if (length(control_overrides) > 0) {
    ctl <- utils::modifyList(ctl, control_overrides)
  }

  est <- estimate_motion(run, control = ctl)

  truth_correction <- -injected_motion
  moved_idx <- which(rowSums(abs(injected_motion)) > 0)
  err <- est$par[moved_idx, , drop = FALSE] - truth_correction[moved_idx, , drop = FALSE]

  trans_active <- any(abs(injected_motion[, 1:3]) > 0)
  rot_active <- any(abs(injected_motion[, 4:6]) > 0)

  trans_rmse <- if (isTRUE(trans_active)) sqrt(mean(err[, 1:3]^2)) else NA_real_
  rot_rmse <- if (isTRUE(rot_active)) sqrt(mean(err[, 4:6]^2)) else NA_real_

  improve_ratio <- NA_real_
  if (length(improve_frames) > 0) {
    corrected <- apply_motion(run, est$par, interp = "lanczos")
    base_run <- replicate_run(base, nrow(injected_motion))
    mse_before <- mean((run[, , , improve_frames, drop = FALSE] - base_run[, , , improve_frames, drop = FALSE])^2)
    mse_after <- mean((corrected[, , , improve_frames, drop = FALSE] - base_run[, , , improve_frames, drop = FALSE])^2)
    improve_ratio <- mse_after / mse_before
  }

  list(
    trans_rmse = trans_rmse,
    rot_rmse = rot_rmse,
    improve_ratio = improve_ratio,
    est = est
  )
}

test_that("oracle translation sweep: sparse shifts are recovered across noise levels", {
  base <- make_vignette_phantom()
  nt <- 10

  inj <- matrix(0, nrow = nt, ncol = 6)
  inj[4, 1] <- 1.2
  inj[7, 1] <- 2.0
  inj[9, 1] <- -1.3

  scenarios <- list(
    list(name = "low-noise", noise_sd = 1.0, seed = 810),
    list(name = "mid-noise", noise_sd = 2.0, seed = 811),
    list(name = "high-noise", noise_sd = 3.0, seed = 812)
  )

  for (sc in scenarios) {
    out <- evaluate_sparse_case(
      base = base,
      noise_sd = sc$noise_sd,
      seed = sc$seed,
      injected_motion = inj,
      improve_frames = c(4, 7, 9)
    )

    expect_true(out$trans_rmse < 0.45, paste0(sc$name, ": trans_rmse=", signif(out$trans_rmse, 4)))
    expect_true(out$improve_ratio < 0.40, paste0(sc$name, ": improve_ratio=", signif(out$improve_ratio, 4)))
  }
})

test_that("oracle rotation sweep: sparse rotations are recovered with bounded angular error", {
  base <- make_vignette_phantom()
  nt <- 10
  r <- pi / 180

  inj <- matrix(0, nrow = nt, ncol = 6)
  inj[6, 6] <- 2.0 * r
  inj[8, 4] <- -2.0 * r

  scenarios <- list(
    list(name = "low-noise", noise_sd = 1.0, seed = 812),
    list(name = "high-noise", noise_sd = 3.0, seed = 813)
  )

  for (sc in scenarios) {
    out <- evaluate_sparse_case(
      base = base,
      noise_sd = sc$noise_sd,
      seed = sc$seed,
      injected_motion = inj
    )

    expect_true(out$rot_rmse < 0.05, paste0(sc$name, ": rot_rmse=", signif(out$rot_rmse, 4)))
    moved_idx <- c(6, 8)
    expect_true(
      max(abs(out$est$par[moved_idx, 4:6])) < 0.25,
      paste0(sc$name, ": unstable rotation magnitude")
    )
  }
})

test_that("oracle combined sweep: mixed sparse translation+rotation remains within tolerance", {
  base <- make_vignette_phantom()
  nt <- 10
  r <- pi / 180

  inj <- matrix(0, nrow = nt, ncol = 6)
  inj[4, 1] <- 1.2
  inj[7, 1] <- 2.0
  inj[9, 1] <- -1.3
  inj[6, 6] <- 2.0 * r
  inj[8, 4] <- -2.0 * r

  scenarios <- list(
    list(name = "mid-noise", noise_sd = 2.0, seed = 700),
    list(name = "high-noise", noise_sd = 3.0, seed = 701)
  )

  for (sc in scenarios) {
    out <- evaluate_sparse_case(
      base = base,
      noise_sd = sc$noise_sd,
      seed = sc$seed,
      injected_motion = inj,
      improve_frames = c(4, 7, 9)
    )

    expect_true(out$trans_rmse < 0.40, paste0(sc$name, ": trans_rmse=", signif(out$trans_rmse, 4)))
    expect_true(out$rot_rmse < 0.10, paste0(sc$name, ": rot_rmse=", signif(out$rot_rmse, 4)))
    expect_true(out$improve_ratio < 0.40, paste0(sc$name, ": improve_ratio=", signif(out$improve_ratio, 4)))
  }
})

test_that("ic stencil mode: sparse translation case satisfies correctness contracts", {
  base <- make_vignette_phantom()
  nt <- 10
  inj <- matrix(0, nrow = nt, ncol = 6)
  inj[4, 1] <- 1.2
  inj[7, 1] <- 2.0
  inj[9, 1] <- -1.3

  out <- evaluate_sparse_case(
    base = base,
    noise_sd = 2.0,
    seed = 811,
    injected_motion = inj,
    improve_frames = c(4, 7, 9),
    control_overrides = list(
      use_ic_stencil = TRUE,
      stencil_n = 3000L,
      stencil_bx = 4L,
      stencil_by = 4L,
      stencil_bz = 2L
    )
  )

  expect_true(out$trans_rmse < 0.35, paste0("ic trans_rmse=", signif(out$trans_rmse, 4)))
  expect_true(out$improve_ratio < 0.20, paste0("ic improve_ratio=", signif(out$improve_ratio, 4)))
  expect_true(all(is.finite(out$est$par)))
  expect_true(all(is.finite(out$est$cost_final)))
})
