.deg <- function(x) x * pi / 180

.make_diff_template <- function(nx = 24, ny = 24, nz = 8) {
  grid <- expand.grid(x = seq_len(nx), y = seq_len(ny), z = seq_len(nz))
  d <- sqrt((grid$x - 12.5)^2 + (grid$y - 12.5)^2 + (grid$z - 4.5)^2)
  tpl <- array(0, dim = c(nx, ny, nz))
  tpl[d <= 7] <- 100
  tpl[8:10, 15:17, 3:5] <- 140
  for (i in seq_len(nx)) {
    tpl[i, , ] <- tpl[i, , ] + 0.8 * i
  }
  tpl
}

.build_source_from_pose <- function(tpl, pose6, noise_sd, seed) {
  set.seed(seed)
  src <- apply_motion(array(tpl, dim = c(dim(tpl), 1L)), matrix(pose6, nrow = 1), interp = "linear")[, , , 1]
  src + array(rnorm(length(src), sd = noise_sd), dim = dim(src))
}

.run_volreg_pair <- function(tpl, src, control_overrides = list()) {
  run <- array(0, dim = c(dim(tpl), 2L))
  run[, , , 1] <- tpl
  run[, , , 2] <- src

  ctl <- utils::modifyList(
    volreg_control(),
    list(
      max_iter = 8L,
      n_samples = 3000L,
      n_restarts = 2L,
      use_ic_stencil = TRUE,
      stencil_n = 2200L,
      whiten_modes = TRUE,
      capture_boost = TRUE
    )
  )
  if (length(control_overrides) > 0) {
    ctl <- utils::modifyList(ctl, control_overrides)
  }

  est <- estimate_motion(run, ref = "index", ref_index = 1, control = ctl)
  list(est = est, corrected = apply_motion(run, est$par, interp = "linear")[, , , 2])
}

.run_rn_pair <- function(tpl, src, estimate_only = FALSE) {
  reg <- RNiftyReg::niftyreg.linear(
    source = src,
    target = tpl,
    scope = "rigid",
    symmetric = FALSE,
    estimateOnly = isTRUE(estimate_only),
    nLevels = 3L,
    maxIterations = 5L,
    useBlockPercentage = 100L,
    interpolation = 1L,
    threads = 1L
  )
  list(
    image = if (isTRUE(estimate_only)) NULL else as.array(reg$image),
    affine = unclass(RNiftyReg::asAffine(reg$forwardTransforms[[1]]))
  )
}

.compute_metrics <- function(tpl, src, vol_corr, rn_corr) {
  mask <- is.finite(tpl) & is.finite(src) & is.finite(vol_corr) & is.finite(rn_corr)

  mse_before <- mean((src[mask] - tpl[mask])^2)
  mse_vol <- mean((vol_corr[mask] - tpl[mask])^2)
  mse_rn <- mean((rn_corr[mask] - tpl[mask])^2)

  corr_before <- stats::cor(src[mask], tpl[mask])
  corr_vol <- stats::cor(vol_corr[mask], tpl[mask])
  corr_rn <- stats::cor(rn_corr[mask], tpl[mask])

  list(
    mask_frac = mean(mask),
    mse_before = mse_before,
    mse_vol = mse_vol,
    mse_rn = mse_rn,
    corr_before = corr_before,
    corr_vol = corr_vol,
    corr_rn = corr_rn
  )
}

.run_dual_oracle_from_source <- function(tpl, src, control_overrides = list()) {
  vol <- .run_volreg_pair(tpl, src, control_overrides = control_overrides)
  rn <- .run_rn_pair(tpl, src, estimate_only = FALSE)
  list(
    vol = vol,
    rn = rn,
    metrics = .compute_metrics(tpl, src, vol$corrected, rn$image)
  )
}

test_that("RNiftyReg rigid registration is near identity on identical volumes", {
  skip_if_not_installed("RNiftyReg")

  set.seed(7001)
  x <- .make_diff_template()
  x <- x + array(rnorm(length(x), sd = 0.2), dim = dim(x))

  rn <- .run_rn_pair(tpl = x, src = x, estimate_only = TRUE)
  aff <- rn$affine
  R <- aff[1:3, 1:3, drop = FALSE]

  expect_equal(dim(aff), c(4, 4))
  expect_true(max(abs(aff - diag(4))) < 1e-3)
  expect_true(max(abs(crossprod(R) - diag(3))) < 1e-3)
  expect_true(all(is.finite(aff)))
})

test_that("Differential grid: volregger and RNiftyReg improve across rigid mismatch scenarios", {
  skip_if_not_installed("RNiftyReg")

  tpl <- .make_diff_template()

  scenarios <- list(
    list(name = "mild-mixed", pose = c(0.8, -0.4, 0.2, .deg(1.0), .deg(-0.5), .deg(0.7)), noise_sd = 1.0, seed = 7101),
    list(name = "translation-heavy", pose = c(-1.4, 1.2, -0.6, .deg(0.7), .deg(-0.5), .deg(0.8)), noise_sd = 1.3, seed = 7102),
    list(name = "rotation-heavy", pose = c(0.5, 0.3, -0.4, .deg(2.0), .deg(-1.7), .deg(1.6)), noise_sd = 1.6, seed = 7103)
  )

  vol_vs_rn <- numeric(length(scenarios))

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    src <- .build_source_from_pose(tpl, pose6 = sc$pose, noise_sd = sc$noise_sd, seed = sc$seed)
    out <- .run_dual_oracle_from_source(tpl, src)
    m <- out$metrics

    expect_true(m$mask_frac > 0.3, paste0(sc$name, ": mask_frac=", signif(m$mask_frac, 4)))
    expect_true(m$mse_vol < m$mse_before * 0.995, paste0(sc$name, ": vol mse ratio=", signif(m$mse_vol / m$mse_before, 4)))
    expect_true(m$mse_rn < m$mse_before * 0.70, paste0(sc$name, ": rn mse ratio=", signif(m$mse_rn / m$mse_before, 4)))
    expect_true(m$corr_vol > m$corr_before, paste0(sc$name, ": vol corr gain=", signif(m$corr_vol - m$corr_before, 4)))
    expect_true(m$corr_rn > m$corr_before + 0.01, paste0(sc$name, ": rn corr gain=", signif(m$corr_rn - m$corr_before, 4)))

    vol_vs_rn[i] <- m$mse_vol / max(m$mse_rn, 1e-10)
  }

  expect_true(median(vol_vs_rn) < 4.5, paste0("median vol/rn mse=", signif(median(vol_vs_rn), 4)))
  expect_true(max(vol_vs_rn) < 8.0, paste0("max vol/rn mse=", signif(max(vol_vs_rn), 4)))
})

test_that("Metamorphic: affine intensity transform preserves oracle-relative improvement ratios", {
  skip_if_not_installed("RNiftyReg")

  tpl <- .make_diff_template()
  pose <- c(1.0, -0.6, 0.4, .deg(1.4), .deg(-1.0), .deg(1.1))
  src <- .build_source_from_pose(tpl, pose6 = pose, noise_sd = 1.2, seed = 7201)

  base <- .run_dual_oracle_from_source(tpl, src)

  gain <- 1.7
  offset <- 25
  tpl_aff <- gain * tpl + offset
  src_aff <- gain * src + offset
  aff <- .run_dual_oracle_from_source(tpl_aff, src_aff)

  vol_ratio_base <- base$metrics$mse_vol / base$metrics$mse_before
  rn_ratio_base <- base$metrics$mse_rn / base$metrics$mse_before
  vol_ratio_aff <- aff$metrics$mse_vol / aff$metrics$mse_before
  rn_ratio_aff <- aff$metrics$mse_rn / aff$metrics$mse_before

  expect_true(vol_ratio_base < 1.0)
  expect_true(rn_ratio_base < 1.0)
  expect_true(vol_ratio_aff < 1.0)
  expect_true(rn_ratio_aff < 1.0)
  expect_true(abs(vol_ratio_base - vol_ratio_aff) < 0.10, paste0("vol ratio delta=", signif(abs(vol_ratio_base - vol_ratio_aff), 4)))
  expect_true(abs(rn_ratio_base - rn_ratio_aff) < 0.06, paste0("rn ratio delta=", signif(abs(rn_ratio_base - rn_ratio_aff), 4)))
})

test_that("Property sweep: seeded random rigid mismatches stay finite and mostly improve", {
  skip_if_not_installed("RNiftyReg")

  tpl <- .make_diff_template()
  set.seed(7301)

  n_cases <- 6L
  poses <- cbind(
    runif(n_cases, -1.4, 1.4),
    runif(n_cases, -1.4, 1.4),
    runif(n_cases, -0.9, 0.9),
    .deg(runif(n_cases, -1.8, 1.8)),
    .deg(runif(n_cases, -1.8, 1.8)),
    .deg(runif(n_cases, -1.8, 1.8))
  )
  noise_levels <- runif(n_cases, 0.9, 1.8)

  vol_improved <- logical(n_cases)
  rn_improved <- logical(n_cases)
  vol_vs_rn <- numeric(n_cases)

  for (i in seq_len(n_cases)) {
    src <- .build_source_from_pose(
      tpl = tpl,
      pose6 = poses[i, ],
      noise_sd = noise_levels[i],
      seed = 7301 + i
    )
    out <- .run_dual_oracle_from_source(tpl, src)
    m <- out$metrics

    expect_true(all(is.finite(unlist(m))), paste0("case=", i, " produced non-finite metrics"))

    vol_improved[i] <- m$mse_vol < m$mse_before
    rn_improved[i] <- m$mse_rn < m$mse_before
    vol_vs_rn[i] <- m$mse_vol / max(m$mse_rn, 1e-10)
  }

  expect_true(mean(rn_improved) >= 5 / 6, paste0("rn success=", sum(rn_improved), "/", n_cases))
  expect_true(mean(vol_improved) >= 4 / 6, paste0("vol success=", sum(vol_improved), "/", n_cases))
  expect_true(median(vol_vs_rn) < 5.0, paste0("median vol/rn mse=", signif(median(vol_vs_rn), 4)))
})
