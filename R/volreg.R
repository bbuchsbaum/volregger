#' Run retrospective motion correction
#'
#' @param x 4D image-like input.
#' @param mask Optional 3D mask.
#' @param ref Reference strategy.
#' @param ref_index Explicit reference index for `ref='index'`.
#' @param engine Motion engine.
#' @param slice_times Optional slice timing offsets.
#' @param mb_groups Optional multiband group IDs.
#' @param resample If TRUE, apply estimated motion to generate corrected output.
#' @param interp_final Final interpolation method. When `NULL`, the active
#'   profile selects the final interpolator, falling back to `"bspline"`.
#' @param pad_final Boundary handling for final resampling.
#' @param zpad_final Number of symmetric zero-valued voxels to add before final
#'   resampling and cropping back to the original grid.
#' @param qc If TRUE, compute QC metrics.
#' @param control Control list.
#' @param profile Optional named estimator profile. Bare `volreg()` uses
#'   `"fast_fmri"` when `control` is not supplied; pass `profile = NULL`
#'   with an explicit `control` list for fully manual settings.
#' @param verbose Verbose logging.
#' @return A `volregger_fit` object.
#' @export
volreg <- function(
    x,
    mask = NULL,
    ref = c("middle", "robust_mean", "index"),
    ref_index = NULL,
    engine = c("rigid_robust", "rigid_spline"),
    slice_times = NULL,
    mb_groups = NULL,
    resample = TRUE,
    interp_final = NULL,
    pad_final = c("clamp", "zero"),
    zpad_final = 0L,
    qc = TRUE,
    control = volreg_control(),
    profile = NULL,
    verbose = FALSE) {
  control_missing <- missing(control)
  engine_missing <- missing(engine)
  ref <- match.arg(ref)
  pad_final <- match.arg(pad_final)
  zpad_final <- as.integer(zpad_final)
  if (length(zpad_final) != 1L || is.na(zpad_final) || zpad_final < 0L) {
    stop("zpad_final must be a non-negative integer", call. = FALSE)
  }
  if (is.null(profile) && control_missing) {
    profile <- "fast_fmri"
  }
  prof <- .apply_volreg_profile(control, profile)
  control <- prof$control
  valid_interp_final <- c("bspline", "lanczos", "lanczos5", "lanczos7", "quintic", "heptic", "linear")
  if (is.null(interp_final)) {
    interp_final <- control$interp_final %||% "bspline"
  }
  interp_final <- match.arg(interp_final, valid_interp_final)
  if (!is.null(prof$engine) && engine_missing) {
    engine <- prof$engine
  }
  engine <- match.arg(engine, c("rigid_robust", "rigid_spline"))

  inp <- .as_volreg_input(
    x = x,
    mask = mask,
    slice_times = slice_times,
    mb_groups = mb_groups
  )

  ctl <- utils::modifyList(control, list(engine = engine, verbose = isTRUE(verbose), tr = inp$tr))
  ctl$ref_mode <- ref

  est <- cpp_volreg_estimate(
    data = inp$data,
    dim = inp$dim,
    pixdim = inp$pixdim,
    mask = inp$mask,
    slice_times = inp$slice_times,
    mb_groups = inp$mb_groups,
    ref_index = .resolve_ref_index(ref, ref_index, inp$dim[4]),
    control = ctl
  )
  shrink <- .apply_low_motion_pose_shrink(est$par, ctl, est$mats)
  if (isTRUE(shrink$applied)) {
    est$par <- shrink$par
    est$mats <- shrink$mats
  }
  est$pose_scale_applied <- shrink$scale_applied
  est$pose_scale_fd_p95 <- shrink$fd_p95

  pose_spline <- NULL
  par_apply <- est$par
  if (identical(engine, "rigid_spline")) {
    pose_spline <- cpp_fit_pose_spline(
      par = est$par,
      tr = inp$tr,
      slice_times = inp$slice_times,
      mb_groups = inp$mb_groups,
      control = ctl
    )
    par_apply <- pose_spline$par_smooth
  }

  corrected <- NULL
  output_smooth_alpha <- 0
  if (isTRUE(resample)) {
    arr_corr <- cpp_volreg_apply(
      data = inp$data,
      dim = inp$dim,
      pixdim = inp$pixdim,
      par = par_apply,
      slice_times = if (identical(engine, "rigid_spline")) inp$slice_times else NULL,
      mb_groups = if (identical(engine, "rigid_spline")) inp$mb_groups else NULL,
      interp = interp_final,
      pad_mode = if (identical(pad_final, "zero")) 1L else 0L,
      zpad = zpad_final,
      nthreads = as.integer(ctl$nthreads)
    )
    arr_corr <- .apply_low_motion_output_smooth(arr_corr, ctl, shrink$fd_p95)
    output_smooth_alpha <- attr(arr_corr, "output_smooth_alpha", exact = TRUE) %||% 0
    corrected <- .rebuild_like(x, arr_corr)
  }

  qc_out <- NULL
  if (isTRUE(qc)) {
    qc_out <- motion_qc(
      x = x,
      motion = est$par,
      corrected = corrected,
      mask = est$mask,
      cost_init = est$cost_init,
      cost_final = est$cost_final
    )
    if (!is.null(pose_spline)) {
      qc_out$packet_correction_mag <- rowSums(abs(pose_spline$par_smooth - est$par))
    }
  }

  structure(
    list(
      corrected = corrected,
      motion = est$par,
      mats = est$mats,
      pose_spline = pose_spline,
      template = .rebuild_3d_like(x, est$template),
      mask = .rebuild_3d_like(x, est$mask),
      qc = qc_out,
      pose_scale_applied = est$pose_scale_applied,
      pose_scale_fd_p95 = est$pose_scale_fd_p95,
      output_smooth_alpha = output_smooth_alpha,
      control = ctl
    ),
    class = "volregger_fit"
  )
}

.apply_low_motion_pose_shrink <- function(par, control, mats = NULL) {
  scale <- as.numeric(control$low_motion_pose_scale %||% 1)
  fd_p95_max <- as.numeric(control$low_motion_fd_p95_max %||% 0)
  use <- isTRUE(control$use_low_motion_pose_shrink)
  if (is.null(mats)) {
    mats <- motion_matrices(par)
    pose <- par
  } else {
    mats <- motion_matrices(mats)
    pose <- .motion_mats_to_pose(mats)
  }
  fd <- .framewise_disp_from_mats(mats)
  fd_p95 <- unname(stats::quantile(fd, probs = 0.95, na.rm = TRUE))
  out <- list(
    par = par,
    mats = mats,
    applied = FALSE,
    scale_applied = 1,
    smooth_alpha_applied = 0,
    fd_p95 = fd_p95
  )
  low_motion <- is.finite(fd_p95) && is.finite(fd_p95_max) && fd_p95_max > 0 && fd_p95 <= fd_p95_max
  if (!low_motion) {
    return(out)
  }
  smooth_alpha <- as.numeric(control$low_motion_pose_smooth_alpha %||% 0)
  if (isTRUE(control$use_low_motion_pose_smooth) && is.finite(smooth_alpha) && smooth_alpha > 0) {
    smooth_alpha <- min(1, max(0, smooth_alpha))
    pose <- .smooth_pose_neighbors(pose, smooth_alpha)
    out$par <- pose
    out$mats <- motion_matrices(pose)
    out$applied <- TRUE
    out$smooth_alpha_applied <- smooth_alpha
  }
  if (use && is.finite(scale) && scale > 0 && scale < 1) {
    out$par <- pose * scale
    out$mats <- motion_matrices(out$par)
    out$applied <- TRUE
    out$scale_applied <- scale
  }
  out
}

.smooth_pose_neighbors <- function(pose, alpha) {
  pose <- as.matrix(pose)
  if (nrow(pose) < 3L || alpha <= 0) {
    return(pose)
  }
  out <- pose
  for (j in seq_len(ncol(pose))) {
    for (t in 2:(nrow(pose) - 1L)) {
      out[t, j] <- (1 - alpha) * pose[t, j] + alpha * 0.5 * (pose[t - 1L, j] + pose[t + 1L, j])
    }
  }
  out
}

.apply_low_motion_output_smooth <- function(x, control, fd_p95) {
  alpha <- as.numeric(control$low_motion_output_smooth_alpha %||% 0)
  fd_p95_max <- as.numeric(control$low_motion_fd_p95_max %||% 0)
  if (!isTRUE(control$use_low_motion_output_smooth) ||
      !is.finite(alpha) || alpha <= 0 ||
      !is.finite(fd_p95) || !is.finite(fd_p95_max) || fd_p95 > fd_p95_max) {
    attr(x, "output_smooth_alpha") <- 0
    return(x)
  }
  alpha <- min(1, max(0, alpha))
  nt <- dim(x)[4]
  if (nt < 3L) {
    attr(x, "output_smooth_alpha") <- 0
    return(x)
  }
  out <- x
  for (t in 2:(nt - 1L)) {
    out[, , , t] <- (1 - alpha) * x[, , , t] + alpha * 0.5 * (x[, , , t - 1L] + x[, , , t + 1L])
  }
  attr(out, "output_smooth_alpha") <- alpha
  out
}

.framewise_disp_from_mats <- function(mats, radius = 50) {
  mats <- motion_matrices(mats)
  nt <- dim(mats)[3]
  if (nt <= 1L) {
    return(rep(0, nt))
  }
  out <- numeric(nt)
  out[-1] <- fd_from_transform(
    mats[, , 2:nt, drop = FALSE],
    radius = radius,
    reference = mats[, , 1:(nt - 1L), drop = FALSE]
  )
  out
}

.motion_mats_to_pose <- function(mats) {
  mats <- motion_matrices(mats)
  nt <- dim(mats)[3]
  out <- matrix(0, nrow = nt, ncol = 6)
  for (i in seq_len(nt)) {
    M <- mats[, , i]
    R <- M[1:3, 1:3, drop = FALSE]
    s <- svd(R)
    R <- s$u %*% t(s$v)
    if (det(R) < 0) {
      s$u[, 3] <- -s$u[, 3]
      R <- s$u %*% t(s$v)
    }
    ry <- asin(max(-1, min(1, -R[3, 1])))
    cy <- cos(ry)
    if (abs(cy) > 1e-8) {
      rx <- atan2(R[3, 2], R[3, 3])
      rz <- atan2(R[2, 1], R[1, 1])
    } else {
      rx <- atan2(-R[2, 3], R[2, 2])
      rz <- 0
    }
    out[i, ] <- c(M[1:3, 4], rx, ry, rz)
  }
  out
}
