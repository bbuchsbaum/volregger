`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

load_volregger_benchmark_env <- function() {
  if (!requireNamespace("RNifti", quietly = TRUE)) {
    stop("Need RNifti installed to run real-data benchmarks", call. = FALSE)
  }

  if (file.exists("DESCRIPTION")) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("Need devtools to load local package from source", call. = FALSE)
    }
    devtools::load_all(".", quiet = TRUE)
  } else if (!requireNamespace("volregger", quietly = TRUE)) {
    stop("Need installed volregger package when running outside source tree", call. = FALSE)
  } else {
    library(volregger)
  }

  invisible(TRUE)
}

parse_benchmark_args <- function(args, defaults = list(), flags = character()) {
  out <- defaults
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }

    key <- sub("^--", "", arg)
    value <- NULL
    has_inline_value <- FALSE
    if (grepl("=", key, fixed = TRUE)) {
      bits <- strsplit(key, "=", fixed = TRUE)[[1]]
      key <- bits[[1]]
      value <- paste(bits[-1], collapse = "=")
      has_inline_value <- TRUE
    }

    field <- gsub("-", "_", key, fixed = TRUE)
    if (field %in% flags) {
      out[[field]] <- TRUE
      i <- i + 1L
      next
    }

    if (!has_inline_value) {
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        stop("Option requires a value: --", key, call. = FALSE)
      }
      value <- args[[i + 1L]]
      i <- i + 1L
    }

    out[[field]] <- value
    i <- i + 1L
  }

  out
}

resolve_benchmark_nthreads <- function(value = "auto", max_auto = 4L) {
  if (is.null(value)) {
    value <- "auto"
  }
  if (identical(tolower(as.character(value)), "auto")) {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores) || cores < 1L) {
      cores <- parallel::detectCores(logical = TRUE)
    }
    if (is.na(cores) || cores < 1L) {
      cores <- 1L
    }
    return(as.integer(max(1L, min(as.integer(max_auto), as.integer(cores)))))
  }

  n <- as.integer(value)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("nthreads must be a positive integer or 'auto'", call. = FALSE)
  }
  n
}

benchmark_scalar_last <- function(value, default = NA) {
  if (is.null(value) || !length(value)) {
    return(default)
  }
  value[[length(value)]]
}

parse_frame_limit <- function(value, total_frames) {
  if (is.null(value)) {
    return(total_frames)
  }
  if (identical(value, "all")) {
    return(total_frames)
  }
  n <- as.integer(value)
  if (!is.finite(n) || n < 1L) {
    stop("Frame count must be a positive integer or 'all'", call. = FALSE)
  }
  min(n, total_frames)
}

load_bold_nifti <- function(path, max_frames = NULL) {
  img <- RNifti::readNifti(path)
  arr <- as.array(img)
  nt <- dim(arr)[4]
  n_keep <- parse_frame_limit(max_frames %||% nt, nt)
  arr <- arr[, , , seq_len(n_keep), drop = FALSE]
  storage.mode(arr) <- "double"
  attr(arr, "pixdim") <- as.numeric(RNifti::pixdim(img)[seq_len(3L)])
  attr(arr, "tr") <- as.numeric(RNifti::pixdim(img)[4])
  attr(arr, "affine") <- RNifti::xform(img)
  list(data = arr, template = img, nt = n_keep)
}

write_bold_nifti <- function(x, path, template) {
  RNifti::writeNifti(x, file = path, template = template, compression = 6L)
  invisible(path)
}

derive_brain_mask <- function(x, quantile_prob = 0.35) {
  ref <- apply(x, c(1L, 2L, 3L), mean)
  positive <- ref[is.finite(ref) & ref > 0]
  thr <- if (length(positive)) {
    unname(stats::quantile(positive, probs = quantile_prob, na.rm = TRUE))
  } else {
    0
  }
  ref > thr
}

select_reference_frame <- function(x, mask) {
  ref <- apply(x, c(1L, 2L, 3L), mean)
  score <- vapply(seq_len(dim(x)[4]), function(i) {
    vol <- x[, , , i]
    if (!any(mask)) {
      return(NA_real_)
    }
    stats::cor(vol[mask], ref[mask])
  }, numeric(1))
  which.max(score)
}

native_global_signal <- function(x, mask) {
  vapply(seq_len(dim(x)[4]), function(i) mean(x[, , , i][mask], na.rm = TRUE), numeric(1))
}

normalized_signal <- function(x) {
  x <- as.numeric(x)
  x <- x - mean(x, na.rm = TRUE)
  scale <- max(abs(x), na.rm = TRUE)
  if (!is.finite(scale) || scale < 1e-8) {
    return(rep(0, length(x)))
  }
  x / scale
}

make_truth_pose <- function(name, nt) {
  tt <- seq(0, 1, length.out = nt)
  deg <- function(x) x * pi / 180

  params <- switch(
    name,
    low_motion = list(
      trans = c(0.9, 0.6, 0.4),
      rot = deg(c(0.9, 0.7, 0.8)),
      drift_amp = 0,
      bias_amp = 0
    ),
    moderate_motion = list(
      trans = c(1.8, 1.2, 0.9),
      rot = deg(c(1.8, 1.5, 1.7)),
      drift_amp = 0,
      bias_amp = 0
    ),
    hard_motion = list(
      trans = c(3.5, 2.5, 1.8),
      rot = deg(c(3.8, 3.2, 3.5)),
      drift_amp = 0,
      bias_amp = 0
    ),
    hard_motion_plus_nuisance = list(
      trans = c(3.5, 2.5, 1.8),
      rot = deg(c(3.8, 3.2, 3.5)),
      drift_amp = 0.06,
      bias_amp = 0.12
    ),
    stop("Unknown scenario: ", name, call. = FALSE)
  )

  pose <- cbind(
    params$trans[1] * sin(2 * pi * tt),
    params$trans[2] * cos(2 * pi * tt + 0.35),
    params$trans[3] * sin(4 * pi * tt + 0.15),
    params$rot[1] * sin(2 * pi * tt),
    params$rot[2] * cos(2 * pi * tt + 0.45),
    params$rot[3] * sin(3 * pi * tt + 0.10)
  )
  pose <- sweep(pose, 2, pose[1, ], FUN = "-")
  params$pose <- pose
  params
}

.make_bias_field <- function(dim3) {
  nx <- dim3[1]
  ny <- dim3[2]
  nz <- dim3[3]

  gx <- array(rep(seq(-1, 1, length.out = nx), each = ny * nz), dim = dim3)
  gy <- array(rep(rep(seq(-1, 1, length.out = ny), each = nz), times = nx), dim = dim3)
  gz <- array(rep(seq(-1, 1, length.out = nz), times = nx * ny), dim = dim3)
  0.55 * gx - 0.35 * gy + 0.20 * gz
}

replicate_reference_run <- function(reference_volume, gains, bias_amp = 0) {
  dim3 <- dim(reference_volume)
  nt <- length(gains)
  out <- array(0, dim = c(dim3, nt))
  bias_field <- .make_bias_field(dim3)
  bias_scale <- sin(seq(0, 2 * pi, length.out = nt))
  for (i in seq_len(nt)) {
    bias <- 1 + bias_amp * bias_field * bias_scale[i]
    out[, , , i] <- reference_volume * gains[i] * bias
  }
  out
}

make_truth_scenario <- function(native_data, mask, scenario_name, truth_frames) {
  nt <- min(truth_frames, dim(native_data)[4])
  ref_idx <- select_reference_frame(native_data[, , , seq_len(nt), drop = FALSE], mask)
  ref_vol <- native_data[, , , ref_idx]
  native_signal <- native_global_signal(native_data[, , , seq_len(nt), drop = FALSE], mask)
  params <- make_truth_pose(scenario_name, nt)
  gains <- 1 + params$drift_amp * normalized_signal(native_signal)

  clean <- replicate_reference_run(ref_vol, gains = gains, bias_amp = params$bias_amp)
  attr(clean, "pixdim") <- attr(native_data, "pixdim")
  attr(clean, "tr") <- attr(native_data, "tr")
  attr(clean, "affine") <- attr(native_data, "affine")

  moved <- apply_motion(clean, params$pose, interp = "lanczos")
  attr(moved, "pixdim") <- attr(native_data, "pixdim")
  attr(moved, "tr") <- attr(native_data, "tr")
  attr(moved, "affine") <- attr(native_data, "affine")

  list(
    input = moved,
    clean = clean,
    truth_pose = params$pose,
    truth_mats = invert_motion(params$pose),
    reference_index = ref_idx,
    drift_amp = params$drift_amp,
    bias_amp = params$bias_amp
  )
}

motion_matrix_colnames <- function() {
  c("frame", sprintf("m%d%d", rep(seq_len(4), each = 4), rep(seq_len(4), times = 4)))
}

write_motion_matrix_csv <- function(mats, path) {
  nt <- dim(mats)[3]
  out <- matrix(NA_real_, nrow = nt, ncol = 17L)
  out[, 1] <- seq_len(nt)
  for (i in seq_len(nt)) {
    out[i, -1] <- as.numeric(mats[, , i])
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- motion_matrix_colnames()
  utils::write.csv(out, path, row.names = FALSE)
  invisible(path)
}

read_motion_matrix_csv <- function(path) {
  dat <- utils::read.csv(path, check.names = FALSE)
  nt <- nrow(dat)
  mats <- array(NA_real_, dim = c(4, 4, nt))
  for (i in seq_len(nt)) {
    mats[, , i] <- matrix(as.numeric(dat[i, -1, drop = TRUE]), nrow = 4, byrow = FALSE)
  }
  mats
}

write_pose_csv <- function(pose, path) {
  out <- data.frame(
    frame = seq_len(nrow(pose)),
    tx = pose[, 1],
    ty = pose[, 2],
    tz = pose[, 3],
    rx = pose[, 4],
    ry = pose[, 5],
    rz = pose[, 6]
  )
  utils::write.csv(out, path, row.names = FALSE)
  invisible(path)
}

project_rotation <- function(R) {
  s <- svd(R)
  R_proj <- s$u %*% t(s$v)
  if (det(R_proj) < 0) {
    s$u[, 3] <- -s$u[, 3]
    R_proj <- s$u %*% t(s$v)
  }
  R_proj
}

motion_mats_to_pose <- function(mats) {
  nt <- dim(mats)[3]
  out <- matrix(0, nrow = nt, ncol = 6)
  for (i in seq_len(nt)) {
    M <- mats[, , i]
    R <- project_rotation(M[1:3, 1:3, drop = FALSE])
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

compute_tsnr <- function(x, mask) {
  means <- apply(x, c(1L, 2L, 3L), mean)
  sds <- apply(x, c(1L, 2L, 3L), stats::sd)
  tsnr <- means / pmax(sds, 1e-6)
  stats::median(tsnr[mask], na.rm = TRUE)
}

compute_blur_proxy <- function(x, mask) {
  ref <- apply(x, c(1L, 2L, 3L), mean)
  dx <- ref[-1, , , drop = FALSE] - ref[-dim(ref)[1], , , drop = FALSE]
  dy <- ref[, -1, , drop = FALSE] - ref[, -dim(ref)[2], , drop = FALSE]
  dz <- ref[, , -1, drop = FALSE] - ref[, , -dim(ref)[3], drop = FALSE]

  mx <- mask[-1, , , drop = FALSE] & mask[-dim(mask)[1], , , drop = FALSE]
  my <- mask[, -1, , drop = FALSE] & mask[, -dim(mask)[2], , drop = FALSE]
  mz <- mask[, , -1, drop = FALSE] & mask[, , -dim(mask)[3], drop = FALSE]

  num <- mean(abs(dx[mx]), na.rm = TRUE) +
    mean(abs(dy[my]), na.rm = TRUE) +
    mean(abs(dz[mz]), na.rm = TRUE)
  as.numeric(num / 3)
}

compute_generic_fail_rate <- function(mats) {
  nt <- dim(mats)[3]
  bad <- logical(nt)
  for (i in seq_len(nt)) {
    M <- mats[, , i]
    R <- M[1:3, 1:3, drop = FALSE]
    bad[i] <- any(!is.finite(M)) || abs(det(project_rotation(R)) - 1) > 0.05
  }
  mean(bad)
}

masked_frame_correlation <- function(a, b, mask) {
  nt <- dim(a)[4]
  vals <- vapply(seq_len(nt), function(i) {
    stats::cor(a[, , , i][mask], b[, , , i][mask])
  }, numeric(1))
  stats::median(vals, na.rm = TRUE)
}

masked_frame_mse <- function(a, b, mask) {
  nt <- dim(a)[4]
  vals <- vapply(seq_len(nt), function(i) {
    mean((a[, , , i][mask] - b[, , , i][mask])^2, na.rm = TRUE)
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

masked_frame_mae <- function(a, b, mask) {
  nt <- dim(a)[4]
  vals <- vapply(seq_len(nt), function(i) {
    mean(abs(a[, , , i][mask] - b[, , , i][mask]), na.rm = TRUE)
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

masked_values_4d <- function(x, mask) {
  dim3 <- dim(mask)
  nxyz <- prod(dim3)
  nt <- dim(x)[4]
  idx <- which(as.array(mask))
  matrix(as.numeric(x), nrow = nxyz, ncol = nt)[idx, , drop = FALSE]
}

masked_frame_nrmse <- function(a, b, mask) {
  mse <- masked_frame_mse(a, b, mask)
  truth_vals <- as.numeric(masked_values_4d(b, mask))
  scale <- stats::sd(truth_vals, na.rm = TRUE)
  if (!is.finite(scale) || scale < 1e-8) {
    scale <- mean(abs(truth_vals), na.rm = TRUE)
  }
  sqrt(mse) / max(scale, 1e-8)
}

masked_residual_dvars <- function(a, b, mask) {
  residual <- a - b
  vals <- dvars(residual, mask = mask, robust = FALSE)
  mean(vals, na.rm = TRUE)
}

masked_residual_dvars_ratio <- function(a, b, mask) {
  resid <- masked_residual_dvars(a, b, mask)
  truth_vals <- as.numeric(masked_values_4d(b, mask))
  scale <- stats::sd(truth_vals, na.rm = TRUE)
  if (!is.finite(scale) || scale < 1e-8) {
    scale <- mean(abs(truth_vals), na.rm = TRUE)
  }
  resid / max(scale, 1e-8)
}

masked_temporal_corr <- function(a, b, mask) {
  amat <- masked_values_4d(a, mask)
  bmat <- masked_values_4d(b, mask)
  if (!nrow(amat)) {
    return(NA_real_)
  }
  vals <- vapply(seq_len(nrow(amat)), function(i) {
    if (stats::sd(amat[i, ], na.rm = TRUE) <= 1e-8 || stats::sd(bmat[i, ], na.rm = TRUE) <= 1e-8) {
      return(NA_real_)
    }
    stats::cor(amat[i, ], bmat[i, ], use = "complete.obs")
  }, numeric(1))
  stats::median(vals, na.rm = TRUE)
}

masked_edge_ratio_to_truth <- function(corrected, clean, mask) {
  clean_edge <- compute_blur_proxy(clean, mask)
  corrected_edge <- compute_blur_proxy(corrected, mask)
  corrected_edge / max(clean_edge, 1e-8)
}

masked_highfreq_residual_ratio <- function(corrected, clean, mask) {
  residual_edge <- compute_blur_proxy(corrected - clean, mask)
  clean_edge <- compute_blur_proxy(clean, mask)
  residual_edge / max(clean_edge, 1e-8)
}

erode_mask_once <- function(mask) {
  mask <- as.array(mask)
  storage.mode(mask) <- "logical"
  out <- array(FALSE, dim = dim(mask))
  if (any(dim(mask) < 3L)) {
    return(out)
  }
  core <- mask[2:(dim(mask)[1] - 1L), 2:(dim(mask)[2] - 1L), 2:(dim(mask)[3] - 1L), drop = FALSE] &
    mask[1:(dim(mask)[1] - 2L), 2:(dim(mask)[2] - 1L), 2:(dim(mask)[3] - 1L), drop = FALSE] &
    mask[3:dim(mask)[1], 2:(dim(mask)[2] - 1L), 2:(dim(mask)[3] - 1L), drop = FALSE] &
    mask[2:(dim(mask)[1] - 1L), 1:(dim(mask)[2] - 2L), 2:(dim(mask)[3] - 1L), drop = FALSE] &
    mask[2:(dim(mask)[1] - 1L), 3:dim(mask)[2], 2:(dim(mask)[3] - 1L), drop = FALSE] &
    mask[2:(dim(mask)[1] - 1L), 2:(dim(mask)[2] - 1L), 1:(dim(mask)[3] - 2L), drop = FALSE] &
    mask[2:(dim(mask)[1] - 1L), 2:(dim(mask)[2] - 1L), 3:dim(mask)[3], drop = FALSE]
  out[2:(dim(mask)[1] - 1L), 2:(dim(mask)[2] - 1L), 2:(dim(mask)[3] - 1L)] <- core
  out
}

masked_frame_nrmse_or_na <- function(a, b, mask) {
  if (!any(mask, na.rm = TRUE)) {
    return(NA_real_)
  }
  masked_frame_nrmse(a, b, mask)
}

truth_image_fidelity_loss <- function(corr_to_truth, nrmse_to_truth, edge_ratio_to_truth, temporal_corr_to_truth = NA_real_) {
  parts <- c(
    nrmse_to_truth,
    abs(edge_ratio_to_truth - 1),
    max(0, 1 - corr_to_truth)
  )
  if (is.finite(temporal_corr_to_truth)) {
    parts <- c(parts, max(0, 1 - temporal_corr_to_truth))
  }
  sum(parts[is.finite(parts)])
}

motion_matrix_error_summary <- function(motion_mats, truth) {
  mats <- motion_matrices(motion_mats)
  ref <- motion_matrices(truth)
  nt <- dim(mats)[3]
  if (dim(ref)[3] != nt) {
    stop("truth must contain the same number of frames as motion", call. = FALSE)
  }

  pose <- motion_mats_to_pose(mats)
  truth_pose <- motion_mats_to_pose(ref)
  err <- pose - truth_pose

  rel_frob <- numeric(nt)
  rel_rot_deg <- numeric(nt)
  det_error <- numeric(nt)
  for (i in seq_len(nt)) {
    dmat <- mats[, , i] %*% solve(ref[, , i])
    rel_frob[i] <- sqrt(sum((dmat - diag(4))^2))
    R <- dmat[1:3, 1:3, drop = FALSE]
    angle <- acos(max(-1, min(1, (sum(diag(R)) - 1) / 2)))
    rel_rot_deg[i] <- angle * 180 / pi
    det_error[i] <- abs(det(mats[1:3, 1:3, i]) - 1)
  }

  list(
    trans_abs_median = stats::median(sqrt(rowSums(err[, 1:3, drop = FALSE]^2)), na.rm = TRUE),
    trans_abs_p95 = unname(stats::quantile(sqrt(rowSums(err[, 1:3, drop = FALSE]^2)), probs = 0.95, na.rm = TRUE)),
    tx_mae = mean(abs(err[, 1]), na.rm = TRUE),
    ty_mae = mean(abs(err[, 2]), na.rm = TRUE),
    tz_mae = mean(abs(err[, 3]), na.rm = TRUE),
    rot_abs_deg_median = stats::median(sqrt(rowSums(err[, 4:6, drop = FALSE]^2)) * 180 / pi, na.rm = TRUE),
    rot_abs_deg_p95 = unname(stats::quantile(sqrt(rowSums(err[, 4:6, drop = FALSE]^2)) * 180 / pi, probs = 0.95, na.rm = TRUE)),
    rx_mae_deg = mean(abs(err[, 4]) * 180 / pi, na.rm = TRUE),
    ry_mae_deg = mean(abs(err[, 5]) * 180 / pi, na.rm = TRUE),
    rz_mae_deg = mean(abs(err[, 6]) * 180 / pi, na.rm = TRUE),
    rel_frob_median = stats::median(rel_frob, na.rm = TRUE),
    rel_frob_p95 = unname(stats::quantile(rel_frob, probs = 0.95, na.rm = TRUE)),
    rel_rot_deg_median = stats::median(rel_rot_deg, na.rm = TRUE),
    rel_rot_deg_p95 = unname(stats::quantile(rel_rot_deg, probs = 0.95, na.rm = TRUE)),
    affine_det_error_max = max(det_error, na.rm = TRUE)
  )
}

motion_dynamics_error_summary <- function(motion_mats, truth) {
  mats <- motion_matrices(motion_mats)
  ref <- motion_matrices(truth)
  nt <- dim(mats)[3]
  if (dim(ref)[3] != nt) {
    stop("truth must contain the same number of frames as motion", call. = FALSE)
  }

  # Keep fd_corr_to_truth aligned with the FD name: correlate frame-to-frame
  # displacement, while fd_abs_corr_to_truth preserves the absolute-pose trace.
  fd_step_est <- framewise_disp_from_transform(mats)
  fd_step_truth <- framewise_disp_from_transform(ref)
  fd_corr <- if (stats::sd(fd_step_est, na.rm = TRUE) > 1e-8 && stats::sd(fd_step_truth, na.rm = TRUE) > 1e-8) {
    stats::cor(fd_step_est, fd_step_truth, use = "complete.obs")
  } else {
    NA_real_
  }
  fd_abs_est <- fd_from_transform(mats)
  fd_abs_truth <- fd_from_transform(ref)
  fd_abs_corr <- if (stats::sd(fd_abs_est, na.rm = TRUE) > 1e-8 && stats::sd(fd_abs_truth, na.rm = TRUE) > 1e-8) {
    stats::cor(fd_abs_est, fd_abs_truth, use = "complete.obs")
  } else {
    NA_real_
  }

  pose <- motion_mats_to_pose(mats)
  truth_pose <- motion_mats_to_pose(ref)
  if (nt < 2L) {
    pose_jitter_error <- NA_real_
  } else {
    dpose <- diff(pose)
    dtruth <- diff(truth_pose)
    dpose[, 4:6] <- dpose[, 4:6] * 50
    dtruth[, 4:6] <- dtruth[, 4:6] * 50
    pose_jitter_error <- sqrt(mean((dpose - dtruth)^2, na.rm = TRUE))
  }

  list(
    fd_corr_to_truth = fd_corr,
    fd_abs_corr_to_truth = fd_abs_corr,
    pose_jitter_error = pose_jitter_error
  )
}

benchmark_metric_names <- function() {
  c(
    "fit_fail_rate",
    "fd_mean",
    "fd_p95",
    "dvars_before_mean",
    "dvars_before_p95",
    "dvars_after_mean",
    "dvars_after_p95",
    "dvars_ratio",
    "tsnr_before",
    "tsnr_after",
    "blur_proxy_before",
    "blur_proxy_after",
    "tsnr_ratio",
    "blur_proxy_ratio",
    "disp_median",
    "disp_p95",
    "disp_max",
    "fd_error",
    "trans_rmse",
    "rot_rmse",
    "trans_abs_median",
    "trans_abs_p95",
    "tx_mae",
    "ty_mae",
    "tz_mae",
    "rot_abs_deg_median",
    "rot_abs_deg_p95",
    "rx_mae_deg",
    "ry_mae_deg",
    "rz_mae_deg",
    "rel_frob_median",
    "rel_frob_p95",
    "rel_rot_deg_median",
    "rel_rot_deg_p95",
    "affine_det_error_max",
    "fd_corr_to_truth",
    "fd_abs_corr_to_truth",
    "pose_jitter_error",
    "corr_to_truth",
    "mse_to_truth",
    "mae_to_truth",
    "nrmse_to_truth",
    "residual_dvars_to_truth",
    "residual_dvars_ratio_to_truth",
    "temporal_corr_to_truth",
    "edge_ratio_to_truth",
    "edge_error_to_truth",
    "highfreq_residual_ratio",
    "nrmse_inner_to_truth",
    "nrmse_boundary_to_truth",
    "truth_image_fidelity_loss"
  )
}

empty_benchmark_metrics <- function() {
  vals <- stats::setNames(rep(NA_real_, length(benchmark_metric_names())), benchmark_metric_names())
  as.data.frame(as.list(vals), stringsAsFactors = FALSE)
}

benchmark_scorecard_specs <- function() {
  rows <- list(
    runtime = data.frame(
      metric = c("elapsed_sec", "estimate_sec", "temp_bytes_written"),
      family = "runtime",
      role = c("claim_efficiency", "claim_efficiency", "diagnostic_resource"),
      stringsAsFactors = FALSE
    ),
    displacement = data.frame(
      metric = c("disp_median", "disp_p95", "disp_max", "fd_error"),
      family = "truth_displacement",
      role = "claim_truth",
      stringsAsFactors = FALSE
    ),
    pose_claim = data.frame(
      metric = c("trans_rmse", "rot_rmse", "trans_abs_median", "trans_abs_p95",
                 "rot_abs_deg_median", "rot_abs_deg_p95"),
      family = "truth_pose",
      role = "claim_truth",
      stringsAsFactors = FALSE
    ),
    pose_axis = data.frame(
      metric = c("tx_mae", "ty_mae", "tz_mae", "rx_mae_deg", "ry_mae_deg", "rz_mae_deg"),
      family = "truth_pose",
      role = "diagnostic_pose_axis",
      stringsAsFactors = FALSE
    ),
    matrix = data.frame(
      metric = c("rel_frob_median", "rel_frob_p95", "rel_rot_deg_median",
                 "rel_rot_deg_p95", "affine_det_error_max"),
      family = "truth_matrix",
      role = "claim_truth",
      stringsAsFactors = FALSE
    ),
    dynamics = data.frame(
      metric = c("fd_corr_to_truth", "pose_jitter_error"),
      family = "truth_dynamics",
      role = "claim_truth",
      stringsAsFactors = FALSE
    ),
    dynamics_diag = data.frame(
      metric = "fd_abs_corr_to_truth",
      family = "truth_dynamics",
      role = "diagnostic_dynamics",
      stringsAsFactors = FALSE
    ),
    image_claim = data.frame(
      metric = c("mse_to_truth", "mae_to_truth", "nrmse_to_truth",
                 "residual_dvars_ratio_to_truth", "corr_to_truth",
                 "temporal_corr_to_truth", "edge_ratio_to_truth",
                 "edge_error_to_truth", "truth_image_fidelity_loss"),
      family = "truth_image",
      role = "claim_truth",
      stringsAsFactors = FALSE
    ),
    residual_diag = data.frame(
      metric = c("residual_dvars_to_truth", "highfreq_residual_ratio"),
      family = "truth_image",
      role = c("diagnostic_residual", "diagnostic_ringing"),
      stringsAsFactors = FALSE
    ),
    smoothing = data.frame(
      metric = c("dvars_after_mean", "dvars_after_p95", "dvars_ratio",
                 "tsnr_ratio", "blur_proxy_ratio"),
      family = "qc",
      role = "diagnostic_smoothing",
      stringsAsFactors = FALSE
    ),
    region_diag = data.frame(
      metric = c("nrmse_inner_to_truth", "nrmse_boundary_to_truth"),
      family = "truth_image",
      role = "diagnostic_region",
      stringsAsFactors = FALSE
    )
  )

  specs <- do.call(rbind, rows)
  specs$direction <- ifelse(
    specs$metric %in% c("corr_to_truth", "temporal_corr_to_truth", "tsnr_ratio", "fd_corr_to_truth", "fd_abs_corr_to_truth"),
    "higher",
    ifelse(specs$metric == "edge_ratio_to_truth", "target", "lower")
  )
  specs$target <- ifelse(specs$metric == "edge_ratio_to_truth", 1, NA_real_)
  specs$tolerance <- 1e-8
  specs$tolerance[specs$metric == "fd_corr_to_truth"] <- 0.02
  specs$tolerance[specs$metric == "fd_abs_corr_to_truth"] <- 0.02
  specs$tolerance[specs$metric == "residual_dvars_ratio_to_truth"] <- 0.02
  specs <- specs[, c("metric", "direction", "target", "tolerance", "family", "role")]
  specs$claim_metric <- startsWith(specs$role, "claim_")
  specs
}

scorecard_metric_result <- function(primary, baseline, direction, target = NA_real_, tolerance = 1e-8) {
  if (!is.finite(primary) || !is.finite(baseline)) {
    return(list(
      primary_margin = NA_real_,
      primary_ratio = NA_real_,
      primary_wins = NA,
      comparison = "unavailable"
    ))
  }

  if (identical(direction, "lower")) {
    margin <- baseline - primary
    ratio <- primary / max(abs(baseline), tolerance)
    wins <- primary < baseline - tolerance
  } else if (identical(direction, "higher")) {
    margin <- primary - baseline
    ratio <- primary / max(abs(baseline), tolerance)
    wins <- primary > baseline + tolerance
  } else if (identical(direction, "target")) {
    primary_error <- abs(primary - target)
    baseline_error <- abs(baseline - target)
    margin <- baseline_error - primary_error
    ratio <- primary_error / max(baseline_error, tolerance)
    wins <- primary_error < baseline_error - tolerance
  } else {
    stop("Unknown scorecard direction: ", direction, call. = FALSE)
  }

  list(
    primary_margin = margin,
    primary_ratio = ratio,
    primary_wins = wins,
    comparison = if (isTRUE(wins)) "win" else if (abs(margin) <= tolerance) "tie" else "loss"
  )
}

make_benchmark_scorecard <- function(raw, primary_method = "volregger", truth_only = TRUE) {
  specs <- benchmark_scorecard_specs()
  if (isTRUE(truth_only)) {
    raw <- raw[raw$dataset_mode == "truth_injected", , drop = FALSE]
  }
  scenarios <- unique(raw$scenario)
  rows <- list()
  row_idx <- 1L

  for (scenario_name in scenarios) {
    primary <- raw[raw$scenario == scenario_name & raw$method == primary_method & raw$status == "ok", , drop = FALSE]
    if (nrow(primary) != 1L) {
      next
    }
    baselines <- raw[raw$scenario == scenario_name & raw$method != primary_method & raw$status == "ok", , drop = FALSE]
    if (!nrow(baselines)) {
      next
    }

    for (baseline_idx in seq_len(nrow(baselines))) {
      baseline <- baselines[baseline_idx, , drop = FALSE]
      for (spec_idx in seq_len(nrow(specs))) {
        metric <- specs$metric[spec_idx]
        if (!metric %in% names(raw)) {
          next
        }
        result <- scorecard_metric_result(
          primary = primary[[metric]],
          baseline = baseline[[metric]],
          direction = specs$direction[spec_idx],
          target = specs$target[spec_idx],
          tolerance = specs$tolerance[spec_idx]
        )
        rows[[row_idx]] <- data.frame(
          scenario = scenario_name,
          dataset_mode = primary$dataset_mode,
          baseline = baseline$method,
          metric = metric,
          family = specs$family[spec_idx],
          role = specs$role[spec_idx],
          claim_metric = specs$claim_metric[spec_idx],
          direction = specs$direction[spec_idx],
          target = specs$target[spec_idx],
          primary_value = primary[[metric]],
          baseline_value = baseline[[metric]],
          primary_margin = result$primary_margin,
          primary_ratio = result$primary_ratio,
          primary_wins = result$primary_wins,
          comparison = result$comparison,
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }
  }

  if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame()
  }
}

summarize_benchmark_scorecard <- function(scorecard, claim_only = FALSE) {
  if (!nrow(scorecard)) {
    return(data.frame())
  }
  if (isTRUE(claim_only)) {
    scorecard <- scorecard[scorecard$claim_metric %in% TRUE, , drop = FALSE]
    if (!nrow(scorecard)) {
      return(data.frame())
    }
  }
  aggregate(
    x = list(
      metrics_available = scorecard$comparison != "unavailable",
      wins = scorecard$comparison == "win",
      losses = scorecard$comparison == "loss",
      ties = scorecard$comparison == "tie"
    ),
    by = list(
      scenario = scorecard$scenario,
      dataset_mode = scorecard$dataset_mode,
      baseline = scorecard$baseline,
      family = scorecard$family,
      role = scorecard$role,
      claim_metric = scorecard$claim_metric
    ),
    FUN = sum,
    na.rm = TRUE
  )
}

no_regression_gate_specs <- function() {
  data.frame(
    metric = c(
      "disp_p95",
      "fd_error",
      "fd_corr_to_truth",
      "pose_jitter_error",
      "nrmse_to_truth",
      "edge_error_to_truth",
      "truth_image_fidelity_loss",
      "tsnr_ratio"
    ),
    direction = c("lower", "lower", "higher", "lower", "lower", "lower", "lower", "higher"),
    tolerance = c(1e-8, 1e-8, 0.02, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8),
    role = c(
      "truth_displacement",
      "truth_dynamics",
      "truth_dynamics",
      "truth_dynamics",
      "truth_image",
      "truth_image",
      "truth_image",
      "diagnostic_tsnr"
    ),
    stringsAsFactors = FALSE
  )
}

make_no_regression_gate <- function(raw, primary_method = "volregger", truth_only = TRUE) {
  specs <- no_regression_gate_specs()
  if (isTRUE(truth_only)) {
    raw <- raw[raw$dataset_mode == "truth_injected", , drop = FALSE]
  }
  scenarios <- unique(raw$scenario)
  rows <- list()
  row_idx <- 1L

  for (scenario_name in scenarios) {
    primary <- raw[raw$scenario == scenario_name & raw$method == primary_method & raw$status == "ok", , drop = FALSE]
    if (nrow(primary) != 1L) {
      next
    }
    baselines <- raw[raw$scenario == scenario_name & raw$method != primary_method & raw$status == "ok", , drop = FALSE]
    if (!nrow(baselines)) {
      next
    }

    for (baseline_idx in seq_len(nrow(baselines))) {
      baseline <- baselines[baseline_idx, , drop = FALSE]
      for (spec_idx in seq_len(nrow(specs))) {
        metric <- specs$metric[spec_idx]
        if (!metric %in% names(raw)) {
          next
        }
        result <- scorecard_metric_result(
          primary = primary[[metric]],
          baseline = baseline[[metric]],
          direction = specs$direction[spec_idx],
          tolerance = specs$tolerance[spec_idx]
        )
        rows[[row_idx]] <- data.frame(
          scenario = scenario_name,
          dataset_mode = primary$dataset_mode,
          baseline = baseline$method,
          metric = metric,
          role = specs$role[spec_idx],
          direction = specs$direction[spec_idx],
          tolerance = specs$tolerance[spec_idx],
          primary_value = primary[[metric]],
          baseline_value = baseline[[metric]],
          primary_margin = result$primary_margin,
          primary_ratio = result$primary_ratio,
          comparison = result$comparison,
          no_regression_pass = result$comparison %in% c("win", "tie"),
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }
  }

  if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame()
  }
}

summarize_no_regression_gate <- function(gate) {
  if (!nrow(gate)) {
    return(data.frame())
  }
  aggregate(
    x = list(
      metrics_available = gate$comparison != "unavailable",
      wins = gate$comparison == "win",
      regressions = gate$comparison == "loss",
      ties = gate$comparison == "tie",
      passes = gate$no_regression_pass %in% TRUE
    ),
    by = list(
      scenario = gate$scenario,
      dataset_mode = gate$dataset_mode,
      baseline = gate$baseline
    ),
    FUN = sum,
    na.rm = TRUE
  )
}

make_tsnr_context <- function(raw, primary_method = "volregger", truth_only = TRUE) {
  if (isTRUE(truth_only)) {
    raw <- raw[raw$dataset_mode == "truth_injected", , drop = FALSE]
  }
  scenarios <- unique(raw$scenario)
  rows <- list()
  row_idx <- 1L

  metric_value <- function(row, metric) {
    if (!metric %in% names(row)) {
      return(NA_real_)
    }
    row[[metric]]
  }

  for (scenario_name in scenarios) {
    primary <- raw[raw$scenario == scenario_name & raw$method == primary_method & raw$status == "ok", , drop = FALSE]
    if (nrow(primary) != 1L) {
      next
    }
    baselines <- raw[raw$scenario == scenario_name & raw$method != primary_method & raw$status == "ok", , drop = FALSE]
    if (!nrow(baselines)) {
      next
    }

    for (baseline_idx in seq_len(nrow(baselines))) {
      baseline <- baselines[baseline_idx, , drop = FALSE]
      rows[[row_idx]] <- data.frame(
        scenario = scenario_name,
        dataset_mode = primary$dataset_mode,
        baseline = baseline$method,
        tsnr_primary = metric_value(primary, "tsnr_ratio"),
        tsnr_baseline = metric_value(baseline, "tsnr_ratio"),
        tsnr_delta = metric_value(primary, "tsnr_ratio") - metric_value(baseline, "tsnr_ratio"),
        blur_primary = metric_value(primary, "blur_proxy_ratio"),
        blur_baseline = metric_value(baseline, "blur_proxy_ratio"),
        blur_delta = metric_value(primary, "blur_proxy_ratio") - metric_value(baseline, "blur_proxy_ratio"),
        disp_p95_margin = metric_value(baseline, "disp_p95") - metric_value(primary, "disp_p95"),
        fd_error_margin = metric_value(baseline, "fd_error") - metric_value(primary, "fd_error"),
        fd_corr_margin = metric_value(primary, "fd_corr_to_truth") - metric_value(baseline, "fd_corr_to_truth"),
        pose_jitter_margin = metric_value(baseline, "pose_jitter_error") - metric_value(primary, "pose_jitter_error"),
        nrmse_margin = metric_value(baseline, "nrmse_to_truth") - metric_value(primary, "nrmse_to_truth"),
        edge_error_margin = metric_value(baseline, "edge_error_to_truth") - metric_value(primary, "edge_error_to_truth"),
        fidelity_loss_margin = metric_value(baseline, "truth_image_fidelity_loss") - metric_value(primary, "truth_image_fidelity_loss"),
        tsnr_loss = metric_value(primary, "tsnr_ratio") < metric_value(baseline, "tsnr_ratio") - 1e-8,
        truth_fidelity_win = metric_value(primary, "truth_image_fidelity_loss") < metric_value(baseline, "truth_image_fidelity_loss") - 1e-8,
        edge_win = metric_value(primary, "edge_error_to_truth") < metric_value(baseline, "edge_error_to_truth") - 1e-8,
        displacement_win = metric_value(primary, "disp_p95") < metric_value(baseline, "disp_p95") - 1e-8,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
    }
  }

  if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame()
  }
}

summarize_interpolation_guardrail <- function(raw, baseline_interp = "lanczos") {
  required <- c("interp", "scenario", "dataset_mode", "method", "status")
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop("raw is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  metrics <- c(
    "elapsed_sec",
    "estimate_sec",
    "apply_sec",
    "disp_p95",
    "fd_error",
    "fd_corr_to_truth",
    "pose_jitter_error",
    "nrmse_to_truth",
    "edge_error_to_truth",
    "highfreq_residual_ratio",
    "truth_image_fidelity_loss",
    "tsnr_ratio",
    "blur_proxy_ratio"
  )
  metrics <- intersect(metrics, names(raw))
  rows <- raw[
    raw$method == "volregger" &
      raw$status == "ok" &
      raw$dataset_mode == "truth_injected",
    c("interp", "scenario", "dataset_mode", metrics),
    drop = FALSE
  ]
  if (!nrow(rows)) {
    return(data.frame())
  }

  out <- rows
  for (metric in metrics) {
    out[[paste0(metric, "_vs_", baseline_interp)]] <- NA_real_
  }

  for (scenario_name in unique(rows$scenario)) {
    base <- rows[rows$scenario == scenario_name & rows$interp == baseline_interp, , drop = FALSE]
    if (nrow(base) != 1L) {
      next
    }
    idx <- which(out$scenario == scenario_name)
    for (metric in metrics) {
      out[idx, paste0(metric, "_vs_", baseline_interp)] <- out[idx, metric] - base[[metric]]
    }
  }
  out
}

compute_benchmark_metrics <- function(input_run, corrected, mask, motion_mats, truth = NULL, clean = NULL) {
  fd <- fd_from_transform(motion_mats)
  dvars_before <- dvars(input_run, robust = FALSE)
  dvars_after <- dvars(corrected, robust = FALSE)

  out <- list(
    fit_fail_rate = compute_generic_fail_rate(motion_mats),
    fd_mean = mean(fd, na.rm = TRUE),
    fd_p95 = unname(stats::quantile(fd, probs = 0.95, na.rm = TRUE)),
    dvars_before_mean = mean(dvars_before, na.rm = TRUE),
    dvars_before_p95 = unname(stats::quantile(dvars_before, probs = 0.95, na.rm = TRUE)),
    dvars_after_mean = mean(dvars_after, na.rm = TRUE),
    dvars_after_p95 = unname(stats::quantile(dvars_after, probs = 0.95, na.rm = TRUE)),
    dvars_ratio = mean(dvars_after, na.rm = TRUE) / max(mean(dvars_before, na.rm = TRUE), 1e-8),
    tsnr_before = compute_tsnr(input_run, mask),
    tsnr_after = compute_tsnr(corrected, mask),
    blur_proxy_before = compute_blur_proxy(input_run, mask),
    blur_proxy_after = compute_blur_proxy(corrected, mask)
  )

  out$tsnr_ratio <- out$tsnr_after / max(out$tsnr_before, 1e-8)
  out$blur_proxy_ratio <- out$blur_proxy_after / max(out$blur_proxy_before, 1e-8)

  if (!is.null(truth)) {
    truth_disp <- motion_displacement_summary(
      mask = mask,
      motion = motion_mats,
      truth = truth,
      pixdim = attr(input_run, "pixdim")
    )
    est_pose <- motion_mats_to_pose(motion_mats)
    truth_pose <- motion_mats_to_pose(truth)
    err <- est_pose - truth_pose

    out$disp_median <- stats::median(truth_disp$disp_median, na.rm = TRUE)
    out$disp_p95 <- stats::median(truth_disp$disp_p95, na.rm = TRUE)
    out$disp_max <- stats::median(truth_disp$disp_max, na.rm = TRUE)
    out$fd_error <- stats::median(truth_disp$fd_error, na.rm = TRUE)
    out$trans_rmse <- sqrt(mean(err[, 1:3, drop = FALSE]^2, na.rm = TRUE))
    out$rot_rmse <- sqrt(mean(err[, 4:6, drop = FALSE]^2, na.rm = TRUE))

    matrix_err <- motion_matrix_error_summary(motion_mats, truth)
    out <- c(out, matrix_err)
    dynamics_err <- motion_dynamics_error_summary(motion_mats, truth)
    out <- c(out, dynamics_err)
  } else {
    out$disp_median <- NA_real_
    out$disp_p95 <- NA_real_
    out$disp_max <- NA_real_
    out$fd_error <- NA_real_
    out$trans_rmse <- NA_real_
    out$rot_rmse <- NA_real_
    out$trans_abs_median <- NA_real_
    out$trans_abs_p95 <- NA_real_
    out$tx_mae <- NA_real_
    out$ty_mae <- NA_real_
    out$tz_mae <- NA_real_
    out$rot_abs_deg_median <- NA_real_
    out$rot_abs_deg_p95 <- NA_real_
    out$rx_mae_deg <- NA_real_
    out$ry_mae_deg <- NA_real_
    out$rz_mae_deg <- NA_real_
    out$rel_frob_median <- NA_real_
    out$rel_frob_p95 <- NA_real_
    out$rel_rot_deg_median <- NA_real_
    out$rel_rot_deg_p95 <- NA_real_
    out$affine_det_error_max <- NA_real_
    out$fd_corr_to_truth <- NA_real_
    out$fd_abs_corr_to_truth <- NA_real_
    out$pose_jitter_error <- NA_real_
  }

  if (!is.null(clean)) {
    out$corr_to_truth <- masked_frame_correlation(corrected, clean, mask)
    out$mse_to_truth <- masked_frame_mse(corrected, clean, mask)
    out$mae_to_truth <- masked_frame_mae(corrected, clean, mask)
    out$nrmse_to_truth <- masked_frame_nrmse(corrected, clean, mask)
    out$residual_dvars_to_truth <- masked_residual_dvars(corrected, clean, mask)
    out$residual_dvars_ratio_to_truth <- masked_residual_dvars_ratio(corrected, clean, mask)
    out$temporal_corr_to_truth <- masked_temporal_corr(corrected, clean, mask)
    out$edge_ratio_to_truth <- masked_edge_ratio_to_truth(corrected, clean, mask)
    out$edge_error_to_truth <- abs(out$edge_ratio_to_truth - 1)
    out$highfreq_residual_ratio <- masked_highfreq_residual_ratio(corrected, clean, mask)
    inner_mask <- erode_mask_once(mask)
    boundary_mask <- mask & !inner_mask
    out$nrmse_inner_to_truth <- masked_frame_nrmse_or_na(corrected, clean, inner_mask)
    out$nrmse_boundary_to_truth <- masked_frame_nrmse_or_na(corrected, clean, boundary_mask)
    out$truth_image_fidelity_loss <- truth_image_fidelity_loss(
      corr_to_truth = out$corr_to_truth,
      nrmse_to_truth = out$nrmse_to_truth,
      edge_ratio_to_truth = out$edge_ratio_to_truth,
      temporal_corr_to_truth = out$temporal_corr_to_truth
    )
  } else {
    out$corr_to_truth <- NA_real_
    out$mse_to_truth <- NA_real_
    out$mae_to_truth <- NA_real_
    out$nrmse_to_truth <- NA_real_
    out$residual_dvars_to_truth <- NA_real_
    out$residual_dvars_ratio_to_truth <- NA_real_
    out$temporal_corr_to_truth <- NA_real_
    out$edge_ratio_to_truth <- NA_real_
    out$edge_error_to_truth <- NA_real_
    out$highfreq_residual_ratio <- NA_real_
    out$nrmse_inner_to_truth <- NA_real_
    out$nrmse_boundary_to_truth <- NA_real_
    out$truth_image_fidelity_loss <- NA_real_
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

mask_world_coords <- function(mask, pixdim) {
  mask <- as.array(mask)
  storage.mode(mask) <- "logical"
  ijk <- which(mask, arr.ind = TRUE)
  if (!nrow(ijk)) {
    stop("mask contains no TRUE voxels", call. = FALSE)
  }
  dim3 <- dim(mask)
  center <- 0.5 * (dim3 - 1)
  cbind(
    (ijk[, 1] - 1 - center[1]) * pixdim[1],
    (ijk[, 2] - 1 - center[2]) * pixdim[2],
    (ijk[, 3] - 1 - center[3]) * pixdim[3]
  )
}

apply_motion_coords <- function(mat, coords) {
  sweep(coords %*% t(mat[1:3, 1:3, drop = FALSE]), 2L, mat[1:3, 4], FUN = "+")
}

framewise_disp_from_transform <- function(motion, radius = 50) {
  mats <- motion_matrices(motion)
  nt <- dim(mats)[3]
  if (nt <= 1L) {
    return(rep(0, nt))
  }
  out <- numeric(nt)
  out[1] <- 0
  out[-1] <- fd_from_transform(
    mats[, , 2:nt, drop = FALSE],
    radius = radius,
    reference = mats[, , 1:(nt - 1L), drop = FALSE]
  )
  out
}

pairwise_solution_comparison <- function(
    mask,
    motion_a,
    motion_b,
    pixdim,
    ref_index = 1L,
    corrected_a = NULL,
    corrected_b = NULL,
    radius = 50) {
  mats_a <- motion_matrices(motion_a)
  mats_b <- motion_matrices(motion_b)
  nt <- dim(mats_a)[3]
  if (dim(mats_b)[3] != nt) {
    stop("motion_a and motion_b must have same frame count", call. = FALSE)
  }
  if (ref_index < 1L || ref_index > nt) {
    stop("ref_index out of range", call. = FALSE)
  }

  coords <- mask_world_coords(mask, pixdim = pixdim)
  ref_a <- apply_motion_coords(mats_a[, , ref_index], coords)
  ref_b <- apply_motion_coords(mats_b[, , ref_index], coords)

  frames <- vector("list", nt)
  for (t in seq_len(nt)) {
    rel_a <- apply_motion_coords(mats_a[, , t], coords) - ref_a
    rel_b <- apply_motion_coords(mats_b[, , t], coords) - ref_b
    disp <- sqrt(rowSums((rel_a - rel_b)^2))
    frames[[t]] <- data.frame(
      frame = t,
      rel_disp_median = stats::median(disp, na.rm = TRUE),
      rel_disp_p95 = unname(stats::quantile(disp, probs = 0.95, na.rm = TRUE)),
      rel_disp_max = max(disp, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  frames <- do.call(rbind, frames)

  fd_a <- framewise_disp_from_transform(mats_a, radius = radius)
  fd_b <- framewise_disp_from_transform(mats_b, radius = radius)
  frames$fd_step_a <- fd_a
  frames$fd_step_b <- fd_b
  frames$fd_step_absdiff <- abs(fd_a - fd_b)

  summary <- data.frame(
    anchor_frame = ref_index,
    rel_disp_median = stats::median(frames$rel_disp_median, na.rm = TRUE),
    rel_disp_p95 = stats::median(frames$rel_disp_p95, na.rm = TRUE),
    rel_disp_max = stats::median(frames$rel_disp_max, na.rm = TRUE),
    fd_step_mean_absdiff = mean(frames$fd_step_absdiff, na.rm = TRUE),
    fd_step_p95_absdiff = unname(stats::quantile(frames$fd_step_absdiff, probs = 0.95, na.rm = TRUE)),
    fd_step_corr = if (sum(is.finite(fd_a) & is.finite(fd_b)) >= 3L) stats::cor(fd_a, fd_b) else NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(corrected_a) && !is.null(corrected_b)) {
    summary$corrected_corr <- masked_frame_correlation(corrected_a, corrected_b, mask)
    summary$corrected_mse <- masked_frame_mse(corrected_a, corrected_b, mask)
  } else {
    summary$corrected_corr <- NA_real_
    summary$corrected_mse <- NA_real_
  }

  list(summary = summary, frames = frames)
}

dir_size_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) {
    return(0)
  }
  sum(file.info(files)$size, na.rm = TRUE)
}

run_profiled_command <- function(command, args, stdout_path, stderr_path) {
  if (!nzchar(command)) {
    writeLines("Command is empty or not found on PATH", con = stderr_path)
    return(list(status = 127L, elapsed_sec = NA_real_, peak_rss_mb = NA_real_))
  }

  start <- proc.time()[["elapsed"]]
  status <- tryCatch(
    system2(command, args = args, stdout = stdout_path, stderr = stderr_path),
    error = function(e) {
      writeLines(conditionMessage(e), con = stderr_path)
      127L
    }
  )
  elapsed <- proc.time()[["elapsed"]] - start

  list(
    status = as.integer(status),
    elapsed_sec = as.numeric(elapsed),
    peak_rss_mb = NA_real_
  )
}

time_expr <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  list(value = value, elapsed_sec = as.numeric(elapsed))
}

default_stage_times <- function(total_sec = NA_real_) {
  data.frame(
    estimate_sec = NA_real_,
    apply_sec = NA_real_,
    report_sec = 0,
    total_sec = as.numeric(total_sec),
    stringsAsFactors = FALSE
  )
}

write_stage_times_csv <- function(stage_times, path) {
  required <- c("estimate_sec", "apply_sec", "report_sec", "total_sec")
  miss <- setdiff(required, names(stage_times))
  if (length(miss)) {
    stop("stage_times is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  utils::write.csv(stage_times[, required, drop = FALSE], path, row.names = FALSE)
  invisible(path)
}

read_stage_times_csv <- function(path, total_sec = NA_real_) {
  if (!file.exists(path)) {
    return(default_stage_times(total_sec = total_sec))
  }
  dat <- utils::read.csv(path, check.names = FALSE)
  if (!nrow(dat)) {
    return(default_stage_times(total_sec = total_sec))
  }
  out <- default_stage_times(total_sec = total_sec)
  for (nm in intersect(names(out), names(dat))) {
    out[[nm]] <- as.numeric(dat[[nm]][1])
  }
  if (!is.finite(out$total_sec) && is.finite(total_sec)) {
    out$total_sec <- as.numeric(total_sec)
  }
  out
}

method_error_message <- function(stderr_path) {
  helper_error <- file.path(dirname(stderr_path), "error.txt")
  if (file.exists(helper_error)) {
    lines <- readLines(helper_error, warn = FALSE)
    if (length(lines)) {
      return(paste(lines, collapse = "\n"))
    }
  }
  if (!file.exists(stderr_path)) {
    return(NA_character_)
  }
  lines <- readLines(stderr_path, warn = FALSE)
  if (!length(lines)) {
    return(NA_character_)
  }
  paste(utils::tail(lines, n = min(10L, length(lines))), collapse = "\n")
}

write_mask_nifti <- function(mask, path, template) {
  RNifti::writeNifti(mask * 1L, file = path, template = template, compression = 6L)
  invisible(path)
}

default_realdata_input <- function() {
  candidate <- file.path("data-raw", "sub-01_task-freerecall_bold.nii.gz")
  if (file.exists(candidate)) candidate else stop("No default real-data NIfTI found", call. = FALSE)
}
