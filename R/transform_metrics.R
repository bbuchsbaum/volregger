#' Convert motion representations to homogeneous transform matrices
#'
#' @param motion Motion represented as one of:
#'   a length-6 pose vector, an `nt x 6` pose matrix, a `4 x 4` matrix,
#'   a `4 x 4 x nt` array, or a `16 x nt` matrix as returned by `estimate_motion()`.
#' @return A `4 x 4 x nt` array of homogeneous transforms mapping moving
#'   world coordinates into fixed world coordinates.
#' @export
motion_matrices <- function(motion) {
  .motion_matrices_array(motion)
}

#' Invert motion transforms
#'
#' @param motion Motion represented as one of:
#'   a length-6 pose vector, an `nt x 6` pose matrix, a `4 x 4` matrix,
#'   a `4 x 4 x nt` array, or a `16 x nt` matrix as returned by `estimate_motion()`.
#' @return A `4 x 4 x nt` array containing the inverse transform for each frame.
#' @export
invert_motion <- function(motion) {
  mats <- .motion_matrices_array(motion)
  nt <- dim(mats)[3]
  out <- array(0, dim = dim(mats))
  for (t in seq_len(nt)) {
    out[, , t] <- solve(mats[, , t, drop = FALSE][, , 1])
  }
  out
}

#' Compute transform displacement on a radius-50 mm head surrogate
#'
#' @param motion Motion represented as one of:
#'   a length-6 pose vector, an `nt x 6` pose matrix, a `4 x 4` matrix,
#'   a `4 x 4 x nt` array, or a `16 x nt` matrix as returned by `estimate_motion()`.
#' @param radius Sphere radius in millimeters used for displacement evaluation.
#' @param reference Optional reference transform(s). If omitted, displacements
#'   are computed against the identity transform.
#' @return Numeric vector of mean point displacements in millimeters.
#' @export
fd_from_transform <- function(motion, radius = 50, reference = NULL) {
  mats <- .motion_matrices_array(motion)
  ref <- .reference_matrices(reference, n_frames = dim(mats)[3])

  fd_coords <- as.matrix(expand.grid(c(radius, -radius), c(radius, -radius), c(radius, -radius)))
  out <- numeric(dim(mats)[3])

  for (t in seq_len(dim(mats)[3])) {
    mapped <- .apply_transform(mats[, , t], fd_coords)
    target <- if (is.null(ref)) fd_coords else .apply_transform(ref[, , t], fd_coords)
    out[t] <- mean(sqrt(rowSums((mapped - target)^2)))
  }

  out
}

#' Compute masked transform displacements
#'
#' @param mask 3D logical mask defining evaluation voxels.
#' @param motion A single transform represented as a length-6 pose vector,
#'   a `4 x 4` matrix, a `4 x 4 x 1` array, or a `16 x 1` matrix.
#' @param reference Optional single reference transform. If omitted, the
#'   identity transform is used.
#' @param pixdim Length-3 voxel size vector in millimeters.
#' @return Numeric vector of Euclidean displacements in millimeters for voxels
#'   inside `mask`.
#' @export
displacements_within_mask <- function(mask, motion, reference = NULL, pixdim = c(1, 1, 1)) {
  coords <- .mask_world_coords(mask, pixdim = pixdim)
  mat <- .single_motion_matrix(motion, arg = "motion")
  ref <- if (is.null(reference)) NULL else .single_motion_matrix(reference, arg = "reference")

  mapped <- .apply_transform(mat, coords)
  target <- if (is.null(ref)) coords else .apply_transform(ref, coords)
  sqrt(rowSums((mapped - target)^2))
}

#' Summarize brain-mask displacement against ground truth motion
#'
#' @param mask 3D logical mask defining evaluation voxels.
#' @param motion Estimated motion represented as one of:
#'   an `nt x 6` pose matrix, a `4 x 4 x nt` array, or a `16 x nt` matrix.
#' @param truth Reference or ground-truth motion in the same accepted formats.
#' @param pixdim Length-3 voxel size vector in millimeters.
#' @param radius Sphere radius in millimeters for transform-space FD summaries.
#' @return A data frame with per-frame transform metrics:
#'   `disp_median`, `disp_p95`, `disp_max`, `fd_est`, `fd_truth`, and `fd_error`.
#' @export
motion_displacement_summary <- function(mask, motion, truth, pixdim = c(1, 1, 1), radius = 50) {
  mats <- .motion_matrices_array(motion)
  ref <- .motion_matrices_array(truth)

  nt <- dim(mats)[3]
  if (dim(ref)[3] != nt) {
    stop("truth must contain the same number of frames as motion", call. = FALSE)
  }

  out <- data.frame(
    frame = seq_len(nt),
    disp_median = numeric(nt),
    disp_p95 = numeric(nt),
    disp_max = numeric(nt),
    fd_est = fd_from_transform(mats, radius = radius),
    fd_truth = fd_from_transform(ref, radius = radius),
    fd_error = fd_from_transform(mats, radius = radius, reference = ref)
  )

  for (t in seq_len(nt)) {
    disp <- displacements_within_mask(
      mask = mask,
      motion = mats[, , t],
      reference = ref[, , t],
      pixdim = pixdim
    )
    out$disp_median[t] <- stats::median(disp)
    out$disp_p95[t] <- unname(stats::quantile(disp, probs = 0.95, names = FALSE))
    out$disp_max[t] <- max(disp)
  }

  out
}

.motion_matrices_array <- function(motion) {
  if (is.numeric(motion) && is.null(dim(motion))) {
    if (length(motion) != 6L) {
      stop("numeric motion vector must have length 6", call. = FALSE)
    }
    motion <- matrix(motion, nrow = 1)
  }

  if (is.matrix(motion) && identical(dim(motion), c(4L, 4L))) {
    out <- array(as.numeric(motion), dim = c(4L, 4L, 1L))
    return(out)
  }

  if (is.array(motion) && length(dim(motion)) == 3L && identical(dim(motion)[1:2], c(4L, 4L))) {
    return(array(as.numeric(motion), dim = dim(motion)))
  }

  if (is.matrix(motion) && ncol(motion) == 6L) {
    nt <- nrow(motion)
    out <- array(0, dim = c(4L, 4L, nt))
    for (t in seq_len(nt)) {
      out[, , t] <- .pose_to_matrix(motion[t, ])
    }
    return(out)
  }

  if (is.matrix(motion) && nrow(motion) == 16L) {
    nt <- ncol(motion)
    out <- array(0, dim = c(4L, 4L, nt))
    for (t in seq_len(nt)) {
      out[, , t] <- array(motion[, t], dim = c(4L, 4L))
    }
    return(out)
  }

  stop(
    "motion must be a length-6 vector, nt x 6 matrix, 4 x 4 matrix, 4 x 4 x nt array, or 16 x nt matrix",
    call. = FALSE
  )
}

.single_motion_matrix <- function(motion, arg = "motion") {
  mats <- .motion_matrices_array(motion)
  if (dim(mats)[3] != 1L) {
    stop(arg, " must describe exactly one transform", call. = FALSE)
  }
  mats[, , 1, drop = FALSE][, , 1]
}

.reference_matrices <- function(reference, n_frames) {
  if (is.null(reference)) {
    return(NULL)
  }

  ref <- .motion_matrices_array(reference)
  n_ref <- dim(ref)[3]
  if (!(n_ref == 1L || n_ref == n_frames)) {
    stop("reference must contain either one transform or the same number of frames as motion", call. = FALSE)
  }

  if (n_ref == 1L && n_frames > 1L) {
    out <- array(0, dim = c(4L, 4L, n_frames))
    for (t in seq_len(n_frames)) {
      out[, , t] <- ref[, , 1]
    }
    return(out)
  }

  ref
}

.pose_to_matrix <- function(pose) {
  stopifnot(length(pose) == 6L)
  R <- .rotation_matrix_zyx(pose[4], pose[5], pose[6])
  out <- diag(4)
  out[1:3, 1:3] <- R
  out[1:3, 4] <- pose[1:3]
  out
}

.rotation_matrix_zyx <- function(rx, ry, rz) {
  cx <- cos(rx)
  sx <- sin(rx)
  cy <- cos(ry)
  sy <- sin(ry)
  cz <- cos(rz)
  sz <- sin(rz)

  matrix(
    c(
      cz * cy,
      cz * sy * sx - sz * cx,
      cz * sy * cx + sz * sx,
      sz * cy,
      sz * sy * sx + cz * cx,
      sz * sy * cx - cz * sx,
      -sy,
      cy * sx,
      cy * cx
    ),
    nrow = 3,
    byrow = TRUE
  )
}

.mask_world_coords <- function(mask, pixdim = c(1, 1, 1)) {
  mask <- as.array(mask)
  storage.mode(mask) <- "logical"
  if (length(dim(mask)) != 3L) {
    stop("mask must be 3D", call. = FALSE)
  }

  pixdim <- as.numeric(pixdim)
  if (length(pixdim) < 3L || any(!is.finite(pixdim[1:3])) || any(pixdim[1:3] <= 0)) {
    stop("pixdim must contain three positive finite values", call. = FALSE)
  }
  pixdim <- pixdim[seq_len(3L)]

  ijk <- which(mask, arr.ind = TRUE)
  if (!nrow(ijk)) {
    stop("mask contains no TRUE voxels", call. = FALSE)
  }

  d <- dim(mask)
  cx <- 0.5 * (d[1] - 1)
  cy <- 0.5 * (d[2] - 1)
  cz <- 0.5 * (d[3] - 1)

  cbind(
    (ijk[, 1] - 1 - cx) * pixdim[1],
    (ijk[, 2] - 1 - cy) * pixdim[2],
    (ijk[, 3] - 1 - cz) * pixdim[3]
  )
}

.apply_transform <- function(mat, coords) {
  sweep(coords %*% t(mat[1:3, 1:3, drop = FALSE]), 2L, mat[1:3, 4], FUN = "+")
}
