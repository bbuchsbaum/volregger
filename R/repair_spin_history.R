#' Optional spin-history intensity repair
#'
#' @param x Corrected 4D image-like input.
#' @param motion Motion matrix (`nt x 6`).
#' @param mask Optional 3D logical mask.
#' @param method Repair method.
#' @param z_threshold Robust z-threshold for residual spike detection.
#' @param smooth_window Odd integer smoothing window for temporal gain smoothing.
#' @return A list with corrected output, gain map, and detection metadata.
#' @export
repair_spin_history <- function(
    x,
    motion,
    mask = NULL,
    method = c("slice_gain"),
    z_threshold = 3,
    smooth_window = 5L) {
  method <- match.arg(method)
  motion <- .require_motion_matrix(motion)
  arr <- as.array(x)
  if (length(dim(arr)) != 4L) {
    stop("x must be 4D", call. = FALSE)
  }

  nx <- dim(arr)[1]
  ny <- dim(arr)[2]
  nz <- dim(arr)[3]
  nt <- dim(arr)[4]

  if (nrow(motion) != nt) {
    stop("nrow(motion) must equal ntime(x)", call. = FALSE)
  }

  m <- array(TRUE, dim = c(nx, ny, nz))
  if (!is.null(mask)) {
    m <- as.array(mask)
    if (!identical(dim(m), c(nx, ny, nz))) {
      stop("mask dimensions must match x spatial dimensions", call. = FALSE)
    }
    m <- m > 0
  }

  if (!identical(method, "slice_gain")) {
    stop("Only slice_gain is currently implemented", call. = FALSE)
  }

  template <- apply(arr, c(1, 2, 3), stats::median, na.rm = TRUE)
  slice_resid <- matrix(0, nrow = nz, ncol = nt)

  for (t in seq_len(nt)) {
    for (k in seq_len(nz)) {
      mk <- m[, , k]
      if (!any(mk)) {
        next
      }
      r <- arr[, , k, t][mk] - template[, , k][mk]
      slice_resid[k, t] <- stats::median(r, na.rm = TRUE)
    }
  }

  tz <- motion[, 3]
  through_plane <- c(0, abs(diff(tz)))
  tp_thr <- stats::quantile(through_plane, probs = 0.75, na.rm = TRUE)

  gains <- matrix(1, nrow = nz, ncol = nt)
  spike_mask <- matrix(FALSE, nrow = nz, ncol = nt)

  for (k in seq_len(nz)) {
    r <- slice_resid[k, ]
    med <- stats::median(r, na.rm = TRUE)
    madv <- stats::mad(r, center = med, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(madv) || madv < 1e-8) {
      next
    }

    z <- (r - med) / madv
    spike <- abs(z) > z_threshold & through_plane >= tp_thr
    spike_mask[k, ] <- spike

    # conservative multiplicative gain around 1 from residual magnitude
    denom <- stats::median(abs(template[, , k][m[, , k]]), na.rm = TRUE)
    if (!is.finite(denom) || denom < 1e-8) {
      denom <- 1
    }
    raw_gain <- 1 + (r / denom)

    w <- as.integer(max(3L, smooth_window))
    if (w %% 2L == 0L) {
      w <- w + 1L
    }
    ker <- rep(1 / w, w)
    smooth_gain <- as.numeric(stats::filter(raw_gain, filter = ker, sides = 2))
    smooth_gain[is.na(smooth_gain)] <- raw_gain[is.na(smooth_gain)]

    # apply only on detected high-confidence spikes
    gains[k, spike] <- smooth_gain[spike]
  }

  out <- arr
  for (t in seq_len(nt)) {
    for (k in seq_len(nz)) {
      g <- gains[k, t]
      if (!is.finite(g) || abs(g - 1) < 1e-6) {
        next
      }
      if (!spike_mask[k, t]) {
        next
      }
      out[, , k, t] <- out[, , k, t] / g
    }
  }

  list(
    corrected = .rebuild_like(x, out),
    gains = gains,
    slice_residual = slice_resid,
    spike_mask = spike_mask,
    through_plane = through_plane
  )
}
