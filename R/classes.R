#' @export
print.volregger_fit <- function(x, ...) {
  nt <- nrow(x$motion)
  has_corr <- !is.null(x$corrected)
  engine <- if (!is.null(x$control$engine)) x$control$engine else NA_character_
  fd_mean <- if (!is.null(x$qc$fd)) mean(x$qc$fd, na.rm = TRUE) else NA_real_
  fit_failures <- if (!is.null(x$qc$fit_failure)) sum(x$qc$fit_failure, na.rm = TRUE) else 0L
  cat("<volregger_fit>\n", sep = "")
  cat("  frames:", nt, "\n")
  cat("  corrected:", has_corr, "\n")
  cat("  engine:", engine, "\n")
  if (is.finite(fd_mean)) {
    cat("  fd_mean:", signif(fd_mean, 4), "\n")
  }
  cat("  fit_failures:", fit_failures, "\n")
  invisible(x)
}

#' @export
summary.volregger_fit <- function(object, ...) {
  qc <- object$qc %||% list()
  fd <- qc$fd
  dvars <- qc$dvars
  robust_dvars <- qc$robust_dvars
  out <- list(
    n_frames = nrow(object$motion),
    engine = object$control$engine %||% NA_character_,
    motion_abs_max = apply(abs(object$motion), 2, max),
    fd_mean = if (is.null(fd)) NA_real_ else mean(fd, na.rm = TRUE),
    fd_p95 = if (is.null(fd)) NA_real_ else unname(stats::quantile(fd, probs = 0.95, na.rm = TRUE)),
    fd_max = if (is.null(fd)) NA_real_ else max(fd, na.rm = TRUE),
    dvars_mean = if (is.null(dvars)) NA_real_ else mean(dvars, na.rm = TRUE),
    dvars_p95 = if (is.null(dvars)) NA_real_ else unname(stats::quantile(dvars, probs = 0.95, na.rm = TRUE)),
    robust_dvars_mean = if (is.null(robust_dvars)) NA_real_ else mean(robust_dvars, na.rm = TRUE),
    robust_dvars_p95 = if (is.null(robust_dvars)) NA_real_ else unname(stats::quantile(robust_dvars, probs = 0.95, na.rm = TRUE)),
    motion_spikes = if (is.null(qc$motion_spike)) 0L else sum(qc$motion_spike, na.rm = TRUE),
    fit_failures = if (is.null(qc$fit_failure)) 0L else sum(qc$fit_failure, na.rm = TRUE),
    censor_suggested = if (is.null(qc$censor_suggest)) 0L else sum(qc$censor_suggest, na.rm = TRUE),
    cost_drop_mean = if (is.null(qc$cost_drop)) NA_real_ else mean(qc$cost_drop, na.rm = TRUE),
    packet_correction_mag_mean = if (is.null(qc$packet_correction_mag)) NA_real_ else mean(qc$packet_correction_mag, na.rm = TRUE),
    packet_correction_mag_max = if (is.null(qc$packet_correction_mag)) NA_real_ else max(qc$packet_correction_mag, na.rm = TRUE)
  )
  class(out) <- "summary.volregger_fit"
  out
}

#' @export
print.summary.volregger_fit <- function(x, ...) {
  cat("<summary.volregger_fit>\n", sep = "")
  cat("  frames:", x$n_frames, "\n")
  cat("  engine:", x$engine, "\n")
  if (is.finite(x$fd_mean)) {
    cat("  fd_mean:", signif(x$fd_mean, 4), "\n")
    cat("  fd_p95:", signif(x$fd_p95, 4), "\n")
    cat("  fd_max:", signif(x$fd_max, 4), "\n")
  }
  if (is.finite(x$dvars_mean)) {
    cat("  dvars_mean:", signif(x$dvars_mean, 4), "\n")
    cat("  dvars_p95:", signif(x$dvars_p95, 4), "\n")
  }
  cat("  motion_spikes:", x$motion_spikes, "\n")
  cat("  fit_failures:", x$fit_failures, "\n")
  cat("  censor_suggested:", x$censor_suggested, "\n")
  invisible(x)
}

#' @export
coef.volregger_fit <- function(object, ...) {
  object$motion
}

#' @export
as.data.frame.volregger_fit <- function(x, ...) {
  out <- data.frame(
    frame = seq_len(nrow(x$motion)),
    tx = x$motion[, 1],
    ty = x$motion[, 2],
    tz = x$motion[, 3],
    rx = x$motion[, 4],
    ry = x$motion[, 5],
    rz = x$motion[, 6]
  )

  if (!is.null(x$qc$fd)) {
    out$fd <- x$qc$fd
  }
  if (!is.null(x$qc$dvars)) {
    out$dvars <- x$qc$dvars
  }
  if (!is.null(x$qc$robust_dvars)) {
    out$robust_dvars <- x$qc$robust_dvars
  }
  if (!is.null(x$qc$motion_spike)) {
    out$motion_spike <- as.logical(x$qc$motion_spike)
  }
  if (!is.null(x$qc$fit_failure)) {
    out$fit_failure <- as.logical(x$qc$fit_failure)
  }
  if (!is.null(x$qc$censor_suggest)) {
    out$censor_suggest <- as.logical(x$qc$censor_suggest)
  }
  if (!is.null(x$qc$cost_drop)) {
    out$cost_drop <- x$qc$cost_drop
  }
  if (!is.null(x$qc$packet_correction_mag)) {
    out$packet_correction_mag <- x$qc$packet_correction_mag
  }

  out
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
