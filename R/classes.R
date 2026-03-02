#' @export
print.volregger_fit <- function(x, ...) {
  nt <- nrow(x$motion)
  has_corr <- !is.null(x$corrected)
  cat("<volregger_fit>\n", sep = "")
  cat("  frames:", nt, "\n")
  cat("  corrected:", has_corr, "\n")
  invisible(x)
}

#' @export
summary.volregger_fit <- function(object, ...) {
  fd <- NULL
  if (!is.null(object$qc$fd)) {
    fd <- object$qc$fd
  }
  list(
    n_frames = nrow(object$motion),
    motion_abs_max = apply(abs(object$motion), 2, max),
    fd_mean = if (is.null(fd)) NA_real_ else mean(fd, na.rm = TRUE)
  )
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

  out
}
