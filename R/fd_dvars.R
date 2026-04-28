#' Compute framewise displacement from motion parameters
#'
#' @param motion nt x 6 motion matrix.
#' @param radius Head radius in millimeters.
#' @return Numeric vector of FD values.
#' @examples
#' motion <- matrix(0, nrow = 5, ncol = 6)
#' motion[3, 1] <- 0.5  # 0.5 mm shift in x at frame 3
#' framewise_disp(motion)
#' @export
framewise_disp <- function(motion, radius = 50) {
  motion <- .require_motion_matrix(motion)
  nt <- nrow(motion)
  if (nt < 1L) {
    return(numeric())
  }

  fd <- numeric(nt)
  if (nt == 1L) {
    return(fd)
  }

  d <- abs(motion[2:nt, , drop = FALSE] - motion[1:(nt - 1L), , drop = FALSE])
  fd[2:nt] <- rowSums(d[, 1:3, drop = FALSE]) + radius * rowSums(d[, 4:6, drop = FALSE])
  fd
}

#' Compute DVARS from a 4D run
#'
#' @param x 4D image-like input.
#' @param mask Optional 3D mask.
#' @param robust If TRUE, return robust DVARS estimate.
#' @return Numeric vector of DVARS values.
#' @examples
#' x <- volregger_example()
#' dvars(x)
#' @export
dvars <- function(x, mask = NULL, robust = TRUE) {
  arr <- as.array(x)
  if (length(dim(arr)) != 4L) {
    stop("x must be 4D", call. = FALSE)
  }

  m <- rep(TRUE, prod(dim(arr)[1:3]))
  if (!is.null(mask)) {
    m <- as.logical(as.array(mask))
    if (length(m) != prod(dim(arr)[1:3])) {
      stop("mask must match x spatial dimensions", call. = FALSE)
    }
  }

  nt <- dim(arr)[4]
  out <- rep(NA_real_, nt)
  if (nt < 2L) {
    return(out)
  }

  for (t in 2:nt) {
    prev <- as.numeric(arr[, , , t - 1L])[m]
    curr <- as.numeric(arr[, , , t])[m]
    delta <- curr - prev
    out[t] <- sqrt(mean(delta * delta))
  }

  if (isTRUE(robust)) {
    med <- stats::median(out, na.rm = TRUE)
    out[is.finite(out)] <- pmin(out[is.finite(out)], 3 * med)
  }

  out
}
