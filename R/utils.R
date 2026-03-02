.resolve_ref_index <- function(ref, ref_index, nt) {
  if (isTRUE(nt < 1L)) {
    stop("Input has no timepoints", call. = FALSE)
  }

  if (identical(ref, "index")) {
    if (is.null(ref_index) || length(ref_index) != 1L) {
      stop("ref_index must be length 1 when ref='index'", call. = FALSE)
    }
    idx <- as.integer(ref_index)
    if (idx < 1L || idx > nt) {
      stop("ref_index out of range", call. = FALSE)
    }
    return(idx)
  }

  if (identical(ref, "middle") || identical(ref, "robust_mean")) {
    return(as.integer((nt + 1L) %/% 2L))
  }

  stop("Unknown ref mode", call. = FALSE)
}

.require_motion_matrix <- function(motion) {
  if (!is.matrix(motion) || ncol(motion) != 6L) {
    stop("motion must be an nt x 6 matrix", call. = FALSE)
  }
  storage.mode(motion) <- "double"
  motion
}
