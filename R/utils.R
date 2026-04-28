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

#' Tiny synthetic 4D fixture used in examples
#'
#' A deterministic 4D array with a small in-plane shift inserted at frame 4,
#' suitable for fast roxygen examples. Not for benchmarking.
#'
#' @param dim Integer vector of length 4 giving `c(nx, ny, nz, nt)`.
#' @param seed Random seed used to fill the volume.
#' @return A 4D numeric array of dimension `dim`.
#' @examples
#' x <- volregger_example()
#' dim(x)
#' @export
volregger_example <- function(dim = c(8L, 8L, 4L, 6L), seed = 1L) {
  dim <- as.integer(dim)
  stopifnot(length(dim) == 4L, all(dim > 0L))
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  base <- array(stats::rnorm(prod(dim[1:3])) * 0.1 + 1, dim = dim[1:3])
  ix <- seq_len(dim[1]); iy <- seq_len(dim[2]); iz <- seq_len(dim[3])
  blob <- outer(outer(exp(-((ix - dim[1] / 2)^2) / 4),
                     exp(-((iy - dim[2] / 2)^2) / 4)),
                exp(-((iz - dim[3] / 2)^2) / 4))
  base <- base + blob
  out <- array(0, dim = dim)
  for (t in seq_len(dim[4])) {
    if (t == 4L && dim[1] >= 2L) {
      shifted <- base
      shifted[-1, , ] <- base[-dim[1], , ]
      shifted[1, , ] <- base[1, , ]
      out[, , , t] <- shifted
    } else {
      out[, , , t] <- base
    }
  }
  out
}
