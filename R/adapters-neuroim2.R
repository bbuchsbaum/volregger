.as_volreg_input <- function(x, mask = NULL, slice_times = NULL, mb_groups = NULL) {
  arr <- as.array(x)
  if (length(dim(arr)) != 4L) {
    stop("x must be 4D", call. = FALSE)
  }
  if (any(!is.finite(arr))) {
    stop("x contains non-finite values", call. = FALSE)
  }
  storage.mode(arr) <- "double"

  pixdim <- attr(x, "pixdim")
  if (is.null(pixdim) || length(pixdim) < 3L) {
    pixdim <- c(1, 1, 1)
  } else {
    pixdim <- as.numeric(pixdim[seq_len(3L)])
  }

  affine <- attr(x, "affine")
  if (is.null(affine)) {
    affine <- diag(4)
  }

  tr <- attr(x, "tr")
  if (is.null(tr)) {
    tr <- 1
  }

  mask_arr <- NULL
  if (!is.null(mask)) {
    mask_arr <- as.array(mask)
    if (!identical(dim(mask_arr), dim(arr)[1:3])) {
      stop("mask must match spatial dimensions of x", call. = FALSE)
    }
    if (any(!is.finite(mask_arr))) {
      stop("mask contains non-finite values", call. = FALSE)
    }
    storage.mode(mask_arr) <- "logical"
  }

  list(
    data = arr,
    dim = as.integer(dim(arr)),
    pixdim = as.numeric(pixdim),
    affine = affine,
    tr = as.numeric(tr),
    mask = mask_arr,
    slice_times = if (is.null(slice_times)) NULL else as.numeric(slice_times),
    mb_groups = if (is.null(mb_groups)) NULL else as.integer(mb_groups)
  )
}

.rebuild_like <- function(x, arr4d) {
  out <- array(as.numeric(arr4d), dim = dim(as.array(x)))
  attr(out, "pixdim") <- attr(x, "pixdim")
  attr(out, "affine") <- attr(x, "affine")
  attr(out, "tr") <- attr(x, "tr")
  out
}

.rebuild_3d_like <- function(x, arr3d) {
  spatial_dim <- dim(as.array(x))[1:3]
  out <- array(arr3d, dim = spatial_dim)
  attr(out, "pixdim") <- attr(x, "pixdim")
  attr(out, "affine") <- attr(x, "affine")
  out
}
