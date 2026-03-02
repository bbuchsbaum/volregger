#' Apply motion parameters to a 4D run
#'
#' @param x 4D image-like input.
#' @param motion nt x 6 motion matrix.
#' @param interp Final interpolation method.
#' @param slice_times Optional slice timing offsets.
#' @param mb_groups Optional multiband group IDs.
#' @param nthreads Number of threads.
#' @return Corrected 4D array-like object.
#' @export
apply_motion <- function(
    x,
    motion,
    interp = c("bspline", "lanczos", "linear"),
    slice_times = NULL,
    mb_groups = NULL,
    nthreads = 1L) {
  interp <- match.arg(interp)
  motion <- .require_motion_matrix(motion)
  inp <- .as_volreg_input(x, mask = NULL, slice_times = slice_times, mb_groups = mb_groups)

  arr_corr <- cpp_volreg_apply(
    data = inp$data,
    dim = inp$dim,
    pixdim = inp$pixdim,
    par = motion,
    slice_times = inp$slice_times,
    mb_groups = inp$mb_groups,
    interp = interp,
    pad_mode = 0L,
    nthreads = as.integer(nthreads)
  )

  .rebuild_like(x, arr_corr)
}
