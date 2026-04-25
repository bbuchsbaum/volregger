#' Apply motion parameters to a 4D run
#'
#' @param x 4D image-like input.
#' @param motion nt x 6 motion matrix.
#' @param interp Final interpolation method.
#' @param pad Boundary handling for samples outside the input field of view.
#' @param zpad Number of symmetric zero-valued voxels to add before resampling
#'   and cropping back to the original grid.
#' @param slice_times Optional slice timing offsets.
#' @param mb_groups Optional multiband group IDs.
#' @param nthreads Number of threads.
#' @return Corrected 4D array-like object.
#' @export
apply_motion <- function(
    x,
    motion,
    interp = c("bspline", "lanczos", "lanczos5", "lanczos7", "quintic", "heptic", "linear"),
    pad = c("clamp", "zero"),
    zpad = 0L,
    slice_times = NULL,
    mb_groups = NULL,
    nthreads = 1L) {
  interp <- match.arg(interp)
  pad <- match.arg(pad)
  zpad <- as.integer(zpad)
  if (length(zpad) != 1L || is.na(zpad) || zpad < 0L) {
    stop("zpad must be a non-negative integer", call. = FALSE)
  }
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
    pad_mode = if (identical(pad, "zero")) 1L else 0L,
    zpad = zpad,
    nthreads = as.integer(nthreads)
  )

  .rebuild_like(x, arr_corr)
}
