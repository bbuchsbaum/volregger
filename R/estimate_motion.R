#' Estimate rigid motion parameters from a 4D run
#'
#' @param x 4D image-like input.
#' @param mask Optional 3D mask.
#' @param ref Reference strategy.
#' @param ref_index Explicit reference index when `ref = "index"`.
#' @param engine Motion engine.
#' @param slice_times Optional slice timing offsets.
#' @param mb_groups Optional multiband group IDs.
#' @param control Control list.
#' @param verbose Verbose logging.
#' @return A list containing estimated motion and fit diagnostics.
#' @export
estimate_motion <- function(
    x,
    mask = NULL,
    ref = c("middle", "robust_mean", "index"),
    ref_index = NULL,
    engine = c("rigid_robust", "rigid_spline"),
    slice_times = NULL,
    mb_groups = NULL,
    control = volreg_control(),
    verbose = FALSE) {
  ref <- match.arg(ref)
  engine <- match.arg(engine)

  inp <- .as_volreg_input(
    x = x,
    mask = mask,
    slice_times = slice_times,
    mb_groups = mb_groups
  )

  ctl <- utils::modifyList(control, list(engine = engine, verbose = isTRUE(verbose), tr = inp$tr))

  cpp_volreg_estimate(
    data = inp$data,
    dim = inp$dim,
    pixdim = inp$pixdim,
    mask = inp$mask,
    slice_times = inp$slice_times,
    mb_groups = inp$mb_groups,
    ref_index = .resolve_ref_index(ref, ref_index, inp$dim[4]),
    control = ctl
  )
}
