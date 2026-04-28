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
#' @param profile Optional named estimator profile. Bare `estimate_motion()`
#'   uses `"fast_fmri"` when `control` is not supplied.
#' @param verbose Verbose logging.
#' @return A list containing estimated motion and fit diagnostics.
#' @examples
#' x <- volregger_example()
#' est <- estimate_motion(x, profile = "fast_native")
#' dim(est$par)  # nt x 6
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
    profile = NULL,
    verbose = FALSE) {
  control_missing <- missing(control)
  engine_missing <- missing(engine)
  ref <- match.arg(ref)
  if (is.null(profile) && control_missing) {
    profile <- "fast_fmri"
  }
  prof <- .apply_volreg_profile(control, profile, control_user = !control_missing)
  control <- prof$control
  volreg_validate_control(control, strict = TRUE)
  if (!is.null(prof$engine) && engine_missing) {
    engine <- prof$engine
  }
  engine <- match.arg(engine, c("rigid_robust", "rigid_spline"))

  inp <- .as_volreg_input(
    x = x,
    mask = mask,
    slice_times = slice_times,
    mb_groups = mb_groups
  )

  ctl <- utils::modifyList(control, list(engine = engine, verbose = isTRUE(verbose), tr = inp$tr))
  ctl$ref_mode <- ref

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
