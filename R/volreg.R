#' Run retrospective motion correction
#'
#' @param x 4D image-like input.
#' @param mask Optional 3D mask.
#' @param ref Reference strategy.
#' @param ref_index Explicit reference index for `ref='index'`.
#' @param engine Motion engine.
#' @param slice_times Optional slice timing offsets.
#' @param mb_groups Optional multiband group IDs.
#' @param resample If TRUE, apply estimated motion to generate corrected output.
#' @param interp_final Final interpolation method.
#' @param pad_final Boundary handling for final resampling.
#' @param zpad_final Number of symmetric zero-valued voxels to add before final
#'   resampling and cropping back to the original grid.
#' @param qc If TRUE, compute QC metrics.
#' @param control Control list.
#' @param verbose Verbose logging.
#' @return A `volregger_fit` object.
#' @export
volreg <- function(
    x,
    mask = NULL,
    ref = c("middle", "robust_mean", "index"),
    ref_index = NULL,
    engine = c("rigid_robust", "rigid_spline"),
    slice_times = NULL,
    mb_groups = NULL,
    resample = TRUE,
    interp_final = c("bspline", "lanczos", "linear"),
    pad_final = c("clamp", "zero"),
    zpad_final = 0L,
    qc = TRUE,
    control = volreg_control(),
    verbose = FALSE) {
  ref <- match.arg(ref)
  engine <- match.arg(engine)
  interp_final <- match.arg(interp_final)
  pad_final <- match.arg(pad_final)
  zpad_final <- as.integer(zpad_final)
  if (length(zpad_final) != 1L || is.na(zpad_final) || zpad_final < 0L) {
    stop("zpad_final must be a non-negative integer", call. = FALSE)
  }

  inp <- .as_volreg_input(
    x = x,
    mask = mask,
    slice_times = slice_times,
    mb_groups = mb_groups
  )

  ctl <- utils::modifyList(control, list(engine = engine, verbose = isTRUE(verbose), tr = inp$tr))
  ctl$ref_mode <- ref

  est <- cpp_volreg_estimate(
    data = inp$data,
    dim = inp$dim,
    pixdim = inp$pixdim,
    mask = inp$mask,
    slice_times = inp$slice_times,
    mb_groups = inp$mb_groups,
    ref_index = .resolve_ref_index(ref, ref_index, inp$dim[4]),
    control = ctl
  )

  pose_spline <- NULL
  par_apply <- est$par
  if (identical(engine, "rigid_spline")) {
    pose_spline <- cpp_fit_pose_spline(
      par = est$par,
      tr = inp$tr,
      slice_times = inp$slice_times,
      mb_groups = inp$mb_groups,
      control = ctl
    )
    par_apply <- pose_spline$par_smooth
  }

  corrected <- NULL
  if (isTRUE(resample)) {
    arr_corr <- cpp_volreg_apply(
      data = inp$data,
      dim = inp$dim,
      pixdim = inp$pixdim,
      par = par_apply,
      slice_times = if (identical(engine, "rigid_spline")) inp$slice_times else NULL,
      mb_groups = if (identical(engine, "rigid_spline")) inp$mb_groups else NULL,
      interp = interp_final,
      pad_mode = if (identical(pad_final, "zero")) 1L else 0L,
      zpad = zpad_final,
      nthreads = as.integer(ctl$nthreads)
    )
    corrected <- .rebuild_like(x, arr_corr)
  }

  qc_out <- NULL
  if (isTRUE(qc)) {
    qc_out <- motion_qc(
      x = x,
      motion = est$par,
      corrected = corrected,
      mask = est$mask,
      cost_init = est$cost_init,
      cost_final = est$cost_final
    )
    if (!is.null(pose_spline)) {
      qc_out$packet_correction_mag <- rowSums(abs(pose_spline$par_smooth - est$par))
    }
  }

  structure(
    list(
      corrected = corrected,
      motion = est$par,
      mats = est$mats,
      pose_spline = pose_spline,
      template = .rebuild_3d_like(x, est$template),
      mask = .rebuild_3d_like(x, est$mask),
      qc = qc_out,
      control = ctl
    ),
    class = "volregger_fit"
  )
}
