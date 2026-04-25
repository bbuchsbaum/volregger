#!/usr/bin/env Rscript

source(file.path("inst", "benchmarks", "realdata_benchmark_lib.R"), local = TRUE)
load_volregger_benchmark_env()

defaults <- list(
  input = NULL,
  mask = NULL,
  out_dir = NULL,
  profile = "ic_nuisance",
  ref = "robust_mean",
  interp = NULL,
  pad = "clamp",
  zpad = "0",
  huber_k = NULL,
  max_iter = NULL,
  edge_exclude_frac = NULL,
  nthreads = "4"
)
opts <- parse_benchmark_args(
  commandArgs(trailingOnly = TRUE),
  defaults = defaults,
  flags = c("with_report")
)
nthreads <- resolve_benchmark_nthreads(opts$nthreads)

if (is.null(opts$input) || is.null(opts$out_dir)) {
  stop("--input and --out-dir are required", call. = FALSE)
}
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

img <- RNifti::readNifti(opts$input)
x <- as.array(img)
storage.mode(x) <- "double"
attr(x, "pixdim") <- as.numeric(RNifti::pixdim(img)[seq_len(3L)])
attr(x, "tr") <- as.numeric(RNifti::pixdim(img)[4])
attr(x, "affine") <- RNifti::xform(img)

mask <- NULL
if (!is.null(opts$mask)) {
  mask <- as.array(RNifti::readNifti(opts$mask)) > 0
}

profile_overrides <- list(nthreads = nthreads)
if (!is.null(opts$huber_k)) {
  profile_overrides$huber_k <- as.numeric(opts$huber_k)
}
if (!is.null(opts$max_iter)) {
  max_iter <- as.integer(opts$max_iter)
  if (length(max_iter) != 1L || is.na(max_iter) || max_iter < 1L) {
    stop("--max-iter must be a positive integer", call. = FALSE)
  }
  profile_overrides$max_iter <- max_iter
}
if (!is.null(opts$edge_exclude_frac)) {
  profile_overrides$edge_exclude_frac <- as.numeric(opts$edge_exclude_frac)
}
profile <- volreg_profile(
  opts$profile,
  overrides = profile_overrides
)
interp_final <- opts$interp %||% profile$control$interp_final %||% "bspline"
if (!interp_final %in% c("bspline", "lanczos", "lanczos5", "lanczos7", "linear")) {
  stop("--interp must be one of: bspline, lanczos, lanczos5, lanczos7, linear", call. = FALSE)
}
if (!opts$pad %in% c("clamp", "zero")) {
  stop("--pad must be one of: clamp, zero", call. = FALSE)
}
zpad_final <- as.integer(opts$zpad)
if (length(zpad_final) != 1L || is.na(zpad_final) || zpad_final < 0L) {
  stop("--zpad must be a non-negative integer", call. = FALSE)
}

estimate_timing <- time_expr(
  estimate_motion(
    x = x,
    mask = mask,
    ref = opts$ref,
    engine = profile$engine,
    control = profile$control
  )
)
est <- estimate_timing$value

pose_spline <- NULL
par_apply <- est$par
if (identical(profile$engine, "rigid_spline")) {
  spline_timing <- time_expr(
    cpp_fit_pose_spline(
      par = est$par,
      tr = attr(x, "tr"),
      slice_times = NULL,
      mb_groups = NULL,
      control = profile$control
    )
  )
  pose_spline <- spline_timing$value
  par_apply <- pose_spline$par_smooth
  estimate_timing$elapsed_sec <- estimate_timing$elapsed_sec + spline_timing$elapsed_sec
}

apply_timing <- time_expr(
  apply_motion(
    x = x,
    motion = par_apply,
    interp = interp_final,
    pad = opts$pad,
    zpad = zpad_final,
    slice_times = if (identical(profile$engine, "rigid_spline")) NULL else NULL,
    mb_groups = if (identical(profile$engine, "rigid_spline")) NULL else NULL,
    nthreads = nthreads
  )
)
corrected <- apply_timing$value

mats <- motion_matrices(est$mats)
pose <- motion_mats_to_pose(mats)

write_bold_nifti(corrected, path = file.path(opts$out_dir, "corrected.nii.gz"), template = img)
write_motion_matrix_csv(mats, file.path(opts$out_dir, "motion_mats.csv"))
write_pose_csv(pose, file.path(opts$out_dir, "motion_pose.csv"))

report_sec <- 0
if (isTRUE(opts$with_report)) {
  report_timing <- time_expr({
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
    fit <- structure(
      list(
        corrected = corrected,
        motion = est$par,
        mats = est$mats,
        pose_spline = pose_spline,
        template = est$template,
        mask = est$mask,
        qc = qc_out,
        control = profile$control
      ),
      class = "volregger_fit"
    )
    report <- volreg_report(fit)
    write_volreg_report(
      report = report,
      out_dir = file.path(opts$out_dir, "report"),
      prefix = "volregger",
      include_plot = TRUE
    )
  })
  report_sec <- report_timing$elapsed_sec
}

write_stage_times_csv(
  data.frame(
    estimate_sec = estimate_timing$elapsed_sec,
    apply_sec = apply_timing$elapsed_sec,
    report_sec = report_sec,
    total_sec = estimate_timing$elapsed_sec + apply_timing$elapsed_sec + report_sec,
    stringsAsFactors = FALSE
  ),
  file.path(opts$out_dir, "stage_times.csv")
)

write.csv(
  data.frame(
    method = "volregger",
    profile = profile$name,
    ref = opts$ref,
    engine = profile$engine,
    interp = interp_final,
    pad = opts$pad,
    zpad = zpad_final,
    huber_k = benchmark_scalar_last(profile$control$huber_k, NA_real_),
    max_iter = benchmark_scalar_last(profile$control$max_iter, NA_integer_),
    edge_exclude_frac = benchmark_scalar_last(profile$control$edge_exclude_frac, 0),
    nthreads = nthreads,
    components = paste(profile$components, collapse = "+"),
    with_report = isTRUE(opts$with_report),
    stringsAsFactors = FALSE
  ),
  file.path(opts$out_dir, "run_metadata.csv"),
  row.names = FALSE
)
