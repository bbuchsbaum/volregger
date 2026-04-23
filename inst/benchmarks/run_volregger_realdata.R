#!/usr/bin/env Rscript

source(file.path("inst", "benchmarks", "realdata_benchmark_lib.R"), local = TRUE)
load_volregger_benchmark_env()

defaults <- list(
  input = NULL,
  mask = NULL,
  out_dir = NULL,
  profile = "ic_nuisance",
  ref = "robust_mean",
  nthreads = "1"
)
opts <- parse_benchmark_args(commandArgs(trailingOnly = TRUE), defaults = defaults)

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

profile <- volreg_profile(
  opts$profile,
  overrides = list(nthreads = as.integer(opts$nthreads))
)

fit <- volreg(
  x = x,
  mask = mask,
  ref = opts$ref,
  engine = profile$engine,
  control = profile$control,
  resample = TRUE,
  qc = TRUE
)

mats <- motion_matrices(fit$mats)
pose <- motion_mats_to_pose(mats)

write_bold_nifti(fit$corrected, path = file.path(opts$out_dir, "corrected.nii.gz"), template = img)
write_motion_matrix_csv(mats, file.path(opts$out_dir, "motion_mats.csv"))
write_pose_csv(pose, file.path(opts$out_dir, "motion_pose.csv"))
write.csv(
  data.frame(
    method = "volregger",
    profile = profile$name,
    ref = opts$ref,
    engine = profile$engine,
    components = paste(profile$components, collapse = "+"),
    stringsAsFactors = FALSE
  ),
  file.path(opts$out_dir, "run_metadata.csv"),
  row.names = FALSE
)

report <- volreg_report(fit)
write_volreg_report(
  report = report,
  out_dir = file.path(opts$out_dir, "report"),
  prefix = "volregger",
  include_plot = TRUE
)
