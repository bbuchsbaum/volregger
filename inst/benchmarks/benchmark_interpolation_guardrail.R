#!/usr/bin/env Rscript

source(file.path("inst", "benchmarks", "realdata_benchmark_lib.R"), local = TRUE)

defaults <- list(
  input = NULL,
  out_dir = file.path("inst", "benchmarks", "interpolation_guardrail"),
  native_frames = "8",
  truth_frames = "24",
  scenarios = "low_motion,hard_motion_plus_nuisance",
  interps = "linear,bspline,lanczos,lanczos5,lanczos7",
  baseline_interp = "lanczos",
  volreg_profile_native = "fast_fmri",
  volreg_profile_truth = "ic_nuisance_adaptive_tuned",
  volreg_pad = "clamp",
  volreg_zpad = "0",
  nthreads = "4"
)
opts <- parse_benchmark_args(commandArgs(trailingOnly = TRUE), defaults = defaults)

interps <- trimws(strsplit(opts$interps, ",", fixed = TRUE)[[1]])
interps <- interps[nzchar(interps)]
valid_interps <- c("linear", "bspline", "lanczos", "lanczos5", "lanczos7")
bad <- setdiff(interps, valid_interps)
if (length(bad)) {
  stop("--interps contains unsupported values: ", paste(bad, collapse = ", "), call. = FALSE)
}
if (!opts$baseline_interp %in% interps) {
  stop("--baseline-interp must be one of --interps", call. = FALSE)
}

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
raw_rows <- list()

for (interp in interps) {
  run_dir <- file.path(opts$out_dir, paste0("interp_", interp))
  args <- c(
    "inst/benchmarks/benchmark_vs_nifreeze_realdata.R",
    "--out-dir", run_dir,
    "--native-frames", opts$native_frames,
    "--truth-frames", opts$truth_frames,
    "--scenarios", opts$scenarios,
    "--volreg-profile-native", opts$volreg_profile_native,
    "--volreg-profile-truth", opts$volreg_profile_truth,
    "--volreg-interp", interp,
    "--volreg-pad", opts$volreg_pad,
    "--volreg-zpad", opts$volreg_zpad,
    "--nthreads", opts$nthreads,
    "--skip-nifreeze",
    "--skip-afni"
  )
  if (!is.null(opts$input)) {
    args <- c(args, "--input", opts$input)
  }
  status <- system2(file.path(R.home("bin"), "Rscript"), args = args)
  if (!identical(status, 0L)) {
    stop("interpolation guardrail run failed for interp=", interp, call. = FALSE)
  }
  raw_path <- file.path(run_dir, "benchmark_vs_nifreeze_realdata_raw.csv")
  raw <- utils::read.csv(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
  raw$interp <- interp
  raw$guardrail_run_dir <- run_dir
  raw_rows[[length(raw_rows) + 1L]] <- raw
}

combined <- do.call(rbind, raw_rows)
summary <- summarize_interpolation_guardrail(combined, baseline_interp = opts$baseline_interp)

raw_csv <- file.path(opts$out_dir, "interpolation_guardrail_raw.csv")
summary_csv <- file.path(opts$out_dir, "interpolation_guardrail_summary.csv")
utils::write.csv(combined, raw_csv, row.names = FALSE)
utils::write.csv(summary, summary_csv, row.names = FALSE)

cat("Saved:\n")
cat(" -", raw_csv, "\n")
cat(" -", summary_csv, "\n")
