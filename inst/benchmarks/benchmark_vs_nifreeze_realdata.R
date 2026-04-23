#!/usr/bin/env Rscript

source(file.path("inst", "benchmarks", "realdata_benchmark_lib.R"), local = TRUE)
load_volregger_benchmark_env()

defaults <- list(
  input = NULL,
  out_dir = file.path("inst", "benchmarks", "realdata_results"),
  native_frames = "120",
  truth_frames = "120",
  scenarios = "low_motion,moderate_motion,hard_motion,hard_motion_plus_nuisance",
  nthreads = "1",
  python = Sys.which("python"),
  volreg_profile = "ic_nuisance",
  nifreeze_model = "mean",
  seed = "4401"
)
opts <- parse_benchmark_args(
  commandArgs(trailingOnly = TRUE),
  defaults = defaults,
  flags = c("skip_nifreeze")
)

input_path <- opts$input %||% default_realdata_input()
out_dir <- opts$out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

native_loaded <- load_bold_nifti(input_path, max_frames = opts$native_frames)
truth_loaded <- load_bold_nifti(input_path, max_frames = opts$truth_frames)

native_mask <- derive_brain_mask(native_loaded$data)
truth_mask <- derive_brain_mask(truth_loaded$data)

native_dir <- file.path(out_dir, "scenarios", "native")
dir.create(native_dir, recursive = TRUE, showWarnings = FALSE)
native_input_path <- file.path(native_dir, "input.nii.gz")
native_mask_path <- file.path(native_dir, "mask.nii.gz")
write_bold_nifti(native_loaded$data, native_input_path, native_loaded$template)
write_mask_nifti(native_mask, native_mask_path, native_loaded$template)

scenario_names <- strsplit(opts$scenarios, ",", fixed = TRUE)[[1]]
scenario_names <- trimws(scenario_names)

scenario_objects <- vector("list", length(scenario_names))
names(scenario_objects) <- scenario_names
manifest_rows <- list(
  data.frame(
    scenario = "native",
    dataset_mode = "native",
    nt = dim(native_loaded$data)[4],
    drift_amp = NA_real_,
    bias_amp = NA_real_,
    reference_index = NA_integer_,
    stringsAsFactors = FALSE
  )
)

for (name in scenario_names) {
  sc <- make_truth_scenario(truth_loaded$data, truth_mask, name, truth_frames = dim(truth_loaded$data)[4])
  scenario_dir <- file.path(out_dir, "scenarios", name)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  input_file <- file.path(scenario_dir, "input.nii.gz")
  clean_file <- file.path(scenario_dir, "clean.nii.gz")
  mask_file <- file.path(scenario_dir, "mask.nii.gz")
  truth_pose_file <- file.path(scenario_dir, "truth_pose.csv")
  truth_mats_file <- file.path(scenario_dir, "truth_mats.csv")

  write_bold_nifti(sc$input, input_file, truth_loaded$template)
  write_bold_nifti(sc$clean, clean_file, truth_loaded$template)
  write_mask_nifti(truth_mask, mask_file, truth_loaded$template)
  write_pose_csv(sc$truth_pose, truth_pose_file)
  write_motion_matrix_csv(sc$truth_mats, truth_mats_file)

  scenario_objects[[name]] <- c(sc, list(input_path = input_file, clean_path = clean_file, mask_path = mask_file))
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    scenario = name,
    dataset_mode = "truth_injected",
    nt = dim(sc$input)[4],
    drift_amp = sc$drift_amp,
    bias_amp = sc$bias_amp,
    reference_index = sc$reference_index,
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
utils::write.csv(manifest, file.path(out_dir, "scenario_manifest.csv"), row.names = FALSE)

probe_nifreeze <- function() {
  if (isTRUE(opts$skip_nifreeze)) {
    return(list(available = FALSE, message = "Skipped by user request"))
  }
  probe <- run_profiled_command(
    command = opts$python,
    args = c("inst/benchmarks/run_nifreeze_realdata.py", "--probe"),
    stdout_path = file.path(out_dir, "nifreeze_probe_stdout.txt"),
    stderr_path = file.path(out_dir, "nifreeze_probe_stderr.txt")
  )
  stdout_lines <- if (file.exists(file.path(out_dir, "nifreeze_probe_stdout.txt"))) {
    readLines(file.path(out_dir, "nifreeze_probe_stdout.txt"), warn = FALSE)
  } else {
    character()
  }
  msg <- if (length(stdout_lines)) paste(stdout_lines, collapse = "\n") else method_error_message(file.path(out_dir, "nifreeze_probe_stderr.txt"))
  list(available = identical(probe$status, 0L), message = msg %||% NA_character_)
}

nifreeze_probe <- probe_nifreeze()

scenario_catalog <- c("native", scenario_names)
methods <- c("volregger", "nifreeze")
rows <- vector("list", length(scenario_catalog) * length(methods))
idx <- 1L

run_one <- function(method, scenario_name) {
  scenario_dir <- file.path(out_dir, "scenarios", scenario_name)
  method_dir <- file.path(out_dir, "runs", scenario_name, method)
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

  if (identical(scenario_name, "native")) {
    input_run <- native_loaded$data
    mask <- native_mask
    input_path <- native_input_path
    mask_path <- native_mask_path
    truth <- NULL
    clean <- NULL
    dataset_mode <- "native"
  } else {
    sc <- scenario_objects[[scenario_name]]
    input_run <- sc$input
    mask <- truth_mask
    input_path <- sc$input_path
    mask_path <- sc$mask_path
    truth <- sc$truth_mats
    clean <- sc$clean
    dataset_mode <- "truth_injected"
  }

  stdout_path <- file.path(method_dir, "stdout.txt")
  stderr_path <- file.path(method_dir, "stderr.txt")

  if (identical(method, "volregger")) {
    prof <- run_profiled_command(
      command = file.path(R.home("bin"), "Rscript"),
      args = c(
        "inst/benchmarks/run_volregger_realdata.R",
        "--input", input_path,
        "--mask", mask_path,
        "--out-dir", method_dir,
        "--profile", opts$volreg_profile,
        "--nthreads", opts$nthreads
      ),
      stdout_path = stdout_path,
      stderr_path = stderr_path
    )
    meta <- data.frame(
      method = method,
      method_model = opts$volreg_profile,
      stringsAsFactors = FALSE
    )
  } else {
    if (!isTRUE(nifreeze_probe$available)) {
      return(data.frame(
        method = method,
        dataset_mode = dataset_mode,
        scenario = scenario_name,
        method_model = opts$nifreeze_model,
        elapsed_sec = NA_real_,
        peak_rss_mb = NA_real_,
        temp_bytes_written = NA_real_,
        status = "unavailable",
        error = nifreeze_probe$message,
        fit_fail_rate = NA_real_,
        fd_mean = NA_real_,
        fd_p95 = NA_real_,
        dvars_before_mean = NA_real_,
        dvars_before_p95 = NA_real_,
        dvars_after_mean = NA_real_,
        dvars_after_p95 = NA_real_,
        dvars_ratio = NA_real_,
        tsnr_before = NA_real_,
        tsnr_after = NA_real_,
        tsnr_ratio = NA_real_,
        blur_proxy_before = NA_real_,
        blur_proxy_after = NA_real_,
        blur_proxy_ratio = NA_real_,
        disp_median = NA_real_,
        disp_p95 = NA_real_,
        disp_max = NA_real_,
        fd_error = NA_real_,
        trans_rmse = NA_real_,
        rot_rmse = NA_real_,
        corr_to_truth = NA_real_,
        mse_to_truth = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    prof <- run_profiled_command(
      command = opts$python,
      args = c(
        "inst/benchmarks/run_nifreeze_realdata.py",
        "--input", input_path,
        "--mask", mask_path,
        "--out-dir", method_dir,
        "--model", opts$nifreeze_model,
        "--nthreads", opts$nthreads,
        "--seed", opts$seed
      ),
      stdout_path = stdout_path,
      stderr_path = stderr_path
    )
    meta <- data.frame(
      method = method,
      method_model = opts$nifreeze_model,
      stringsAsFactors = FALSE
    )
  }

  if (!identical(prof$status, 0L)) {
    return(data.frame(
      method = method,
      dataset_mode = dataset_mode,
      scenario = scenario_name,
      method_model = meta$method_model,
      elapsed_sec = prof$elapsed_sec,
      peak_rss_mb = prof$peak_rss_mb,
      temp_bytes_written = dir_size_bytes(method_dir),
      status = "error",
      error = method_error_message(stderr_path),
      fit_fail_rate = NA_real_,
      fd_mean = NA_real_,
      fd_p95 = NA_real_,
      dvars_before_mean = NA_real_,
      dvars_before_p95 = NA_real_,
      dvars_after_mean = NA_real_,
      dvars_after_p95 = NA_real_,
      dvars_ratio = NA_real_,
      tsnr_before = NA_real_,
      tsnr_after = NA_real_,
      tsnr_ratio = NA_real_,
      blur_proxy_before = NA_real_,
      blur_proxy_after = NA_real_,
      blur_proxy_ratio = NA_real_,
      disp_median = NA_real_,
      disp_p95 = NA_real_,
      disp_max = NA_real_,
      fd_error = NA_real_,
      trans_rmse = NA_real_,
      rot_rmse = NA_real_,
      corr_to_truth = NA_real_,
      mse_to_truth = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  corrected <- as.array(RNifti::readNifti(file.path(method_dir, "corrected.nii.gz")))
  storage.mode(corrected) <- "double"
  attr(corrected, "pixdim") <- attr(input_run, "pixdim")
  attr(corrected, "tr") <- attr(input_run, "tr")
  attr(corrected, "affine") <- attr(input_run, "affine")
  mats <- read_motion_matrix_csv(file.path(method_dir, "motion_mats.csv"))
  metrics <- compute_benchmark_metrics(
    input_run = input_run,
    corrected = corrected,
    mask = mask,
    motion_mats = mats,
    truth = truth,
    clean = clean
  )

  cbind(
    data.frame(
      method = method,
      dataset_mode = dataset_mode,
      scenario = scenario_name,
      method_model = meta$method_model,
      elapsed_sec = prof$elapsed_sec,
      peak_rss_mb = prof$peak_rss_mb,
      temp_bytes_written = dir_size_bytes(method_dir),
      status = "ok",
      error = NA_character_,
      stringsAsFactors = FALSE
    ),
    metrics
  )
}

for (sc in scenario_catalog) {
  for (method in methods) {
    rows[[idx]] <- run_one(method, sc)
    idx <- idx + 1L
  }
}

raw <- do.call(rbind, rows)
raw_csv <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_raw.csv")
utils::write.csv(raw, raw_csv, row.names = FALSE)

paired <- list()
pair_idx <- 1L
for (sc in unique(raw$scenario)) {
  a <- raw[raw$scenario == sc & raw$method == "volregger" & raw$status == "ok", ]
  b <- raw[raw$scenario == sc & raw$method == "nifreeze" & raw$status == "ok", ]
  if (!nrow(a) || !nrow(b)) {
    next
  }
  paired[[pair_idx]] <- data.frame(
    scenario = sc,
    dataset_mode = a$dataset_mode[1],
    elapsed_ratio = a$elapsed_sec / b$elapsed_sec,
    temp_bytes_ratio = a$temp_bytes_written / max(b$temp_bytes_written, 1e-8),
    disp_median_ratio = a$disp_median / max(b$disp_median, 1e-8),
    disp_p95_ratio = a$disp_p95 / max(b$disp_p95, 1e-8),
    fd_error_ratio = a$fd_error / max(b$fd_error, 1e-8),
    tsnr_ratio_ratio = a$tsnr_ratio / max(b$tsnr_ratio, 1e-8),
    blur_ratio_ratio = a$blur_proxy_ratio / max(b$blur_proxy_ratio, 1e-8),
    stringsAsFactors = FALSE
  )
  pair_idx <- pair_idx + 1L
}
paired <- if (length(paired)) do.call(rbind, paired) else data.frame()
paired_csv <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_paired.csv")
utils::write.csv(paired, paired_csv, row.names = FALSE)

cat("Saved:\n")
cat(" -", raw_csv, "\n")
cat(" -", paired_csv, "\n")
cat(" -", file.path(out_dir, "scenario_manifest.csv"), "\n\n")
print(raw[, c("method", "dataset_mode", "scenario", "status", "elapsed_sec", "disp_p95", "fd_error", "tsnr_ratio", "blur_proxy_ratio")])
