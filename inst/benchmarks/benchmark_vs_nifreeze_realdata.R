#!/usr/bin/env Rscript

source(file.path("inst", "benchmarks", "realdata_benchmark_lib.R"), local = TRUE)
load_volregger_benchmark_env()

defaults <- list(
  input = NULL,
  out_dir = file.path("inst", "benchmarks", "realdata_results"),
  native_frames = "120",
  truth_frames = "120",
  scenarios = "low_motion,moderate_motion,hard_motion,hard_motion_plus_nuisance",
  nthreads = "4",
  python = Sys.which("python"),
  nifreeze_backend = "docker",
  nifreeze_docker_image = "volregger-nifreeze-benchmark:latest",
  nifreeze_docker_runner = file.path("inst", "benchmarks", "run_nifreeze_realdata_docker.sh"),
  volreg_profile = NULL,
  volreg_profile_native = "fast_fmri",
  volreg_profile_truth = "ic_nuisance_adaptive_tuned",
  volreg_interp = NULL,
  volreg_pad = "clamp",
  volreg_zpad = "0",
  volreg_huber_k = NULL,
  volreg_max_iter = NULL,
  volreg_edge_exclude_frac = NULL,
  volreg_ref = "robust_mean",
  nifreeze_model = "mean",
  afni_image = "afni/afni_make_build",
  afni_platform = "linux/amd64",
  seed = "4401"
)
opts <- parse_benchmark_args(
  commandArgs(trailingOnly = TRUE),
  defaults = defaults,
  flags = c("skip_nifreeze", "skip_afni", "native_only")
)
nthreads <- resolve_benchmark_nthreads(opts$nthreads)

input_path <- opts$input %||% default_realdata_input()
out_dir <- opts$out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

native_loaded <- load_bold_nifti(input_path, max_frames = opts$native_frames)
native_mask <- derive_brain_mask(native_loaded$data)

native_dir <- file.path(out_dir, "scenarios", "native")
dir.create(native_dir, recursive = TRUE, showWarnings = FALSE)
native_input_path <- file.path(native_dir, "input.nii.gz")
native_mask_path <- file.path(native_dir, "mask.nii.gz")
write_bold_nifti(native_loaded$data, native_input_path, native_loaded$template)
write_mask_nifti(native_mask, native_mask_path, native_loaded$template)
native_base_index <- select_reference_frame(native_loaded$data, native_mask) - 1L

scenario_names <- strsplit(opts$scenarios, ",", fixed = TRUE)[[1]]
scenario_names <- trimws(scenario_names)
scenario_names <- scenario_names[nzchar(scenario_names)]
if (isTRUE(opts$native_only)) {
  scenario_names <- character()
}

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

truth_loaded <- NULL
truth_mask <- NULL
if (length(scenario_names)) {
  truth_loaded <- load_bold_nifti(input_path, max_frames = opts$truth_frames)
  truth_mask <- derive_brain_mask(truth_loaded$data)

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
}

manifest <- do.call(rbind, manifest_rows)
utils::write.csv(manifest, file.path(out_dir, "scenario_manifest.csv"), row.names = FALSE)

probe_nifreeze <- function() {
  if (isTRUE(opts$skip_nifreeze)) {
    return(list(available = FALSE, message = "Skipped by user request"))
  }
  if (identical(opts$nifreeze_backend, "docker")) {
    probe <- run_profiled_command(
      command = "bash",
      args = c(
        opts$nifreeze_docker_runner,
        "--image", opts$nifreeze_docker_image,
        "--probe"
      ),
      stdout_path = file.path(out_dir, "nifreeze_probe_stdout.txt"),
      stderr_path = file.path(out_dir, "nifreeze_probe_stderr.txt")
    )
  } else if (identical(opts$nifreeze_backend, "local")) {
    probe <- run_profiled_command(
      command = opts$python,
      args = c("inst/benchmarks/run_nifreeze_realdata.py", "--probe"),
      stdout_path = file.path(out_dir, "nifreeze_probe_stdout.txt"),
      stderr_path = file.path(out_dir, "nifreeze_probe_stderr.txt")
    )
  } else {
    stop("Unknown nifreeze backend: ", opts$nifreeze_backend, call. = FALSE)
  }
  stdout_lines <- if (file.exists(file.path(out_dir, "nifreeze_probe_stdout.txt"))) {
    readLines(file.path(out_dir, "nifreeze_probe_stdout.txt"), warn = FALSE)
  } else {
    character()
  }
  msg <- if (length(stdout_lines)) paste(stdout_lines, collapse = "\n") else method_error_message(file.path(out_dir, "nifreeze_probe_stderr.txt"))
  list(available = identical(probe$status, 0L), message = msg %||% NA_character_)
}

nifreeze_probe <- probe_nifreeze()

probe_afni <- function() {
  if (isTRUE(opts$skip_afni)) {
    return(list(available = FALSE, message = "Skipped by user request"))
  }
  probe <- run_profiled_command(
    command = file.path(R.home("bin"), "Rscript"),
    args = c(
      "inst/benchmarks/run_afni_3dvolreg_realdata.R",
      "--probe",
      "--image", opts$afni_image,
      "--platform", opts$afni_platform
    ),
    stdout_path = file.path(out_dir, "afni_probe_stdout.txt"),
    stderr_path = file.path(out_dir, "afni_probe_stderr.txt")
  )
  stdout_lines <- if (file.exists(file.path(out_dir, "afni_probe_stdout.txt"))) {
    readLines(file.path(out_dir, "afni_probe_stdout.txt"), warn = FALSE)
  } else {
    character()
  }
  msg <- if (length(stdout_lines)) paste(stdout_lines, collapse = "\n") else method_error_message(file.path(out_dir, "afni_probe_stderr.txt"))
  list(available = identical(probe$status, 0L), message = msg %||% NA_character_)
}

afni_probe <- probe_afni()

scenario_catalog <- c("native", scenario_names)
methods <- c("volregger", "nifreeze", "afni_3dvolreg")
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
    base_index <- native_base_index
    dataset_mode <- "native"
  } else {
    sc <- scenario_objects[[scenario_name]]
    input_run <- sc$input
    mask <- truth_mask
    input_path <- sc$input_path
    mask_path <- sc$mask_path
    truth <- sc$truth_mats
    clean <- sc$clean
    base_index <- sc$reference_index - 1L
    dataset_mode <- "truth_injected"
  }

  stdout_path <- file.path(method_dir, "stdout.txt")
  stderr_path <- file.path(method_dir, "stderr.txt")
  stage_path <- file.path(method_dir, "stage_times.csv")

  if (identical(method, "volregger")) {
    volreg_profile_name <- opts$volreg_profile %||% if (identical(dataset_mode, "native")) {
      opts$volreg_profile_native
    } else {
      opts$volreg_profile_truth
    }
    volreg_profile_spec <- volreg_profile(volreg_profile_name)
    volreg_interp <- opts$volreg_interp %||% volreg_profile_spec$control$interp_final %||% "bspline"
    if (!opts$volreg_pad %in% c("clamp", "zero")) {
      stop("--volreg-pad must be one of: clamp, zero", call. = FALSE)
    }
    volreg_zpad <- as.integer(opts$volreg_zpad)
    if (length(volreg_zpad) != 1L || is.na(volreg_zpad) || volreg_zpad < 0L) {
      stop("--volreg-zpad must be a non-negative integer", call. = FALSE)
    }
    prof <- run_profiled_command(
      command = file.path(R.home("bin"), "Rscript"),
      args = c(
        "inst/benchmarks/run_volregger_realdata.R",
        "--input", input_path,
        "--mask", mask_path,
        "--out-dir", method_dir,
        "--profile", volreg_profile_name,
        "--ref", opts$volreg_ref,
        "--interp", volreg_interp,
        "--pad", opts$volreg_pad,
        "--zpad", as.character(volreg_zpad),
        if (!is.null(opts$volreg_huber_k)) c("--huber-k", opts$volreg_huber_k),
        if (!is.null(opts$volreg_max_iter)) c("--max-iter", opts$volreg_max_iter),
        if (!is.null(opts$volreg_edge_exclude_frac)) c("--edge-exclude-frac", opts$volreg_edge_exclude_frac),
        "--nthreads", as.character(nthreads)
      ),
      stdout_path = stdout_path,
      stderr_path = stderr_path
    )
    run_meta <- NULL
    run_meta_path <- file.path(method_dir, "run_metadata.csv")
    if (file.exists(run_meta_path)) {
      run_meta <- utils::read.csv(run_meta_path, stringsAsFactors = FALSE)
      if (nrow(run_meta)) {
        run_meta <- run_meta[1, , drop = FALSE]
      } else {
        run_meta <- NULL
      }
    }
    actual_huber <- opts$volreg_huber_k %||% "default"
    actual_iter <- opts$volreg_max_iter %||% "default"
    actual_edge <- opts$volreg_edge_exclude_frac %||% 0
    if (!is.null(run_meta)) {
      if ("huber_k" %in% names(run_meta) && is.finite(run_meta$huber_k[1])) {
        actual_huber <- run_meta$huber_k[1]
      }
      if ("max_iter" %in% names(run_meta) && is.finite(run_meta$max_iter[1])) {
        actual_iter <- run_meta$max_iter[1]
      }
      if ("edge_exclude_frac" %in% names(run_meta) && is.finite(run_meta$edge_exclude_frac[1])) {
        actual_edge <- run_meta$edge_exclude_frac[1]
      }
    }
    meta <- data.frame(
      method = method,
      method_model = paste(
        volreg_profile_name,
        opts$volreg_ref,
        volreg_interp,
        opts$volreg_pad,
        paste0("zpad", volreg_zpad),
        paste0("huber", actual_huber),
        paste0("iter", actual_iter),
        paste0("edge", actual_edge),
        paste0("threads", nthreads),
        sep = ":"
      ),
      stringsAsFactors = FALSE
    )
  } else if (identical(method, "nifreeze")) {
    if (!isTRUE(nifreeze_probe$available)) {
      return(cbind(
        data.frame(
          method = method,
          dataset_mode = dataset_mode,
          scenario = scenario_name,
          method_model = opts$nifreeze_model,
          elapsed_sec = NA_real_,
          estimate_sec = NA_real_,
          apply_sec = NA_real_,
          report_sec = NA_real_,
          peak_rss_mb = NA_real_,
          temp_bytes_written = NA_real_,
          status = "unavailable",
          error = nifreeze_probe$message,
          stringsAsFactors = FALSE
        ),
        empty_benchmark_metrics()
      ))
    }

    if (identical(opts$nifreeze_backend, "docker")) {
      prof <- run_profiled_command(
        command = "bash",
        args = c(
          opts$nifreeze_docker_runner,
          "--image", opts$nifreeze_docker_image,
          "--input", input_path,
          "--mask", mask_path,
          "--out-dir", method_dir,
          "--model", opts$nifreeze_model,
          "--nthreads", as.character(nthreads),
          "--seed", opts$seed
        ),
        stdout_path = stdout_path,
        stderr_path = stderr_path
      )
    } else if (identical(opts$nifreeze_backend, "local")) {
      prof <- run_profiled_command(
        command = opts$python,
        args = c(
          "inst/benchmarks/run_nifreeze_realdata.py",
          "--input", input_path,
          "--mask", mask_path,
          "--out-dir", method_dir,
          "--model", opts$nifreeze_model,
          "--nthreads", as.character(nthreads),
          "--seed", opts$seed
        ),
        stdout_path = stdout_path,
        stderr_path = stderr_path
      )
    } else {
      stop("Unknown nifreeze backend: ", opts$nifreeze_backend, call. = FALSE)
    }
    meta <- data.frame(
      method = method,
      method_model = paste(opts$nifreeze_backend, opts$nifreeze_model, paste0("threads", nthreads), sep = ":"),
      stringsAsFactors = FALSE
    )
  } else if (identical(method, "afni_3dvolreg")) {
    if (!isTRUE(afni_probe$available)) {
      return(cbind(
        data.frame(
          method = method,
          dataset_mode = dataset_mode,
          scenario = scenario_name,
          method_model = paste(opts$afni_image, base_index, sep = ":base="),
          elapsed_sec = NA_real_,
          estimate_sec = NA_real_,
          apply_sec = NA_real_,
          report_sec = NA_real_,
          peak_rss_mb = NA_real_,
          temp_bytes_written = NA_real_,
          status = "unavailable",
          error = afni_probe$message,
          stringsAsFactors = FALSE
        ),
        empty_benchmark_metrics()
      ))
    }
    prof <- run_profiled_command(
      command = file.path(R.home("bin"), "Rscript"),
      args = c(
        "inst/benchmarks/run_afni_3dvolreg_realdata.R",
        "--input", input_path,
        "--out-dir", method_dir,
        "--base-index", as.character(base_index),
        "--image", opts$afni_image,
        "--platform", opts$afni_platform
      ),
      stdout_path = stdout_path,
      stderr_path = stderr_path
    )
    meta <- data.frame(
      method = method,
      method_model = paste(opts$afni_image, base_index, sep = ":base="),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown method: ", method, call. = FALSE)
  }

  if (!identical(prof$status, 0L)) {
    stage_times <- read_stage_times_csv(stage_path, total_sec = prof$elapsed_sec)
    return(cbind(
      data.frame(
        method = method,
        dataset_mode = dataset_mode,
        scenario = scenario_name,
        method_model = meta$method_model,
        elapsed_sec = prof$elapsed_sec,
        estimate_sec = stage_times$estimate_sec,
        apply_sec = stage_times$apply_sec,
        report_sec = stage_times$report_sec,
        peak_rss_mb = prof$peak_rss_mb,
        temp_bytes_written = dir_size_bytes(method_dir),
        status = "error",
        error = method_error_message(stderr_path),
        stringsAsFactors = FALSE
      ),
      empty_benchmark_metrics()
    ))
  }

  corrected <- as.array(RNifti::readNifti(file.path(method_dir, "corrected.nii.gz")))
  storage.mode(corrected) <- "double"
  attr(corrected, "pixdim") <- attr(input_run, "pixdim")
  attr(corrected, "tr") <- attr(input_run, "tr")
  attr(corrected, "affine") <- attr(input_run, "affine")
  mats <- read_motion_matrix_csv(file.path(method_dir, "motion_mats.csv"))
  stage_times <- read_stage_times_csv(stage_path, total_sec = prof$elapsed_sec)
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
      estimate_sec = stage_times$estimate_sec,
      apply_sec = stage_times$apply_sec,
      report_sec = stage_times$report_sec,
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

scorecard <- make_benchmark_scorecard(raw)
scorecard_summary <- summarize_benchmark_scorecard(scorecard)
scorecard_claim_summary <- summarize_benchmark_scorecard(scorecard, claim_only = TRUE)
scorecard_csv <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_scorecard.csv")
scorecard_summary_csv <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_scorecard_summary.csv")
scorecard_claim_summary_csv <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_claim_summary.csv")
utils::write.csv(scorecard, scorecard_csv, row.names = FALSE)
utils::write.csv(scorecard_summary, scorecard_summary_csv, row.names = FALSE)
utils::write.csv(scorecard_claim_summary, scorecard_claim_summary_csv, row.names = FALSE)

safe_ratio <- function(x, y) {
  if (!is.finite(x) || !is.finite(y)) {
    return(NA_real_)
  }
  x / max(y, 1e-8)
}

safe_delta <- function(x, y) {
  if (!is.finite(x) || !is.finite(y)) {
    return(NA_real_)
  }
  x - y
}

safe_unit_error_ratio <- function(x, y) {
  if (!is.finite(x) || !is.finite(y)) {
    return(NA_real_)
  }
  abs(x - 1) / max(abs(y - 1), 1e-8)
}

paired <- list()
pair_idx <- 1L
for (sc in unique(raw$scenario)) {
  a <- raw[raw$scenario == sc & raw$method == "volregger" & raw$status == "ok", ]
  if (!nrow(a)) {
    next
  }
  for (baseline in setdiff(unique(raw$method), "volregger")) {
    b <- raw[raw$scenario == sc & raw$method == baseline & raw$status == "ok", ]
    if (!nrow(b)) {
      next
    }
    paired[[pair_idx]] <- data.frame(
      scenario = sc,
      dataset_mode = a$dataset_mode[1],
      baseline = baseline,
      elapsed_ratio = safe_ratio(a$elapsed_sec, b$elapsed_sec),
      estimate_ratio = safe_ratio(a$estimate_sec, b$estimate_sec),
      apply_ratio = safe_ratio(a$apply_sec, b$apply_sec),
      report_ratio = safe_ratio(a$report_sec, b$report_sec),
      temp_bytes_ratio = safe_ratio(a$temp_bytes_written, b$temp_bytes_written),
      disp_median_ratio = safe_ratio(a$disp_median, b$disp_median),
      disp_p95_ratio = safe_ratio(a$disp_p95, b$disp_p95),
      fd_error_ratio = safe_ratio(a$fd_error, b$fd_error),
      trans_abs_p95_ratio = safe_ratio(a$trans_abs_p95, b$trans_abs_p95),
      rot_abs_deg_p95_ratio = safe_ratio(a$rot_abs_deg_p95, b$rot_abs_deg_p95),
      rel_frob_p95_ratio = safe_ratio(a$rel_frob_p95, b$rel_frob_p95),
      rel_rot_deg_p95_ratio = safe_ratio(a$rel_rot_deg_p95, b$rel_rot_deg_p95),
      mse_to_truth_ratio = safe_ratio(a$mse_to_truth, b$mse_to_truth),
      mae_to_truth_ratio = safe_ratio(a$mae_to_truth, b$mae_to_truth),
      nrmse_to_truth_ratio = safe_ratio(a$nrmse_to_truth, b$nrmse_to_truth),
      residual_dvars_to_truth_ratio = safe_ratio(a$residual_dvars_to_truth, b$residual_dvars_to_truth),
      tsnr_ratio_ratio = safe_ratio(a$tsnr_ratio, b$tsnr_ratio),
      blur_ratio_ratio = safe_ratio(a$blur_proxy_ratio, b$blur_proxy_ratio),
      corr_to_truth_delta = safe_delta(a$corr_to_truth, b$corr_to_truth),
      temporal_corr_to_truth_delta = safe_delta(a$temporal_corr_to_truth, b$temporal_corr_to_truth),
      edge_ratio_to_truth_error_ratio = safe_unit_error_ratio(a$edge_ratio_to_truth, b$edge_ratio_to_truth),
      stringsAsFactors = FALSE
    )
    pair_idx <- pair_idx + 1L
  }
}
paired <- if (length(paired)) do.call(rbind, paired) else data.frame()
paired_csv <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_paired.csv")
utils::write.csv(paired, paired_csv, row.names = FALSE)

solution_rows <- list()
solution_frames <- list()
sol_idx <- 1L
frame_idx <- 1L
for (sc in unique(raw$scenario)) {
  ok_rows <- raw[raw$scenario == sc & raw$status == "ok", ]
  if (nrow(ok_rows) < 2L) {
    next
  }
  scenario_mask <- if (identical(sc, "native")) native_mask else truth_mask
  pixdim <- if (identical(sc, "native")) attr(native_loaded$data, "pixdim") else attr(truth_loaded$data, "pixdim")
  ref_index <- if (identical(sc, "native")) (native_base_index + 1L) else scenario_objects[[sc]]$reference_index

  for (i in seq_len(nrow(ok_rows) - 1L)) {
    for (j in (i + 1L):nrow(ok_rows)) {
      method_a <- ok_rows$method[i]
      method_b <- ok_rows$method[j]
      dir_a <- file.path(out_dir, "runs", sc, method_a)
      dir_b <- file.path(out_dir, "runs", sc, method_b)
      mats_a <- read_motion_matrix_csv(file.path(dir_a, "motion_mats.csv"))
      mats_b <- read_motion_matrix_csv(file.path(dir_b, "motion_mats.csv"))
      corrected_a <- as.array(RNifti::readNifti(file.path(dir_a, "corrected.nii.gz")))
      corrected_b <- as.array(RNifti::readNifti(file.path(dir_b, "corrected.nii.gz")))
      storage.mode(corrected_a) <- "double"
      storage.mode(corrected_b) <- "double"

      cmp <- pairwise_solution_comparison(
        mask = scenario_mask,
        motion_a = mats_a,
        motion_b = mats_b,
        pixdim = pixdim,
        ref_index = ref_index,
        corrected_a = corrected_a,
        corrected_b = corrected_b
      )
      solution_rows[[sol_idx]] <- cbind(
        data.frame(
          scenario = sc,
          dataset_mode = ok_rows$dataset_mode[i],
          method_a = method_a,
          method_b = method_b,
          stringsAsFactors = FALSE
        ),
        cmp$summary
      )
      solution_frames[[frame_idx]] <- cbind(
        data.frame(
          scenario = sc,
          dataset_mode = ok_rows$dataset_mode[i],
          method_a = method_a,
          method_b = method_b,
          stringsAsFactors = FALSE
        ),
        cmp$frames
      )
      sol_idx <- sol_idx + 1L
      frame_idx <- frame_idx + 1L
    }
  }
}
solution_summary <- if (length(solution_rows)) do.call(rbind, solution_rows) else data.frame()
solution_frames <- if (length(solution_frames)) do.call(rbind, solution_frames) else data.frame()
solution_summary_csv <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_solution_summary.csv")
solution_frames_csv <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_solution_frames.csv")
utils::write.csv(solution_summary, solution_summary_csv, row.names = FALSE)
utils::write.csv(solution_frames, solution_frames_csv, row.names = FALSE)

cat("Saved:\n")
cat(" -", raw_csv, "\n")
cat(" -", scorecard_csv, "\n")
cat(" -", scorecard_summary_csv, "\n")
cat(" -", scorecard_claim_summary_csv, "\n")
cat(" -", paired_csv, "\n")
cat(" -", solution_summary_csv, "\n")
cat(" -", solution_frames_csv, "\n")
cat(" -", file.path(out_dir, "scenario_manifest.csv"), "\n\n")
print(raw[, c(
  "method",
  "dataset_mode",
  "scenario",
  "status",
  "elapsed_sec",
  "estimate_sec",
  "apply_sec",
  "report_sec",
  "disp_p95",
  "fd_error",
  "trans_abs_p95",
  "rot_abs_deg_p95",
  "mae_to_truth",
  "nrmse_to_truth",
  "temporal_corr_to_truth",
  "truth_image_fidelity_loss",
  "tsnr_ratio",
  "blur_proxy_ratio"
)])

if (nrow(scorecard_summary)) {
  cat("\nScorecard summary:\n")
  print(scorecard_summary)
}

if (nrow(scorecard_claim_summary)) {
  cat("\nClaim-metric summary:\n")
  print(scorecard_claim_summary)
}
