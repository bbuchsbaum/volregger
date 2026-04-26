#!/usr/bin/env Rscript

source(file.path("inst", "benchmarks", "realdata_benchmark_lib.R"), local = TRUE)
load_volregger_benchmark_env()

defaults <- list(
  comparison_dir = NULL,
  out_dir = NULL,
  scenarios = "low_motion,hard_motion_plus_nuisance",
  volreg_interps = "heptic",
  afni_finals = "quintic,wsinc5",
  afni_image = "afni/afni_make_build",
  afni_platform = "linux/amd64"
)
opts <- parse_benchmark_args(
  commandArgs(trailingOnly = TRUE),
  defaults = defaults,
  flags = c("include_afni_apply")
)

if (is.null(opts$comparison_dir)) {
  stop("--comparison-dir is required", call. = FALSE)
}
comparison_dir <- normalizePath(opts$comparison_dir, winslash = "/", mustWork = TRUE)
out_dir <- opts$out_dir %||% file.path(comparison_dir, "resampler_swap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario_manifest_path <- file.path(comparison_dir, "scenario_manifest.csv")
if (!file.exists(scenario_manifest_path)) {
  stop("Missing scenario manifest under comparison dir: ", scenario_manifest_path, call. = FALSE)
}
manifest <- utils::read.csv(scenario_manifest_path, stringsAsFactors = FALSE)
scenario_names <- trimws(strsplit(opts$scenarios, ",", fixed = TRUE)[[1]])
scenario_names <- scenario_names[nzchar(scenario_names)]
scenario_names <- intersect(scenario_names, manifest$scenario[manifest$dataset_mode == "truth_injected"])
if (!length(scenario_names)) {
  stop("--scenarios did not match any truth-injected scenarios in the comparison dir", call. = FALSE)
}

volreg_interps <- trimws(strsplit(opts$volreg_interps, ",", fixed = TRUE)[[1]])
volreg_interps <- volreg_interps[nzchar(volreg_interps)]
valid_volreg_interps <- c("linear", "bspline", "lanczos", "lanczos5", "lanczos7", "quintic", "heptic")
bad <- setdiff(volreg_interps, valid_volreg_interps)
if (length(bad)) {
  stop("--volreg-interps contains unsupported values: ", paste(bad, collapse = ", "), call. = FALSE)
}

afni_finals <- trimws(strsplit(opts$afni_finals, ",", fixed = TRUE)[[1]])
afni_finals <- afni_finals[nzchar(afni_finals)]

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
docker_args <- function(image, platform, forwarded) {
  c(
    "run", "--rm",
    "--platform", platform,
    "-v", paste0(repo_root, ":/work"),
    "-v", "/tmp:/tmp",
    "-w", "/work",
    image,
    forwarded
  )
}

docker_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (startsWith(path, "/private/tmp/")) {
    return(sub("^/private/tmp", "/tmp", path))
  }
  path
}

internal_basis_affine <- function(img) {
  aff <- RNifti::xform(img)
  dim3 <- dim(img)[seq_len(3L)]
  pixdim <- as.numeric(RNifti::pixdim(img)[seq_len(3L)])
  center <- 0.5 * (dim3 - 1)
  Dinv <- diag(c(1 / pixdim, 1))
  C <- diag(4)
  C[1:3, 4] <- center
  aff %*% C %*% Dinv
}

internal_to_afni_aff12 <- function(mats, img) {
  H <- internal_basis_affine(img)
  Hinv <- solve(H)
  lps_to_ras <- diag(c(-1, -1, 1, 1))
  nt <- dim(mats)[3]
  out <- array(0, dim = c(4, 4, nt))
  for (i in seq_len(nt)) {
    source_to_base_lps <- lps_to_ras %*% H %*% mats[, , i] %*% Hinv %*% lps_to_ras
    out[, , i] <- solve(source_to_base_lps)
  }
  out
}

write_afni_aff12 <- function(aff12, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  nt <- dim(aff12)[3]
  for (i in seq_len(nt)) {
    vals <- as.numeric(t(aff12[1:3, 1:4, i]))
    writeLines(paste(format(vals, digits = 12, scientific = FALSE), collapse = " "), con)
  }
  invisible(path)
}

apply_afni_matrix <- function(input_path, mats, template_img, out_path, final = "quintic", platform = opts$afni_platform) {
  aff12_path <- sub("\\.nii(\\.gz)?$", ".aff12.1D", out_path)
  write_afni_aff12(internal_to_afni_aff12(mats, template_img), aff12_path)
  args <- docker_args(
    image = opts$afni_image,
    platform = platform,
    forwarded = c(
      "3dAllineate",
      "-source", docker_path(input_path),
      "-master", docker_path(input_path),
      "-final", final,
      "-float",
      "-1Dmatrix_apply", docker_path(aff12_path),
      "-prefix", docker_path(out_path)
    )
  )
  start <- proc.time()[["elapsed"]]
  status <- system2("docker", args = args)
  elapsed <- proc.time()[["elapsed"]] - start
  list(status = as.integer(status), elapsed_sec = elapsed, aff12_path = aff12_path)
}

load_corrected <- function(path, template_run) {
  corrected <- as.array(RNifti::readNifti(path))
  storage.mode(corrected) <- "double"
  attr(corrected, "pixdim") <- attr(template_run, "pixdim")
  attr(corrected, "tr") <- attr(template_run, "tr")
  attr(corrected, "affine") <- attr(template_run, "affine")
  corrected
}

make_row <- function(method, method_model, scenario_name, corrected, mats, input_run, mask, truth, clean,
                     elapsed_sec = NA_real_, apply_sec = NA_real_, status = "ok", error = NA_character_) {
  cbind(
    data.frame(
      method = method,
      dataset_mode = "truth_injected",
      scenario = scenario_name,
      method_model = method_model,
      elapsed_sec = elapsed_sec,
      estimate_sec = NA_real_,
      apply_sec = apply_sec,
      report_sec = 0,
      peak_rss_mb = NA_real_,
      temp_bytes_written = NA_real_,
      status = status,
      error = error,
      stringsAsFactors = FALSE
    ),
    if (identical(status, "ok")) {
      compute_benchmark_metrics(input_run, corrected, mask, mats, truth = truth, clean = clean)
    } else {
      empty_benchmark_metrics()
    }
  )
}

rows <- list()
idx <- 1L

for (scenario_name in scenario_names) {
  scenario_dir <- file.path(comparison_dir, "scenarios", scenario_name)
  input_path <- file.path(scenario_dir, "input.nii.gz")
  clean_path <- file.path(scenario_dir, "clean.nii.gz")
  mask_path <- file.path(scenario_dir, "mask.nii.gz")
  truth_path <- file.path(scenario_dir, "truth_mats.csv")
  if (!all(file.exists(input_path, clean_path, mask_path, truth_path))) {
    stop("Scenario directory is missing required files: ", scenario_dir, call. = FALSE)
  }

  input_loaded <- load_bold_nifti(input_path, max_frames = "all")
  clean_loaded <- load_bold_nifti(clean_path, max_frames = "all")
  input_run <- input_loaded$data
  clean <- clean_loaded$data
  mask <- as.array(RNifti::readNifti(mask_path)) > 0
  truth <- read_motion_matrix_csv(truth_path)
  template_img <- RNifti::readNifti(input_path)

  volreg_mats <- read_motion_matrix_csv(file.path(comparison_dir, "runs", scenario_name, "volregger", "motion_mats.csv"))
  afni_mats <- read_motion_matrix_csv(file.path(comparison_dir, "runs", scenario_name, "afni_3dvolreg", "motion_mats.csv"))

  for (method in c("volregger", "afni_3dvolreg")) {
    corrected_path <- file.path(comparison_dir, "runs", scenario_name, method, "corrected.nii.gz")
    if (!file.exists(corrected_path)) {
      next
    }
    mats <- if (identical(method, "volregger")) volreg_mats else afni_mats
    corrected <- load_corrected(corrected_path, input_run)
    rows[[idx]] <- make_row(
      method = method,
      method_model = paste0(method, ":native_resampler"),
      scenario_name = scenario_name,
      corrected = corrected,
      mats = mats,
      input_run = input_run,
      mask = mask,
      truth = truth,
      clean = clean
    )
    idx <- idx + 1L
  }

  for (source_name in c("volregger", "afni_3dvolreg")) {
    mats <- if (identical(source_name, "volregger")) volreg_mats else afni_mats
    pose <- motion_mats_to_pose(mats)
    for (interp in volreg_interps) {
      start <- proc.time()[["elapsed"]]
      corrected <- apply_motion(input_run, pose, interp = interp)
      apply_sec <- proc.time()[["elapsed"]] - start
      rows[[idx]] <- make_row(
        method = paste(source_name, "volreg_resampler", sep = "_x_"),
        method_model = paste(source_name, "volregger_apply", interp, sep = ":"),
        scenario_name = scenario_name,
        corrected = corrected,
        mats = mats,
        input_run = input_run,
        mask = mask,
        truth = truth,
        clean = clean,
        elapsed_sec = apply_sec,
        apply_sec = apply_sec
      )
      idx <- idx + 1L
    }
  }

  if (isTRUE(opts$include_afni_apply)) {
    for (source_name in c("volregger", "afni_3dvolreg")) {
      mats <- if (identical(source_name, "volregger")) volreg_mats else afni_mats
      for (final in afni_finals) {
        out_path <- file.path(out_dir, paste0(paste(scenario_name, source_name, paste0("afni_", final), sep = "_"), ".nii.gz"))
        prof <- apply_afni_matrix(input_path, mats, template_img, out_path, final = final)
        if (!identical(prof$status, 0L)) {
          rows[[idx]] <- make_row(
            method = paste(source_name, "afni_resampler", sep = "_x_"),
            method_model = paste(source_name, "afni_3dAllineate", final, sep = ":"),
            scenario_name = scenario_name,
            corrected = input_run,
            mats = mats,
            input_run = input_run,
            mask = mask,
            truth = truth,
            clean = clean,
            elapsed_sec = prof$elapsed_sec,
            apply_sec = prof$elapsed_sec,
            status = "error",
            error = paste("3dAllineate exited", prof$status)
          )
        } else {
          corrected <- load_corrected(out_path, input_run)
          rows[[idx]] <- make_row(
            method = paste(source_name, "afni_resampler", sep = "_x_"),
            method_model = paste(source_name, "afni_3dAllineate", final, sep = ":"),
            scenario_name = scenario_name,
            corrected = corrected,
            mats = mats,
            input_run = input_run,
            mask = mask,
            truth = truth,
            clean = clean,
            elapsed_sec = prof$elapsed_sec,
            apply_sec = prof$elapsed_sec
          )
        }
        idx <- idx + 1L
      }
    }
  }
}

raw <- do.call(rbind, rows)
summary <- summarize_resampler_swap(raw)

raw_csv <- file.path(out_dir, "resampler_swap_raw.csv")
summary_csv <- file.path(out_dir, "resampler_swap_summary.csv")
utils::write.csv(raw, raw_csv, row.names = FALSE)
utils::write.csv(summary, summary_csv, row.names = FALSE)

cat("Saved:\n")
cat(" -", raw_csv, "\n")
cat(" -", summary_csv, "\n")
