#!/usr/bin/env Rscript

source(file.path("inst", "benchmarks", "realdata_benchmark_lib.R"), local = TRUE)
load_volregger_benchmark_env()

defaults <- list(
  input = NULL,
  out_dir = NULL,
  base_index = "0",
  image = "afni/afni_make_build",
  platform = "linux/amd64"
)
opts <- parse_benchmark_args(commandArgs(trailingOnly = TRUE), defaults = defaults, flags = "probe")

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

run_docker <- function(args) {
  suppressWarnings(system2("docker", args = args, stdout = TRUE, stderr = TRUE))
}

probe_afni <- function(image, platform) {
  out <- tryCatch(
    run_docker(docker_args(image, platform, c("3dvolreg", "-help"))),
    warning = function(w) character(),
    error = function(e) structure(character(), status = 1L, msg = conditionMessage(e))
  )
  status <- attr(out, "status", exact = TRUE) %||% 0L
  msg <- if (length(out)) paste(utils::head(out, 5L), collapse = "\n") else "No AFNI output"
  list(ok = identical(as.integer(status), 0L), message = msg)
}

read_afni_aff12 <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
  mats <- lapply(lines, function(line) {
    vals <- scan(text = line, quiet = TRUE)
    if (length(vals) != 12L) {
      stop("Expected 12 values per AFNI aff12 row", call. = FALSE)
    }
    M <- diag(4)
    M[1:3, 1:4] <- matrix(vals, nrow = 3, byrow = TRUE)
    M
  })
  simplify2array(mats)
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

convert_afni_mats <- function(aff12, img) {
  H <- internal_basis_affine(img)
  Hinv <- solve(H)
  lps_to_ras <- diag(c(-1, -1, 1, 1))
  nt <- dim(aff12)[3]
  out <- array(0, dim = c(4, 4, nt))

  for (i in seq_len(nt)) {
    base_to_input_lps <- aff12[, , i]
    input_to_base_lps <- solve(base_to_input_lps)
    out[, , i] <- Hinv %*% lps_to_ras %*% input_to_base_lps %*% lps_to_ras %*% H
  }

  out
}

if (isTRUE(opts$probe)) {
  probe <- probe_afni(opts$image, opts$platform)
  cat(if (probe$ok) "ok" else probe$message, "\n", sep = "")
  quit(status = if (probe$ok) 0L else 2L)
}

if (is.null(opts$input) || is.null(opts$out_dir)) {
  stop("--input and --out-dir are required", call. = FALSE)
}

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

corrected_path <- file.path(opts$out_dir, "corrected.nii.gz")
dfile_path <- file.path(opts$out_dir, "motion.1D")
aff12_path <- file.path(opts$out_dir, "motion.aff12.1D")

afni_args <- docker_args(
  image = opts$image,
  platform = opts$platform,
  forwarded = c(
    "3dvolreg",
    "-base", as.character(as.integer(opts$base_index)),
    "-1Dfile", dfile_path,
    "-1Dmatrix_save", aff12_path,
    "-prefix", corrected_path,
    opts$input
  )
)
status <- system2("docker", args = afni_args)
if (!identical(as.integer(status), 0L)) {
  quit(status = as.integer(status))
}

img <- RNifti::readNifti(opts$input)
aff12 <- read_afni_aff12(aff12_path)
motion_mats <- convert_afni_mats(aff12, img)
write_motion_matrix_csv(motion_mats, file.path(opts$out_dir, "motion_mats.csv"))

write.csv(
  data.frame(
    method = "afni_3dvolreg",
    image = opts$image,
    platform = opts$platform,
    base_index = as.integer(opts$base_index),
    stringsAsFactors = FALSE
  ),
  file.path(opts$out_dir, "run_metadata.csv"),
  row.names = FALSE
)
