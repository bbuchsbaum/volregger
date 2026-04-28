# I/O helpers for the volregger CLI.
# NIfTI handling uses neuroim2 lazily so non-NIfTI commands stay dependency-free.

.cli_require_neuroim2 <- function() {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    .cli_abort(
      "this command requires the 'neuroim2' package; install it with install.packages('neuroim2')",
      class = "volregger_cli_domain_error"
    )
  }
}

.cli_is_nifti_path <- function(path) {
  grepl("\\.nii(\\.gz)?$", path, ignore.case = TRUE)
}

.cli_read_4d <- function(path, label = "input") {
  if (!file.exists(path)) {
    .cli_abort(
      sprintf("%s file does not exist: %s", label, path),
      class = "volregger_cli_domain_error"
    )
  }
  if (.cli_is_nifti_path(path)) {
    .cli_require_neuroim2()
    v <- tryCatch(
      neuroim2::read_vec(path),
      error = function(e) {
        .cli_abort(
          sprintf("failed to read %s '%s': %s", label, path, conditionMessage(e)),
          class = "volregger_cli_domain_error"
        )
      }
    )
    arr <- as.array(v)
    if (length(dim(arr)) != 4L) {
      .cli_abort(
        sprintf("%s '%s' is %dD; a 4D NIfTI is required", label, path, length(dim(arr))),
        class = "volregger_cli_domain_error"
      )
    }
    spc <- neuroim2::spacing(v)
    attr(arr, "pixdim") <- as.numeric(spc[seq_len(min(3L, length(spc)))])
    attr(arr, "affine") <- neuroim2::trans(v)
    attr(arr, "tr") <- 1
    attr(arr, ".neuroim2_template") <- v
    return(arr)
  }
  .cli_read_rds(path, label)
}

.cli_read_3d <- function(path, label = "mask") {
  if (!file.exists(path)) {
    .cli_abort(
      sprintf("%s file does not exist: %s", label, path),
      class = "volregger_cli_domain_error"
    )
  }
  if (.cli_is_nifti_path(path)) {
    .cli_require_neuroim2()
    v <- tryCatch(
      neuroim2::read_vol(path),
      error = function(e) {
        .cli_abort(
          sprintf("failed to read %s '%s': %s", label, path, conditionMessage(e)),
          class = "volregger_cli_domain_error"
        )
      }
    )
    arr <- as.array(v)
    if (length(dim(arr)) != 3L) {
      .cli_abort(
        sprintf("%s '%s' is %dD; a 3D mask is required", label, path, length(dim(arr))),
        class = "volregger_cli_domain_error"
      )
    }
    attr(arr, "pixdim") <- as.numeric(neuroim2::spacing(v))[seq_len(3L)]
    return(arr)
  }
  .cli_read_rds(path, label)
}

.cli_read_rds <- function(path, label) {
  if (!file.exists(path)) {
    .cli_abort(
      sprintf("%s file does not exist: %s", label, path),
      class = "volregger_cli_domain_error"
    )
  }
  tryCatch(
    readRDS(path),
    error = function(e) {
      .cli_abort(
        sprintf("failed to read %s '%s': %s", label, path, conditionMessage(e)),
        class = "volregger_cli_domain_error"
      )
    }
  )
}

.cli_write_nifti <- function(arr4d, path, template) {
  .cli_require_neuroim2()
  if (is.null(template)) {
    .cli_abort(
      "writing NIfTI requires a NIfTI template; pass --input as a .nii or .nii.gz file",
      class = "volregger_cli_domain_error"
    )
  }
  spc <- neuroim2::space(template)
  v <- neuroim2::NeuroVec(as.vector(as.array(arr4d)), spc)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  neuroim2::write_vec(v, path)
  invisible(path)
}

.cli_read_numeric_vector <- function(path, label, integer = FALSE) {
  if (!file.exists(path)) {
    .cli_abort(
      sprintf("%s file does not exist: %s", label, path),
      class = "volregger_cli_domain_error"
    )
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    val <- .cli_read_rds(path, label)
  } else {
    val <- tryCatch(
      scan(path, what = if (integer) integer() else numeric(), quiet = TRUE),
      error = function(e) {
        .cli_abort(
          sprintf("failed to read %s '%s': %s", label, path, conditionMessage(e)),
          class = "volregger_cli_domain_error"
        )
      }
    )
  }
  val <- as.vector(val)
  if (integer) val <- as.integer(val) else val <- as.numeric(val)
  if (any(is.na(val))) {
    .cli_abort(
      sprintf("%s '%s' contains NA values", label, path),
      class = "volregger_cli_domain_error"
    )
  }
  val
}

.cli_write_motion_tsv <- function(path, motion) {
  motion <- as.matrix(motion)
  if (ncol(motion) != 6L) {
    .cli_abort(
      sprintf("motion matrix has %d columns; expected 6", ncol(motion)),
      class = "volregger_cli_domain_error"
    )
  }
  colnames(motion) <- c("tx", "ty", "tz", "rx", "ry", "rz")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    motion,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  invisible(path)
}
