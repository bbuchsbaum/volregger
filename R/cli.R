#' Run the volregger command line interface
#'
#' Programmatic entrypoint for the `volregger` shell wrapper. Tests and
#' embedding callers should invoke this directly with an explicit `args`
#' vector; the wrapper script handles `commandArgs()` parsing and exit
#' status propagation.
#'
#' Available subcommands: `estimate`, `apply`, `run`, `report`. Run
#' `volregger help <command>` for per-command documentation.
#'
#' @param args Character vector of command line arguments.
#'
#' @return Integer status code: `0` on success, `1` on a domain failure
#'   (e.g. invalid input file), `2` on a usage or runtime error.
#' @examples
#' volregger_cli("--version")
#' \dontrun{
#' volregger_cli(c("estimate", "--input", "run.nii.gz", "--out-fit", "fit.rds"))
#' }
#' @export
volregger_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args)) {
    cat(.volregger_usage(), sep = "\n")
    return(0L)
  }

  if (args[[1L]] %in% c("--version", "-V")) {
    cat(.volregger_version_string(), "\n", sep = "")
    return(0L)
  }

  if (args[[1L]] %in% c("--help", "-h")) {
    cat(.volregger_usage(), sep = "\n")
    return(0L)
  }

  if (identical(args[[1L]], "help")) {
    rest <- args[-1L]
    if (!length(rest)) {
      cat(.volregger_usage(), sep = "\n")
      return(0L)
    }
    if (!(rest[[1L]] %in% .volregger_commands())) {
      message("volregger: no help topic for '", rest[[1L]], "'")
      cat(.volregger_usage(), sep = "\n", file = stderr())
      return(2L)
    }
    cat(.volregger_usage(rest[[1L]]), sep = "\n")
    return(0L)
  }

  command <- args[[1L]]
  rest <- args[-1L]

  tryCatch(
    {
      handler <- switch(
        command,
        estimate = .volregger_cli_estimate,
        apply    = .volregger_cli_apply,
        run      = .volregger_cli_run,
        report   = .volregger_cli_report,
        NULL
      )
      if (is.null(handler)) {
        .cli_abort(
          sprintf("unknown command '%s' (try 'volregger help')", command),
          class = "volregger_cli_usage_error"
        )
      }
      handler(rest)
    },
    volregger_cli_usage_error = function(e) {
      message("volregger: ", conditionMessage(e))
      cat(.volregger_usage(command), sep = "\n", file = stderr())
      2L
    },
    volregger_cli_domain_error = function(e) {
      message("volregger: ", conditionMessage(e))
      1L
    },
    error = function(e) {
      message("volregger: ", conditionMessage(e))
      2L
    }
  )
}

#' Install the packaged volregger command wrapper
#'
#' @param dest_dir Destination directory. Created if it does not exist.
#' @param overwrite If `TRUE`, overwrite existing files in `dest_dir`.
#' @param commands Character vector of wrapper names to install.
#' @param quiet If `TRUE`, suppress informational messages.
#'
#' @return Invisibly returns the installed file paths.
#' @examples
#' \dontrun{
#' install_cli("~/.local/bin", overwrite = TRUE)
#' }
#' @export
install_cli <- function(dest_dir = "~/.local/bin",
                        overwrite = FALSE,
                        commands = "volregger",
                        quiet = FALSE) {
  exec_dir <- .find_exec_dir()
  dest_dir <- path.expand(dest_dir)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  src_paths <- file.path(exec_dir, commands)
  if (any(!file.exists(src_paths))) {
    missing <- commands[!file.exists(src_paths)]
    stop("missing wrapper(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  dest_paths <- file.path(dest_dir, basename(src_paths))
  if (!isTRUE(overwrite) && any(file.exists(dest_paths))) {
    existing <- basename(dest_paths[file.exists(dest_paths)])
    stop("refusing to overwrite existing wrapper(s): ", paste(existing, collapse = ", "), call. = FALSE)
  }

  ok <- file.copy(src_paths, dest_paths, overwrite = overwrite, copy.mode = TRUE)
  if (!all(ok)) {
    stop("failed to copy one or more CLI wrappers", call. = FALSE)
  }
  try(Sys.chmod(dest_paths, mode = "755"), silent = TRUE)

  if (!isTRUE(quiet)) {
    message("Installed: ", paste(dest_paths, collapse = ", "))
    if (!.dir_on_path(dest_dir)) {
      message(
        "Note: '", dest_dir, "' is not on $PATH. ",
        "Add `export PATH=\"", dest_dir, ":$PATH\"` to your shell profile."
      )
    }
  }

  invisible(dest_paths)
}

.volregger_commands <- function() c("estimate", "apply", "run", "report")

# ---------------------------------------------------------------------------
# Command: estimate
# ---------------------------------------------------------------------------

.volregger_cli_estimate <- function(args) {
  spec <- list(
    string = c("input", "mask", "engine", "profile", "ref", "ref-index",
               "slice-times", "mb-groups", "threads",
               "out-fit", "out-motion"),
    toggle = list(
      qc      = c("qc", "no-qc"),
      verbose = c("verbose", "no-verbose"),
      json    = c("json", "no-json"),
      quiet   = c("quiet", "no-quiet")
    ),
    defaults = list(
      ref = "middle", engine = "rigid_robust", qc = TRUE,
      verbose = FALSE, json = FALSE, quiet = FALSE
    ),
    required = c("input", "out-fit"),
    help_command = "estimate"
  )
  opts <- .cli_parse(args, spec)
  if (isTRUE(opts$help)) {
    cat(.volregger_usage("estimate"), sep = "\n")
    return(0L)
  }

  x <- .cli_read_4d(opts$input, "input")
  mask <- if (is.null(opts$mask)) NULL else .cli_read_3d(opts$mask, "mask")
  slice_times <- if (is.null(opts[["slice-times"]])) NULL else
    .cli_read_numeric_vector(opts[["slice-times"]], "slice-times")
  mb_groups <- if (is.null(opts[["mb-groups"]])) NULL else
    .cli_read_numeric_vector(opts[["mb-groups"]], "mb-groups", integer = TRUE)

  control <- volreg_control()
  if (!is.null(opts$threads)) {
    control$nthreads <- .cli_parse_positive_int(opts$threads, "--threads")
  }

  ref_index <- if (is.null(opts[["ref-index"]])) NULL else
    .cli_parse_positive_int(opts[["ref-index"]], "--ref-index")

  fit <- volreg(
    x,
    mask = mask,
    ref = opts$ref,
    ref_index = ref_index,
    engine = opts$engine,
    slice_times = slice_times,
    mb_groups = mb_groups,
    resample = FALSE,
    qc = isTRUE(opts$qc),
    control = control,
    profile = opts$profile,
    verbose = isTRUE(opts$verbose)
  )

  out_fit <- opts[["out-fit"]]
  dir.create(dirname(out_fit), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, out_fit)

  written <- c(fit = out_fit)
  if (!is.null(opts[["out-motion"]])) {
    motion_path <- .cli_write_motion_tsv(opts[["out-motion"]], fit$motion)
    written <- c(written, motion = motion_path)
  }

  .cli_emit_paths(written, json = isTRUE(opts$json), quiet = isTRUE(opts$quiet))
  0L
}

# ---------------------------------------------------------------------------
# Command: apply
# ---------------------------------------------------------------------------

.volregger_cli_apply <- function(args) {
  spec <- list(
    string = c("input", "fit", "out", "interp", "pad", "zpad", "threads"),
    toggle = list(
      json  = c("json", "no-json"),
      quiet = c("quiet", "no-quiet")
    ),
    defaults = list(
      interp = "bspline", pad = "clamp", zpad = "0",
      json = FALSE, quiet = FALSE
    ),
    required = c("input", "fit", "out"),
    help_command = "apply"
  )
  opts <- .cli_parse(args, spec)
  if (isTRUE(opts$help)) {
    cat(.volregger_usage("apply"), sep = "\n")
    return(0L)
  }

  x <- .cli_read_4d(opts$input, "input")
  fit <- .cli_read_rds(opts$fit, "fit")
  if (!inherits(fit, "volregger_fit")) {
    .cli_abort("fit file must contain a 'volregger_fit' object",
               class = "volregger_cli_domain_error")
  }

  threads <- if (is.null(opts$threads)) 1L else
    .cli_parse_positive_int(opts$threads, "--threads")
  zpad <- .cli_parse_nonneg_int(opts$zpad, "--zpad")

  corrected <- apply_motion(
    x,
    motion = fit$motion,
    interp = opts$interp,
    pad = opts$pad,
    zpad = zpad,
    nthreads = threads
  )

  template <- attr(x, ".neuroim2_template")
  if (.cli_is_nifti_path(opts$out)) {
    .cli_write_nifti(corrected, opts$out, template)
  } else {
    dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
    saveRDS(corrected, opts$out)
  }

  .cli_emit_paths(c(corrected = opts$out), json = isTRUE(opts$json),
                  quiet = isTRUE(opts$quiet))
  0L
}

# ---------------------------------------------------------------------------
# Command: run (estimate + apply + report)
# ---------------------------------------------------------------------------

.volregger_cli_run <- function(args) {
  spec <- list(
    string = c("input", "mask", "out-dir", "prefix",
               "engine", "profile", "ref", "ref-index",
               "slice-times", "mb-groups",
               "interp", "pad", "zpad", "threads"),
    toggle = list(
      qc        = c("qc", "no-qc"),
      plot      = c("plot", "no-plot"),
      corrected = c("corrected", "no-corrected"),
      report    = c("report", "no-report"),
      verbose   = c("verbose", "no-verbose"),
      json      = c("json", "no-json"),
      quiet     = c("quiet", "no-quiet")
    ),
    defaults = list(
      ref = "middle", engine = "rigid_robust",
      interp = "bspline", pad = "clamp", zpad = "0",
      qc = TRUE, plot = TRUE, corrected = TRUE, report = TRUE,
      verbose = FALSE, json = FALSE, quiet = FALSE
    ),
    required = c("input", "out-dir"),
    help_command = "run"
  )
  opts <- .cli_parse(args, spec)
  if (isTRUE(opts$help)) {
    cat(.volregger_usage("run"), sep = "\n")
    return(0L)
  }

  out_dir <- opts[["out-dir"]]
  prefix <- opts$prefix %||% .cli_default_prefix(opts$input)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  x <- .cli_read_4d(opts$input, "input")
  mask <- if (is.null(opts$mask)) NULL else .cli_read_3d(opts$mask, "mask")
  slice_times <- if (is.null(opts[["slice-times"]])) NULL else
    .cli_read_numeric_vector(opts[["slice-times"]], "slice-times")
  mb_groups <- if (is.null(opts[["mb-groups"]])) NULL else
    .cli_read_numeric_vector(opts[["mb-groups"]], "mb-groups", integer = TRUE)

  control <- volreg_control()
  if (!is.null(opts$threads)) {
    control$nthreads <- .cli_parse_positive_int(opts$threads, "--threads")
  }
  ref_index <- if (is.null(opts[["ref-index"]])) NULL else
    .cli_parse_positive_int(opts[["ref-index"]], "--ref-index")
  zpad <- .cli_parse_nonneg_int(opts$zpad, "--zpad")

  fit <- volreg(
    x,
    mask = mask,
    ref = opts$ref,
    ref_index = ref_index,
    engine = opts$engine,
    slice_times = slice_times,
    mb_groups = mb_groups,
    resample = isTRUE(opts$corrected),
    interp_final = opts$interp,
    pad_final = opts$pad,
    zpad_final = zpad,
    qc = isTRUE(opts$qc),
    control = control,
    profile = opts$profile,
    verbose = isTRUE(opts$verbose)
  )

  written <- c(
    fit    = file.path(out_dir, paste0(prefix, "_fit.rds")),
    motion = file.path(out_dir, paste0(prefix, "_motion.tsv"))
  )
  saveRDS(fit, written[["fit"]])
  .cli_write_motion_tsv(written[["motion"]], fit$motion)

  if (isTRUE(opts$corrected) && !is.null(fit$corrected)) {
    template <- attr(x, ".neuroim2_template")
    has_nifti_input <- !is.null(template)
    corrected_path <- file.path(
      out_dir,
      paste0(prefix, "_corrected", if (has_nifti_input) ".nii.gz" else ".rds")
    )
    if (has_nifti_input) {
      .cli_write_nifti(fit$corrected, corrected_path, template)
    } else {
      saveRDS(fit$corrected, corrected_path)
    }
    written <- c(written, corrected = corrected_path)
  }

  if (isTRUE(opts$report)) {
    report <- volreg_report(fit, mask = mask)
    paths <- write_volreg_report(
      report = report,
      out_dir = out_dir,
      prefix = prefix,
      include_plot = isTRUE(opts$plot)
    )
    written <- c(written, paths)
  }

  .cli_emit_paths(written, json = isTRUE(opts$json), quiet = isTRUE(opts$quiet))
  0L
}

# ---------------------------------------------------------------------------
# Command: report (existing)
# ---------------------------------------------------------------------------

.volregger_cli_report <- function(args) {
  spec <- list(
    string = c("fit", "out-dir", "prefix", "truth", "mask"),
    toggle = list(
      plot  = c("plot", "no-plot"),
      json  = c("json", "no-json"),
      quiet = c("quiet", "no-quiet")
    ),
    defaults = list(
      "out-dir" = ".", prefix = "volregger",
      plot = TRUE, json = FALSE, quiet = FALSE
    ),
    required = c("fit"),
    help_command = "report"
  )
  opts <- .cli_parse(args, spec)
  if (isTRUE(opts$help)) {
    cat(.volregger_usage("report"), sep = "\n")
    return(0L)
  }

  fit <- .cli_read_rds(opts$fit, "fit")
  if (!inherits(fit, "volregger_fit")) {
    .cli_abort("fit file must contain a 'volregger_fit' object",
               class = "volregger_cli_domain_error")
  }

  truth <- if (is.null(opts$truth)) NULL else .cli_read_rds(opts$truth, "truth")
  mask <- if (is.null(opts$mask)) NULL else .cli_read_3d(opts$mask, "mask")

  report <- volreg_report(fit, truth = truth, mask = mask)
  paths <- write_volreg_report(
    report = report,
    out_dir = opts[["out-dir"]],
    prefix = opts$prefix,
    include_plot = isTRUE(opts$plot)
  )

  .cli_emit_paths(paths, json = isTRUE(opts$json), quiet = isTRUE(opts$quiet))
  0L
}

# ---------------------------------------------------------------------------
# Generic parser
# ---------------------------------------------------------------------------

.cli_parse <- function(args, spec) {
  opts <- list()
  for (key in names(spec$defaults)) opts[[key]] <- spec$defaults[[key]]
  opts$json <- opts$json %||% FALSE
  opts$quiet <- opts$quiet %||% FALSE
  opts$help <- FALSE

  toggle_key_index <- list()
  for (canonical in names(spec$toggle)) {
    keys <- spec$toggle[[canonical]]
    on_key <- keys[[1L]]
    off_key <- if (length(keys) >= 2L) keys[[2L]] else NA_character_
    toggle_key_index[[on_key]] <- list(canonical = canonical, value = TRUE)
    if (!is.na(off_key)) {
      toggle_key_index[[off_key]] <- list(canonical = canonical, value = FALSE)
    }
  }

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("--help", "-h")) {
      opts$help <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(arg, "--")) {
      .cli_abort(
        sprintf("unexpected positional argument '%s'", arg),
        class = "volregger_cli_usage_error"
      )
    }

    if (grepl("=", arg, fixed = TRUE)) {
      parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
      key <- parts[[1L]]
      value <- paste(parts[-1L], collapse = "=")
      has_inline_value <- TRUE
    } else {
      key <- sub("^--", "", arg)
      value <- NULL
      has_inline_value <- FALSE
    }

    if (!is.null(toggle_key_index[[key]])) {
      if (has_inline_value) {
        .cli_abort(
          sprintf("flag '--%s' does not take a value", key),
          class = "volregger_cli_usage_error"
        )
      }
      tk <- toggle_key_index[[key]]
      opts[[tk$canonical]] <- tk$value
      i <- i + 1L
      next
    }

    if (!(key %in% spec$string)) {
      .cli_abort(
        sprintf("unknown option '--%s'", key),
        class = "volregger_cli_usage_error"
      )
    }

    if (!has_inline_value) {
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        .cli_abort(
          sprintf("option '--%s' requires a value", key),
          class = "volregger_cli_usage_error"
        )
      }
      value <- args[[i + 1L]]
      i <- i + 1L
    }

    opts[[key]] <- value
    i <- i + 1L
  }

  if (!isTRUE(opts$help)) {
    for (req in spec$required) {
      if (is.null(opts[[req]]) || identical(opts[[req]], "")) {
        .cli_abort(
          sprintf("option '--%s' is required", req),
          class = "volregger_cli_usage_error"
        )
      }
    }
  }

  opts
}

.cli_parse_positive_int <- function(value, label) {
  n <- suppressWarnings(as.integer(value))
  if (length(n) != 1L || is.na(n) || n <= 0L) {
    .cli_abort(
      sprintf("%s must be a positive integer (got '%s')", label, value),
      class = "volregger_cli_usage_error"
    )
  }
  n
}

.cli_parse_nonneg_int <- function(value, label) {
  n <- suppressWarnings(as.integer(value))
  if (length(n) != 1L || is.na(n) || n < 0L) {
    .cli_abort(
      sprintf("%s must be a non-negative integer (got '%s')", label, value),
      class = "volregger_cli_usage_error"
    )
  }
  n
}

.cli_default_prefix <- function(input_path) {
  base <- basename(input_path)
  base <- sub("\\.nii\\.gz$", "", base, ignore.case = TRUE)
  base <- sub("\\.nii$", "", base, ignore.case = TRUE)
  base <- sub("\\.rds$", "", base, ignore.case = TRUE)
  if (!nzchar(base)) "volregger" else base
}

.cli_emit_paths <- function(paths, json = FALSE, quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))
  if (isTRUE(json)) {
    cat(.json_object(as.list(paths)), "\n", sep = "")
  } else {
    cat(paste(names(paths), unname(paths), sep = "="), sep = "\n")
    cat("\n")
  }
  invisible(NULL)
}

.find_exec_dir <- function() {
  exec_dir <- system.file("exec", package = "volregger")
  if (nzchar(exec_dir)) {
    return(exec_dir)
  }

  local_exec <- normalizePath(file.path(getwd(), "exec"), mustWork = FALSE)
  if (dir.exists(local_exec)) {
    return(local_exec)
  }

  stop(
    "could not locate exec/; call install_cli() from an installed package or the source checkout root",
    call. = FALSE
  )
}

.dir_on_path <- function(dir) {
  path_env <- Sys.getenv("PATH", unset = "")
  if (!nzchar(path_env)) return(FALSE)
  parts <- strsplit(path_env, .Platform$path.sep, fixed = TRUE)[[1L]]
  target <- normalizePath(dir, mustWork = FALSE)
  any(vapply(parts, function(p) {
    isTRUE(normalizePath(p, mustWork = FALSE) == target)
  }, logical(1L)))
}

.volregger_version_string <- function() {
  ver <- tryCatch(as.character(utils::packageVersion("volregger")), error = function(e) "unknown")
  paste0("volregger ", ver)
}

# ---------------------------------------------------------------------------
# Help text
# ---------------------------------------------------------------------------

.volregger_usage <- function(command = NULL) {
  if (length(command) && identical(command[[1L]], "estimate")) {
    return(c(
      "Usage:",
      "  volregger estimate --input INPUT --out-fit FIT.rds [options]",
      "",
      "Estimate rigid-body motion parameters for a 4D run.",
      "",
      "Required:",
      "  --input FILE        4D input (.nii / .nii.gz / .rds)",
      "  --out-fit FILE      where to save the volregger_fit RDS",
      "",
      "Options:",
      "  --mask FILE         3D mask (.nii / .nii.gz / .rds)",
      "  --engine NAME       rigid_robust (default) or rigid_spline",
      "  --profile NAME      preset name passed to volreg() (e.g. fast_fmri)",
      "  --ref MODE          middle (default), robust_mean, or index",
      "  --ref-index N       1-based reference frame index (when --ref index)",
      "  --slice-times FILE  per-slice timing offsets (.rds or whitespace text)",
      "  --mb-groups FILE    per-slice multiband group ids (.rds or text)",
      "  --threads N         control$nthreads (default: 1)",
      "  --qc, --no-qc       compute QC metrics (default: --qc)",
      "  --out-motion FILE   also write motion params as a 6-column TSV",
      "  --verbose           verbose volreg() logging",
      "  --quiet             suppress stdout summary",
      "  --json              emit artifact paths as a JSON object",
      "  -h, --help          show this help",
      "",
      "Examples:",
      "  volregger estimate --input run.nii.gz --out-fit fit.rds",
      "  volregger estimate --input run.nii.gz --mask mask.nii.gz \\",
      "                     --engine rigid_spline --slice-times slice_times.txt \\",
      "                     --threads 4 --out-fit fit.rds --out-motion motion.tsv"
    ))
  }

  if (length(command) && identical(command[[1L]], "apply")) {
    return(c(
      "Usage:",
      "  volregger apply --input INPUT --fit FIT.rds --out OUT [options]",
      "",
      "Apply a saved fit to a 4D run, writing the corrected volume.",
      "",
      "Required:",
      "  --input FILE        4D input (.nii / .nii.gz / .rds)",
      "  --fit FILE          saved volregger_fit RDS (from `volregger estimate`)",
      "  --out FILE          output path (.nii / .nii.gz writes a NIfTI;",
      "                      any other extension writes an RDS array)",
      "",
      "Options:",
      "  --interp METHOD     bspline (default), lanczos, lanczos5, lanczos7, linear",
      "  --pad MODE          clamp (default) or zero",
      "  --zpad N            symmetric zero-padding before resampling (default: 0)",
      "  --threads N         number of threads for apply (default: 1)",
      "  --quiet             suppress stdout summary",
      "  --json              emit artifact paths as a JSON object",
      "  -h, --help          show this help",
      "",
      "Examples:",
      "  volregger apply --input run.nii.gz --fit fit.rds --out corrected.nii.gz",
      "  volregger apply --input run.nii.gz --fit fit.rds --out corrected.nii.gz \\",
      "                  --interp lanczos --threads 4"
    ))
  }

  if (length(command) && identical(command[[1L]], "run")) {
    return(c(
      "Usage:",
      "  volregger run --input INPUT --out-dir DIR [options]",
      "",
      "End-to-end pipeline: estimate motion, write corrected NIfTI,",
      "save fit, motion TSV, and a structured report bundle.",
      "",
      "Required:",
      "  --input FILE        4D input (.nii / .nii.gz / .rds)",
      "  --out-dir DIR       output directory (created if missing)",
      "",
      "Options:",
      "  --prefix NAME       output filename prefix (default: from --input)",
      "  --mask FILE         3D mask (.nii / .nii.gz / .rds)",
      "  --engine NAME       rigid_robust (default) or rigid_spline",
      "  --profile NAME      preset name passed to volreg()",
      "  --ref MODE          middle (default), robust_mean, or index",
      "  --ref-index N       1-based reference index (when --ref index)",
      "  --slice-times FILE  per-slice timing offsets",
      "  --mb-groups FILE    per-slice multiband group ids",
      "  --interp METHOD     bspline (default), lanczos, lanczos5, lanczos7, linear",
      "  --pad MODE          clamp (default) or zero",
      "  --zpad N            symmetric zero-padding (default: 0)",
      "  --threads N         worker threads (default: 1)",
      "  --qc, --no-qc       QC metrics (default: --qc)",
      "  --plot, --no-plot   overview PNG (default: --plot)",
      "  --report, --no-report   write report bundle (default: --report)",
      "  --corrected, --no-corrected   write corrected output (default: --corrected)",
      "  --verbose           verbose logging",
      "  --json              emit artifact paths as a JSON object",
      "  --quiet             suppress stdout summary",
      "  -h, --help          show this help",
      "",
      "Outputs (under --out-dir, named with --prefix):",
      "  <prefix>_corrected.nii.gz  resampled 4D run (if NIfTI input)",
      "  <prefix>_fit.rds           saved volregger_fit",
      "  <prefix>_motion.tsv        6-column motion params",
      "  <prefix>_summary.csv       single-row run summary",
      "  <prefix>_frames.csv        per-frame metrics",
      "  <prefix>_matrices.csv      per-frame 4x4 motion matrices",
      "  <prefix>_summary.txt       printed report",
      "  <prefix>_report.rds        report object",
      "  <prefix>_overview.png      QC overview (unless --no-plot)",
      "",
      "Examples:",
      "  volregger run --input run.nii.gz --out-dir out/",
      "  volregger run --input run.nii.gz --mask brain.nii.gz \\",
      "                --out-dir out/ --prefix sub01 --threads 4 --json"
    ))
  }

  if (length(command) && identical(command[[1L]], "report")) {
    return(c(
      "Usage:",
      "  volregger report --fit FIT.rds [options]",
      "",
      "Build a structured report bundle from a saved volregger_fit object.",
      "Writes summary.csv, frames.csv, matrices.csv, summary.txt, report.rds,",
      "and (unless --no-plot) overview.png under --out-dir.",
      "",
      "Required:",
      "  --fit FILE          saved volregger_fit RDS",
      "",
      "Options:",
      "  --out-dir DIR       output directory (default: .)",
      "  --prefix NAME       output filename prefix (default: volregger)",
      "  --truth FILE        ground-truth motion RDS (4x4xT or pose matrix)",
      "  --mask FILE         3D mask (.nii / .nii.gz / .rds); required for",
      "                      truth-based metrics if the fit has no embedded mask",
      "  --plot, --no-plot   write/skip the overview PNG (default: --plot)",
      "  --json              emit artifact paths as a single JSON object",
      "  --quiet             suppress stdout (still writes files)",
      "  -h, --help          show this help",
      "",
      "Examples:",
      "  volregger report --fit fit.rds --out-dir report/",
      "  volregger report --fit fit.rds --truth truth.rds --mask mask.nii.gz \\",
      "                   --prefix subject01 --json"
    ))
  }

  c(
    "Usage:",
    "  volregger <command> [options]",
    "  volregger help [<command>]",
    "  volregger --version",
    "",
    "Commands:",
    "  estimate  Estimate motion parameters for a 4D run",
    "  apply     Apply a saved fit to a 4D run, writing corrected output",
    "  run       End-to-end pipeline: estimate + apply + report",
    "  report    Build a report bundle from a saved fit",
    "",
    "Options:",
    "  -h, --help     Show this help",
    "  -V, --version  Print package version",
    "",
    "Exit status:",
    "  0  success",
    "  1  domain failure (e.g. invalid input file)",
    "  2  usage or runtime error",
    "",
    "NIfTI I/O (.nii / .nii.gz) requires the 'neuroim2' package.",
    "See 'volregger help <command>' for command-specific help."
  )
}

# ---------------------------------------------------------------------------
# JSON + helpers
# ---------------------------------------------------------------------------

.json_object <- function(x) {
  if (!length(x)) return("{}")
  parts <- vapply(seq_along(x), function(i) {
    sprintf('"%s":"%s"', .json_escape(names(x)[[i]]), .json_escape(as.character(x[[i]])))
  }, character(1L))
  paste0("{", paste(parts, collapse = ","), "}")
}

.json_escape <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  s
}

.cli_abort <- function(message, class = "volregger_cli_domain_error") {
  err <- structure(
    list(message = message, call = NULL),
    class = c(class, "error", "condition")
  )
  stop(err)
}
