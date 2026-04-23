#' Run the volregger command line interface
#'
#' @param args Character vector of command line arguments.
#'
#' @return Integer status code.
#' @export
volregger_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || args[[1]] %in% c("help", "--help", "-h")) {
    cat(.volregger_usage(args[-1]), sep = "\n")
    return(0L)
  }

  command <- args[[1]]
  rest <- args[-1]

  tryCatch(
    {
      if (identical(command, "report")) {
        return(.volregger_cli_report(rest))
      }
      .cli_abort(
        sprintf("unknown command '%s'", command),
        class = "volregger_cli_usage_error"
      )
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
#' @param dest_dir Destination directory.
#' @param overwrite If `TRUE`, overwrite existing files.
#' @param commands Character vector of wrapper names to install.
#'
#' @return Invisibly returns the installed file paths.
#' @export
install_cli <- function(dest_dir = "~/.local/bin", overwrite = FALSE, commands = "volregger") {
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

  invisible(dest_paths)
}

.volregger_cli_report <- function(args) {
  opts <- .parse_report_args(args)
  if (isTRUE(opts$help)) {
    cat(.volregger_usage("report"), sep = "\n")
    return(0L)
  }

  fit <- .read_cli_rds(opts$fit, "fit")
  if (!inherits(fit, "volregger_fit")) {
    .cli_abort("fit file must contain a 'volregger_fit' object", class = "volregger_cli_domain_error")
  }

  truth <- if (is.null(opts$truth)) NULL else .read_cli_rds(opts$truth, "truth")
  mask <- if (is.null(opts$mask)) NULL else .read_cli_rds(opts$mask, "mask")

  report <- volreg_report(fit, truth = truth, mask = mask)
  paths <- write_volreg_report(
    report = report,
    out_dir = opts$out_dir,
    prefix = opts$prefix,
    include_plot = opts$include_plot
  )

  cat(
    paste(names(paths), unname(paths), sep = "="),
    sep = "\n"
  )
  cat("\n")

  0L
}

.parse_report_args <- function(args) {
  opts <- list(
    fit = NULL,
    out_dir = ".",
    prefix = "volregger",
    truth = NULL,
    mask = NULL,
    include_plot = TRUE,
    help = FALSE
  )

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
      parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      key <- parts[[1]]
      value <- paste(parts[-1], collapse = "=")
      has_inline_value <- TRUE
    } else {
      key <- sub("^--", "", arg)
      value <- NULL
      has_inline_value <- FALSE
    }

    if (key %in% c("plot", "no-plot")) {
      opts$include_plot <- identical(key, "plot")
      i <- i + 1L
      next
    }

    if (!(key %in% c("fit", "out-dir", "prefix", "truth", "mask"))) {
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

    field <- sub("-", "_", key, fixed = TRUE)
    opts[[field]] <- value
    i <- i + 1L
  }

  if (!isTRUE(opts$help) && is.null(opts$fit)) {
    .cli_abort("option '--fit' is required", class = "volregger_cli_usage_error")
  }

  opts
}

.read_cli_rds <- function(path, label) {
  if (!file.exists(path)) {
    .cli_abort(
      sprintf("%s file does not exist: %s", label, path),
      class = "volregger_cli_domain_error"
    )
  }
  readRDS(path)
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

.volregger_usage <- function(command = NULL) {
  if (length(command) && identical(command[[1]], "report")) {
    return(c(
      "Usage:",
      "  volregger report --fit FIT.rds [--out-dir DIR] [--prefix NAME] [--truth TRUTH.rds] [--mask MASK.rds] [--plot|--no-plot]",
      "",
      "Description:",
      "  Build a structured report bundle from a saved volregger_fit object."
    ))
  }

  c(
    "Usage:",
    "  volregger help",
    "  volregger report --fit FIT.rds [--out-dir DIR] [--prefix NAME] [--truth TRUTH.rds] [--mask MASK.rds] [--plot|--no-plot]"
  )
}

.cli_abort <- function(message, class = "volregger_cli_domain_error") {
  err <- structure(
    list(message = message, call = NULL),
    class = c(class, "error", "condition")
  )
  stop(err)
}
