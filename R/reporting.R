#' Build a structured motion-correction report
#'
#' @param fit A `volregger_fit` object.
#' @param truth Optional ground-truth motion in any format accepted by
#'   [motion_matrices()].
#' @param mask Optional 3D mask used for truth-based displacement summaries.
#'   Defaults to `fit$mask`.
#' @param pixdim Optional voxel sizes in millimeters. Defaults to
#'   `attr(mask, "pixdim")` when available, otherwise `c(1, 1, 1)`.
#'
#' @return A `volregger_report` object with summary rows, per-frame metrics,
#'   motion matrices, and component labels.
#' @export
volreg_report <- function(fit, truth = NULL, mask = NULL, pixdim = NULL) {
  if (!inherits(fit, "volregger_fit")) {
    stop("fit must inherit from 'volregger_fit'", call. = FALSE)
  }

  frame_metrics <- as.data.frame(fit)
  mats <- motion_matrices(fit$mats)
  report_summary <- .volreg_report_summary(summary(fit), fit$control)

  if (!is.null(truth)) {
    if (is.null(mask)) {
      mask <- fit$mask
    }
    if (is.null(mask)) {
      stop("mask must be supplied when truth metrics are requested", call. = FALSE)
    }
    if (is.null(pixdim)) {
      pixdim <- attr(mask, "pixdim") %||% c(1, 1, 1)
    }

    truth_metrics <- motion_displacement_summary(
      mask = mask,
      motion = mats,
      truth = truth,
      pixdim = pixdim
    )
    frame_metrics <- merge(frame_metrics, truth_metrics, by = "frame", all.x = TRUE, sort = FALSE)
    report_summary <- utils::modifyList(report_summary, .truth_summary_row(truth_metrics))
  }

  structure(
    list(
      summary = as.data.frame(report_summary, stringsAsFactors = FALSE),
      frames = frame_metrics,
      matrices = .motion_matrix_data_frame(mats),
      components = .control_components(fit$control),
      engine = fit$control$engine %||% NA_character_
    ),
    class = "volregger_report"
  )
}

#' Write a structured motion-correction report
#'
#' @param report A `volregger_report` object.
#' @param out_dir Output directory.
#' @param prefix Output filename prefix.
#' @param include_plot If `TRUE`, write an overview PNG.
#' @param width Plot width in pixels.
#' @param height Plot height in pixels.
#' @param res Plot DPI.
#'
#' @return Invisibly returns a named character vector of written files.
#' @export
write_volreg_report <- function(
    report,
    out_dir,
    prefix = "volregger",
    include_plot = TRUE,
    width = 1200,
    height = 900,
    res = 120) {
  if (!inherits(report, "volregger_report")) {
    stop("report must inherit from 'volregger_report'", call. = FALSE)
  }

  out_dir <- path.expand(out_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- c(
    summary = file.path(out_dir, paste0(prefix, "_summary.csv")),
    frames = file.path(out_dir, paste0(prefix, "_frames.csv")),
    matrices = file.path(out_dir, paste0(prefix, "_matrices.csv")),
    text = file.path(out_dir, paste0(prefix, "_summary.txt")),
    rds = file.path(out_dir, paste0(prefix, "_report.rds"))
  )

  utils::write.csv(report$summary, paths[["summary"]], row.names = FALSE)
  utils::write.csv(report$frames, paths[["frames"]], row.names = FALSE)
  utils::write.csv(report$matrices, paths[["matrices"]], row.names = FALSE)
  writeLines(capture.output(print(report)), con = paths[["text"]])
  saveRDS(report, paths[["rds"]])

  if (isTRUE(include_plot)) {
    plot_path <- file.path(out_dir, paste0(prefix, "_overview.png"))
    grDevices::png(filename = plot_path, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(report)
    paths <- c(paths, overview = plot_path)
  }

  invisible(paths)
}

#' @export
print.volregger_report <- function(x, ...) {
  cat("<volregger_report>\n", sep = "")
  cat("  engine:", x$engine, "\n")
  cat("  components:", paste(x$components, collapse = ", "), "\n")
  cat("  frames:", nrow(x$frames), "\n")
  cat("  truth_metrics:", all(c("disp_median", "disp_p95", "disp_max") %in% names(x$frames)), "\n")
  invisible(x)
}

#' @export
plot.volregger_report <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = c(2, 2))

  graphics::matplot(
    x$frames[, c("tx", "ty", "tz"), drop = FALSE],
    type = "l",
    lty = 1,
    xlab = "Frame",
    ylab = "Translation (mm)",
    main = "Translations"
  )

  graphics::matplot(
    x$frames[, c("rx", "ry", "rz"), drop = FALSE],
    type = "l",
    lty = 1,
    xlab = "Frame",
    ylab = "Rotation (rad)",
    main = "Rotations"
  )

  if (all(c("disp_median", "disp_p95", "disp_max") %in% names(x$frames))) {
    graphics::matplot(
      x$frames[, c("disp_median", "disp_p95", "disp_max"), drop = FALSE],
      type = "l",
      lty = 1,
      col = c("black", "steelblue", "firebrick"),
      xlab = "Frame",
      ylab = "Displacement (mm)",
      main = "Mask Displacement"
    )
    graphics::legend(
      "topright",
      legend = c("median", "p95", "max"),
      col = c("black", "steelblue", "firebrick"),
      lty = 1,
      bty = "n"
    )
  } else if ("fd" %in% names(x$frames)) {
    graphics::plot(
      x$frames$fd,
      type = "l",
      xlab = "Frame",
      ylab = "FD",
      main = "Framewise Displacement"
    )
    graphics::abline(h = 0.5, lty = 2, col = "firebrick")
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No displacement metrics available")
  }

  .plot_qc_panel(x$frames)

  invisible(x)
}

#' @export
as.data.frame.volregger_report <- function(x, ...) {
  x$frames
}

.plot_qc_panel <- function(frame_metrics) {
  if ("fd" %in% names(frame_metrics)) {
    graphics::plot(
      frame_metrics$fd,
      type = "l",
      xlab = "Frame",
      ylab = "FD",
      main = "QC Overview"
    )
    if ("fd_error" %in% names(frame_metrics)) {
      graphics::lines(frame_metrics$fd_error, col = "steelblue")
      graphics::legend(
        "topright",
        legend = c("fd", "fd_error"),
        col = c("black", "steelblue"),
        lty = 1,
        bty = "n"
      )
    }
    if ("motion_spike" %in% names(frame_metrics)) {
      spike_idx <- which(as.logical(frame_metrics$motion_spike))
      if (length(spike_idx)) {
        graphics::points(spike_idx, frame_metrics$fd[spike_idx], col = "firebrick", pch = 19)
      }
    }
    return(invisible(NULL))
  }

  if ("dvars" %in% names(frame_metrics) || "cost_drop" %in% names(frame_metrics)) {
    y <- if ("cost_drop" %in% names(frame_metrics)) frame_metrics$cost_drop else frame_metrics$dvars
    ylab <- if ("cost_drop" %in% names(frame_metrics)) "Cost drop" else "DVARS"
    graphics::plot(y, type = "l", xlab = "Frame", ylab = ylab, main = "QC Overview")
    if ("dvars" %in% names(frame_metrics) && "cost_drop" %in% names(frame_metrics)) {
      graphics::par(new = TRUE)
      graphics::plot(
        frame_metrics$dvars,
        type = "l",
        axes = FALSE,
        xlab = "",
        ylab = "",
        col = "steelblue"
      )
      graphics::axis(4, col.axis = "steelblue")
      graphics::mtext("DVARS", side = 4, line = 2.5, col = "steelblue")
    }
    return(invisible(NULL))
  }

  graphics::plot.new()
  graphics::text(0.5, 0.5, "No QC metrics available")
  invisible(NULL)
}

.volreg_report_summary <- function(fit_summary, control) {
  motion_abs_max <- fit_summary$motion_abs_max
  data.frame(
    engine = fit_summary$engine,
    components = paste(.control_components(control), collapse = ","),
    n_frames = fit_summary$n_frames,
    tx_abs_max = unname(motion_abs_max[1]),
    ty_abs_max = unname(motion_abs_max[2]),
    tz_abs_max = unname(motion_abs_max[3]),
    rx_abs_max = unname(motion_abs_max[4]),
    ry_abs_max = unname(motion_abs_max[5]),
    rz_abs_max = unname(motion_abs_max[6]),
    fd_mean = fit_summary$fd_mean,
    fd_p95 = fit_summary$fd_p95,
    fd_max = fit_summary$fd_max,
    dvars_mean = fit_summary$dvars_mean,
    dvars_p95 = fit_summary$dvars_p95,
    robust_dvars_mean = fit_summary$robust_dvars_mean,
    robust_dvars_p95 = fit_summary$robust_dvars_p95,
    motion_spikes = fit_summary$motion_spikes,
    fit_failures = fit_summary$fit_failures,
    censor_suggested = fit_summary$censor_suggested,
    cost_drop_mean = fit_summary$cost_drop_mean,
    packet_correction_mag_mean = fit_summary$packet_correction_mag_mean,
    packet_correction_mag_max = fit_summary$packet_correction_mag_max,
    stringsAsFactors = FALSE
  )
}

.truth_summary_row <- function(truth_metrics) {
  data.frame(
    truth_disp_median_median = stats::median(truth_metrics$disp_median, na.rm = TRUE),
    truth_disp_p95_median = stats::median(truth_metrics$disp_p95, na.rm = TRUE),
    truth_disp_max_max = max(truth_metrics$disp_max, na.rm = TRUE),
    truth_fd_error_median = stats::median(truth_metrics$fd_error, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

.control_components <- function(control) {
  components <- "rigid"
  if (isTRUE(control$robust_template)) {
    components <- c(components, "robust_template")
  }
  if (isTRUE(control$use_overlap_taper)) {
    components <- c(components, "overlap_taper")
  }
  if (isTRUE(control$capture_boost)) {
    components <- c(components, "capture_boost")
  }
  if (isTRUE(control$use_ic_stencil)) {
    components <- c(components, "ic_stencil")
  }
  if (isTRUE(control$whiten_modes)) {
    components <- c(components, "whitening")
  }
  if (isTRUE(control$use_nuisance_marginalization)) {
    components <- c(components, "nuisance")
  }
  if (isTRUE(control$nuisance_include_template)) {
    components <- c(components, "template_basis")
  }
  if (isTRUE(control$use_slice_spline) || identical(control$engine, "rigid_spline")) {
    components <- c(components, "slice_spline")
  }
  unique(components)
}

.motion_matrix_data_frame <- function(mats) {
  nt <- dim(mats)[3]
  out <- matrix(NA_real_, nrow = nt, ncol = 16L)
  for (i in seq_len(nt)) {
    out[i, ] <- as.numeric(mats[, , i])
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- sprintf(
    "m%d%d",
    rep(seq_len(4), each = 4),
    rep(seq_len(4), times = 4)
  )
  out$frame <- seq_len(nt)
  out[, c("frame", setdiff(names(out), "frame")), drop = FALSE]
}
