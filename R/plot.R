#' @export
plot.volregger_fit <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = c(2, 2))

  graphics::matplot(
    x$motion[, 1:3, drop = FALSE],
    type = "l",
    lty = 1,
    xlab = "Frame",
    ylab = "Translation (mm)",
    main = "Translations"
  )

  graphics::matplot(
    x$motion[, 4:6, drop = FALSE],
    type = "l",
    lty = 1,
    xlab = "Frame",
    ylab = "Rotation (rad)",
    main = "Rotations"
  )

  if (!is.null(x$qc$fd)) {
    graphics::plot(
      x$qc$fd,
      type = "l",
      xlab = "Frame",
      ylab = "FD",
      main = "Framewise displacement"
    )
    graphics::abline(h = 0.5, lty = 2, col = "firebrick")
    if (!is.null(x$qc$motion_spike)) {
      spike_idx <- which(as.logical(x$qc$motion_spike))
      if (length(spike_idx)) {
        graphics::points(spike_idx, x$qc$fd[spike_idx], col = "firebrick", pch = 19)
      }
    }
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No FD available")
  }

  if (!is.null(x$qc$cost_drop) || !is.null(x$qc$dvars)) {
    y <- if (!is.null(x$qc$cost_drop)) x$qc$cost_drop else x$qc$dvars
    ylab <- if (!is.null(x$qc$cost_drop)) "Cost drop" else "DVARS"
    main <- if (!is.null(x$qc$cost_drop)) "Fit diagnostics" else "DVARS"
    graphics::plot(y, type = "l", xlab = "Frame", ylab = ylab, main = main)
    if (!is.null(x$qc$dvars) && !is.null(x$qc$cost_drop)) {
      graphics::par(new = TRUE)
      graphics::plot(
        x$qc$dvars,
        type = "l",
        axes = FALSE,
        xlab = "",
        ylab = "",
        col = "steelblue"
      )
      graphics::axis(4, col.axis = "steelblue")
      graphics::mtext("DVARS", side = 4, line = 2.5, col = "steelblue")
    }
    if (!is.null(x$qc$fit_failure)) {
      fail_idx <- which(as.logical(x$qc$fit_failure))
      if (length(fail_idx)) {
        y_fail <- y[fail_idx]
        graphics::points(fail_idx, y_fail, col = "firebrick", pch = 19)
      }
    }
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No QC diagnostics available")
  }

  invisible(x)
}
