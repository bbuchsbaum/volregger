#' @export
plot.volregger_fit <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = c(2, 1))

  graphics::matplot(
    x$motion,
    type = "l",
    lty = 1,
    xlab = "Frame",
    ylab = "Motion parameters",
    main = "volregger motion traces"
  )

  if (!is.null(x$qc$fd)) {
    graphics::plot(
      x$qc$fd,
      type = "l",
      xlab = "Frame",
      ylab = "FD",
      main = "Framewise displacement"
    )
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No FD available")
  }

  invisible(x)
}
