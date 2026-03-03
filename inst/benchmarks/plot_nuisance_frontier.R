#!/usr/bin/env Rscript

in_csv <- file.path("inst", "benchmarks", "nuisance_basis_frontier_summary.csv")
if (!file.exists(in_csv)) {
  stop("Missing input CSV: ", in_csv, "\nRun benchmark_nuisance_basis_frontier.R first.", call. = FALSE)
}

df <- read.csv(in_csv, stringsAsFactors = FALSE)
out_dir <- file.path("inst", "benchmarks", "figures")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

scenarios <- unique(df$scenario)
methods <- c("off", "small", "full")
cols <- c(off = "#4D4D4D", small = "#1b9e77", full = "#d95f02")
pchs <- c(off = 15, small = 16, full = 17)

plot_metric <- function(metric_col, ylab, outfile) {
  png(outfile, width = 1500, height = 520)
  op <- par(mfrow = c(1, length(scenarios)), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 1, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  yvals <- df[[metric_col]]
  yr <- range(c(0.75, 1.25, yvals), finite = TRUE)
  yr[1] <- min(yr[1], 0.75)
  yr[2] <- max(yr[2], 1.25)

  for (sc in scenarios) {
    d <- df[df$scenario == sc, ]
    d <- d[match(methods, d$method), ]

    plot(
      d$elapsed_ratio_median,
      d[[metric_col]],
      xlim = c(0.95, max(1.55, max(d$elapsed_ratio_median, na.rm = TRUE) + 0.05)),
      ylim = yr,
      pch = pchs[d$method],
      col = cols[d$method],
      cex = 1.3,
      xlab = "Runtime ratio vs off",
      ylab = ylab,
      main = sc
    )
    abline(v = 1, h = 1, lty = 2, col = "gray60")
    text(
      d$elapsed_ratio_median,
      d[[metric_col]],
      labels = d$method,
      pos = 3,
      cex = 0.9,
      col = cols[d$method]
    )

    # Draw improvement arrows from off -> small/full.
    off <- d[d$method == "off", ]
    for (m in c("small", "full")) {
      mm <- d[d$method == m, ]
      arrows(
        x0 = off$elapsed_ratio_median, y0 = off[[metric_col]],
        x1 = mm$elapsed_ratio_median, y1 = mm[[metric_col]],
        length = 0.08, col = cols[m], lwd = 1.5
      )
    }
  }
  mtext("Nuisance basis frontier (lower y is better)", outer = TRUE, cex = 1.15, font = 2)
}

fig1 <- file.path(out_dir, "nuisance_frontier_runtime_vs_disp50.png")
fig2 <- file.path(out_dir, "nuisance_frontier_runtime_vs_rot.png")

plot_metric("disp50_ratio_median", "Median disp50 ratio vs off", fig1)
plot_metric("rot_ratio_median", "Rotation RMSE ratio vs off", fig2)

cat("Saved figures:\n")
cat(" -", fig1, "\n")
cat(" -", fig2, "\n")
