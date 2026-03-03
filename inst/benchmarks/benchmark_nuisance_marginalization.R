#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("Need devtools to load local package from source", call. = FALSE)
    }
    devtools::load_all(".", quiet = TRUE)
  } else {
    if (!requireNamespace("volregger", quietly = TRUE)) {
      stop("Need installed volregger package when running outside source tree", call. = FALSE)
    }
    library(volregger)
  }
})

helper_path <- file.path("tests", "testthat", "helper-continuous-sim.R")
if (!file.exists(helper_path)) {
  stop("Missing helper file: ", helper_path, call. = FALSE)
}
source(helper_path, local = TRUE)

deg <- function(x) x * pi / 180

pose_metrics <- function(est, truth_pose) {
  err <- est$par + truth_pose
  disp50 <- sqrt(rowSums(err[, 1:3, drop = FALSE]^2) +
    (50 * sqrt(rowSums(err[, 4:6, drop = FALSE]^2)))^2)
  c(
    trans_rmse = sqrt(mean(err[, 1:3]^2)),
    rot_rmse = sqrt(mean(err[, 4:6]^2)),
    disp50_median = median(disp50),
    fail_rate = mean(est$flags != 0),
    mean_iters = mean(est$iters)
  )
}

run_one <- function(seed, scenario, use_nuisance) {
  subject <- make_continuous_subject(seed = 6200 + seed)
  nt <- 12
  tt <- seq(0, 1, length.out = nt)

  poses <- cbind(
    2.0 * sin(2 * pi * tt),
    1.1 * cos(2 * pi * tt),
    0.8 * sin(3 * pi * tt),
    deg(2.5) * sin(2 * pi * tt),
    deg(2.0) * cos(2 * pi * tt),
    deg(2.2) * sin(3 * pi * tt)
  )
  poses <- sweep(poses, 2, poses[1, ], FUN = "-")

  sim <- simulate_continuous_run(
    subject = subject,
    nx = 16,
    ny = 16,
    nz = 8,
    nt = nt,
    poses = poses,
    noise_sd = scenario$noise_sd,
    drift_amp = scenario$drift_amp,
    bias_amp = scenario$bias_amp,
    supersample = 1L,
    seed = 7200 + 100L * seed + as.integer(1000 * scenario$bias_amp)
  )

  ctl <- utils::modifyList(
    volreg_control(),
    list(
      max_iter = 8L,
      n_samples = 2200L,
      n_restarts = 2L,
      use_ic_stencil = TRUE,
      stencil_n = 1800L,
      whiten_modes = TRUE,
      capture_boost = TRUE,
      capture_trans_mm = 10,
      capture_rot_deg = 10,
      capture_topk = 8L,
      use_nuisance_marginalization = isTRUE(use_nuisance),
      nuisance_include_template = FALSE,
      nuisance_orthogonalize = TRUE,
      nuisance_ridge = 1e-4
    )
  )

  elapsed <- system.time({
    est <- estimate_motion(sim$run, ref = "index", ref_index = 1, control = ctl)
  })[["elapsed"]]
  met <- pose_metrics(est, poses)

  data.frame(
    scenario = scenario$name,
    seed = seed,
    nuisance = if (use_nuisance) "on" else "off",
    elapsed_sec = unname(elapsed),
    trans_rmse = met[["trans_rmse"]],
    rot_rmse = met[["rot_rmse"]],
    disp50_median = met[["disp50_median"]],
    fail_rate = met[["fail_rate"]],
    mean_iters = met[["mean_iters"]],
    stringsAsFactors = FALSE
  )
}

benchmark_nuisance <- function(seeds = 1:8) {
  scenarios <- list(
    list(name = "low_nuisance", noise_sd = 1.8, drift_amp = 0.02, bias_amp = 0.04),
    list(name = "moderate_nuisance", noise_sd = 2.0, drift_amp = 0.04, bias_amp = 0.08),
    list(name = "strong_nuisance", noise_sd = 2.2, drift_amp = 0.06, bias_amp = 0.12)
  )

  rows <- vector("list", length(scenarios) * length(seeds) * 2L)
  idx <- 1L
  for (sc in scenarios) {
    for (sd in seeds) {
      rows[[idx]] <- run_one(sd, sc, FALSE); idx <- idx + 1L
      rows[[idx]] <- run_one(sd, sc, TRUE); idx <- idx + 1L
    }
  }
  raw <- do.call(rbind, rows)

  # Paired ratios: on/off for each scenario and seed.
  keys <- unique(raw[, c("scenario", "seed")])
  paired <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    sc <- keys$scenario[i]
    sd <- keys$seed[i]
    a <- raw[raw$scenario == sc & raw$seed == sd & raw$nuisance == "off", ]
    b <- raw[raw$scenario == sc & raw$seed == sd & raw$nuisance == "on", ]
    paired[[i]] <- data.frame(
      scenario = sc,
      seed = sd,
      elapsed_ratio = b$elapsed_sec / a$elapsed_sec,
      trans_ratio = b$trans_rmse / a$trans_rmse,
      rot_ratio = b$rot_rmse / a$rot_rmse,
      disp50_ratio = b$disp50_median / a$disp50_median,
      fail_delta = b$fail_rate - a$fail_rate,
      iter_ratio = b$mean_iters / a$mean_iters,
      stringsAsFactors = FALSE
    )
  }
  paired <- do.call(rbind, paired)

  summary <- do.call(
    rbind,
    lapply(split(paired, paired$scenario), function(df) {
      data.frame(
        scenario = df$scenario[1],
        n_pairs = nrow(df),
        elapsed_ratio_median = median(df$elapsed_ratio),
        trans_ratio_median = median(df$trans_ratio),
        rot_ratio_median = median(df$rot_ratio),
        disp50_ratio_median = median(df$disp50_ratio),
        fail_delta_median = median(df$fail_delta),
        iter_ratio_median = median(df$iter_ratio),
        stringsAsFactors = FALSE
      )
    })
  )

  list(raw = raw, paired = paired, summary = summary)
}

res <- benchmark_nuisance(seeds = 1:8)

raw_csv <- file.path("inst", "benchmarks", "nuisance_marginalization_raw.csv")
paired_csv <- file.path("inst", "benchmarks", "nuisance_marginalization_paired.csv")
summary_csv <- file.path("inst", "benchmarks", "nuisance_marginalization_summary.csv")

write.csv(res$raw, raw_csv, row.names = FALSE)
write.csv(res$paired, paired_csv, row.names = FALSE)
write.csv(res$summary, summary_csv, row.names = FALSE)

cat("Saved:\n")
cat(" -", raw_csv, "\n")
cat(" -", paired_csv, "\n")
cat(" -", summary_csv, "\n\n")
print(res$summary)
