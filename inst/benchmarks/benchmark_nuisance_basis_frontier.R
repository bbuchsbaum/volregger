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

run_one <- function(seed, scenario, method) {
  subject <- make_continuous_subject(seed = 7000 + seed)
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
    seed = 8200 + 100L * seed + as.integer(1000 * scenario$bias_amp)
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
      use_nuisance_marginalization = isTRUE(method$use_nuisance),
      nuisance_include_template = isTRUE(method$include_template),
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
    method = method$name,
    nuisance = method$use_nuisance,
    include_template = method$include_template,
    elapsed_sec = unname(elapsed),
    trans_rmse = met[["trans_rmse"]],
    rot_rmse = met[["rot_rmse"]],
    disp50_median = met[["disp50_median"]],
    fail_rate = met[["fail_rate"]],
    mean_iters = met[["mean_iters"]],
    stringsAsFactors = FALSE
  )
}

benchmark_frontier <- function(seeds = 1:8) {
  scenarios <- list(
    list(name = "low_nuisance", noise_sd = 1.8, drift_amp = 0.02, bias_amp = 0.04),
    list(name = "moderate_nuisance", noise_sd = 2.0, drift_amp = 0.04, bias_amp = 0.08),
    list(name = "strong_nuisance", noise_sd = 2.2, drift_amp = 0.06, bias_amp = 0.12)
  )

  methods <- list(
    list(name = "off", use_nuisance = FALSE, include_template = FALSE),
    list(name = "small", use_nuisance = TRUE, include_template = FALSE),
    list(name = "full", use_nuisance = TRUE, include_template = TRUE)
  )

  rows <- vector("list", length(scenarios) * length(seeds) * length(methods))
  idx <- 1L
  for (sc in scenarios) {
    for (sd in seeds) {
      for (md in methods) {
        rows[[idx]] <- run_one(sd, sc, md)
        idx <- idx + 1L
      }
    }
  }
  raw <- do.call(rbind, rows)

  # Paired deltas/ratios against "off" baseline for same scenario+seed.
  keys <- unique(raw[, c("scenario", "seed")])
  paired <- list()
  pidx <- 1L
  for (i in seq_len(nrow(keys))) {
    sc <- keys$scenario[i]
    sd <- keys$seed[i]
    base <- raw[raw$scenario == sc & raw$seed == sd & raw$method == "off", ]
    for (m in c("small", "full")) {
      alt <- raw[raw$scenario == sc & raw$seed == sd & raw$method == m, ]
      paired[[pidx]] <- data.frame(
        scenario = sc,
        seed = sd,
        method = m,
        elapsed_ratio = alt$elapsed_sec / base$elapsed_sec,
        trans_ratio = alt$trans_rmse / base$trans_rmse,
        rot_ratio = alt$rot_rmse / base$rot_rmse,
        disp50_ratio = alt$disp50_median / base$disp50_median,
        fail_delta = alt$fail_rate - base$fail_rate,
        iter_ratio = alt$mean_iters / base$mean_iters,
        stringsAsFactors = FALSE
      )
      pidx <- pidx + 1L
    }
  }
  paired <- do.call(rbind, paired)

  # Summary per scenario x method.
  summary_rows <- list()
  sidx <- 1L
  for (sc in unique(raw$scenario)) {
    for (m in c("off", "small", "full")) {
      r <- raw[raw$scenario == sc & raw$method == m, ]
      if (m == "off") {
        summary_rows[[sidx]] <- data.frame(
          scenario = sc,
          method = m,
          n_pairs = nrow(r),
          elapsed_sec_median = median(r$elapsed_sec),
          trans_rmse_median = median(r$trans_rmse),
          rot_rmse_median = median(r$rot_rmse),
          disp50_median = median(r$disp50_median),
          fail_rate_median = median(r$fail_rate),
          mean_iters_median = median(r$mean_iters),
          elapsed_ratio_median = 1.0,
          trans_ratio_median = 1.0,
          rot_ratio_median = 1.0,
          disp50_ratio_median = 1.0,
          fail_delta_median = 0.0,
          iter_ratio_median = 1.0,
          stringsAsFactors = FALSE
        )
      } else {
        p <- paired[paired$scenario == sc & paired$method == m, ]
        summary_rows[[sidx]] <- data.frame(
          scenario = sc,
          method = m,
          n_pairs = nrow(p),
          elapsed_sec_median = median(r$elapsed_sec),
          trans_rmse_median = median(r$trans_rmse),
          rot_rmse_median = median(r$rot_rmse),
          disp50_median = median(r$disp50_median),
          fail_rate_median = median(r$fail_rate),
          mean_iters_median = median(r$mean_iters),
          elapsed_ratio_median = median(p$elapsed_ratio),
          trans_ratio_median = median(p$trans_ratio),
          rot_ratio_median = median(p$rot_ratio),
          disp50_ratio_median = median(p$disp50_ratio),
          fail_delta_median = median(p$fail_delta),
          iter_ratio_median = median(p$iter_ratio),
          stringsAsFactors = FALSE
        )
      }
      sidx <- sidx + 1L
    }
  }
  summary <- do.call(rbind, summary_rows)

  list(raw = raw, paired = paired, summary = summary)
}

res <- benchmark_frontier(seeds = 1:8)

raw_csv <- file.path("inst", "benchmarks", "nuisance_basis_frontier_raw.csv")
paired_csv <- file.path("inst", "benchmarks", "nuisance_basis_frontier_paired.csv")
summary_csv <- file.path("inst", "benchmarks", "nuisance_basis_frontier_summary.csv")

write.csv(res$raw, raw_csv, row.names = FALSE)
write.csv(res$paired, paired_csv, row.names = FALSE)
write.csv(res$summary, summary_csv, row.names = FALSE)

cat("Saved:\n")
cat(" -", raw_csv, "\n")
cat(" -", paired_csv, "\n")
cat(" -", summary_csv, "\n\n")
print(res$summary)
