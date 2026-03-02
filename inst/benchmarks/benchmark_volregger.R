#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(volregger))

set.seed(20260302)
fixture_dim <- c(24, 24, 10, 20)
x <- array(rnorm(prod(fixture_dim)), dim = fixture_dim)
ctl <- utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 3000L, n_restarts = 1L))

bench <- system.time({
  fit <- volreg(x, resample = TRUE, qc = TRUE, control = ctl)
})

cat("volregger benchmark\n")
cat("fixture_dim:", paste(fixture_dim, collapse = "x"), "\n")
cat("elapsed_sec:", unname(bench[["elapsed"]]), "\n")
cat("user_sec:", unname(bench[["user.self"]]), "\n")
cat("sys_sec:", unname(bench[["sys.self"]]), "\n")
cat("motion_rows:", nrow(fit$motion), "\n")
cat("motion_cols:", ncol(fit$motion), "\n")
cat("qc_fields:", paste(sort(names(fit$qc)), collapse = ","), "\n")
cat("fit_bytes:", as.numeric(object.size(fit)), "\n")
cat("session_info:\n")
print(utils::sessionInfo())
