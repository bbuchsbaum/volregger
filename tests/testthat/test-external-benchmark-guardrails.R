docker_benchmark_ready <- function() {
  docker <- Sys.which("docker")
  if (!nzchar(docker)) {
    return(FALSE)
  }
  status <- suppressWarnings(system2(docker, "info", stdout = FALSE, stderr = FALSE))
  identical(as.integer(status), 0L)
}

run_realdata_benchmark <- function(args, out_dir) {
  repo_root <- if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::pkg_path()
  } else {
    testthat::test_path("..", "..")
  }
  input_path <- file.path(repo_root, "data-raw", "sub-01_task-freerecall_bold.nii.gz")
  old <- setwd(repo_root)
  on.exit(setwd(old), add = TRUE)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("inst/benchmarks/benchmark_vs_nifreeze_realdata.R", "--input", input_path, args, "--out-dir", out_dir),
    stdout = FALSE,
    stderr = FALSE
  )
  expect_identical(as.integer(status), 0L)
  utils::read.csv(file.path(out_dir, "benchmark_vs_nifreeze_realdata_raw.csv"), check.names = FALSE)
}

test_that("external native benchmark keeps volregger closer to AFNI than nifreeze", {
  skip_on_cran()
  skip_if_not(identical(tolower(Sys.getenv("VOLREGGER_RUN_EXTERNAL_BENCH")), "true"))
  skip_if_not(docker_benchmark_ready())
  skip_if_not(requireNamespace("pkgload", quietly = TRUE))
  skip_if_not(file.exists(file.path(pkgload::pkg_path(), "data-raw", "sub-01_task-freerecall_bold.nii.gz")))

  out_dir <- file.path(tempdir(), "volregger_native_guardrail")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  raw <- run_realdata_benchmark(
    c("--native-only", "--native-frames", "8", "--nthreads", "2"),
    out_dir = out_dir
  )
  native <- raw[raw$scenario == "native" & raw$status == "ok", ]
  if (!all(c("volregger", "nifreeze", "afni_3dvolreg") %in% native$method)) {
    skip("native benchmark backends unavailable in this environment")
  }

  soln_path <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_solution_summary.csv")
  skip_if_not(file.exists(soln_path))
  skip_if_not(file.info(soln_path)$size > 16)
  soln <- utils::read.csv(
    soln_path,
    check.names = FALSE
  )
  vol <- native[native$method == "volregger", ]
  nif <- native[native$method == "nifreeze", ]
  afni <- native[native$method == "afni_3dvolreg", ]

  expect_true(all(c("estimate_sec", "apply_sec", "report_sec") %in% names(raw)))
  expect_equal(vol$report_sec, 0, tolerance = 1e-8)
  expect_lt(vol$elapsed_sec, nif$elapsed_sec)
  expect_gt(vol$tsnr_ratio, nif$tsnr_ratio)

  cmp_va <- soln[soln$scenario == "native" & soln$method_a == "volregger" & soln$method_b == "afni_3dvolreg", ]
  cmp_vn <- soln[soln$scenario == "native" & soln$method_a == "volregger" & soln$method_b == "nifreeze", ]

  expect_lt(cmp_va$rel_disp_p95, cmp_vn$rel_disp_p95)
  expect_gt(cmp_va$corrected_corr, cmp_vn$corrected_corr)
})

test_that("external truth benchmark keeps volregger ahead of nifreeze on low-motion truth", {
  skip_on_cran()
  skip_if_not(identical(tolower(Sys.getenv("VOLREGGER_RUN_EXTERNAL_BENCH")), "true"))
  skip_if_not(docker_benchmark_ready())
  skip_if_not(requireNamespace("pkgload", quietly = TRUE))
  skip_if_not(file.exists(file.path(pkgload::pkg_path(), "data-raw", "sub-01_task-freerecall_bold.nii.gz")))

  out_dir <- file.path(tempdir(), "volregger_truth_guardrail")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  raw <- run_realdata_benchmark(
    c("--native-frames", "4", "--truth-frames", "8", "--scenarios", "low_motion", "--skip-afni", "--nthreads", "2"),
    out_dir = out_dir
  )

  truth_rows <- raw[raw$scenario == "low_motion" & raw$status == "ok", ]
  vol <- truth_rows[truth_rows$method == "volregger", ]
  nif <- truth_rows[truth_rows$method == "nifreeze", ]
  if (nrow(vol) != 1 || nrow(nif) != 1) {
    skip("truth benchmark backends unavailable in this environment")
  }

  expect_match(vol$method_model, "^fast_fmri:")
  expect_lt(vol$disp_p95, nif$disp_p95)
  expect_lt(vol$fd_error, nif$fd_error)

  scorecard_path <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_scorecard.csv")
  skip_if_not(file.exists(scorecard_path))
  scorecard <- utils::read.csv(scorecard_path, check.names = FALSE)
  high_value <- scorecard[
    scorecard$scenario == "low_motion" &
      scorecard$baseline == "nifreeze" &
      scorecard$claim_metric &
      scorecard$family %in% c("truth_displacement", "truth_matrix"),
  ]
  expect_gt(nrow(high_value), 0)
  expect_false(any(high_value$comparison == "loss", na.rm = TRUE))
})

test_that("external truth benchmark claim metrics hold across varied scenarios", {
  skip_on_cran()
  skip_if_not(identical(tolower(Sys.getenv("VOLREGGER_RUN_EXTERNAL_BENCH")), "true"))
  skip_if_not(docker_benchmark_ready())
  skip_if_not(requireNamespace("pkgload", quietly = TRUE))
  skip_if_not(file.exists(file.path(pkgload::pkg_path(), "data-raw", "sub-01_task-freerecall_bold.nii.gz")))

  out_dir <- file.path(tempdir(), "volregger_varied_truth_guardrail")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- c("low_motion", "moderate_motion", "hard_motion_plus_nuisance")

  raw <- run_realdata_benchmark(
    c(
      "--native-frames", "4",
      "--truth-frames", "6",
      "--scenarios", paste(scenarios, collapse = ","),
      "--nthreads", "2"
    ),
    out_dir = out_dir
  )

  truth_rows <- raw[raw$scenario %in% scenarios & raw$status == "ok", ]
  if (!all(scenarios %in% truth_rows$scenario) ||
      !all(c("volregger", "nifreeze") %in% truth_rows$method)) {
    skip("varied truth benchmark backends unavailable in this environment")
  }

  claim_path <- file.path(out_dir, "benchmark_vs_nifreeze_realdata_claim_summary.csv")
  skip_if_not(file.exists(claim_path))
  claim_summary <- utils::read.csv(claim_path, check.names = FALSE)
  rows <- claim_summary[
    claim_summary$scenario %in% scenarios &
      claim_summary$baseline == "nifreeze" &
      claim_summary$claim_metric,
  ]

  expect_equal(sort(unique(rows$scenario)), sort(scenarios))
  expect_false(any(rows$losses > 0, na.rm = TRUE))
  expect_true(all(rows$wins > 0 | rows$ties > 0))

  afni_rows <- claim_summary[
    claim_summary$scenario %in% scenarios &
      claim_summary$baseline == "afni_3dvolreg" &
      claim_summary$claim_metric &
      claim_summary$family %in% c("truth_displacement", "truth_pose", "truth_matrix"),
  ]
  if (!nrow(afni_rows)) {
    skip("AFNI benchmark backend unavailable in this environment")
  }

  expect_equal(sort(unique(afni_rows$scenario)), sort(scenarios))
  expect_false(any(afni_rows$losses > 0, na.rm = TRUE))
  expect_true(all(afni_rows$wins > 0 | afni_rows$ties > 0))
})
