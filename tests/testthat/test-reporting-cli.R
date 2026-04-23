test_that("volreg_report merges frame metrics and truth displacement summaries", {
  set.seed(301)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 5), dim = c(8, 8, 4, 5))
  mask <- array(TRUE, dim = c(8, 8, 4))
  truth <- array(diag(4), dim = c(4, 4, 5))

  fit <- volreg(
    x,
    resample = FALSE,
    qc = TRUE,
    control = utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 1800L))
  )
  report <- volreg_report(fit, truth = truth, mask = mask)

  expect_s3_class(report, "volregger_report")
  expect_true(all(c("disp_median", "disp_p95", "disp_max", "fd_error") %in% names(report$frames)))
  expect_true(all(c("engine", "components", "truth_disp_median_median") %in% names(report$summary)))
  expect_equal(nrow(report$matrices), 5)
  expect_equal(report$matrices$frame, seq_len(5))
  expect_match(paste(report$components, collapse = ","), "capture_boost")
})

test_that("write_volreg_report writes the expected artifact bundle", {
  set.seed(302)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  fit <- volreg(
    x,
    resample = FALSE,
    qc = TRUE,
    control = utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 1800L))
  )
  report <- volreg_report(fit)

  out_dir <- tempfile("volregger-report-")
  paths <- write_volreg_report(report, out_dir = out_dir, prefix = "smoke")

  expect_true(all(file.exists(unname(paths))))
  expect_true(all(c("summary", "frames", "matrices", "text", "rds", "overview") %in% names(paths)))
})

test_that("volregger_cli report command writes report outputs", {
  set.seed(303)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  fit <- volreg(
    x,
    resample = FALSE,
    qc = TRUE,
    control = utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 1800L))
  )
  truth <- array(diag(4), dim = c(4, 4, 4))
  mask <- array(TRUE, dim = c(8, 8, 4))

  fit_file <- tempfile("volregger-fit-", fileext = ".rds")
  truth_file <- tempfile("volregger-truth-", fileext = ".rds")
  mask_file <- tempfile("volregger-mask-", fileext = ".rds")
  out_dir <- tempfile("volregger-cli-")

  saveRDS(fit, fit_file)
  saveRDS(truth, truth_file)
  saveRDS(mask, mask_file)

  status <- suppressMessages(
    capture.output(
      volregger_cli(c(
        "report",
        "--fit", fit_file,
        "--truth", truth_file,
        "--mask", mask_file,
        "--out-dir", out_dir,
        "--prefix", "cli-smoke"
      ))
    )
  )

  expect_true(any(grepl("cli-smoke_summary.csv", status, fixed = TRUE)))
  expect_true(file.exists(file.path(out_dir, "cli-smoke_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "cli-smoke_frames.csv")))
  expect_true(file.exists(file.path(out_dir, "cli-smoke_matrices.csv")))
  expect_true(file.exists(file.path(out_dir, "cli-smoke_overview.png")))
})

test_that("volregger_cli returns stable help and usage statuses", {
  help_lines <- capture.output(status <- volregger_cli(character()))
  expect_equal(status, 0L)
  expect_equal(help_lines[[1]], "Usage:")
  expect_equal(suppressMessages(volregger_cli(c("report", "--help"))), 0L)
  expect_equal(suppressMessages(volregger_cli(c("report"))), 2L)
  expect_equal(suppressMessages(volregger_cli(c("report", "--unknown"))), 2L)
})

test_that("install_cli copies the public wrapper", {
  dest_dir <- tempfile("volregger-bin-")
  dir.create(dest_dir)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(normalizePath(file.path(testthat::test_path(), "..", "..")))

  paths <- install_cli(dest_dir = dest_dir, overwrite = TRUE)

  expect_length(paths, 1L)
  expect_true(file.exists(paths[[1]]))
  expect_equal(basename(paths[[1]]), "volregger")
})
