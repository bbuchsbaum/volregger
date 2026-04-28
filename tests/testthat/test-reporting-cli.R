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

test_that("volregger_cli --json emits a parseable single-line JSON object", {
  set.seed(304)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  fit <- volreg(
    x,
    resample = FALSE,
    qc = TRUE,
    control = utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 1800L))
  )

  fit_file <- tempfile("volregger-fit-", fileext = ".rds")
  out_dir <- tempfile("volregger-cli-json-")
  saveRDS(fit, fit_file)

  out <- capture.output(
    status <- volregger_cli(c(
      "report",
      "--fit", fit_file,
      "--out-dir", out_dir,
      "--prefix", "j",
      "--no-plot",
      "--json"
    ))
  )
  expect_equal(status, 0L)

  json_lines <- out[nzchar(out)]
  expect_length(json_lines, 1L)
  expect_true(startsWith(json_lines[[1L]], "{"))
  expect_true(endsWith(json_lines[[1L]], "}"))
  expect_true(grepl('"summary":"', json_lines, fixed = TRUE))
  expect_true(grepl('"frames":"', json_lines, fixed = TRUE))
  expect_true(grepl('"rds":"', json_lines, fixed = TRUE))
  expect_false(grepl("overview", json_lines, fixed = TRUE))
})

test_that("volregger_cli --quiet suppresses stdout but still writes files", {
  set.seed(305)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  fit <- volreg(
    x,
    resample = FALSE,
    qc = TRUE,
    control = utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 1800L))
  )

  fit_file <- tempfile("volregger-fit-", fileext = ".rds")
  out_dir <- tempfile("volregger-cli-quiet-")
  saveRDS(fit, fit_file)

  out <- capture.output(
    status <- volregger_cli(c(
      "report",
      "--fit", fit_file,
      "--out-dir", out_dir,
      "--prefix", "q",
      "--no-plot",
      "--quiet"
    ))
  )
  expect_equal(status, 0L)
  expect_equal(length(out[nzchar(out)]), 0L)
  expect_true(file.exists(file.path(out_dir, "q_summary.csv")))
})

test_that("volregger_cli returns stable help, version, and usage statuses", {
  help_lines <- capture.output(status <- volregger_cli(character()))
  expect_equal(status, 0L)
  expect_equal(help_lines[[1L]], "Usage:")

  expect_equal(suppressMessages(volregger_cli(c("report", "--help"))), 0L)
  expect_equal(suppressMessages(volregger_cli(c("report"))), 2L)
  expect_equal(suppressMessages(volregger_cli(c("report", "--unknown"))), 2L)

  expect_equal(suppressMessages(volregger_cli(c("--help"))), 0L)
  expect_equal(suppressMessages(volregger_cli(c("-h"))), 0L)

  ver_out <- capture.output(ver_status <- volregger_cli("--version"))
  expect_equal(ver_status, 0L)
  expect_match(ver_out[[1L]], "^volregger ")

  help_cmd_out <- capture.output(help_cmd_status <- volregger_cli(c("help", "report")))
  expect_equal(help_cmd_status, 0L)
  expect_true(any(grepl("Build a structured report bundle", help_cmd_out, fixed = TRUE)))
  expect_true(any(grepl("Examples:", help_cmd_out, fixed = TRUE)))

  expect_equal(suppressMessages(volregger_cli(c("help", "garbage"))), 2L)
})

test_that("volregger_cli parser handles --key=value, missing values, and unknown flags", {
  fit_file <- tempfile("volregger-fit-", fileext = ".rds")
  saveRDS(list(), fit_file)

  expect_equal(suppressMessages(volregger_cli(c("report", "--fit"))), 2L)
  expect_equal(suppressMessages(volregger_cli(c("report", "--out-dir"))), 2L)
  expect_equal(suppressMessages(volregger_cli(c("report", "positional"))), 2L)
  expect_equal(suppressMessages(volregger_cli(c("report", "--plot=yes", "--fit", fit_file))), 2L)

  expect_equal(suppressMessages(volregger_cli(c("report", "--fit", fit_file))), 1L)
})

test_that("install_cli copies the public wrapper", {
  dest_dir <- tempfile("volregger-bin-")
  dir.create(dest_dir)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(normalizePath(file.path(testthat::test_path(), "..", "..")))

  paths <- suppressMessages(install_cli(dest_dir = dest_dir, overwrite = TRUE, quiet = TRUE))

  expect_length(paths, 1L)
  expect_true(file.exists(paths[[1L]]))
  expect_equal(basename(paths[[1L]]), "volregger")
})

test_that("install_cli refuses overwrite and warns about PATH", {
  dest_dir <- tempfile("volregger-bin-")
  dir.create(dest_dir)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(normalizePath(file.path(testthat::test_path(), "..", "..")))

  suppressMessages(install_cli(dest_dir = dest_dir, overwrite = TRUE, quiet = TRUE))
  expect_error(
    install_cli(dest_dir = dest_dir, overwrite = FALSE, quiet = TRUE),
    "refusing to overwrite"
  )

  msgs <- capture.output(
    install_cli(dest_dir = dest_dir, overwrite = TRUE, quiet = FALSE),
    type = "message"
  )
  expect_true(any(grepl("Installed:", msgs, fixed = TRUE)))
  if (!any(strsplit(Sys.getenv("PATH", unset = ""), .Platform$path.sep, fixed = TRUE)[[1L]] ==
           normalizePath(dest_dir, mustWork = FALSE))) {
    expect_true(any(grepl("not on $PATH", msgs, fixed = TRUE)))
  }
})

test_that("volregger_cli estimate writes fit + motion TSV from RDS input", {
  set.seed(401)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  attr(x, "pixdim") <- c(2, 2, 3)
  attr(x, "tr") <- 1
  input <- tempfile("volregger-input-", fileext = ".rds")
  saveRDS(x, input)

  out_fit <- tempfile("volregger-fit-", fileext = ".rds")
  out_motion <- tempfile("volregger-motion-", fileext = ".tsv")

  status <- suppressMessages(capture.output(
    s <- volregger_cli(c(
      "estimate",
      "--input", input,
      "--out-fit", out_fit,
      "--out-motion", out_motion,
      "--profile", "fast_fmri",
      "--quiet"
    ))
  ))
  expect_equal(s, 0L)
  expect_true(file.exists(out_fit))
  fit <- readRDS(out_fit)
  expect_s3_class(fit, "volregger_fit")
  expect_true(is.null(fit$corrected))

  expect_true(file.exists(out_motion))
  motion <- utils::read.table(out_motion, header = TRUE, sep = "\t")
  expect_equal(colnames(motion), c("tx", "ty", "tz", "rx", "ry", "rz"))
  expect_equal(nrow(motion), 4L)
})

test_that("volregger_cli apply round-trips a corrected RDS array", {
  set.seed(402)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  attr(x, "pixdim") <- c(2, 2, 3)
  input <- tempfile("volregger-input-", fileext = ".rds")
  saveRDS(x, input)

  fit <- volreg(
    x, resample = FALSE, qc = TRUE,
    control = utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 1800L))
  )
  fit_path <- tempfile("volregger-fit-", fileext = ".rds")
  saveRDS(fit, fit_path)

  out <- tempfile("volregger-corrected-", fileext = ".rds")
  status <- suppressMessages(capture.output(
    s <- volregger_cli(c(
      "apply",
      "--input", input,
      "--fit", fit_path,
      "--out", out,
      "--quiet"
    ))
  ))
  expect_equal(s, 0L)
  expect_true(file.exists(out))
  corrected <- readRDS(out)
  expect_equal(dim(corrected), dim(x))
})

test_that("volregger_cli run produces fit, motion, and report bundle", {
  set.seed(403)
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  x <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  attr(x, "pixdim") <- c(2, 2, 3)
  input <- tempfile("volregger-run-", fileext = ".rds")
  saveRDS(x, input)

  out_dir <- tempfile("volregger-run-out-")

  out <- capture.output(
    s <- volregger_cli(c(
      "run",
      "--input", input,
      "--out-dir", out_dir,
      "--prefix", "r1",
      "--profile", "fast_fmri",
      "--no-plot",
      "--json"
    ))
  )
  expect_equal(s, 0L)

  json_lines <- out[nzchar(out)]
  expect_length(json_lines, 1L)
  expect_true(startsWith(json_lines[[1L]], "{"))

  expect_true(file.exists(file.path(out_dir, "r1_fit.rds")))
  expect_true(file.exists(file.path(out_dir, "r1_motion.tsv")))
  expect_true(file.exists(file.path(out_dir, "r1_corrected.rds")))
  expect_true(file.exists(file.path(out_dir, "r1_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "r1_frames.csv")))
  expect_true(file.exists(file.path(out_dir, "r1_matrices.csv")))
  expect_false(file.exists(file.path(out_dir, "r1_overview.png")))
})

test_that("volregger_cli run handles NIfTI round-trip via neuroim2", {
  skip_on_cran()
  skip_if_not_installed("neuroim2")
  suppressPackageStartupMessages(library(neuroim2))

  set.seed(404)
  spc <- neuroim2::NeuroSpace(c(8L, 8L, 4L, 4L), spacing = c(2, 2, 3))
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  data4d <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  v <- neuroim2::NeuroVec(data4d, spc)

  input <- tempfile("volregger-nii-", fileext = ".nii.gz")
  neuroim2::write_vec(v, input)

  out_dir <- tempfile("volregger-nii-out-")

  status <- suppressMessages(capture.output(
    s <- volregger_cli(c(
      "run",
      "--input", input,
      "--out-dir", out_dir,
      "--prefix", "n1",
      "--profile", "fast_fmri",
      "--no-plot",
      "--quiet"
    ))
  ))
  expect_equal(s, 0L)

  corrected_nii <- file.path(out_dir, "n1_corrected.nii.gz")
  expect_true(file.exists(corrected_nii))
  expect_true(file.exists(file.path(out_dir, "n1_fit.rds")))
  expect_true(file.exists(file.path(out_dir, "n1_motion.tsv")))

  v2 <- neuroim2::read_vec(corrected_nii)
  expect_equal(as.integer(dim(v2)), as.integer(dim(v)))
  expect_equal(neuroim2::spacing(v2), neuroim2::spacing(v))
})

test_that("volregger_cli apply writes a NIfTI when input is NIfTI", {
  skip_on_cran()
  skip_if_not_installed("neuroim2")
  suppressPackageStartupMessages(library(neuroim2))

  set.seed(405)
  spc <- neuroim2::NeuroSpace(c(8L, 8L, 4L, 4L), spacing = c(2, 2, 3))
  x0 <- array(rnorm(8 * 8 * 4), dim = c(8, 8, 4))
  data4d <- array(rep(as.numeric(x0), 4), dim = c(8, 8, 4, 4))
  v <- neuroim2::NeuroVec(data4d, spc)

  input <- tempfile("volregger-nii-apply-", fileext = ".nii.gz")
  neuroim2::write_vec(v, input)

  fit <- volreg(
    data4d, resample = FALSE, qc = TRUE,
    control = utils::modifyList(volreg_control(), list(max_iter = 4L, n_samples = 1800L))
  )
  fit_path <- tempfile("volregger-fit-", fileext = ".rds")
  saveRDS(fit, fit_path)

  out_path <- tempfile("volregger-corrected-", fileext = ".nii.gz")
  status <- suppressMessages(capture.output(
    s <- volregger_cli(c(
      "apply",
      "--input", input,
      "--fit", fit_path,
      "--out", out_path,
      "--quiet"
    ))
  ))
  expect_equal(s, 0L)
  expect_true(file.exists(out_path))
  v2 <- neuroim2::read_vec(out_path)
  expect_equal(as.integer(dim(v2)), as.integer(dim(v)))
})

test_that("volregger_cli new commands expose help and required-arg errors", {
  for (cmd in c("estimate", "apply", "run")) {
    expect_equal(suppressMessages(volregger_cli(c(cmd, "--help"))), 0L)
    expect_equal(suppressMessages(volregger_cli(c(cmd))), 2L)
  }
  help_out <- capture.output(s <- volregger_cli(c("help", "run")))
  expect_equal(s, 0L)
  expect_true(any(grepl("End-to-end pipeline", help_out, fixed = TRUE)))
})

test_that("volregger_cli estimate/apply require neuroim2 for NIfTI inputs", {
  if (requireNamespace("neuroim2", quietly = TRUE)) {
    succeed("neuroim2 is installed; skipping unavailability test")
    return()
  }
  fake <- tempfile("volregger-fake-", fileext = ".nii.gz")
  writeLines("not a real nifti", fake)
  expect_equal(
    suppressMessages(volregger_cli(c(
      "estimate", "--input", fake, "--out-fit", tempfile(fileext = ".rds")
    ))),
    1L
  )
})

test_that("exec/volregger wrapper script runs end-to-end via Rscript", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.which("Rscript")), "Rscript not available")

  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
  wrapper <- file.path(repo_root, "exec", "volregger")
  skip_if_not(file.exists(wrapper), "wrapper script not present")

  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    args = c(shQuote(wrapper), "--version"),
    stdout = TRUE,
    stderr = TRUE
  ))
  expect_true(any(grepl("^volregger ", out)))

  help_out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    args = c(shQuote(wrapper), "help"),
    stdout = TRUE,
    stderr = TRUE
  ))
  expect_true(any(grepl("Usage:", help_out, fixed = TRUE)))
  expect_true(any(grepl("report", help_out, fixed = TRUE)))
})
