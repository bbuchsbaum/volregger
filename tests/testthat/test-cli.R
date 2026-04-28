.cli_fixture <- function(seed = 1L, dim = c(8L, 8L, 4L, 5L)) {
  set.seed(seed)
  array(stats::rnorm(prod(dim)), dim = dim)
}

.cli_capture <- function(args) {
  status <- NULL
  out <- capture.output(
    msgs <- capture.output(status <- volregger_cli(args), type = "message")
  )
  list(status = status, stdout = out, stderr = msgs)
}

test_that("no args prints usage and returns 0", {
  r <- .cli_capture(character(0))
  expect_equal(r$status, 0L)
  expect_true(any(grepl("Usage:", r$stdout, fixed = TRUE)))
})

test_that("--help, -h, --version return 0", {
  for (flag in c("--help", "-h")) {
    r <- .cli_capture(flag)
    expect_equal(r$status, 0L)
    expect_true(any(grepl("Usage:", r$stdout, fixed = TRUE)))
  }
  r <- .cli_capture("--version")
  expect_equal(r$status, 0L)
  expect_true(any(grepl("^volregger ", r$stdout)))
})

test_that("help <command> prints command-specific usage", {
  for (cmd in c("estimate", "apply", "run", "report")) {
    r <- .cli_capture(c("help", cmd))
    expect_equal(r$status, 0L)
    expect_true(any(grepl("Usage:", r$stdout, fixed = TRUE)))
  }
})

test_that("unknown command returns 2", {
  r <- .cli_capture("nosuchcmd")
  expect_equal(r$status, 2L)
})

test_that("estimate happy path works (RDS in/out, no QC)", {
  in_path <- tempfile(fileext = ".rds")
  out_fit <- tempfile(fileext = ".rds")
  saveRDS(.cli_fixture(), in_path)
  r <- .cli_capture(c("estimate", "--input", in_path, "--out-fit", out_fit,
                      "--no-qc", "--profile", "fast_native"))
  expect_equal(r$status, 0L)
  expect_true(file.exists(out_fit))
  fit <- readRDS(out_fit)
  expect_true(inherits(fit, "volregger_fit"))
})

test_that("estimate --json emits a JSON object on stdout", {
  in_path <- tempfile(fileext = ".rds")
  out_fit <- tempfile(fileext = ".rds")
  saveRDS(.cli_fixture(), in_path)
  r <- .cli_capture(c("estimate", "--input", in_path, "--out-fit", out_fit,
                      "--no-qc", "--profile", "fast_native", "--json"))
  expect_equal(r$status, 0L)
  joined <- paste(r$stdout, collapse = "")
  expect_match(joined, "^\\{.*\\}$")
  expect_match(joined, "\"fit\":")
})

test_that("estimate without --json emits key=value lines", {
  in_path <- tempfile(fileext = ".rds")
  out_fit <- tempfile(fileext = ".rds")
  saveRDS(.cli_fixture(), in_path)
  r <- .cli_capture(c("estimate", "--input", in_path, "--out-fit", out_fit,
                      "--no-qc", "--profile", "fast_native"))
  expect_equal(r$status, 0L)
  expect_true(any(grepl("^fit=", r$stdout)))
})

test_that("estimate --quiet suppresses stdout but still writes the fit", {
  in_path <- tempfile(fileext = ".rds")
  out_fit <- tempfile(fileext = ".rds")
  saveRDS(.cli_fixture(), in_path)
  r <- .cli_capture(c("estimate", "--input", in_path, "--out-fit", out_fit,
                      "--no-qc", "--profile", "fast_native", "--quiet"))
  expect_equal(r$status, 0L)
  joined <- paste(r$stdout, collapse = "")
  expect_false(grepl("fit=", joined, fixed = TRUE))
  expect_true(file.exists(out_fit))
})

test_that("apply --json emits JSON object", {
  in_path <- tempfile(fileext = ".rds")
  out_fit <- tempfile(fileext = ".rds")
  saveRDS(.cli_fixture(), in_path)
  expect_equal(
    .cli_capture(c("estimate", "--input", in_path, "--out-fit", out_fit,
                   "--no-qc", "--profile", "fast_native", "--quiet"))$status,
    0L
  )
  out_corr <- tempfile(fileext = ".rds")
  r <- .cli_capture(c("apply", "--input", in_path, "--fit", out_fit,
                      "--out", out_corr, "--json"))
  expect_equal(r$status, 0L)
  expect_match(paste(r$stdout, collapse = ""), "\"corrected\":")
  expect_true(file.exists(out_corr))
})

test_that("inline --key=value parsing works", {
  in_path <- tempfile(fileext = ".rds")
  out_fit <- tempfile(fileext = ".rds")
  saveRDS(.cli_fixture(), in_path)
  r <- .cli_capture(c("estimate",
                      paste0("--input=", in_path),
                      paste0("--out-fit=", out_fit),
                      "--no-qc", "--profile=fast_native"))
  expect_equal(r$status, 0L)
  expect_true(file.exists(out_fit))
})

test_that("missing required option returns 2", {
  r <- .cli_capture(c("estimate", "--input", "x.rds"))
  expect_equal(r$status, 2L)
  expect_true(any(grepl("--out-fit", r$stderr)))
})

test_that("missing input file returns 1 (domain error)", {
  r <- .cli_capture(c("estimate", "--input", "/no/such/file.rds",
                      "--out-fit", tempfile(fileext = ".rds")))
  expect_equal(r$status, 1L)
})

test_that("unknown option returns 2", {
  r <- .cli_capture(c("estimate", "--input", "x.rds",
                      "--out-fit", "y.rds", "--bogus", "z"))
  expect_equal(r$status, 2L)
})

test_that("apply with non-volregger_fit RDS returns 1", {
  in_path <- tempfile(fileext = ".rds"); saveRDS(.cli_fixture(), in_path)
  fit_path <- tempfile(fileext = ".rds"); saveRDS("not a fit", fit_path)
  r <- .cli_capture(c("apply", "--input", in_path, "--fit", fit_path,
                      "--out", tempfile(fileext = ".rds")))
  expect_equal(r$status, 1L)
})

test_that("--ref index with out-of-range --ref-index returns non-zero", {
  # ref-index parses cleanly, but the index is rejected at runtime by
  # .resolve_ref_index() — that's a generic error, not a CLI domain error,
  # so the wrapper returns 2 ('runtime error'). The contract here is that
  # the run does not silently succeed.
  in_path <- tempfile(fileext = ".rds"); saveRDS(.cli_fixture(), in_path)
  r <- .cli_capture(c("estimate", "--input", in_path,
                      "--out-fit", tempfile(fileext = ".rds"),
                      "--ref", "index", "--ref-index", "999",
                      "--no-qc", "--profile", "fast_native"))
  expect_gt(r$status, 0L)
})

test_that("run --no-corrected --no-report --no-plot still writes fit and motion", {
  in_path <- tempfile(fileext = ".rds"); saveRDS(.cli_fixture(), in_path)
  out_dir <- tempfile()
  r <- .cli_capture(c("run", "--input", in_path, "--out-dir", out_dir,
                      "--prefix", "smoke",
                      "--no-corrected", "--no-report", "--no-plot",
                      "--no-qc", "--profile", "fast_native"))
  expect_equal(r$status, 0L)
  expect_true(file.exists(file.path(out_dir, "smoke_fit.rds")))
  expect_true(file.exists(file.path(out_dir, "smoke_motion.tsv")))
})

test_that("report --json from a saved fit emits JSON paths", {
  in_path <- tempfile(fileext = ".rds"); saveRDS(.cli_fixture(), in_path)
  fit_path <- tempfile(fileext = ".rds")
  expect_equal(
    .cli_capture(c("estimate", "--input", in_path, "--out-fit", fit_path,
                   "--profile", "fast_native", "--quiet"))$status,
    0L
  )
  out_dir <- tempfile()
  r <- .cli_capture(c("report", "--fit", fit_path, "--out-dir", out_dir,
                      "--prefix", "rep", "--no-plot", "--json"))
  expect_equal(r$status, 0L)
  expect_match(paste(r$stdout, collapse = ""), "^\\{.*\\}$")
})
