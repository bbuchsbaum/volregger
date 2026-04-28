test_that("volreg_control() default round-trips through validation", {
  expect_silent(volreg_validate_control(volreg_control()))
  expect_invisible(volreg_validate_control(volreg_control()))
})

test_that("volreg_control() respects per-field overrides", {
  ctl <- volreg_control(huber_k = 0.05, max_iter = c(60L, 40L, 30L), nthreads = 4L)
  expect_equal(ctl$huber_k, 0.05)
  expect_equal(ctl$max_iter, c(60L, 40L, 30L))
  expect_equal(ctl$nthreads, 4L)
  expect_silent(volreg_validate_control(ctl))
})

test_that("every named profile produces a control list that validates", {
  profile_names <- eval(formals(volreg_profile)$name)
  for (nm in profile_names) {
    prof <- volreg_profile(nm)
    expect_silent(volreg_validate_control(prof$control))
  }
})

test_that("validator rejects bad scalar types and ranges", {
  bad_cases <- list(
    list(huber_k = -1),
    list(huber_k = c(1, 2)),
    list(lambda0 = -0.1),
    list(tol_step = "x"),
    list(n_restarts = 1.5),
    list(capture_trans_mm = 0),
    list(capture_topk = 0L),
    list(stencil_n = -1L),
    list(stencil_gamma = 0),
    list(temporal_regularization_alpha = 1.5),
    list(low_motion_pose_scale = 0),
    list(low_motion_pose_scale = 1.5),
    list(edge_exclude_frac = 1),
    list(edge_exclude_frac = -0.1),
    list(nthreads = 0L),
    list(nthreads = -2L),
    list(tr = 0),
    list(engine = "rigid_other"),
    list(se3_update_order = "middle"),
    list(dense_interp = "foo"),
    list(adaptive_huber = "yes"),
    list(use_pyramid = c(TRUE, FALSE))
  )
  for (override in bad_cases) {
    ctl <- utils::modifyList(volreg_control(), override)
    expect_error(volreg_validate_control(ctl, strict = TRUE),
                 "invalid volreg control")
  }
})

test_that("validator catches use_pyramid with mismatched max_iter length", {
  ctl <- volreg_control(use_pyramid = TRUE,
                        downsample = c(4L, 2L, 1L),
                        max_iter = c(30L, 20L))
  expect_error(volreg_validate_control(ctl), "max_iter length")
})

test_that("validator tolerates extra fields set by profiles", {
  ctl <- volreg_control()
  ctl$interp_final <- "heptic"
  ctl$profile <- "fast_fmri"
  ctl$ref_mode <- "middle"
  expect_silent(volreg_validate_control(ctl))
})

test_that("strict = FALSE returns issue strings instead of stopping", {
  ctl <- utils::modifyList(volreg_control(), list(huber_k = -1, nthreads = 0L))
  issues <- volreg_validate_control(ctl, strict = FALSE)
  expect_type(issues, "character")
  expect_length(issues, 2L)
  expect_true(any(grepl("huber_k", issues)))
  expect_true(any(grepl("nthreads", issues)))
})

test_that("non-list control is rejected", {
  expect_error(volreg_validate_control(NULL), "named list")
  expect_error(volreg_validate_control(list(1, 2, 3)), "named list")
})

test_that("volreg() surfaces validator errors before reaching C++", {
  x <- array(stats::rnorm(8 * 8 * 4 * 5), dim = c(8, 8, 4, 5))
  ctl <- volreg_control(huber_k = -1)
  expect_error(
    volreg(x, control = ctl, profile = NULL, qc = FALSE, resample = FALSE),
    "huber_k"
  )
})
