test_that("profile silently used when control is at defaults", {
  expect_silent(.apply_volreg_profile(volreg_control(), "fast_fmri", control_user = FALSE))
  expect_silent(.apply_volreg_profile(volreg_control(), "fast_fmri", control_user = TRUE))
})

test_that("warning is emitted only when user-supplied control conflicts with profile fields", {
  ctl <- volreg_control()
  ctl$huber_k <- 0.05  # explicit user override that fast_fmri does NOT touch
  expect_silent(
    .apply_volreg_profile(ctl, "fast_fmri", control_user = TRUE)
  )

  ctl$max_iter <- c(60L, 40L, 30L)  # field that fast_fmri sets to a different value
  expect_warning(
    .apply_volreg_profile(ctl, "fast_fmri", control_user = TRUE),
    "max_iter"
  )
})

test_that("warning lists every conflicting field", {
  ctl <- volreg_control()
  ctl$max_iter <- c(60L, 40L, 30L)
  ctl$stencil_n <- 9999L
  w <- tryCatch(
    {
      .apply_volreg_profile(ctl, "fast_fmri", control_user = TRUE)
      NULL
    },
    warning = function(w) conditionMessage(w)
  )
  expect_match(w, "max_iter")
  expect_match(w, "stencil_n")
})

test_that("warning is suppressed by options(volregger.silent_profile = TRUE)", {
  ctl <- volreg_control()
  ctl$max_iter <- c(60L, 40L, 30L)
  withr::with_options(list(volregger.silent_profile = TRUE), {
    expect_silent(.apply_volreg_profile(ctl, "fast_fmri", control_user = TRUE))
  })
})

test_that("warning does not fire when control_user = FALSE (auto-applied profile)", {
  ctl <- volreg_control()
  ctl$max_iter <- c(60L, 40L, 30L)
  expect_silent(.apply_volreg_profile(ctl, "fast_fmri", control_user = FALSE))
})

test_that("volreg() passes control_user correctly: warning for explicit control + profile", {
  skip_if_not(exists("volreg"))
  x <- array(stats::rnorm(8 * 8 * 4 * 5), dim = c(8, 8, 4, 5))
  ctl <- volreg_control()
  ctl$max_iter <- c(60L, 40L, 30L)
  expect_warning(
    fit <- volreg(x, control = ctl, profile = "fast_native", qc = FALSE, resample = FALSE),
    "max_iter"
  )
  expect_true(inherits(fit, "volregger_fit"))
})

test_that("volreg() default call (no control, no profile) is silent", {
  x <- array(stats::rnorm(8 * 8 * 4 * 5), dim = c(8, 8, 4, 5))
  expect_silent(
    fit <- volreg(x, qc = FALSE, resample = FALSE)
  )
})
