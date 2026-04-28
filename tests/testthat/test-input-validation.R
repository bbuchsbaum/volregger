test_that("volreg() rejects non-4D input through .as_volreg_input", {
  expect_error(volreg(array(0, dim = c(8, 8, 4))), regexp = ".*")
  expect_error(volreg(matrix(0, 4, 4)), regexp = ".*")
})

test_that("volreg() rejects mask whose spatial dims do not match x", {
  x <- volregger_example()
  bad_mask <- array(TRUE, dim = c(4, 4, 4))
  expect_error(volreg(x, mask = bad_mask, qc = FALSE, resample = FALSE),
               regexp = ".*")
})

test_that("volreg(ref = 'index') without ref_index errors before C++", {
  x <- volregger_example()
  expect_error(
    volreg(x, ref = "index", qc = FALSE, resample = FALSE),
    "ref_index"
  )
})

test_that("volreg(ref_index out of range) errors", {
  x <- volregger_example()
  expect_error(
    volreg(x, ref = "index", ref_index = 999L,
           qc = FALSE, resample = FALSE),
    "ref_index"
  )
})

test_that("zpad_final must be non-negative integer", {
  x <- volregger_example()
  expect_error(volreg(x, zpad_final = -1L, qc = FALSE, resample = FALSE),
               "zpad_final")
  expect_error(volreg(x, zpad_final = NA_integer_, qc = FALSE, resample = FALSE),
               "zpad_final")
})

test_that("apply_motion() requires nt x 6 motion matrix", {
  x <- volregger_example()
  expect_error(apply_motion(x, motion = matrix(0, nrow = dim(x)[4], ncol = 5)),
               "nt x 6 matrix")
  expect_error(apply_motion(x, motion = "not a matrix"), "nt x 6 matrix")
})

test_that("apply_motion() rejects negative zpad", {
  x <- volregger_example()
  motion <- matrix(0, nrow = dim(x)[4], ncol = 6)
  expect_error(apply_motion(x, motion, zpad = -1L), "zpad")
})

test_that("framewise_disp() requires nt x 6 motion matrix", {
  expect_error(framewise_disp(matrix(0, 5, 5)), "nt x 6 matrix")
  expect_error(framewise_disp(c(0, 0, 0, 0, 0, 0)), "nt x 6 matrix")
})

test_that("dvars() rejects non-4D input", {
  expect_error(dvars(array(0, dim = c(8, 8, 4))), "must be 4D")
  expect_error(dvars(matrix(0, 4, 4)), "must be 4D")
})

test_that("dvars() rejects mask with wrong spatial size", {
  x <- volregger_example()
  expect_error(dvars(x, mask = array(TRUE, dim = c(4, 4, 4))),
               "spatial dimensions")
})

test_that("repair_spin_history() rejects non-4D input", {
  expect_error(repair_spin_history(array(0, dim = c(8, 8, 4)),
                                   matrix(0, 1, 6)),
               "must be 4D")
})

test_that("repair_spin_history() rejects motion length mismatch", {
  x <- volregger_example()
  bad_motion <- matrix(0, nrow = dim(x)[4] + 1L, ncol = 6)
  expect_error(repair_spin_history(x, bad_motion), "ntime")
})

test_that("repair_spin_history() rejects mask shape mismatch", {
  x <- volregger_example()
  motion <- matrix(0, nrow = dim(x)[4], ncol = 6)
  expect_error(
    repair_spin_history(x, motion, mask = array(TRUE, dim = c(4, 4, 4))),
    "dimensions"
  )
})

test_that("motion_matrices() rejects unsupported shapes", {
  expect_error(motion_matrices(c(1, 2, 3)), "length 6")
  expect_error(motion_matrices(matrix(0, 5, 5)), "must be a length-6")
})

test_that("displacements_within_mask() rejects non-3D mask", {
  expect_error(displacements_within_mask(array(TRUE, dim = c(4, 4)),
                                         c(0, 0, 0, 0, 0, 0)),
               "must be 3D")
})

test_that("displacements_within_mask() rejects empty mask", {
  expect_error(displacements_within_mask(array(FALSE, dim = c(4, 4, 2)),
                                         c(0, 0, 0, 0, 0, 0)),
               "no TRUE voxels")
})

test_that("displacements_within_mask() rejects bad pixdim", {
  m <- array(TRUE, dim = c(4, 4, 2))
  expect_error(displacements_within_mask(m, c(0, 0, 0, 0, 0, 0),
                                         pixdim = c(-1, 1, 1)),
               "positive finite")
})

test_that("motion_displacement_summary() rejects truth/motion length mismatch", {
  m <- array(TRUE, dim = c(4, 4, 2))
  motion <- matrix(0, nrow = 5, ncol = 6)
  truth <- matrix(0, nrow = 4, ncol = 6)
  expect_error(motion_displacement_summary(m, motion, truth),
               "same number of frames")
})
