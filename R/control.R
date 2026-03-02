#' Default control parameters for volregger
#'
#' @return A named list of control parameters.
#' @export
volreg_control <- function() {
  list(
    downsample = c(4L, 2L, 1L),
    max_iter = c(40L, 30L, 20L),
    n_samples = c(12000L, 24000L, 40000L),
    huber_k = 1.5,
    lambda0 = 1e-2,
    tol_step = 1e-5,
    tol_cost = 1e-6,
    n_restarts = 0L,
    robust_template = TRUE,
    use_overlap_taper = TRUE,
    use_slice_spline = FALSE,
    nthreads = 1L,
    tr = 1,
    engine = "rigid_robust",
    verbose = FALSE
  )
}
