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
    adaptive_huber = FALSE,
    adaptive_huber_guard = FALSE,
    lambda0 = 1e-2,
    tol_step = 1e-5,
    tol_cost = 1e-6,
    n_restarts = 0L,
    robust_template = TRUE,
    use_overlap_taper = TRUE,
    edge_exclude_frac = 0,
    edge_exclude_vox = 0L,
    capture_boost = TRUE,
    capture_trans_mm = 8,
    capture_rot_deg = 8,
    capture_topk = 6L,
    use_ic_stencil = FALSE,
    stencil_n = 4000L,
    stencil_bx = 4L,
    stencil_by = 4L,
    stencil_bz = 2L,
    stencil_gamma = 0.5,
    use_dense_polish = FALSE,
    polish_iter = 5L,
    polish_n_samples = 0L,
    whiten_modes = FALSE,
    whiten_reg = 1e-6,
    use_nuisance_marginalization = FALSE,
    nuisance_include_template = FALSE,
    nuisance_orthogonalize = TRUE,
    nuisance_ridge = 1e-4,
    nuisance_orth_reg = 1e-6,
    use_slice_spline = FALSE,
    nthreads = 1L,
    tr = 1,
    engine = "rigid_robust",
    verbose = FALSE
  )
}
