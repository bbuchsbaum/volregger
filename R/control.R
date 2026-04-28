#' Control parameters for the volregger optimizer
#'
#' Returns a list of all control fields used by [volreg()] and
#' [estimate_motion()]. Every field has a documented default below; pass
#' overrides as named arguments (e.g. `volreg_control(huber_k = 0.05)`)
#' or via `utils::modifyList()`. Use [volreg_validate_control()] to check
#' that a custom list still satisfies the field-level invariants.
#'
#' @section Pyramid (coarse-to-fine):
#'
#' @param downsample Integer vector of pyramid downsampling factors, finest
#'   level last. Each entry must be `>= 1L`. Default `c(4L, 2L, 1L)`.
#' @param use_pyramid Logical scalar. When `TRUE`, the optimizer iterates
#'   over the full `downsample` schedule; when `FALSE`, only the finest
#'   level is used. Default `FALSE`.
#' @param max_iter Integer vector of maximum LM iterations per pyramid
#'   level. Length must equal `length(downsample)` when `use_pyramid =
#'   TRUE`. A scalar is broadcast. Default `c(30L, 20L, 12L)`.
#' @param n_samples Integer vector of voxel-sample counts per pyramid level.
#'   Same length rule as `max_iter`. Default `c(12000L, 24000L, 40000L)`.
#'
#' @section Optimization:
#'
#' @param huber_k Positive scalar. Huber-loss break point in robust-residual
#'   units; lower values down-weight outliers more aggressively. Default
#'   `1.5`.
#' @param adaptive_huber Logical. If `TRUE`, `huber_k` is rescaled by the
#'   per-iteration robust residual sigma. Default `FALSE`.
#' @param adaptive_huber_guard Logical. If `TRUE`, guards adaptive scaling
#'   from collapsing to zero on near-perfect frames. Default `FALSE`.
#' @param lambda0 Non-negative scalar. Initial Levenberg-Marquardt damping.
#'   Default `1e-2`.
#' @param tol_step Non-negative scalar. Convergence threshold on the
#'   pose-step norm. Default `1e-5`.
#' @param tol_cost Non-negative scalar. Convergence threshold on the
#'   relative cost decrease between iterations. Default `1e-6`.
#' @param n_restarts Non-negative integer. Number of additional
#'   capture-boost restart seeds tried when the primary fit is poor.
#'   Default `0L`.
#' @param robust_template Logical. If `TRUE`, the per-voxel template is a
#'   trimmed mean across frames; if `FALSE`, a plain mean. Default `TRUE`.
#' @param template_refresh_valid_only Logical. If `TRUE`, the template
#'   refresh step only averages over frames whose fits passed validity
#'   checks. Default `TRUE`.
#'
#' @section Coverage / mask:
#'
#' @param use_overlap_taper Logical. If `TRUE`, the residual is tapered
#'   along the boundary of the voxel-overlap region between reference
#'   and warped frame. Default `TRUE`.
#' @param edge_exclude_frac Scalar in `[0, 1)`. Fraction of voxels at each
#'   edge to drop from the optimization mask. Default `0.05`.
#' @param edge_exclude_vox Non-negative integer. Absolute number of edge
#'   voxels to drop in addition to `edge_exclude_frac`. Default `0L`.
#'
#' @section Capture boost (initial pose seeding):
#'
#' @param capture_boost Logical. If `TRUE`, seed each frame with a small
#'   grid of translation/rotation candidates and pick the best initialization
#'   before LM iteration. Default `TRUE`.
#' @param capture_trans_mm Positive scalar. Half-width (mm) of the
#'   translation grid evaluated by capture boost. Default `8`.
#' @param capture_rot_deg Positive scalar. Half-width (deg) of the rotation
#'   grid evaluated by capture boost. Default `8`.
#' @param capture_topk Positive integer. Number of best-scoring grid
#'   seeds carried forward into LM iteration. Default `4L`.
#'
#' @section IC stencil sampling:
#'
#' @param use_ic_stencil Logical. If `TRUE`, sample voxels through a
#'   structured information-content stencil rather than uniformly. Default
#'   `TRUE`.
#' @param stencil_n Positive integer. Number of stencil sample points per
#'   frame. Default `1000L`.
#' @param stencil_bx,stencil_by,stencil_bz Positive integers. Bin counts
#'   along each axis used to bucket voxels by information content.
#'   Defaults `4L`, `4L`, `2L`.
#' @param stencil_gamma Positive scalar. Tempering exponent applied to
#'   per-bin importance weights. Default `0.5`.
#'
#' @section Polish (post-LM refinement):
#'
#' @param use_ic_polish Logical. If `TRUE`, run an extra IC-stencil polish
#'   pass after LM. Default `FALSE`.
#' @param use_dense_polish Logical. If `TRUE`, run a dense (full-voxel)
#'   polish pass after LM. Default `FALSE`.
#' @param polish_iter Non-negative integer. Iterations per polish pass.
#'   Default `5L`.
#' @param polish_n_samples Non-negative integer. Voxel sample count for
#'   dense polish (`0L` lets the optimizer choose). Default `0L`.
#' @param dense_interp,dense_polish_interp Strings. Interpolation methods
#'   used during dense sampling and polish, respectively. Must be one of
#'   `"linear"`, `"bspline"`, `"lanczos"`, `"lanczos5"`, `"lanczos7"`,
#'   `"quintic"`, or `"heptic"`. Defaults `"linear"`, `"linear"`.
#'
#' @section Temporal regularization:
#'
#' @param use_temporal_regularization Logical. If `TRUE`, enable the
#'   temporal-smoothing pose update with cost-improvement guard.
#'   Default `FALSE`.
#' @param temporal_regularization_iter Positive integer. Number of
#'   smoothing passes applied. Default `1L`.
#' @param temporal_regularization_alpha Scalar in `[0, 1]`. Mixing weight
#'   between original and neighbor-averaged pose. Default `0.25`.
#' @param temporal_regularization_cost_tol Non-negative scalar. Maximum
#'   tolerated cost increase per smoothing pass. Default `0.002`.
#' @param temporal_regularization_fd_max Positive scalar. FD threshold
#'   above which the smoothing guard kicks in. Default `0.8`.
#'
#' @section SE(3) updates:
#'
#' @param use_se3_updates Logical. If `TRUE`, apply pose updates in
#'   SE(3) (composition) rather than additive parameter space.
#'   Default `FALSE`.
#' @param se3_update_order String, one of `"left"` or `"right"`. Side on
#'   which the SE(3) increment is composed. Default `"left"`.
#'
#' @section Temporal normal prior:
#'
#' @param use_temporal_normal_prior Logical. If `TRUE`, add a Gaussian
#'   neighbor-pose prior to the cost. Default `FALSE`.
#' @param temporal_normal_prior_weight Non-negative scalar. Prior weight.
#'   Default `0`.
#' @param temporal_normal_prior_rot_radius Positive scalar. Equivalent
#'   radius (mm) used to scale rotation contributions to the prior.
#'   Default `50`.
#'
#' @section Low-motion stabilization:
#'
#' @param use_low_motion_pose_shrink Logical. If `TRUE`, shrink pose
#'   parameters toward zero on runs whose `fd_p95` is below
#'   `low_motion_fd_p95_max`. Default `FALSE`.
#' @param low_motion_pose_scale Scalar in `(0, 1]`. Multiplicative shrink
#'   factor. Default `0.95`.
#' @param low_motion_fd_p95_max Non-negative scalar. FD-p95 ceiling for
#'   the low-motion regime. Default `0.8`.
#' @param use_low_motion_pose_smooth Logical. If `TRUE`, neighbor-smooth
#'   the pose trace under the low-motion regime. Default `FALSE`.
#' @param low_motion_pose_smooth_alpha Scalar in `[0, 1]`. Smoothing
#'   weight. Default `0.3`.
#' @param low_motion_pose_smooth_iter Positive integer. Number of
#'   neighbor-smoothing passes when low-motion pose smoothing is enabled.
#'   Default `1L`.
#' @param use_low_motion_output_smooth Logical. If `TRUE`, neighbor-smooth
#'   the resampled output under the low-motion regime. Default `FALSE`.
#' @param low_motion_output_smooth_alpha Scalar in `[0, 1]`. Output-smooth
#'   weight. Default `0.02`.
#' @param use_low_motion_spatial_smooth Logical. If `TRUE`, apply a
#'   low-gradient, edge-preserving spatial smooth to the resampled output
#'   under the low-motion regime. Default `FALSE`.
#' @param low_motion_spatial_smooth_alpha Scalar in `[0, 1]`. Spatial
#'   smoothing weight. Default `0.08`.
#' @param low_motion_spatial_smooth_quantile Scalar in `[0, 1]`. Quantile of
#'   mean-image local-gradient voxels eligible for spatial smoothing. Default
#'   `0.5`.
#'
#' @section Whitening:
#'
#' @param whiten_modes Logical. If `TRUE`, whiten residuals against the
#'   leading template modes. Default `TRUE`.
#' @param whiten_reg Positive scalar. Ridge added to the whitening
#'   covariance for numerical stability. Default `1e-6`.
#'
#' @section Nuisance marginalization:
#'
#' @param use_nuisance_marginalization Logical. If `TRUE`, project a
#'   nuisance basis out of the residual at each iteration. Default
#'   `FALSE`.
#' @param nuisance_include_template Logical. If `TRUE`, include template
#'   modes in the nuisance basis. Default `FALSE`.
#' @param nuisance_orthogonalize Logical. If `TRUE`, orthogonalize the
#'   nuisance basis before regression. Default `TRUE`.
#' @param nuisance_ridge Positive scalar. Ridge added to the nuisance
#'   normal equations. Default `1e-4`.
#' @param nuisance_orth_reg Positive scalar. Ridge applied during
#'   orthogonalization. Default `1e-6`.
#'
#' @section Engine, threading, misc:
#'
#' @param use_slice_spline Logical. If `TRUE`, treat the run as
#'   slice/packet-timed and apply the spline-smoothed engine. Set
#'   automatically by `engine = "rigid_spline"`. Default `FALSE`.
#' @param parallel_frames Logical. If `TRUE`, fit independent frames in
#'   parallel using up to `nthreads` workers. Default `FALSE`.
#' @param nthreads Positive integer. Worker thread count. Default `1L`.
#' @param tr Positive scalar. Repetition time (s); set automatically from
#'   input metadata when available. Default `1`.
#' @param engine String, one of `"rigid_robust"` or `"rigid_spline"`.
#'   Engine identifier. Default `"rigid_robust"`.
#' @param verbose Logical. If `TRUE`, surface per-iteration diagnostics
#'   from the C++ core. Default `FALSE`.
#'
#' @return A named list of every control field shown above, in the order
#'   listed.
#'
#' @seealso [volreg_validate_control()] for invariant checks; [volreg_profile()]
#'   for named presets.
#' @examples
#' ctl <- volreg_control(huber_k = 0.05, nthreads = 2L)
#' ctl$huber_k
#' @export
volreg_control <- function(
    downsample = c(4L, 2L, 1L),
    use_pyramid = FALSE,
    max_iter = c(30L, 20L, 12L),
    n_samples = c(12000L, 24000L, 40000L),
    huber_k = 1.5,
    adaptive_huber = FALSE,
    adaptive_huber_guard = FALSE,
    lambda0 = 1e-2,
    tol_step = 1e-5,
    tol_cost = 1e-6,
    n_restarts = 0L,
    robust_template = TRUE,
    template_refresh_valid_only = TRUE,
    use_overlap_taper = TRUE,
    edge_exclude_frac = 0.05,
    edge_exclude_vox = 0L,
    capture_boost = TRUE,
    capture_trans_mm = 8,
    capture_rot_deg = 8,
    capture_topk = 4L,
    use_ic_stencil = TRUE,
    stencil_n = 1000L,
    stencil_bx = 4L,
    stencil_by = 4L,
    stencil_bz = 2L,
    stencil_gamma = 0.5,
    use_ic_polish = FALSE,
    use_dense_polish = FALSE,
    polish_iter = 5L,
    polish_n_samples = 0L,
    dense_interp = "linear",
    dense_polish_interp = "linear",
    use_temporal_regularization = FALSE,
    temporal_regularization_iter = 1L,
    temporal_regularization_alpha = 0.25,
    temporal_regularization_cost_tol = 0.002,
    temporal_regularization_fd_max = 0.8,
    use_se3_updates = FALSE,
    se3_update_order = "left",
    use_temporal_normal_prior = FALSE,
    temporal_normal_prior_weight = 0,
    temporal_normal_prior_rot_radius = 50,
    use_low_motion_pose_shrink = FALSE,
    low_motion_pose_scale = 0.95,
    low_motion_fd_p95_max = 0.8,
    use_low_motion_pose_smooth = FALSE,
    low_motion_pose_smooth_alpha = 0.3,
    low_motion_pose_smooth_iter = 1L,
    use_low_motion_output_smooth = FALSE,
    low_motion_output_smooth_alpha = 0.02,
    use_low_motion_spatial_smooth = FALSE,
    low_motion_spatial_smooth_alpha = 0.08,
    low_motion_spatial_smooth_quantile = 0.5,
    parallel_frames = FALSE,
    whiten_modes = TRUE,
    whiten_reg = 1e-6,
    use_nuisance_marginalization = FALSE,
    nuisance_include_template = FALSE,
    nuisance_orthogonalize = TRUE,
    nuisance_ridge = 1e-4,
    nuisance_orth_reg = 1e-6,
    use_slice_spline = FALSE,
    nthreads = 1L,
    tr = 1,
    engine = c("rigid_robust", "rigid_spline"),
    verbose = FALSE) {
  if (length(engine) > 1L) engine <- engine[[1L]]
  list(
    downsample = downsample,
    use_pyramid = use_pyramid,
    max_iter = max_iter,
    n_samples = n_samples,
    huber_k = huber_k,
    adaptive_huber = adaptive_huber,
    adaptive_huber_guard = adaptive_huber_guard,
    lambda0 = lambda0,
    tol_step = tol_step,
    tol_cost = tol_cost,
    n_restarts = n_restarts,
    robust_template = robust_template,
    template_refresh_valid_only = template_refresh_valid_only,
    use_overlap_taper = use_overlap_taper,
    edge_exclude_frac = edge_exclude_frac,
    edge_exclude_vox = edge_exclude_vox,
    capture_boost = capture_boost,
    capture_trans_mm = capture_trans_mm,
    capture_rot_deg = capture_rot_deg,
    capture_topk = capture_topk,
    use_ic_stencil = use_ic_stencil,
    stencil_n = stencil_n,
    stencil_bx = stencil_bx,
    stencil_by = stencil_by,
    stencil_bz = stencil_bz,
    stencil_gamma = stencil_gamma,
    use_ic_polish = use_ic_polish,
    use_dense_polish = use_dense_polish,
    polish_iter = polish_iter,
    polish_n_samples = polish_n_samples,
    dense_interp = dense_interp,
    dense_polish_interp = dense_polish_interp,
    use_temporal_regularization = use_temporal_regularization,
    temporal_regularization_iter = temporal_regularization_iter,
    temporal_regularization_alpha = temporal_regularization_alpha,
    temporal_regularization_cost_tol = temporal_regularization_cost_tol,
    temporal_regularization_fd_max = temporal_regularization_fd_max,
    use_se3_updates = use_se3_updates,
    se3_update_order = se3_update_order,
    use_temporal_normal_prior = use_temporal_normal_prior,
    temporal_normal_prior_weight = temporal_normal_prior_weight,
    temporal_normal_prior_rot_radius = temporal_normal_prior_rot_radius,
    use_low_motion_pose_shrink = use_low_motion_pose_shrink,
    low_motion_pose_scale = low_motion_pose_scale,
    low_motion_fd_p95_max = low_motion_fd_p95_max,
    use_low_motion_pose_smooth = use_low_motion_pose_smooth,
    low_motion_pose_smooth_alpha = low_motion_pose_smooth_alpha,
    low_motion_pose_smooth_iter = low_motion_pose_smooth_iter,
    use_low_motion_output_smooth = use_low_motion_output_smooth,
    low_motion_output_smooth_alpha = low_motion_output_smooth_alpha,
    use_low_motion_spatial_smooth = use_low_motion_spatial_smooth,
    low_motion_spatial_smooth_alpha = low_motion_spatial_smooth_alpha,
    low_motion_spatial_smooth_quantile = low_motion_spatial_smooth_quantile,
    parallel_frames = parallel_frames,
    whiten_modes = whiten_modes,
    whiten_reg = whiten_reg,
    use_nuisance_marginalization = use_nuisance_marginalization,
    nuisance_include_template = nuisance_include_template,
    nuisance_orthogonalize = nuisance_orthogonalize,
    nuisance_ridge = nuisance_ridge,
    nuisance_orth_reg = nuisance_orth_reg,
    use_slice_spline = use_slice_spline,
    nthreads = nthreads,
    tr = tr,
    engine = engine,
    verbose = verbose
  )
}

#' Validate a volregger control list
#'
#' Checks that every recognized field of `control` has the right type and
#' lies in the allowed range, and that cross-field invariants hold.
#' Unknown fields (e.g. profile-specific extras like `interp_final`,
#' `profile`, `ref_mode`) are tolerated but not checked.
#'
#' @param control A named list, typically from [volreg_control()] or
#'   [volreg_profile()]`$control`.
#' @param strict Logical. If `TRUE` (the default), invalid fields cause
#'   `stop()`; if `FALSE`, returns the diagnostic messages without raising.
#'
#' @return Invisibly, the validated `control` list. With `strict = FALSE`,
#'   a character vector of issue strings (empty when valid).
#' @examples
#' volreg_validate_control(volreg_control())
#' bad <- volreg_control(huber_k = -1)
#' volreg_validate_control(bad, strict = FALSE)
#' @export
volreg_validate_control <- function(control, strict = TRUE) {
  if (!is.list(control) || is.null(names(control))) {
    msg <- "control must be a named list"
    if (isTRUE(strict)) stop(msg, call. = FALSE) else return(msg)
  }
  issues <- character(0)
  add <- function(...) issues <<- c(issues, sprintf(...))

  is_scalar <- function(x) length(x) == 1L && !is.na(x)
  is_pos_int <- function(x) is_scalar(x) && is.finite(x) && x == round(x) && x >= 1L
  is_nn_int  <- function(x) is_scalar(x) && is.finite(x) && x == round(x) && x >= 0L
  is_pos_num <- function(x) is_scalar(x) && is.numeric(x) && is.finite(x) && x > 0
  is_nn_num  <- function(x) is_scalar(x) && is.numeric(x) && is.finite(x) && x >= 0
  is_in_unit <- function(x) is_scalar(x) && is.numeric(x) && is.finite(x) && x >= 0 && x <= 1
  is_lgl     <- function(x) is_scalar(x) && is.logical(x)

  has <- function(k) k %in% names(control)
  v <- function(k) control[[k]]

  # Pyramid
  if (has("downsample")) {
    d <- v("downsample")
    if (!is.numeric(d) || length(d) < 1L || any(!is.finite(d)) || any(d < 1) || any(d != round(d))) {
      add("downsample must be an integer vector with all entries >= 1")
    }
  }
  if (has("use_pyramid") && !is_lgl(v("use_pyramid"))) add("use_pyramid must be a single TRUE/FALSE")
  for (k in c("max_iter", "n_samples")) {
    if (has(k)) {
      m <- v(k)
      if (!is.numeric(m) || length(m) < 1L || any(!is.finite(m)) || any(m < 0) || any(m != round(m))) {
        add("%s must be a non-negative integer vector", k)
      }
    }
  }
  if (isTRUE(v("use_pyramid")) && has("downsample") && has("max_iter")) {
    if (length(v("max_iter")) > 1L && length(v("max_iter")) != length(v("downsample"))) {
      add("max_iter length (%d) must equal downsample length (%d) when use_pyramid = TRUE",
          length(v("max_iter")), length(v("downsample")))
    }
  }

  # Optimization
  if (has("huber_k") && !is_pos_num(v("huber_k"))) add("huber_k must be a positive scalar")
  for (k in c("lambda0", "tol_step", "tol_cost", "whiten_reg",
              "nuisance_ridge", "nuisance_orth_reg",
              "stencil_gamma", "temporal_normal_prior_rot_radius",
              "temporal_regularization_fd_max",
              "capture_trans_mm", "capture_rot_deg")) {
    if (has(k) && !is_pos_num(v(k))) add("%s must be a positive scalar", k)
  }
  for (k in c("temporal_regularization_cost_tol", "temporal_normal_prior_weight",
              "low_motion_fd_p95_max")) {
    if (has(k) && !is_nn_num(v(k))) add("%s must be a non-negative scalar", k)
  }
  for (k in c("n_restarts", "edge_exclude_vox", "polish_iter", "polish_n_samples")) {
    if (has(k) && !is_nn_int(v(k))) add("%s must be a non-negative integer scalar", k)
  }
  for (k in c("capture_topk", "stencil_n", "stencil_bx", "stencil_by", "stencil_bz",
              "temporal_regularization_iter", "low_motion_pose_smooth_iter",
              "nthreads")) {
    if (has(k) && !is_pos_int(v(k))) add("%s must be a positive integer scalar", k)
  }
  for (k in c("temporal_regularization_alpha",
              "low_motion_pose_smooth_alpha",
              "low_motion_output_smooth_alpha",
              "low_motion_spatial_smooth_alpha",
              "low_motion_spatial_smooth_quantile")) {
    if (has(k) && !is_in_unit(v(k))) add("%s must be a scalar in [0, 1]", k)
  }
  for (k in c("adaptive_huber", "adaptive_huber_guard", "robust_template",
              "template_refresh_valid_only", "use_overlap_taper",
              "capture_boost", "use_ic_stencil", "use_ic_polish",
              "use_dense_polish", "use_temporal_regularization",
              "use_se3_updates", "use_temporal_normal_prior",
              "use_low_motion_pose_shrink", "use_low_motion_pose_smooth",
              "use_low_motion_output_smooth", "use_low_motion_spatial_smooth",
              "parallel_frames",
              "whiten_modes", "use_nuisance_marginalization",
              "nuisance_include_template", "nuisance_orthogonalize",
              "use_slice_spline", "verbose")) {
    if (has(k) && !is_lgl(v(k))) add("%s must be a single TRUE/FALSE", k)
  }
  if (has("edge_exclude_frac") &&
      !(is_scalar(v("edge_exclude_frac")) &&
        is.numeric(v("edge_exclude_frac")) &&
        v("edge_exclude_frac") >= 0 && v("edge_exclude_frac") < 1)) {
    add("edge_exclude_frac must be a scalar in [0, 1)")
  }
  if (has("low_motion_pose_scale")) {
    s <- v("low_motion_pose_scale")
    if (!is_scalar(s) || !is.numeric(s) || s <= 0 || s > 1) {
      add("low_motion_pose_scale must be a scalar in (0, 1]")
    }
  }
  if (has("tr") && !is_pos_num(v("tr"))) add("tr must be a positive scalar")

  valid_interp <- c("linear", "bspline", "lanczos", "lanczos5", "lanczos7",
                    "quintic", "heptic")
  for (k in c("dense_interp", "dense_polish_interp")) {
    if (has(k)) {
      x <- v(k)
      if (!is_scalar(x) || !is.character(x) || !x %in% valid_interp) {
        add("%s must be one of: %s", k, paste(valid_interp, collapse = ", "))
      }
    }
  }
  if (has("se3_update_order")) {
    x <- v("se3_update_order")
    if (!is_scalar(x) || !is.character(x) || !x %in% c("left", "right")) {
      add("se3_update_order must be 'left' or 'right'")
    }
  }
  if (has("engine")) {
    x <- v("engine")
    if (!is_scalar(x) || !is.character(x) || !x %in% c("rigid_robust", "rigid_spline")) {
      add("engine must be 'rigid_robust' or 'rigid_spline'")
    }
  }

  if (length(issues)) {
    msg <- paste0("invalid volreg control:\n  - ", paste(issues, collapse = "\n  - "))
    if (isTRUE(strict)) stop(msg, call. = FALSE) else return(issues)
  }
  if (isTRUE(strict)) invisible(control) else character(0)
}
