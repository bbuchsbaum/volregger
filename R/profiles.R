#' Named estimator profiles for ablation and benchmark work
#'
#' @param name Profile name.
#' @param control Base control list to modify.
#' @param overrides Named list of control overrides applied last.
#' @return A `volreg_profile` object with fields `name`, `engine`,
#'   `components`, and `control`.
#' @export
volreg_profile <- function(
    name = c(
      "dense_baseline",
      "fast_native",
      "fast_native_parallel",
      "fast_native_pyramid",
      "fast_fmri",
      "ic_stencil",
      "ic_whiten",
      "ic_nuisance",
      "ic_nuisance_adaptive",
      "ic_nuisance_adaptive_tuned",
      "ic_nuisance_adaptive_guarded",
      "ic_nuisance_polish",
      "ic_nuisance_kernel_polish",
      "ic_nuisance_temporal_guard",
      "ic_nuisance_se3",
      "ic_nuisance_temporal_prior",
      "ic_nuisance_low_motion_tsnr",
      "ic_nuisance_full",
      "slice_spline"
    ),
    control = volreg_control(),
    overrides = list()) {
  name <- match.arg(name)

  spec <- switch(
    name,
    dense_baseline = list(
      engine = "rigid_robust",
      components = c("capture_boost", "dense_sampling"),
      control = list(
        use_ic_stencil = FALSE,
        whiten_modes = FALSE,
        use_nuisance_marginalization = FALSE,
        nuisance_include_template = FALSE,
        edge_exclude_frac = 0,
        template_refresh_valid_only = FALSE
      )
    ),
    fast_native = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "fast_native"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = FALSE,
        nuisance_include_template = FALSE,
        stencil_n = 1000L,
        capture_topk = 4L,
        max_iter = c(30L, 20L, 12L)
      )
    ),
    fast_native_parallel = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "fast_native", "parallel_frames"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = FALSE,
        nuisance_include_template = FALSE,
        stencil_n = 1000L,
        capture_topk = 4L,
        max_iter = c(30L, 20L, 12L),
        parallel_frames = TRUE
      )
    ),
    fast_native_pyramid = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "fast_native", "pyramid"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = FALSE,
        nuisance_include_template = FALSE,
        stencil_n = 1000L,
        capture_topk = 4L,
        use_pyramid = TRUE,
        downsample = c(4L, 2L, 1L),
        max_iter = c(4L, 3L, 4L),
        n_samples = c(1000L, 2000L, 8000L)
      )
    ),
    fast_fmri = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "fast_native", "edge_exclude", "valid_template_refresh", "parallel_frames", "heptic_final"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = FALSE,
        nuisance_include_template = FALSE,
        stencil_n = 1000L,
        capture_topk = 4L,
        max_iter = c(30L, 20L, 12L),
        edge_exclude_frac = 0.05,
        template_refresh_valid_only = TRUE,
        parallel_frames = TRUE,
        interp_final = "heptic"
      )
    ),
    ic_stencil = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = FALSE,
        use_nuisance_marginalization = FALSE,
        nuisance_include_template = FALSE
      )
    ),
    ic_whiten = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = FALSE,
        nuisance_include_template = FALSE
      )
    ),
    ic_nuisance = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4
      )
    ),
    ic_nuisance_adaptive = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE
      )
    ),
    ic_nuisance_adaptive_tuned = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber", "huber_0.15", "iter_12", "parallel_frames", "heptic_final"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE,
        template_refresh_valid_only = TRUE,
        huber_k = 0.15,
        max_iter = 12L,
        parallel_frames = TRUE,
        interp_final = "heptic"
      )
    ),
    ic_nuisance_adaptive_guarded = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber", "huber_0.5", "guarded"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE,
        adaptive_huber_guard = TRUE,
        huber_k = 0.5
      )
    ),
    ic_nuisance_polish = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber", "huber_0.15", "iter_12", "parallel_frames", "ic_polish"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE,
        template_refresh_valid_only = TRUE,
        huber_k = 0.15,
        max_iter = 12L,
        parallel_frames = TRUE,
        use_ic_polish = TRUE,
        polish_iter = 3L,
        polish_n_samples = 2000L
      )
    ),
    ic_nuisance_kernel_polish = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber", "huber_0.15", "iter_12", "parallel_frames", "dense_polish", "kernel_matched_polish"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE,
        template_refresh_valid_only = TRUE,
        huber_k = 0.15,
        max_iter = 12L,
        parallel_frames = TRUE,
        use_dense_polish = TRUE,
        polish_iter = 2L,
        polish_n_samples = 1500L,
        dense_polish_interp = "lanczos"
      )
    ),
    ic_nuisance_temporal_guard = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber", "huber_0.15", "iter_12", "parallel_frames", "temporal_guard"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE,
        template_refresh_valid_only = TRUE,
        huber_k = 0.15,
        max_iter = 12L,
        parallel_frames = TRUE,
        use_temporal_regularization = TRUE,
        temporal_regularization_iter = 1L,
        temporal_regularization_alpha = 0.25,
        temporal_regularization_cost_tol = 0.002,
        temporal_regularization_fd_max = 0.8
      )
    ),
    ic_nuisance_se3 = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber", "huber_0.15", "iter_12", "parallel_frames", "se3_updates"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE,
        template_refresh_valid_only = TRUE,
        huber_k = 0.15,
        max_iter = 12L,
        parallel_frames = TRUE,
        use_se3_updates = TRUE,
        se3_update_order = "left"
      )
    ),
    ic_nuisance_temporal_prior = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber", "huber_0.15", "iter_12", "parallel_frames", "temporal_normal_prior"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE,
        template_refresh_valid_only = TRUE,
        huber_k = 0.15,
        max_iter = 12L,
        parallel_frames = TRUE,
        use_temporal_normal_prior = TRUE,
        temporal_normal_prior_weight = 0.02,
        temporal_normal_prior_rot_radius = 50
      )
    ),
    ic_nuisance_low_motion_tsnr = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "adaptive_huber", "huber_0.15", "iter_12", "parallel_frames", "heptic_final", "low_motion_pose_smooth", "low_motion_pose_shrink", "low_motion_output_smooth"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = FALSE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4,
        adaptive_huber = TRUE,
        template_refresh_valid_only = TRUE,
        huber_k = 0.15,
        max_iter = 12L,
        parallel_frames = TRUE,
        interp_final = "heptic",
        use_low_motion_pose_smooth = TRUE,
        low_motion_pose_smooth_alpha = 0.3,
        use_low_motion_pose_shrink = TRUE,
        low_motion_pose_scale = 0.97,
        low_motion_fd_p95_max = 0.8,
        use_low_motion_output_smooth = TRUE,
        low_motion_output_smooth_alpha = 0.02
      )
    ),
    ic_nuisance_full = list(
      engine = "rigid_robust",
      components = c("capture_boost", "ic_stencil", "whiten", "nuisance", "template_basis"),
      control = list(
        use_ic_stencil = TRUE,
        whiten_modes = TRUE,
        use_nuisance_marginalization = TRUE,
        nuisance_include_template = TRUE,
        nuisance_orthogonalize = TRUE,
        nuisance_ridge = 1e-4
      )
    ),
    slice_spline = list(
      engine = "rigid_spline",
      components = c("capture_boost", "slice_spline"),
      control = list(
        use_slice_spline = TRUE
      )
    )
  )

  ctl <- utils::modifyList(control, spec$control)
  if (length(overrides) > 0L) {
    ctl <- utils::modifyList(ctl, overrides)
  }

  structure(
    list(
      name = name,
      engine = spec$engine,
      components = spec$components,
      control = ctl
    ),
    class = "volreg_profile"
  )
}

#' @export
print.volreg_profile <- function(x, ...) {
  cat("<volreg_profile>\n", sep = "")
  cat("  name:", x$name, "\n")
  cat("  engine:", x$engine, "\n")
  cat("  components:", paste(x$components, collapse = ", "), "\n")
  invisible(x)
}

.apply_volreg_profile <- function(control, profile) {
  if (is.null(profile) || identical(profile, FALSE) || identical(profile, "none")) {
    return(list(control = control, engine = NULL, name = NULL))
  }
  prof <- if (inherits(profile, "volreg_profile")) {
    profile
  } else {
    volreg_profile(as.character(profile), control = control)
  }
  prof$control$profile <- prof$name
  list(control = prof$control, engine = prof$engine, name = prof$name)
}
