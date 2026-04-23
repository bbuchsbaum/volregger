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
      "ic_stencil",
      "ic_whiten",
      "ic_nuisance",
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
        nuisance_include_template = FALSE
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
