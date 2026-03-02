#' Compute QC metrics for a motion estimate
#'
#' @param x Original 4D image-like input.
#' @param motion Motion matrix (nt x 6).
#' @param corrected Optional corrected 4D data.
#' @param mask Optional 3D mask.
#' @param cost_init Optional initial costs.
#' @param cost_final Optional final costs.
#' @return A list of QC metrics.
#' @export
motion_qc <- function(
    x,
    motion,
    corrected = NULL,
    mask = NULL,
    cost_init = NULL,
    cost_final = NULL) {
  inp <- .as_volreg_input(x = x, mask = mask)
  motion <- .require_motion_matrix(motion)

  nt <- inp$dim[4]
  if (is.null(cost_init)) {
    cost_init <- rep(0, nt)
  }
  if (is.null(cost_final)) {
    cost_final <- rep(0, nt)
  }

  cpp_volreg_metrics(
    original = inp$data,
    corrected = if (is.null(corrected)) NULL else as.array(corrected),
    dim = inp$dim,
    mask = inp$mask,
    par = motion,
    cost_init = as.numeric(cost_init),
    cost_final = as.numeric(cost_final)
  )
}
