.rotation_matrix_zyx <- function(rx, ry, rz) {
  cx <- cos(rx); sx <- sin(rx)
  cy <- cos(ry); sy <- sin(ry)
  cz <- cos(rz); sz <- sin(rz)

  matrix(
    c(
      cz * cy,
      cz * sy * sx - sz * cx,
      cz * sy * cx + sz * sx,
      sz * cy,
      sz * sy * sx + cz * cx,
      sz * sy * cx - cz * sx,
      -sy,
      cy * sx,
      cy * cx
    ),
    nrow = 3,
    byrow = TRUE
  )
}

.make_subject_params <- function(seed = 1L) {
  set.seed(seed)
  list(
    head_axes = c(22, 26, 20),
    brain_axes = c(18, 21, 16),
    wm_axes = c(15.5, 18, 13.5),
    asym = c(runif(1, -1.5, 1.5), runif(1, -1.0, 1.0), runif(1, -0.8, 0.8)),
    vent1_center = c(2.8, 0.0, 0.0),
    vent2_center = c(-2.3, 0.1, 0.0),
    vent_axes = c(2.4, 1.8, 2.2),
    tex_amp = runif(1, 5, 9),
    tex_basis = {
      nb <- 12
      list(
        ax = runif(nb, 0.03, 0.10),
        ay = runif(nb, 0.03, 0.10),
        az = runif(nb, 0.03, 0.10),
        ph = runif(nb, 0, 2 * pi),
        w = runif(nb, 0.4, 1.0) * sample(c(-1, 1), nb, replace = TRUE)
      )
    }
  )
}

.continuous_texture <- function(x, y, z, tex) {
  out <- 0
  for (m in seq_along(tex$w)) {
    out <- out + tex$w[m] * sin(2 * pi * (tex$ax[m] * x + tex$ay[m] * y + tex$az[m] * z) + tex$ph[m])
  }
  out / (sqrt(length(tex$w)) + 1e-8)
}

.smooth_inside <- function(r, sharpness = 16) {
  1 / (1 + exp(sharpness * (r - 1)))
}

.ellip_r <- function(x, y, z, center, axes) {
  dx <- (x - center[1]) / axes[1]
  dy <- (y - center[2]) / axes[2]
  dz <- (z - center[3]) / axes[3]
  sqrt(dx * dx + dy * dy + dz * dz)
}

make_continuous_subject <- function(seed = 1L) {
  prm <- .make_subject_params(seed)

  eval_fun <- function(x, y, z) {
    head <- .smooth_inside(.ellip_r(x, y, z, center = c(0, 0, 0), axes = prm$head_axes))
    brain <- .smooth_inside(.ellip_r(x, y, z, center = prm$asym, axes = prm$brain_axes))
    wm <- .smooth_inside(.ellip_r(x, y, z, center = prm$asym + c(0.4, -0.2, 0), axes = prm$wm_axes))

    gm <- pmax(brain - wm, 0)
    csf <- pmax(head - brain, 0)

    v1 <- .smooth_inside(.ellip_r(x, y, z, center = prm$vent1_center + prm$asym, axes = prm$vent_axes), sharpness = 20)
    v2 <- .smooth_inside(.ellip_r(x, y, z, center = prm$vent2_center + prm$asym, axes = prm$vent_axes), sharpness = 20)
    ventricles <- pmin(v1 + v2, 1)

    tex <- .continuous_texture(x, y, z, prm$tex_basis)

    # EPI-like synthetic contrast with asymmetry and texture.
    intensity <- 10 * head + 20 * csf + 75 * gm + 100 * wm - 35 * ventricles + prm$tex_amp * tex * (0.3 * gm + 0.6 * wm)
    pmax(intensity, 0)
  }

  list(params = prm, eval = eval_fun)
}

.render_single_voxel <- function(subject, x, y, z, pose) {
  R <- .rotation_matrix_zyx(pose[4], pose[5], pose[6])
  xf <- x - pose[1]
  yf <- y - pose[2]
  zf <- z - pose[3]
  xm <- R[1, 1] * xf + R[2, 1] * yf + R[3, 1] * zf
  ym <- R[1, 2] * xf + R[2, 2] * yf + R[3, 2] * zf
  zm <- R[1, 3] * xf + R[2, 3] * yf + R[3, 3] * zf
  subject$eval(xm, ym, zm)
}

render_continuous_volume <- function(subject,
                                     nx,
                                     ny,
                                     nz,
                                     pixdim = c(2, 2, 2),
                                     pose = rep(0, 6),
                                     slice_poses = NULL,
                                     supersample = 2L,
                                     bias_amp = 0,
                                     gain = 1,
                                     noise_sd = 0,
                                     seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  supersample <- as.integer(max(1L, supersample))
  offs <- if (supersample == 1L) 0 else seq(-0.5 + 0.5 / supersample, 0.5 - 0.5 / supersample, length.out = supersample)

  out <- array(0, dim = c(nx, ny, nz))
  cx <- 0.5 * (nx - 1)
  cy <- 0.5 * (ny - 1)
  cz <- 0.5 * (nz - 1)

  fx <- max(1e-6, cx * pixdim[1])
  fy <- max(1e-6, cy * pixdim[2])
  fz <- max(1e-6, cz * pixdim[3])

  for (k in seq_len(nz)) {
    pk <- if (is.null(slice_poses)) pose else slice_poses[k, ]
    for (j in seq_len(ny)) {
      for (i in seq_len(nx)) {
        acc <- 0
        nacc <- 0
        for (oz in offs) {
          for (oy in offs) {
            for (ox in offs) {
              x <- ((i - 1) - cx + ox) * pixdim[1]
              y <- ((j - 1) - cy + oy) * pixdim[2]
              z <- ((k - 1) - cz + oz) * pixdim[3]
              acc <- acc + .render_single_voxel(subject, x, y, z, pk)
              nacc <- nacc + 1
            }
          }
        }
        val <- acc / nacc
        bx <- (((i - 1) - cx) * pixdim[1]) / fx
        by <- (((j - 1) - cy) * pixdim[2]) / fy
        bz <- (((k - 1) - cz) * pixdim[3]) / fz
        bias <- 1 + bias_amp * (0.55 * bx - 0.35 * by + 0.20 * bz)
        out[i, j, k] <- gain * bias * val
      }
    }
  }

  if (noise_sd > 0) {
    out <- out + array(rnorm(length(out), sd = noise_sd), dim = dim(out))
  }

  out
}

continuous_brain_mask <- function(subject,
                                  nx,
                                  ny,
                                  nz,
                                  pixdim = c(2, 2, 2),
                                  threshold = 15,
                                  supersample = 1L) {
  ref <- render_continuous_volume(
    subject = subject,
    nx = nx,
    ny = ny,
    nz = nz,
    pixdim = pixdim,
    pose = rep(0, 6),
    supersample = supersample,
    bias_amp = 0,
    gain = 1,
    noise_sd = 0
  )
  ref > threshold
}

simulate_continuous_run <- function(subject,
                                    nx = 16,
                                    ny = 16,
                                    nz = 8,
                                    nt = 10,
                                    pixdim = c(2, 2, 2),
                                    poses,
                                    noise_sd = 2,
                                    drift_amp = 0,
                                    bias_amp = 0,
                                    supersample = 2L,
                                    slice_times = NULL,
                                    pose_slopes = NULL,
                                    seed = 1L) {
  stopifnot(is.matrix(poses), nrow(poses) == nt, ncol(poses) == 6)
  if (!is.null(pose_slopes)) {
    stopifnot(is.matrix(pose_slopes), all(dim(pose_slopes) == dim(poses)))
  }

  if (is.null(slice_times)) {
    slice_times <- rep(0, nz)
  }
  slice_times <- as.numeric(slice_times)
  if (length(slice_times) != nz) {
    stop("slice_times must have length nz", call. = FALSE)
  }

  set.seed(seed)
  run <- array(0, dim = c(nx, ny, nz, nt))

  for (t in seq_len(nt)) {
    phase <- if (nt > 1) (t - 1) / (nt - 1) else 0
    gain_t <- 1 + drift_amp * sin(2 * pi * phase)

    slice_poses <- NULL
    if (!is.null(pose_slopes)) {
      slice_poses <- matrix(0, nrow = nz, ncol = 6)
      for (k in seq_len(nz)) {
        slice_poses[k, ] <- poses[t, ] + pose_slopes[t, ] * slice_times[k]
      }
    }

    run[, , , t] <- render_continuous_volume(
      subject = subject,
      nx = nx,
      ny = ny,
      nz = nz,
      pixdim = pixdim,
      pose = poses[t, ],
      slice_poses = slice_poses,
      supersample = supersample,
      bias_amp = bias_amp,
      gain = gain_t,
      noise_sd = noise_sd
    )
  }

  list(
    run = run,
    poses = poses,
    pose_slopes = pose_slopes,
    pixdim = pixdim,
    slice_times = slice_times
  )
}
