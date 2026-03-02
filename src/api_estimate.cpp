#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

inline size_t idx4(int i, int j, int k, int t, int nx, int ny, int nz) {
  return static_cast<size_t>(i + nx * (j + ny * (k + nz * t)));
}

inline int clamp_int(int x, int lo, int hi) {
  return std::max(lo, std::min(hi, x));
}

struct Pose6 {
  double tx;
  double ty;
  double tz;
  double rx;
  double ry;
  double rz;
};

Pose6 make_pose(double tx = 0.0,
                double ty = 0.0,
                double tz = 0.0,
                double rx = 0.0,
                double ry = 0.0,
                double rz = 0.0) {
  return Pose6{tx, ty, tz, rx, ry, rz};
}

Pose6 add_pose(const Pose6& a, const Pose6& b) {
  return make_pose(a.tx + b.tx,
                   a.ty + b.ty,
                   a.tz + b.tz,
                   a.rx + b.rx,
                   a.ry + b.ry,
                   a.rz + b.rz);
}

Pose6 sub_pose(const Pose6& a, const Pose6& b) {
  return make_pose(a.tx - b.tx,
                   a.ty - b.ty,
                   a.tz - b.tz,
                   a.rx - b.rx,
                   a.ry - b.ry,
                   a.rz - b.rz);
}

Pose6 scale_pose(const Pose6& a, double s) {
  return make_pose(a.tx * s, a.ty * s, a.tz * s, a.rx * s, a.ry * s, a.rz * s);
}

NumericVector extract_ref_volume(const NumericVector& data,
                                 int nx,
                                 int ny,
                                 int nz,
                                 int ref_index_1based) {
  NumericVector out(nx * ny * nz);
  int t = std::max(0, ref_index_1based - 1);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        size_t src = idx4(i, j, k, t, nx, ny, nz);
        size_t dst = static_cast<size_t>(i + nx * (j + ny * k));
        out[dst] = data[src];
      }
    }
  }
  out.attr("dim") = IntegerVector::create(nx, ny, nz);
  return out;
}

NumericVector trimmed_mean_template(const NumericVector& data,
                                    int nx,
                                    int ny,
                                    int nz,
                                    int nt,
                                    double trim_frac) {
  int nxyz = nx * ny * nz;
  NumericVector out(nxyz);
  std::vector<double> vals(nt);
  int trim_n = static_cast<int>(std::floor(trim_frac * nt));
  trim_n = std::max(0, std::min(trim_n, (nt - 1) / 2));

  for (int v = 0; v < nxyz; ++v) {
    for (int t = 0; t < nt; ++t) {
      vals[t] = data[static_cast<size_t>(v + nxyz * t)];
    }
    std::sort(vals.begin(), vals.end());
    int lo = trim_n;
    int hi = nt - trim_n;
    double s = 0.0;
    for (int i = lo; i < hi; ++i) {
      s += vals[i];
    }
    out[v] = s / static_cast<double>(hi - lo);
  }

  out.attr("dim") = IntegerVector::create(nx, ny, nz);
  return out;
}

LogicalVector build_mask(const Nullable<LogicalVector>& mask, int nxyz) {
  if (mask.isNotNull()) {
    LogicalVector m(mask);
    if (m.size() != nxyz) {
      stop("mask has incompatible size");
    }
    return m;
  }
  return LogicalVector(nxyz, true);
}

NumericMatrix identity_mats(int nt) {
  NumericMatrix mats(16, nt);
  for (int t = 0; t < nt; ++t) {
    mats(0, t) = 1.0;
    mats(5, t) = 1.0;
    mats(10, t) = 1.0;
    mats(15, t) = 1.0;
  }
  return mats;
}

inline void rotation_matrix(double rx, double ry, double rz, double* R) {
  // R = Rz * Ry * Rx
  double cx = std::cos(rx), sx = std::sin(rx);
  double cy = std::cos(ry), sy = std::sin(ry);
  double cz = std::cos(rz), sz = std::sin(rz);

  R[0] = cz * cy;
  R[1] = cz * sy * sx - sz * cx;
  R[2] = cz * sy * cx + sz * sx;

  R[3] = sz * cy;
  R[4] = sz * sy * sx + cz * cx;
  R[5] = sz * sy * cx - cz * sx;

  R[6] = -sy;
  R[7] = cy * sx;
  R[8] = cy * cx;
}

struct CostResult {
  double cost;
  double overlap;
};

double sample_trilinear_frame(const NumericVector& data,
                              int nx,
                              int ny,
                              int nz,
                              int nxyz,
                              int t,
                              double x,
                              double y,
                              double z) {
  int x0 = static_cast<int>(std::floor(x));
  int y0 = static_cast<int>(std::floor(y));
  int z0 = static_cast<int>(std::floor(z));
  int x1 = x0 + 1;
  int y1 = y0 + 1;
  int z1 = z0 + 1;

  double ax = x - x0;
  double ay = y - y0;
  double az = z - z0;

  x0 = clamp_int(x0, 0, nx - 1);
  y0 = clamp_int(y0, 0, ny - 1);
  z0 = clamp_int(z0, 0, nz - 1);
  x1 = clamp_int(x1, 0, nx - 1);
  y1 = clamp_int(y1, 0, ny - 1);
  z1 = clamp_int(z1, 0, nz - 1);

  size_t off = static_cast<size_t>(nxyz) * t;
  auto at = [&](int i, int j, int k) {
    return data[off + static_cast<size_t>(i + nx * (j + ny * k))];
  };

  double c000 = at(x0, y0, z0);
  double c100 = at(x1, y0, z0);
  double c010 = at(x0, y1, z0);
  double c110 = at(x1, y1, z0);
  double c001 = at(x0, y0, z1);
  double c101 = at(x1, y0, z1);
  double c011 = at(x0, y1, z1);
  double c111 = at(x1, y1, z1);

  double c00 = c000 * (1.0 - ax) + c100 * ax;
  double c10 = c010 * (1.0 - ax) + c110 * ax;
  double c01 = c001 * (1.0 - ax) + c101 * ax;
  double c11 = c011 * (1.0 - ax) + c111 * ax;

  double c0 = c00 * (1.0 - ay) + c10 * ay;
  double c1 = c01 * (1.0 - ay) + c11 * ay;

  return c0 * (1.0 - az) + c1 * az;
}

inline void fixed_to_moving_vox(int i,
                                int j,
                                int k,
                                int nx,
                                int ny,
                                int nz,
                                const NumericVector& pixdim,
                                const Pose6& pose,
                                double* src_xyz,
                                bool* in_bounds) {
  double cx = 0.5 * (nx - 1.0);
  double cy = 0.5 * (ny - 1.0);
  double cz = 0.5 * (nz - 1.0);

  double xw = (i - cx) * pixdim[0];
  double yw = (j - cy) * pixdim[1];
  double zw = (k - cz) * pixdim[2];

  double xf = xw - pose.tx;
  double yf = yw - pose.ty;
  double zf = zw - pose.tz;

  double R[9];
  rotation_matrix(pose.rx, pose.ry, pose.rz, R);

  // Inverse map by R^T
  double xm = R[0] * xf + R[3] * yf + R[6] * zf;
  double ym = R[1] * xf + R[4] * yf + R[7] * zf;
  double zm = R[2] * xf + R[5] * yf + R[8] * zf;

  double sx = xm / pixdim[0] + cx;
  double sy = ym / pixdim[1] + cy;
  double sz = zm / pixdim[2] + cz;

  src_xyz[0] = sx;
  src_xyz[1] = sy;
  src_xyz[2] = sz;

  *in_bounds = (sx >= 0.0 && sx <= (nx - 1.0) && sy >= 0.0 && sy <= (ny - 1.0) && sz >= 0.0 && sz <= (nz - 1.0));
}

CostResult evaluate_cost(const NumericVector& data,
                         const NumericVector& tmpl,
                         int nx,
                         int ny,
                         int nz,
                         int nxyz,
                         int t,
                         const NumericVector& pixdim,
                         const std::vector<int>& sample_points,
                         const Pose6& pose,
                         double huber_k) {
  double ss = 0.0;
  double n = 0.0;
  double inb = 0.0;
  double src[3];

  for (size_t s = 0; s < sample_points.size(); ++s) {
    int v = sample_points[s];
    int i = v % nx;
    int j = (v / nx) % ny;
    int k = v / (nx * ny);

    bool in_bounds = false;
    fixed_to_moving_vox(i, j, k, nx, ny, nz, pixdim, pose, src, &in_bounds);
    if (in_bounds) {
      inb += 1.0;
    }

    double mv = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
    double r = mv - tmpl[v];
    double ar = std::abs(r);

    double rho = 0.0;
    if (ar <= huber_k) {
      rho = r * r;
    } else {
      rho = 2.0 * huber_k * ar - huber_k * huber_k;
    }

    ss += rho;
    n += 1.0;
  }

  if (n <= 0.0) {
    return CostResult{0.0, 0.0};
  }
  return CostResult{ss / n, inb / n};
}

std::vector<int> make_sample_points(const LogicalVector& mask, int nxyz, int target) {
  std::vector<int> all;
  all.reserve(nxyz);
  for (int v = 0; v < nxyz; ++v) {
    if (mask[v] == TRUE) {
      all.push_back(v);
    }
  }
  if (all.empty()) {
    return all;
  }

  if (target <= 0 || static_cast<int>(all.size()) <= target) {
    return all;
  }

  std::vector<int> out;
  out.reserve(target);
  double step = static_cast<double>(all.size()) / static_cast<double>(target);
  for (int m = 0; m < target; ++m) {
    int idx = static_cast<int>(std::floor(m * step));
    idx = clamp_int(idx, 0, static_cast<int>(all.size()) - 1);
    out.push_back(all[idx]);
  }
  return out;
}

struct Vec3 {
  double x;
  double y;
  double z;
};

Vec3 center_of_mass_3d(const NumericVector& vol,
                       int nx,
                       int ny,
                       int nz,
                       const NumericVector& pixdim,
                       const LogicalVector& mask) {
  double sx = 0.0, sy = 0.0, sz = 0.0, wsum = 0.0;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int v = i + nx * (j + ny * k);
        if (mask[v] != TRUE) {
          continue;
        }
        double w = std::abs(vol[v]);
        sx += w * (i * pixdim[0]);
        sy += w * (j * pixdim[1]);
        sz += w * (k * pixdim[2]);
        wsum += w;
      }
    }
  }
  if (wsum <= 0.0) {
    return Vec3{0.5 * (nx - 1) * pixdim[0],
                0.5 * (ny - 1) * pixdim[1],
                0.5 * (nz - 1) * pixdim[2]};
  }
  return Vec3{sx / wsum, sy / wsum, sz / wsum};
}

Vec3 center_of_mass_frame(const NumericVector& data,
                          int nx,
                          int ny,
                          int nz,
                          int nxyz,
                          int t,
                          const NumericVector& pixdim,
                          const LogicalVector& mask) {
  double sx = 0.0, sy = 0.0, sz = 0.0, wsum = 0.0;
  size_t off = static_cast<size_t>(nxyz) * t;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int v = i + nx * (j + ny * k);
        if (mask[v] != TRUE) {
          continue;
        }
        double w = std::abs(data[off + v]);
        sx += w * (i * pixdim[0]);
        sy += w * (j * pixdim[1]);
        sz += w * (k * pixdim[2]);
        wsum += w;
      }
    }
  }
  if (wsum <= 0.0) {
    return Vec3{0.5 * (nx - 1) * pixdim[0],
                0.5 * (ny - 1) * pixdim[1],
                0.5 * (nz - 1) * pixdim[2]};
  }
  return Vec3{sx / wsum, sy / wsum, sz / wsum};
}

struct FrameFit {
  Pose6 pose;
  double cost_init;
  double cost_final;
  double overlap;
  int iters;
  bool restart_used;
  int flags;
};

FrameFit optimize_frame_pose(const NumericVector& data,
                             const NumericVector& tmpl,
                             int nx,
                             int ny,
                             int nz,
                             int nxyz,
                             int t,
                             const NumericVector& pixdim,
                             const std::vector<int>& samples,
                             const Pose6& init_pose,
                             int max_iter,
                             double lambda0,
                             double tol_step,
                             double tol_cost,
                             double huber_k,
                             int n_restarts) {
  std::vector<Pose6> seeds;
  seeds.push_back(init_pose);

  double rmm = 2.0;
  double rdeg = 2.0 * M_PI / 180.0;

  if (n_restarts > 0) {
    std::vector<Pose6> cand = {
        make_pose(init_pose.tx + rmm, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz),
        make_pose(init_pose.tx - rmm, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty + rmm, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty - rmm, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty, init_pose.tz + rmm, init_pose.rx, init_pose.ry, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty, init_pose.tz - rmm, init_pose.rx, init_pose.ry, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx + rdeg, init_pose.ry, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx - rdeg, init_pose.ry, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry + rdeg, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry - rdeg, init_pose.rz),
        make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz + rdeg),
        make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz - rdeg)};

    for (int i = 0; i < n_restarts && i < static_cast<int>(cand.size()); ++i) {
      seeds.push_back(cand[i]);
    }
  }

  Pose6 best_pose = init_pose;
  CostResult best_eval = evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, init_pose, huber_k);
  double best_c_init = best_eval.cost;
  double best_c_final = best_eval.cost;
  int best_iters = 0;
  bool best_restart = false;

  double h[6] = {
      std::max(0.2, 0.25 * pixdim[0]),
      std::max(0.2, 0.25 * pixdim[1]),
      std::max(0.2, 0.25 * pixdim[2]),
      0.01,
      0.01,
      0.01};

  for (size_t seed_idx = 0; seed_idx < seeds.size(); ++seed_idx) {
    Pose6 cur = seeds[seed_idx];
    double lambda = lambda0;
    CostResult eval0 = evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, cur, huber_k);
    double c_init = eval0.cost;
    double c_prev = c_init;
    int used_iters = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
      used_iters = iter + 1;
      double grad[6] = {0, 0, 0, 0, 0, 0};
      double hdiag[6] = {0, 0, 0, 0, 0, 0};

      auto set_param = [&](Pose6 p, int d, double v) {
        if (d == 0) p.tx = v;
        if (d == 1) p.ty = v;
        if (d == 2) p.tz = v;
        if (d == 3) p.rx = v;
        if (d == 4) p.ry = v;
        if (d == 5) p.rz = v;
        return p;
      };
      auto get_param = [&](const Pose6& p, int d) {
        if (d == 0) return p.tx;
        if (d == 1) return p.ty;
        if (d == 2) return p.tz;
        if (d == 3) return p.rx;
        if (d == 4) return p.ry;
        return p.rz;
      };

      for (int d = 0; d < 6; ++d) {
        double v = get_param(cur, d);
        Pose6 p_plus = set_param(cur, d, v + h[d]);
        Pose6 p_minus = set_param(cur, d, v - h[d]);

        double cp = evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, p_plus, huber_k).cost;
        double cm = evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, p_minus, huber_k).cost;

        grad[d] = (cp - cm) / (2.0 * h[d]);
        hdiag[d] = (cp - 2.0 * c_prev + cm) / (h[d] * h[d]);
        if (!R_finite(hdiag[d]) || hdiag[d] < 1e-6) {
          hdiag[d] = 1e-6;
        }
      }

      Pose6 delta = make_pose(
          -grad[0] / (hdiag[0] + lambda),
          -grad[1] / (hdiag[1] + lambda),
          -grad[2] / (hdiag[2] + lambda),
          -grad[3] / (hdiag[3] + lambda),
          -grad[4] / (hdiag[4] + lambda),
          -grad[5] / (hdiag[5] + lambda));

      double dnorm = std::sqrt(delta.tx * delta.tx + delta.ty * delta.ty + delta.tz * delta.tz +
                               delta.rx * delta.rx + delta.ry * delta.ry + delta.rz * delta.rz);

      bool accepted = false;
      double alpha = 1.0;
      double c_new = c_prev;

      for (int ls = 0; ls < 8; ++ls) {
        Pose6 cand = add_pose(cur, scale_pose(delta, alpha));
        c_new = evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, cand, huber_k).cost;
        if (c_new < c_prev) {
          cur = cand;
          accepted = true;
          break;
        }
        alpha *= 0.5;
      }

      if (accepted) {
        lambda *= 0.7;
      } else {
        lambda *= 2.0;
      }

      if (std::abs(c_prev - c_new) < tol_cost || dnorm < tol_step) {
        c_prev = c_new;
        break;
      }
      c_prev = c_new;
    }

    CostResult evalf = evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, cur, huber_k);

    if (evalf.cost < best_c_final) {
      best_pose = cur;
      best_c_init = c_init;
      best_c_final = evalf.cost;
      best_iters = used_iters;
      best_eval = evalf;
      best_restart = (seed_idx > 0);
    }
  }

  int flags = 0;
  if (best_eval.overlap < 0.85) {
    flags |= 1;
  }
  if (best_c_final > best_c_init * 0.999) {
    flags |= 2;
  }

  return FrameFit{best_pose,
                  best_c_init,
                  best_c_final,
                  best_eval.overlap,
                  best_iters,
                  best_restart,
                  flags};
}

NumericVector build_aligned_trimmed_template(const NumericVector& data,
                                             const std::vector<Pose6>& poses,
                                             int nx,
                                             int ny,
                                             int nz,
                                             int nt,
                                             int nxyz,
                                             const NumericVector& pixdim,
                                             double trim_frac) {
  NumericVector out(nxyz);
  std::vector<double> vals(nt);
  int trim_n = static_cast<int>(std::floor(trim_frac * nt));
  trim_n = std::max(0, std::min(trim_n, (nt - 1) / 2));

  double src[3];
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int v = i + nx * (j + ny * k);

        for (int t = 0; t < nt; ++t) {
          bool in_bounds = false;
          fixed_to_moving_vox(i, j, k, nx, ny, nz, pixdim, poses[t], src, &in_bounds);
          vals[t] = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
        }

        std::sort(vals.begin(), vals.end());
        int lo = trim_n;
        int hi = nt - trim_n;
        double s = 0.0;
        for (int tt = lo; tt < hi; ++tt) {
          s += vals[tt];
        }
        out[v] = s / static_cast<double>(hi - lo);
      }
    }
  }

  out.attr("dim") = IntegerVector::create(nx, ny, nz);
  return out;
}

void fill_matrix_from_pose(const Pose6& p, NumericMatrix& mats, int t) {
  double R[9];
  rotation_matrix(p.rx, p.ry, p.rz, R);

  mats(0, t) = R[0];
  mats(1, t) = R[3];
  mats(2, t) = R[6];
  mats(3, t) = 0.0;

  mats(4, t) = R[1];
  mats(5, t) = R[4];
  mats(6, t) = R[7];
  mats(7, t) = 0.0;

  mats(8, t) = R[2];
  mats(9, t) = R[5];
  mats(10, t) = R[8];
  mats(11, t) = 0.0;

  mats(12, t) = p.tx;
  mats(13, t) = p.ty;
  mats(14, t) = p.tz;
  mats(15, t) = 1.0;
}

int get_int_ctrl(const List& control, const std::string& name, int default_value) {
  if (!control.containsElementNamed(name.c_str())) {
    return default_value;
  }
  SEXP v = control[name];
  if (Rf_isInteger(v)) {
    IntegerVector iv(v);
    if (iv.size() > 0) {
      return iv[iv.size() - 1];
    }
  }
  if (Rf_isReal(v)) {
    NumericVector nv(v);
    if (nv.size() > 0) {
      return static_cast<int>(std::round(nv[nv.size() - 1]));
    }
  }
  return default_value;
}

double get_dbl_ctrl(const List& control, const std::string& name, double default_value) {
  if (!control.containsElementNamed(name.c_str())) {
    return default_value;
  }
  return as<double>(control[name]);
}

}  // namespace

// [[Rcpp::export]]
List cpp_volreg_estimate(
    NumericVector data,
    IntegerVector dim,
    NumericVector pixdim,
    Nullable<LogicalVector> mask = R_NilValue,
    Nullable<NumericVector> slice_times = R_NilValue,
    Nullable<IntegerVector> mb_groups = R_NilValue,
    int ref_index = 1,
    List control = List::create()) {
  if (dim.size() != 4) {
    stop("dim must be length 4: nx, ny, nz, nt");
  }
  if (pixdim.size() < 3) {
    stop("pixdim must have at least 3 elements");
  }

  int nx = dim[0];
  int ny = dim[1];
  int nz = dim[2];
  int nt = dim[3];
  int nxyz = nx * ny * nz;

  if (data.size() != static_cast<R_xlen_t>(nxyz) * nt) {
    stop("data size does not match dim");
  }
  if (ref_index < 1 || ref_index > nt) {
    stop("ref_index out of range");
  }

  bool robust_template = true;
  double trim_frac = 0.2;
  if (control.containsElementNamed("robust_template")) {
    robust_template = as<bool>(control["robust_template"]);
  }
  if (control.containsElementNamed("template_trim")) {
    trim_frac = as<double>(control["template_trim"]);
  }

  int max_iter = get_int_ctrl(control, "max_iter", 20);
  int n_samples = get_int_ctrl(control, "n_samples", 20000);
  int n_restarts = get_int_ctrl(control, "n_restarts", 0);
  double huber_k = get_dbl_ctrl(control, "huber_k", 1.5);
  double lambda0 = get_dbl_ctrl(control, "lambda0", 1e-2);
  double tol_step = get_dbl_ctrl(control, "tol_step", 1e-5);
  double tol_cost = get_dbl_ctrl(control, "tol_cost", 1e-6);

  NumericMatrix par(nt, 6);
  NumericMatrix mats = identity_mats(nt);

  NumericVector template_vol = extract_ref_volume(data, nx, ny, nz, ref_index);
  if (robust_template && nt > 1) {
    template_vol = trimmed_mean_template(data, nx, ny, nz, nt, trim_frac);
  }

  LogicalVector fit_mask = build_mask(mask, nxyz);
  std::vector<int> samples = make_sample_points(fit_mask, nxyz, n_samples);
  if (samples.empty()) {
    stop("mask produced zero fitting samples");
  }

  NumericVector cost_init(nt, NA_REAL);
  NumericVector cost_final(nt, NA_REAL);
  IntegerVector iters(nt, 0);
  LogicalVector restart_used(nt, false);
  NumericVector overlap(nt, NA_REAL);
  IntegerVector flags(nt, 0);

  std::vector<Pose6> poses(nt, make_pose());

  auto store_fit = [&](int t, const FrameFit& fit) {
    poses[t] = fit.pose;
    par(t, 0) = fit.pose.tx;
    par(t, 1) = fit.pose.ty;
    par(t, 2) = fit.pose.tz;
    par(t, 3) = fit.pose.rx;
    par(t, 4) = fit.pose.ry;
    par(t, 5) = fit.pose.rz;
    fill_matrix_from_pose(fit.pose, mats, t);

    cost_init[t] = fit.cost_init;
    cost_final[t] = fit.cost_final;
    iters[t] = fit.iters;
    restart_used[t] = fit.restart_used;
    overlap[t] = fit.overlap;
    flags[t] = fit.flags;
  };

  int ref_t = ref_index - 1;
  Vec3 com_ref = center_of_mass_3d(template_vol, nx, ny, nz, pixdim, fit_mask);
  Vec3 com_init_ref = center_of_mass_frame(data, nx, ny, nz, nxyz, ref_t, pixdim, fit_mask);
  Pose6 ref_init = make_pose(com_ref.x - com_init_ref.x,
                             com_ref.y - com_init_ref.y,
                             com_ref.z - com_init_ref.z,
                             0.0,
                             0.0,
                             0.0);
  auto solve_sequence = [&](const Pose6& ref_seed, bool blend_current_pose_seed) {
    store_fit(ref_t,
              optimize_frame_pose(data,
                                  template_vol,
                                  nx,
                                  ny,
                                  nz,
                                  nxyz,
                                  ref_t,
                                  pixdim,
                                  samples,
                                  ref_seed,
                                  max_iter,
                                  lambda0,
                                  tol_step,
                                  tol_cost,
                                  huber_k,
                                  n_restarts));

    // Backward from reference
    for (int t = ref_t - 1; t >= 0; --t) {
      Pose6 init = poses[t + 1];
      if (t + 2 <= ref_t) {
        init = add_pose(poses[t + 1], sub_pose(poses[t + 1], poses[t + 2]));
      }
      if (blend_current_pose_seed) {
        init = scale_pose(add_pose(init, poses[t]), 0.5);
      }

      store_fit(t,
                optimize_frame_pose(data,
                                    template_vol,
                                    nx,
                                    ny,
                                    nz,
                                    nxyz,
                                    t,
                                    pixdim,
                                    samples,
                                    init,
                                    max_iter,
                                    lambda0,
                                    tol_step,
                                    tol_cost,
                                    huber_k,
                                    n_restarts));
    }

    // Forward from reference
    for (int t = ref_t + 1; t < nt; ++t) {
      Pose6 init = poses[t - 1];
      if (t - 2 >= ref_t) {
        init = add_pose(poses[t - 1], sub_pose(poses[t - 1], poses[t - 2]));
      }
      if (blend_current_pose_seed) {
        init = scale_pose(add_pose(init, poses[t]), 0.5);
      }

      store_fit(t,
                optimize_frame_pose(data,
                                    template_vol,
                                    nx,
                                    ny,
                                    nz,
                                    nxyz,
                                    t,
                                    pixdim,
                                    samples,
                                    init,
                                    max_iter,
                                    lambda0,
                                    tol_step,
                                    tol_cost,
                                    huber_k,
                                    n_restarts));
    }
  };

  solve_sequence(ref_init, false);

  bool template_refresh = robust_template && nt > 1;
  if (control.containsElementNamed("template_refresh")) {
    template_refresh = as<bool>(control["template_refresh"]);
  }

  if (template_refresh) {
    template_vol = build_aligned_trimmed_template(
        data, poses, nx, ny, nz, nt, nxyz, pixdim, trim_frac);
    solve_sequence(poses[ref_t], true);
  }

  // Placeholder: optional inputs are intentionally accepted to stabilize API.
  (void) slice_times;
  (void) mb_groups;

  fit_mask.attr("dim") = IntegerVector::create(nx, ny, nz);

  return List::create(_["par"] = par,
                      _["mats"] = mats,
                      _["template"] = template_vol,
                      _["mask"] = fit_mask,
                      _["cost_init"] = cost_init,
                      _["cost_final"] = cost_final,
                      _["iters"] = iters,
                      _["restart_used"] = restart_used,
                      _["overlap"] = overlap,
                      _["flags"] = flags);
}
