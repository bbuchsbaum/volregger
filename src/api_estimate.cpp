#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <thread>
#include <vector>

using namespace Rcpp;

namespace {

inline size_t idx4(int i, int j, int k, int t, int nx, int ny, int nz) {
  return static_cast<size_t>(i + nx * (j + ny * (k + nz * t)));
}

inline int clamp_int(int x, int lo, int hi) {
  return std::max(lo, std::min(hi, x));
}

inline int effective_thread_count(int requested, size_t work_items) {
  if (requested <= 1 || work_items <= 1) {
    return 1;
  }
  unsigned int hw = std::thread::hardware_concurrency();
  int max_threads = hw > 0 ? static_cast<int>(hw) : requested;
  int out = std::max(1, std::min(requested, max_threads));
  if (work_items < static_cast<size_t>(2 * out)) {
    return 1;
  }
  return out;
}

template <typename Fn>
void parallel_for_range(size_t n, int nthreads, Fn fn) {
  int workers = effective_thread_count(nthreads, n);
  if (workers <= 1) {
    fn(0, n);
    return;
  }

  size_t chunk = (n + static_cast<size_t>(workers) - 1) / static_cast<size_t>(workers);
  std::vector<std::thread> pool;
  pool.reserve(static_cast<size_t>(workers));

  for (int tid = 0; tid < workers; ++tid) {
    size_t begin = static_cast<size_t>(tid) * chunk;
    if (begin >= n) {
      break;
    }
    size_t end = std::min(n, begin + chunk);
    pool.emplace_back([&, begin, end]() { fn(begin, end); });
  }

  for (size_t i = 0; i < pool.size(); ++i) {
    pool[i].join();
  }
}

struct Pose6 {
  double tx;
  double ty;
  double tz;
  double rx;
  double ry;
  double rz;
};

double wrap_pi(double a) {
  const double two_pi = 2.0 * M_PI;
  while (a > M_PI) {
    a -= two_pi;
  }
  while (a < -M_PI) {
    a += two_pi;
  }
  return a;
}

Pose6 make_pose(double tx = 0.0,
                double ty = 0.0,
                double tz = 0.0,
                double rx = 0.0,
                double ry = 0.0,
                double rz = 0.0) {
  return Pose6{tx, ty, tz, wrap_pi(rx), wrap_pi(ry), wrap_pi(rz)};
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
                                    double trim_frac,
                                    int nthreads = 1) {
  int nxyz = nx * ny * nz;
  NumericVector out(nxyz);
  int trim_n = static_cast<int>(std::floor(trim_frac * nt));
  trim_n = std::max(0, std::min(trim_n, (nt - 1) / 2));

  parallel_for_range(static_cast<size_t>(nxyz), nthreads, [&](size_t begin, size_t end) {
    std::vector<double> vals(static_cast<size_t>(nt));
    for (size_t vv = begin; vv < end; ++vv) {
      int v = static_cast<int>(vv);
      for (int t = 0; t < nt; ++t) {
        vals[static_cast<size_t>(t)] = data[static_cast<size_t>(v + nxyz * t)];
      }
      std::sort(vals.begin(), vals.end());
      int lo = trim_n;
      int hi = nt - trim_n;
      double s = 0.0;
      for (int i = lo; i < hi; ++i) {
        s += vals[static_cast<size_t>(i)];
      }
      out[v] = s / static_cast<double>(hi - lo);
    }
  });

  out.attr("dim") = IntegerVector::create(nx, ny, nz);
  return out;
}

double quantile_inplace(std::vector<double>& x, double q) {
  if (x.empty()) {
    return 0.0;
  }
  q = std::max(0.0, std::min(1.0, q));
  size_t k = static_cast<size_t>(std::floor(q * static_cast<double>(x.size() - 1)));
  std::nth_element(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(k), x.end());
  return x[k];
}

double robust_sigma_from_residuals(std::vector<double>& residuals) {
  if (residuals.empty()) {
    return 1.0;
  }

  std::vector<double> med_vals = residuals;
  double med = quantile_inplace(med_vals, 0.5);

  std::vector<double> absdev;
  absdev.reserve(residuals.size());
  std::vector<double> absres;
  absres.reserve(residuals.size());
  for (size_t i = 0; i < residuals.size(); ++i) {
    absdev.push_back(std::abs(residuals[i] - med));
    absres.push_back(std::abs(residuals[i]));
  }

  double mad = quantile_inplace(absdev, 0.5);
  double sigma = 1.4826 * mad;
  if (!R_finite(sigma) || sigma <= 1e-8) {
    double mar = quantile_inplace(absres, 0.5);
    sigma = 1.4826 * mar;
  }
  if (!R_finite(sigma) || sigma <= 1e-8) {
    sigma = 1.0;
  }
  return sigma;
}

std::string get_str_ctrl(const List& control,
                         const std::string& name,
                         const std::string& default_value) {
  if (!control.containsElementNamed(name.c_str())) {
    return default_value;
  }
  SEXP v = control[name];
  if (TYPEOF(v) == STRSXP && Rf_length(v) > 0) {
    CharacterVector cv(v);
    return as<std::string>(cv[cv.size() - 1]);
  }
  return default_value;
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

int select_reference_frame_from_template(const NumericVector& data,
                                         const NumericVector& tmpl,
                                         int nt,
                                         int nxyz,
                                         const LogicalVector& mask) {
  double sum_t = 0.0;
  double sum_t2 = 0.0;
  int nmask = 0;
  for (int v = 0; v < nxyz; ++v) {
    if (!mask[v]) {
      continue;
    }
    double tv = tmpl[v];
    if (!R_finite(tv)) {
      continue;
    }
    sum_t += tv;
    sum_t2 += tv * tv;
    ++nmask;
  }
  if (nmask == 0) {
    return std::max(0, (nt - 1) / 2);
  }

  double mean_t = sum_t / static_cast<double>(nmask);
  double var_t = std::max(1e-12, sum_t2 / static_cast<double>(nmask) - mean_t * mean_t);

  int best_t = 0;
  double best_score = -std::numeric_limits<double>::infinity();

  for (int t = 0; t < nt; ++t) {
    double sum_x = 0.0;
    double sum_x2 = 0.0;
    double sum_tx = 0.0;
    int n = 0;
    size_t offset = static_cast<size_t>(nxyz) * static_cast<size_t>(t);
    for (int v = 0; v < nxyz; ++v) {
      if (!mask[v]) {
        continue;
      }
      double tv = tmpl[v];
      double xv = data[offset + static_cast<size_t>(v)];
      if (!R_finite(tv) || !R_finite(xv)) {
        continue;
      }
      sum_x += xv;
      sum_x2 += xv * xv;
      sum_tx += tv * xv;
      ++n;
    }
    if (n == 0) {
      continue;
    }
    double mean_x = sum_x / static_cast<double>(n);
    double var_x = sum_x2 / static_cast<double>(n) - mean_x * mean_x;
    if (var_x <= 1e-12) {
      continue;
    }
    double cov = sum_tx / static_cast<double>(n) - mean_t * mean_x;
    double corr = cov / std::sqrt(var_t * var_x);
    if (corr > best_score) {
      best_score = corr;
      best_t = t;
    }
  }

  return best_t;
}

LogicalVector auto_mask_from_template(const NumericVector& tmpl, int nxyz) {
  std::vector<double> vals;
  vals.reserve(nxyz);
  for (int v = 0; v < nxyz; ++v) {
    double tv = tmpl[v];
    if (R_finite(tv)) {
      vals.push_back(tv);
    }
  }
  if (vals.empty()) {
    return LogicalVector(nxyz, true);
  }

  std::vector<double> vals_med = vals;
  double med = quantile_inplace(vals_med, 0.5);

  std::vector<double> absdev;
  absdev.reserve(vals.size());
  for (size_t i = 0; i < vals.size(); ++i) {
    absdev.push_back(std::abs(vals[i] - med));
  }
  double thr = quantile_inplace(absdev, 0.7);
  if (thr <= 0.0) {
    thr = quantile_inplace(absdev, 0.5);
  }

  LogicalVector m(nxyz, false);
  int kept = 0;
  for (int v = 0; v < nxyz; ++v) {
    bool keep = R_finite(tmpl[v]) && std::abs(tmpl[v] - med) >= thr;
    m[v] = keep;
    if (keep) {
      kept++;
    }
  }

  double frac = static_cast<double>(kept) / static_cast<double>(nxyz);
  if (frac < 0.05 || frac > 0.95) {
    return LogicalVector(nxyz, true);
  }
  return m;
}

LogicalVector exclude_mask_edges(const LogicalVector& mask,
                                 int nx,
                                 int ny,
                                 int nz,
                                 double edge_frac,
                                 int edge_vox) {
  int ex = std::max(edge_vox, static_cast<int>(std::round(edge_frac * nx)));
  int ey = std::max(edge_vox, static_cast<int>(std::round(edge_frac * ny)));
  int ez = std::max(edge_vox, static_cast<int>(std::round(edge_frac * nz)));
  ex = clamp_int(ex, 0, std::max(0, (nx - 1) / 2));
  ey = clamp_int(ey, 0, std::max(0, (ny - 1) / 2));
  ez = clamp_int(ez, 0, std::max(0, (nz - 1) / 2));
  if (ex == 0 && ey == 0 && ez == 0) {
    return mask;
  }

  LogicalVector out = clone(mask);
  int kept = 0;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int v = i + nx * (j + ny * k);
        if (i < ex || i >= nx - ex || j < ey || j >= ny - ey || k < ez || k >= nz - ez) {
          out[v] = false;
        }
        if (out[v] == TRUE) {
          ++kept;
        }
      }
    }
  }

  if (kept == 0) {
    stop("edge exclusion removed all fitting voxels");
  }
  return out;
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

bool solve_linear6(double* A, double* b, double* x);
bool solve_linear_small(double* A, double* b, int n, double* x);

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
                         double huber_k,
                         bool use_overlap_taper) {
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
    double mv = 0.0;
    if (in_bounds) {
      mv = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
      inb += 1.0;
    }
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
    return CostResult{1e12, 0.0};
  }
  double overlap = inb / n;
  double c = ss / n;
  if (use_overlap_taper) {
    double miss = 1.0 - overlap;
    c *= (1.0 + 4.0 * miss * miss);
  }
  if (overlap <= 0.0) {
    c = 1e12;
  }
  return CostResult{c, overlap};
}

double estimate_residual_sigma_dense(const NumericVector& data,
                                     const NumericVector& tmpl,
                                     int nx,
                                     int ny,
                                     int nz,
                                     int nxyz,
                                     int t,
                                     const NumericVector& pixdim,
                                     const std::vector<int>& sample_points,
                                     const Pose6& pose) {
  std::vector<double> residuals;
  residuals.reserve(sample_points.size());
  double src[3];

  for (size_t s = 0; s < sample_points.size(); ++s) {
    int v = sample_points[s];
    int i = v % nx;
    int j = (v / nx) % ny;
    int k = v / (nx * ny);

    bool in_bounds = false;
    fixed_to_moving_vox(i, j, k, nx, ny, nz, pixdim, pose, src, &in_bounds);
    if (!in_bounds) {
      continue;
    }
    double mv = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
    residuals.push_back(mv - tmpl[v]);
  }

  return robust_sigma_from_residuals(residuals);
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

LogicalVector make_downsampled_mask(const LogicalVector& mask,
                                    int nx,
                                    int ny,
                                    int nz,
                                    int downsample) {
  if (downsample <= 1) {
    return clone(mask);
  }
  LogicalVector out(mask.size(), false);
  int kept = 0;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int v = i + nx * (j + ny * k);
        bool keep = mask[v] == TRUE &&
                    (i % downsample == 0) &&
                    (j % downsample == 0) &&
                    (k % downsample == 0);
        out[v] = keep;
        if (keep) {
          ++kept;
        }
      }
    }
  }
  if (kept == 0) {
    return clone(mask);
  }
  return out;
}

struct StencilPoint {
  int v;
  int i;
  int j;
  int k;
  double tval;
  double jcol[6];
  double jwhite[6];
  double bvec[5];
  double w0;
};

struct WhitenMetric {
  bool ok;
  double L[36];
  double Linv[36];
  double LinvT[36];
};

void set_identity6(double* M) {
  for (int i = 0; i < 36; ++i) {
    M[i] = 0.0;
  }
  for (int i = 0; i < 6; ++i) {
    M[i * 6 + i] = 1.0;
  }
}

bool chol6_lower(const double* A, double* L) {
  for (int i = 0; i < 36; ++i) {
    L[i] = 0.0;
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = A[i * 6 + j];
      for (int k = 0; k < j; ++k) {
        s -= L[i * 6 + k] * L[j * 6 + k];
      }
      if (i == j) {
        if (s <= 1e-12 || !R_finite(s)) {
          return false;
        }
        L[i * 6 + j] = std::sqrt(s);
      } else {
        L[i * 6 + j] = s / L[j * 6 + j];
      }
    }
  }
  return true;
}

bool invert_lower6(const double* L, double* Linv) {
  for (int i = 0; i < 36; ++i) {
    Linv[i] = 0.0;
  }
  for (int c = 0; c < 6; ++c) {
    for (int i = 0; i < 6; ++i) {
      double s = (i == c) ? 1.0 : 0.0;
      for (int k = 0; k < i; ++k) {
        s -= L[i * 6 + k] * Linv[k * 6 + c];
      }
      double d = L[i * 6 + i];
      if (std::abs(d) < 1e-12 || !R_finite(d)) {
        return false;
      }
      Linv[i * 6 + c] = s / d;
    }
  }
  return true;
}

WhitenMetric build_whiten_metric(const std::vector<StencilPoint>& stencil, double reg_scale) {
  WhitenMetric wm;
  wm.ok = false;
  set_identity6(wm.L);
  set_identity6(wm.Linv);
  set_identity6(wm.LinvT);

  if (stencil.empty()) {
    return wm;
  }

  double H[36] = {0};
  for (size_t s = 0; s < stencil.size(); ++s) {
    const StencilPoint& p = stencil[s];
    for (int a = 0; a < 6; ++a) {
      for (int b = a; b < 6; ++b) {
        H[a * 6 + b] += p.w0 * p.jcol[a] * p.jcol[b];
      }
    }
  }
  for (int a = 0; a < 6; ++a) {
    for (int b = a + 1; b < 6; ++b) {
      H[b * 6 + a] = H[a * 6 + b];
    }
  }

  double tr = 0.0;
  for (int a = 0; a < 6; ++a) {
    tr += H[a * 6 + a];
  }
  double eps = std::max(1e-8, std::max(0.0, reg_scale) * std::max(1.0, tr / 6.0));
  for (int a = 0; a < 6; ++a) {
    H[a * 6 + a] += eps;
  }

  if (!chol6_lower(H, wm.L)) {
    return wm;
  }
  if (!invert_lower6(wm.L, wm.Linv)) {
    return wm;
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      wm.LinvT[i * 6 + j] = wm.Linv[j * 6 + i];
    }
  }
  wm.ok = true;
  return wm;
}

void apply_whiten_metric_to_stencil(std::vector<StencilPoint>* stencil, const WhitenMetric& wm) {
  if (stencil == nullptr) {
    return;
  }
  for (size_t s = 0; s < stencil->size(); ++s) {
    StencilPoint& p = (*stencil)[s];
    for (int a = 0; a < 6; ++a) {
      double v = 0.0;
      for (int b = 0; b < 6; ++b) {
        v += p.jcol[b] * wm.LinvT[b * 6 + a];
      }
      p.jwhite[a] = v;
    }
  }
}

void copy_raw_jacobian_to_whitened(std::vector<StencilPoint>* stencil) {
  if (stencil == nullptr) {
    return;
  }
  for (size_t s = 0; s < stencil->size(); ++s) {
    for (int a = 0; a < 6; ++a) {
      (*stencil)[s].jwhite[a] = (*stencil)[s].jcol[a];
    }
  }
}

inline double template_at_vox(const NumericVector& tmpl,
                              int nx,
                              int ny,
                              int nz,
                              int i,
                              int j,
                              int k) {
  i = clamp_int(i, 0, nx - 1);
  j = clamp_int(j, 0, ny - 1);
  k = clamp_int(k, 0, nz - 1);
  return tmpl[i + nx * (j + ny * k)];
}

void template_gradient_vox(const NumericVector& tmpl,
                           int nx,
                           int ny,
                           int nz,
                           int i,
                           int j,
                           int k,
                           double* gx,
                           double* gy,
                           double* gz) {
  *gx = 0.5 * (template_at_vox(tmpl, nx, ny, nz, i + 1, j, k) -
               template_at_vox(tmpl, nx, ny, nz, i - 1, j, k));
  *gy = 0.5 * (template_at_vox(tmpl, nx, ny, nz, i, j + 1, k) -
               template_at_vox(tmpl, nx, ny, nz, i, j - 1, k));
  *gz = 0.5 * (template_at_vox(tmpl, nx, ny, nz, i, j, k + 1) -
               template_at_vox(tmpl, nx, ny, nz, i, j, k - 1));
}

std::vector<double> temporal_reliability_map(const NumericVector& data, int nxyz, int nt) {
  std::vector<double> rel(nxyz, 1.0);
  if (nt <= 1) {
    return rel;
  }

  for (int v = 0; v < nxyz; ++v) {
    double m = 0.0;
    for (int t = 0; t < nt; ++t) {
      m += data[static_cast<size_t>(v + nxyz * t)];
    }
    m /= static_cast<double>(nt);

    double ss = 0.0;
    for (int t = 0; t < nt; ++t) {
      double d = data[static_cast<size_t>(v + nxyz * t)] - m;
      ss += d * d;
    }
    double sd = std::sqrt(ss / std::max(1.0, static_cast<double>(nt - 1)));
    rel[v] = 1.0 / (1.0 + sd);
  }

  std::vector<double> rel_copy = rel;
  double hi = quantile_inplace(rel_copy, 0.95);
  if (hi <= 0.0 || !R_finite(hi)) {
    return rel;
  }
  for (int v = 0; v < nxyz; ++v) {
    double r = rel[v] / hi;
    rel[v] = std::max(0.05, std::min(1.0, r));
  }
  return rel;
}

std::vector<StencilPoint> build_ic_stencil(const NumericVector& tmpl,
                                           const NumericVector& data,
                                           int nx,
                                           int ny,
                                           int nz,
                                           int nt,
                                           int nxyz,
                                           const NumericVector& pixdim,
                                           const LogicalVector& mask,
                                           int target_points,
                                           int blocks_x,
                                           int blocks_y,
                                           int blocks_z,
                                           double gamma) {
  struct Candidate {
    int v;
    int i;
    int j;
    int k;
    int bid;
    double gx;
    double gy;
    double gz;
    double score;
  };

  target_points = std::max(256, target_points);
  blocks_x = std::max(1, blocks_x);
  blocks_y = std::max(1, blocks_y);
  blocks_z = std::max(1, blocks_z);
  gamma = std::max(0.0, gamma);

  std::vector<double> rel = temporal_reliability_map(data, nxyz, nt);
  std::vector<Candidate> candidates;
  candidates.reserve(nxyz / 4);

  double cx = 0.5 * (nx - 1.0);
  double cy = 0.5 * (ny - 1.0);
  double cz = 0.5 * (nz - 1.0);

  int nblocks = blocks_x * blocks_y * blocks_z;
  std::vector<std::vector<int>> bins(nblocks);

  for (int k = 1; k < nz - 1; ++k) {
    for (int j = 1; j < ny - 1; ++j) {
      for (int i = 1; i < nx - 1; ++i) {
        int v = i + nx * (j + ny * k);
        if (mask[v] != TRUE) {
          continue;
        }

        double gx = 0.0, gy = 0.0, gz = 0.0;
        template_gradient_vox(tmpl, nx, ny, nz, i, j, k, &gx, &gy, &gz);
        double gmag = std::sqrt(gx * gx + gy * gy + gz * gz);
        if (gmag < 1e-6 || !R_finite(gmag)) {
          continue;
        }

        double xw = (i - cx) * pixdim[0];
        double yw = (j - cy) * pixdim[1];
        double zw = (k - cz) * pixdim[2];

        double c1 = yw * gz - zw * gy;
        double c2 = zw * gx - xw * gz;
        double c3 = xw * gy - yw * gx;
        double leverage = std::sqrt(c1 * c1 + c2 * c2 + c3 * c3) + 0.25 * gmag;

        double score = std::pow(gmag, gamma) * rel[v] * leverage;
        if (!R_finite(score) || score <= 0.0) {
          continue;
        }

        int bi = std::min(blocks_x - 1, (i * blocks_x) / nx);
        int bj = std::min(blocks_y - 1, (j * blocks_y) / ny);
        int bk = std::min(blocks_z - 1, (k * blocks_z) / nz);
        int bid = bi + blocks_x * (bj + blocks_y * bk);

        int id = static_cast<int>(candidates.size());
        candidates.push_back(Candidate{v, i, j, k, bid, gx, gy, gz, score});
        bins[bid].push_back(id);
      }
    }
  }

  if (candidates.empty()) {
    return {};
  }

  auto by_score = [&](int a, int b) { return candidates[a].score > candidates[b].score; };
  for (int b = 0; b < nblocks; ++b) {
    std::sort(bins[b].begin(), bins[b].end(), by_score);
  }

  std::vector<int> selected;
  selected.reserve(target_points);
  std::vector<int> cursor(nblocks, 0);
  int per_block = std::max(1, target_points / nblocks);

  for (int b = 0; b < nblocks; ++b) {
    int take = std::min(per_block, static_cast<int>(bins[b].size()));
    for (int n = 0; n < take; ++n) {
      selected.push_back(bins[b][n]);
    }
    cursor[b] = take;
  }

  if (static_cast<int>(selected.size()) < target_points) {
    std::vector<int> rest;
    for (int b = 0; b < nblocks; ++b) {
      for (int n = cursor[b]; n < static_cast<int>(bins[b].size()); ++n) {
        rest.push_back(bins[b][n]);
      }
    }
    std::sort(rest.begin(), rest.end(), by_score);
    int need = target_points - static_cast<int>(selected.size());
    for (int n = 0; n < need && n < static_cast<int>(rest.size()); ++n) {
      selected.push_back(rest[n]);
    }
  }

  if (selected.empty()) {
    return {};
  }

  if (static_cast<int>(selected.size()) > target_points) {
    std::sort(selected.begin(), selected.end(), by_score);
    selected.resize(target_points);
  }

  double max_score = 0.0;
  for (size_t s = 0; s < selected.size(); ++s) {
    max_score = std::max(max_score, candidates[selected[s]].score);
  }
  if (max_score <= 0.0 || !R_finite(max_score)) {
    max_score = 1.0;
  }

  std::vector<StencilPoint> stencil;
  stencil.reserve(selected.size());

  double invdx = 1.0 / pixdim[0];
  double invdy = 1.0 / pixdim[1];
  double invdz = 1.0 / pixdim[2];

  for (size_t s = 0; s < selected.size(); ++s) {
    const Candidate& c = candidates[selected[s]];
    double xw = (c.i - cx) * pixdim[0];
    double yw = (c.j - cy) * pixdim[1];
    double zw = (c.k - cz) * pixdim[2];

    StencilPoint p;
    p.v = c.v;
    p.i = c.i;
    p.j = c.j;
    p.k = c.k;
    p.tval = tmpl[c.v];
    p.w0 = std::max(0.05, std::min(1.0, c.score / max_score));

    // Translation terms.
    p.jcol[0] = c.gx * (-invdx);
    p.jcol[1] = c.gy * (-invdy);
    p.jcol[2] = c.gz * (-invdz);
    // Rotation terms at identity using d(R^T x)/dtheta = -(e x x).
    p.jcol[3] = c.gy * (zw * invdy) + c.gz * (-yw * invdz);
    p.jcol[4] = c.gx * (-zw * invdx) + c.gz * (xw * invdz);
    p.jcol[5] = c.gx * (yw * invdx) + c.gy * (-xw * invdy);
    for (int a = 0; a < 6; ++a) {
      p.jwhite[a] = p.jcol[a];
    }
    for (int q = 0; q < 5; ++q) {
      p.bvec[q] = 0.0;
    }

    stencil.push_back(p);
  }

  return stencil;
}

CostResult evaluate_cost_stencil(const NumericVector& data,
                                 int nx,
                                 int ny,
                                 int nz,
                                 int nxyz,
                                 int t,
                                 const NumericVector& pixdim,
                                 const std::vector<StencilPoint>& stencil,
                                 const Pose6& pose,
                                 double huber_k,
                                 bool use_overlap_taper,
                                 int nthreads = 1) {
  int workers = effective_thread_count(nthreads, stencil.size());
  std::vector<double> ss_local(static_cast<size_t>(workers), 0.0);
  std::vector<double> wsum_local(static_cast<size_t>(workers), 0.0);
  std::vector<double> inb_local(static_cast<size_t>(workers), 0.0);

  parallel_for_range(stencil.size(), workers, [&](size_t begin, size_t end) {
    size_t tid = begin == end ? 0 : begin / std::max<size_t>(1, (stencil.size() + static_cast<size_t>(workers) - 1) / static_cast<size_t>(workers));
    if (tid >= static_cast<size_t>(workers)) {
      tid = static_cast<size_t>(workers - 1);
    }
    double ss = 0.0;
    double wsum = 0.0;
    double inb = 0.0;
    double src[3];

    for (size_t s = begin; s < end; ++s) {
      const StencilPoint& p = stencil[s];
      bool in_bounds = false;
      fixed_to_moving_vox(p.i, p.j, p.k, nx, ny, nz, pixdim, pose, src, &in_bounds);

      double mv = 0.0;
      if (in_bounds) {
        mv = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
        inb += p.w0;
      }

      double r = mv - p.tval;
      double ar = std::abs(r);
      double rho = (ar <= huber_k) ? (r * r) : (2.0 * huber_k * ar - huber_k * huber_k);

      ss += p.w0 * rho;
      wsum += p.w0;
    }

    ss_local[tid] = ss;
    wsum_local[tid] = wsum;
    inb_local[tid] = inb;
  });

  double ss = 0.0;
  double wsum = 0.0;
  double inb = 0.0;
  for (int tid = 0; tid < workers; ++tid) {
    ss += ss_local[static_cast<size_t>(tid)];
    wsum += wsum_local[static_cast<size_t>(tid)];
    inb += inb_local[static_cast<size_t>(tid)];
  }

  if (wsum <= 0.0) {
    return CostResult{1e12, 0.0};
  }
  double overlap = inb / wsum;
  double c = ss / wsum;
  if (use_overlap_taper) {
    double miss = 1.0 - overlap;
    c *= (1.0 + 4.0 * miss * miss);
  }
  if (overlap <= 0.0) {
    c = 1e12;
  }
  return CostResult{c, overlap};
}

double estimate_residual_sigma_stencil(const NumericVector& data,
                                       int nx,
                                       int ny,
                                       int nz,
                                       int nxyz,
                                       int t,
                                       const NumericVector& pixdim,
                                       const std::vector<StencilPoint>& stencil,
                                       const Pose6& pose) {
  std::vector<double> residuals;
  residuals.reserve(stencil.size());
  double src[3];

  for (size_t s = 0; s < stencil.size(); ++s) {
    const StencilPoint& p = stencil[s];
    bool in_bounds = false;
    fixed_to_moving_vox(p.i, p.j, p.k, nx, ny, nz, pixdim, pose, src, &in_bounds);
    if (!in_bounds) {
      continue;
    }
    double mv = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
    residuals.push_back(mv - p.tval);
  }

  return robust_sigma_from_residuals(residuals);
}

int prepare_nuisance_basis(std::vector<StencilPoint>* stencil,
                           int nx,
                           int ny,
                           int nz,
                           bool include_template,
                           bool orthogonalize,
                           bool use_whiten_jacobian,
                           double orth_reg) {
  if (stencil == nullptr || stencil->empty()) {
    return 0;
  }

  int p = include_template ? 5 : 4;  // [1, tmpl?, x, y, z]
  double cx = 0.5 * (nx - 1.0);
  double cy = 0.5 * (ny - 1.0);
  double cz = 0.5 * (nz - 1.0);
  double sx = std::max(1.0, cx);
  double sy = std::max(1.0, cy);
  double sz = std::max(1.0, cz);

  for (size_t s = 0; s < stencil->size(); ++s) {
    StencilPoint& pt = (*stencil)[s];
    int q = 0;
    pt.bvec[q++] = 1.0;
    if (include_template) {
      pt.bvec[q++] = pt.tval;
    }
    pt.bvec[q++] = (pt.i - cx) / sx;
    pt.bvec[q++] = (pt.j - cy) / sy;
    pt.bvec[q++] = (pt.k - cz) / sz;
    for (; q < 5; ++q) {
      pt.bvec[q] = 0.0;
    }
  }

  // Weighted center/scale for non-constant columns.
  for (int q = 1; q < p; ++q) {
    double wsum = 0.0;
    double mean = 0.0;
    for (size_t s = 0; s < stencil->size(); ++s) {
      const StencilPoint& pt = (*stencil)[s];
      wsum += pt.w0;
      mean += pt.w0 * pt.bvec[q];
    }
    if (wsum <= 0.0) {
      continue;
    }
    mean /= wsum;
    double var = 0.0;
    for (size_t s = 0; s < stencil->size(); ++s) {
      const StencilPoint& pt = (*stencil)[s];
      double d = pt.bvec[q] - mean;
      var += pt.w0 * d * d;
    }
    var /= wsum;
    double sd = std::sqrt(std::max(var, 1e-8));
    for (size_t s = 0; s < stencil->size(); ++s) {
      (*stencil)[s].bvec[q] = ((*stencil)[s].bvec[q] - mean) / sd;
    }
  }

  if (orthogonalize && p > 1) {
    double H[36] = {0};
    for (size_t s = 0; s < stencil->size(); ++s) {
      const StencilPoint& pt = (*stencil)[s];
      const double* jac = use_whiten_jacobian ? pt.jwhite : pt.jcol;
      for (int a = 0; a < 6; ++a) {
        for (int b = a; b < 6; ++b) {
          H[a * 6 + b] += pt.w0 * jac[a] * jac[b];
        }
      }
    }
    for (int a = 0; a < 6; ++a) {
      for (int b = a + 1; b < 6; ++b) {
        H[b * 6 + a] = H[a * 6 + b];
      }
    }

    double tr = 0.0;
    for (int a = 0; a < 6; ++a) {
      tr += H[a * 6 + a];
    }
    double reg = std::max(1e-8, std::max(0.0, orth_reg) * std::max(1.0, tr / 6.0));
    for (int a = 0; a < 6; ++a) {
      H[a * 6 + a] += reg;
    }

    for (int q = 1; q < p; ++q) {
      double h[6] = {0};
      for (size_t s = 0; s < stencil->size(); ++s) {
        const StencilPoint& pt = (*stencil)[s];
        const double* jac = use_whiten_jacobian ? pt.jwhite : pt.jcol;
        for (int a = 0; a < 6; ++a) {
          h[a] += pt.w0 * jac[a] * pt.bvec[q];
        }
      }
      double A[36];
      double rhs[6];
      double alpha[6] = {0};
      for (int a = 0; a < 6; ++a) {
        rhs[a] = h[a];
        for (int b = 0; b < 6; ++b) {
          A[a * 6 + b] = H[a * 6 + b];
        }
      }
      if (!solve_linear6(A, rhs, alpha)) {
        continue;
      }
      for (size_t s = 0; s < stencil->size(); ++s) {
        StencilPoint& pt = (*stencil)[s];
        const double* jac = use_whiten_jacobian ? pt.jwhite : pt.jcol;
        double proj = 0.0;
        for (int a = 0; a < 6; ++a) {
          proj += jac[a] * alpha[a];
        }
        pt.bvec[q] -= proj;
      }
    }

    // Re-center/scale after orthogonalization.
    for (int q = 1; q < p; ++q) {
      double wsum = 0.0;
      double mean = 0.0;
      for (size_t s = 0; s < stencil->size(); ++s) {
        const StencilPoint& pt = (*stencil)[s];
        wsum += pt.w0;
        mean += pt.w0 * pt.bvec[q];
      }
      if (wsum <= 0.0) {
        continue;
      }
      mean /= wsum;
      double var = 0.0;
      for (size_t s = 0; s < stencil->size(); ++s) {
        const StencilPoint& pt = (*stencil)[s];
        double d = pt.bvec[q] - mean;
        var += pt.w0 * d * d;
      }
      var /= wsum;
      double sd = std::sqrt(std::max(var, 1e-8));
      for (size_t s = 0; s < stencil->size(); ++s) {
        (*stencil)[s].bvec[q] = ((*stencil)[s].bvec[q] - mean) / sd;
      }
    }
  }

  return p;
}

bool solve_linear_small(double* A, double* b, int n, double* x) {
  if (n <= 0) {
    return true;
  }
  for (int i = 0; i < n; ++i) {
    int piv = i;
    double apiv = std::abs(A[i * n + i]);
    for (int r = i + 1; r < n; ++r) {
      double av = std::abs(A[r * n + i]);
      if (av > apiv) {
        apiv = av;
        piv = r;
      }
    }
    if (apiv < 1e-12 || !R_finite(apiv)) {
      return false;
    }
    if (piv != i) {
      for (int c = i; c < n; ++c) {
        std::swap(A[i * n + c], A[piv * n + c]);
      }
      std::swap(b[i], b[piv]);
    }

    double diag = A[i * n + i];
    for (int c = i; c < n; ++c) {
      A[i * n + c] /= diag;
    }
    b[i] /= diag;

    for (int r = i + 1; r < n; ++r) {
      double f = A[r * n + i];
      if (f == 0.0) {
        continue;
      }
      for (int c = i; c < n; ++c) {
        A[r * n + c] -= f * A[i * n + c];
      }
      b[r] -= f * b[i];
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    double s = b[i];
    for (int c = i + 1; c < n; ++c) {
      s -= A[i * n + c] * x[c];
    }
    x[i] = s;
  }
  return true;
}

CostResult evaluate_cost_stencil_nuisance(const NumericVector& data,
                                          int nx,
                                          int ny,
                                          int nz,
                                          int nxyz,
                                          int t,
                                          const NumericVector& pixdim,
                                          const std::vector<StencilPoint>& stencil,
                                          const Pose6& pose,
                                          double huber_k,
                                          bool use_overlap_taper,
                                          int nuisance_dim,
                                          double nuisance_ridge,
                                          int nthreads = 1) {
  if (nuisance_dim <= 0) {
    return evaluate_cost_stencil(
        data, nx, ny, nz, nxyz, t, pixdim, stencil, pose, huber_k, use_overlap_taper, nthreads);
  }

  std::vector<double> residual(stencil.size(), 0.0);
  const int p = nuisance_dim;
  int workers = effective_thread_count(nthreads, stencil.size());
  std::vector<double> inb_local(static_cast<size_t>(workers), 0.0);
  std::vector<double> wsum_local(static_cast<size_t>(workers), 0.0);
  std::vector< std::vector<double> > bt_local(static_cast<size_t>(workers), std::vector<double>(5, 0.0));
  std::vector< std::vector<double> > btb_local(static_cast<size_t>(workers), std::vector<double>(25, 0.0));

  parallel_for_range(stencil.size(), workers, [&](size_t begin, size_t end) {
    size_t tid = begin == end ? 0 : begin / std::max<size_t>(1, (stencil.size() + static_cast<size_t>(workers) - 1) / static_cast<size_t>(workers));
    if (tid >= static_cast<size_t>(workers)) {
      tid = static_cast<size_t>(workers - 1);
    }
    double inb = 0.0;
    double wsum = 0.0;
    double src[3];

    for (size_t s = begin; s < end; ++s) {
      const StencilPoint& pt = stencil[s];
      bool in_bounds = false;
      fixed_to_moving_vox(pt.i, pt.j, pt.k, nx, ny, nz, pixdim, pose, src, &in_bounds);
      double mv = 0.0;
      if (in_bounds) {
        mv = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
        inb += pt.w0;
      }
      double r = mv - pt.tval;
      residual[s] = r;
      wsum += pt.w0;
      for (int q = 0; q < p; ++q) {
        double bq = pt.bvec[q];
        bt_local[tid][static_cast<size_t>(q)] += pt.w0 * bq * r;
        for (int qq = q; qq < p; ++qq) {
          btb_local[tid][static_cast<size_t>(q * p + qq)] += pt.w0 * bq * pt.bvec[qq];
        }
      }
    }

    inb_local[tid] = inb;
    wsum_local[tid] = wsum;
  });

  double inb = 0.0;
  double wsum = 0.0;
  double BtB[25] = {0};
  double bt[5] = {0};
  for (int tid = 0; tid < workers; ++tid) {
    inb += inb_local[static_cast<size_t>(tid)];
    wsum += wsum_local[static_cast<size_t>(tid)];
    for (int q = 0; q < p; ++q) {
      bt[q] += bt_local[static_cast<size_t>(tid)][static_cast<size_t>(q)];
      for (int qq = q; qq < p; ++qq) {
        BtB[q * p + qq] += btb_local[static_cast<size_t>(tid)][static_cast<size_t>(q * p + qq)];
      }
    }
  }
  if (wsum <= 0.0) {
    return CostResult{1e12, 0.0};
  }
  for (int q = 0; q < p; ++q) {
    for (int qq = q + 1; qq < p; ++qq) {
      BtB[qq * p + q] = BtB[q * p + qq];
    }
  }

  double tr = 0.0;
  for (int q = 0; q < p; ++q) {
    tr += BtB[q * p + q];
  }
  double eps = std::max(1e-8, std::max(0.0, nuisance_ridge) * std::max(1.0, tr / std::max(1, p)));
  for (int q = 0; q < p; ++q) {
    BtB[q * p + q] += eps;
  }

  double A[25];
  double rhs[5];
  double theta[5] = {0};
  for (int q = 0; q < p; ++q) {
    rhs[q] = bt[q];
    for (int qq = 0; qq < p; ++qq) {
      A[q * p + qq] = BtB[q * p + qq];
    }
  }
  bool solved = solve_linear_small(A, rhs, p, theta);

  std::vector<double> ss_local(static_cast<size_t>(workers), 0.0);
  parallel_for_range(stencil.size(), workers, [&](size_t begin, size_t end) {
    size_t tid = begin == end ? 0 : begin / std::max<size_t>(1, (stencil.size() + static_cast<size_t>(workers) - 1) / static_cast<size_t>(workers));
    if (tid >= static_cast<size_t>(workers)) {
      tid = static_cast<size_t>(workers - 1);
    }
    double ss = 0.0;
    for (size_t s = begin; s < end; ++s) {
      const StencilPoint& pt = stencil[s];
      double r = residual[s];
      if (solved) {
        double proj = 0.0;
        for (int q = 0; q < p; ++q) {
          proj += pt.bvec[q] * theta[q];
        }
        r -= proj;
      }
      double ar = std::abs(r);
      double rho = (ar <= huber_k) ? (r * r) : (2.0 * huber_k * ar - huber_k * huber_k);
      ss += pt.w0 * rho;
    }
    ss_local[tid] = ss;
  });

  double ss = 0.0;
  for (int tid = 0; tid < workers; ++tid) {
    ss += ss_local[static_cast<size_t>(tid)];
  }

  double overlap = inb / wsum;
  double c = ss / wsum;
  if (use_overlap_taper) {
    double miss = 1.0 - overlap;
    c *= (1.0 + 4.0 * miss * miss);
  }
  if (overlap <= 0.0) {
    c = 1e12;
  }
  return CostResult{c, overlap};
}

bool solve_linear6(double* A, double* b, double* x) {
  for (int i = 0; i < 6; ++i) {
    int piv = i;
    double apiv = std::abs(A[i * 6 + i]);
    for (int r = i + 1; r < 6; ++r) {
      double av = std::abs(A[r * 6 + i]);
      if (av > apiv) {
        apiv = av;
        piv = r;
      }
    }
    if (apiv < 1e-10 || !R_finite(apiv)) {
      return false;
    }
    if (piv != i) {
      for (int c = i; c < 6; ++c) {
        std::swap(A[i * 6 + c], A[piv * 6 + c]);
      }
      std::swap(b[i], b[piv]);
    }

    double diag = A[i * 6 + i];
    for (int c = i; c < 6; ++c) {
      A[i * 6 + c] /= diag;
    }
    b[i] /= diag;

    for (int r = i + 1; r < 6; ++r) {
      double f = A[r * 6 + i];
      if (f == 0.0) {
        continue;
      }
      for (int c = i; c < 6; ++c) {
        A[r * 6 + c] -= f * A[i * 6 + c];
      }
      b[r] -= f * b[i];
    }
  }

  for (int i = 5; i >= 0; --i) {
    double s = b[i];
    for (int c = i + 1; c < 6; ++c) {
      s -= A[i * 6 + c] * x[c];
    }
    x[i] = s;
  }
  return true;
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

std::vector<Pose6> build_restart_candidates(const Pose6& init_pose,
                                            bool capture_boost,
                                            double capture_trans_mm,
                                            double capture_rot_deg) {
  std::vector<Pose6> cand;
  cand.reserve(48);

  auto add_xyz = [&](double dmm) {
    cand.push_back(make_pose(init_pose.tx + dmm, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx - dmm, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty + dmm, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty - dmm, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty, init_pose.tz + dmm, init_pose.rx, init_pose.ry, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty, init_pose.tz - dmm, init_pose.rx, init_pose.ry, init_pose.rz));
  };
  auto add_rxyz = [&](double ddeg) {
    double dr = ddeg * M_PI / 180.0;
    cand.push_back(make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx + dr, init_pose.ry, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx - dr, init_pose.ry, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry + dr, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry - dr, init_pose.rz));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz + dr));
    cand.push_back(make_pose(init_pose.tx, init_pose.ty, init_pose.tz, init_pose.rx, init_pose.ry, init_pose.rz - dr));
  };

  cand.push_back(make_pose());  // identity
  add_xyz(2.0);
  add_rxyz(2.0);

  if (capture_boost) {
    double ct = std::max(4.0, capture_trans_mm);
    double cr = std::max(4.0, capture_rot_deg);
    add_xyz(0.5 * ct);
    add_rxyz(0.5 * cr);
    add_xyz(ct);
    add_rxyz(cr);
  }
  return cand;
}

template <typename EvalFn>
std::vector<Pose6> select_best_restart_seeds(const Pose6& init_pose,
                                             const std::vector<Pose6>& candidates,
                                             int n_keep,
                                             EvalFn eval_cost) {
  std::vector<Pose6> out;
  out.push_back(init_pose);
  if (n_keep <= 0 || candidates.empty()) {
    return out;
  }

  struct ScoredPose {
    double cost;
    Pose6 pose;
  };
  std::vector<ScoredPose> scored;
  scored.reserve(candidates.size());
  for (size_t i = 0; i < candidates.size(); ++i) {
    double c = eval_cost(candidates[i]);
    scored.push_back(ScoredPose{c, candidates[i]});
  }
  std::sort(scored.begin(), scored.end(), [](const ScoredPose& a, const ScoredPose& b) { return a.cost < b.cost; });
  int keep = std::min(n_keep, static_cast<int>(scored.size()));
  for (int i = 0; i < keep; ++i) {
    out.push_back(scored[i].pose);
  }
  return out;
}

FrameFit optimize_frame_pose_ic(const NumericVector& data,
                                int nx,
                                int ny,
                                int nz,
                                int nxyz,
                                int t,
                                const NumericVector& pixdim,
                                const std::vector<StencilPoint>& stencil,
                                const Pose6& init_pose,
                                int max_iter,
                                double lambda0,
                                double tol_step,
                                double tol_cost,
                                double huber_k,
                                bool adaptive_huber,
                                int n_restarts,
                                bool use_overlap_taper,
                                bool capture_boost,
                                double capture_trans_mm,
                                double capture_rot_deg,
                                int capture_topk,
                                bool use_whiten_modes,
                                const double* Linv_u2xi,
                                bool use_nuisance_marginalization,
                                int nuisance_dim,
                                double nuisance_ridge,
                                int nthreads) {
  bool use_nuisance = use_nuisance_marginalization && nuisance_dim > 0;
  double huber_eff = huber_k;
  if (adaptive_huber) {
    double sigma = estimate_residual_sigma_stencil(
        data, nx, ny, nz, nxyz, t, pixdim, stencil, init_pose);
    huber_eff = std::max(1e-6, huber_k * sigma);
  }
  auto eval_pose = [&](const Pose6& pose) {
    if (use_nuisance) {
      return evaluate_cost_stencil_nuisance(data,
                                            nx,
                                            ny,
                                            nz,
                                            nxyz,
                                            t,
                                            pixdim,
                                            stencil,
                                            pose,
                                            huber_eff,
                                            use_overlap_taper,
                                            nuisance_dim,
                                            nuisance_ridge,
                                            nthreads);
    }
    return evaluate_cost_stencil(
        data, nx, ny, nz, nxyz, t, pixdim, stencil, pose, huber_eff, use_overlap_taper, nthreads);
  };

  std::vector<Pose6> seeds;
  int n_keep = std::max(n_restarts, capture_boost ? capture_topk : n_restarts);
  if (n_keep > 0 || capture_boost) {
    std::vector<Pose6> cand =
        build_restart_candidates(init_pose, capture_boost, capture_trans_mm, capture_rot_deg);
    seeds = select_best_restart_seeds(
        init_pose,
        cand,
        n_keep,
        [&](const Pose6& p) {
          return eval_pose(p).cost;
        });
  } else {
    seeds.push_back(init_pose);
  }

  Pose6 best_pose = init_pose;
  CostResult best_eval = eval_pose(init_pose);
  double best_c_init = best_eval.cost;
  double best_c_final = best_eval.cost;
  int best_iters = 0;
  bool best_restart = false;
  int stencil_workers = effective_thread_count(nthreads, stencil.size());

  for (size_t seed_idx = 0; seed_idx < seeds.size(); ++seed_idx) {
    Pose6 cur = seeds[seed_idx];
    double lambda = lambda0;
    CostResult eval0 = eval_pose(cur);
    double c_init = eval0.cost;
    double c_prev = c_init;
    int used_iters = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
      used_iters = iter + 1;
      double H[36] = {0};
      double g[6] = {0};
      double C[30] = {0};      // 6 x max 5 nuisance columns
      double BtB[25] = {0};    // max 5 x 5
      double bt[5] = {0};
      double valid_w = 0.0;
      std::vector<double> H_local(static_cast<size_t>(stencil_workers) * 36, 0.0);
      std::vector<double> g_local(static_cast<size_t>(stencil_workers) * 6, 0.0);
      std::vector<double> C_local(static_cast<size_t>(stencil_workers) * 30, 0.0);
      std::vector<double> BtB_local(static_cast<size_t>(stencil_workers) * 25, 0.0);
      std::vector<double> bt_local(static_cast<size_t>(stencil_workers) * 5, 0.0);
      std::vector<double> valid_w_local(static_cast<size_t>(stencil_workers), 0.0);

      parallel_for_range(stencil.size(), stencil_workers, [&](size_t begin, size_t end) {
        size_t tid = begin == end ? 0 : begin / std::max<size_t>(1, (stencil.size() + static_cast<size_t>(stencil_workers) - 1) / static_cast<size_t>(stencil_workers));
        if (tid >= static_cast<size_t>(stencil_workers)) {
          tid = static_cast<size_t>(stencil_workers - 1);
        }
        double* Ht = &H_local[tid * 36];
        double* gt = &g_local[tid * 6];
        double* Ct = &C_local[tid * 30];
        double* BtBt = &BtB_local[tid * 25];
        double* btt = &bt_local[tid * 5];
        double valid_w_thread = 0.0;
        double src[3];

        for (size_t s = begin; s < end; ++s) {
          const StencilPoint& p = stencil[s];
          bool in_bounds = false;
          fixed_to_moving_vox(p.i, p.j, p.k, nx, ny, nz, pixdim, cur, src, &in_bounds);
          if (!in_bounds) {
            continue;
          }

          double mv = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
          double r = mv - p.tval;
          double ar = std::abs(r);
          double wr = (ar <= huber_eff) ? 1.0 : (huber_eff / (ar + 1e-12));
          double wt = p.w0 * wr;
          valid_w_thread += wt;
          const double* jac = use_whiten_modes ? p.jwhite : p.jcol;

          for (int a = 0; a < 6; ++a) {
            double ja = jac[a];
            gt[a] += wt * ja * r;
            for (int b = a; b < 6; ++b) {
              Ht[a * 6 + b] += wt * ja * jac[b];
            }
            if (use_nuisance) {
              for (int q = 0; q < nuisance_dim; ++q) {
                Ct[a * nuisance_dim + q] += wt * ja * p.bvec[q];
              }
            }
          }
          if (use_nuisance) {
            for (int q = 0; q < nuisance_dim; ++q) {
              double bq = p.bvec[q];
              btt[q] += wt * bq * r;
              for (int qq = q; qq < nuisance_dim; ++qq) {
                BtBt[q * nuisance_dim + qq] += wt * bq * p.bvec[qq];
              }
            }
          }
        }

        valid_w_local[tid] = valid_w_thread;
      });

      for (int tid = 0; tid < stencil_workers; ++tid) {
        valid_w += valid_w_local[static_cast<size_t>(tid)];
        const double* Ht = &H_local[static_cast<size_t>(tid) * 36];
        const double* gt = &g_local[static_cast<size_t>(tid) * 6];
        const double* Ct = &C_local[static_cast<size_t>(tid) * 30];
        const double* BtBt = &BtB_local[static_cast<size_t>(tid) * 25];
        const double* btt = &bt_local[static_cast<size_t>(tid) * 5];
        for (int a = 0; a < 6; ++a) {
          g[a] += gt[a];
          for (int b = a; b < 6; ++b) {
            H[a * 6 + b] += Ht[a * 6 + b];
          }
          if (use_nuisance) {
            for (int q = 0; q < nuisance_dim; ++q) {
              C[a * nuisance_dim + q] += Ct[a * nuisance_dim + q];
            }
          }
        }
        if (use_nuisance) {
          for (int q = 0; q < nuisance_dim; ++q) {
            bt[q] += btt[q];
            for (int qq = q; qq < nuisance_dim; ++qq) {
              BtB[q * nuisance_dim + qq] += BtBt[q * nuisance_dim + qq];
            }
          }
        }
      }

      if (valid_w <= 0.0) {
        break;
      }
      for (int a = 0; a < 6; ++a) {
        for (int b = a + 1; b < 6; ++b) {
          H[b * 6 + a] = H[a * 6 + b];
        }
      }

      double Hsolve[36];
      double gsolve[6];
      for (int a = 0; a < 6; ++a) {
        gsolve[a] = g[a];
        for (int b = 0; b < 6; ++b) {
          Hsolve[a * 6 + b] = H[a * 6 + b];
        }
      }

      if (use_nuisance) {
        for (int q = 0; q < nuisance_dim; ++q) {
          for (int qq = q + 1; qq < nuisance_dim; ++qq) {
            BtB[qq * nuisance_dim + q] = BtB[q * nuisance_dim + qq];
          }
        }

        double trb = 0.0;
        for (int q = 0; q < nuisance_dim; ++q) {
          trb += BtB[q * nuisance_dim + q];
        }
        double eps_b =
            std::max(1e-8,
                     std::max(0.0, nuisance_ridge) *
                         std::max(1.0, trb / std::max(1, nuisance_dim)));
        for (int q = 0; q < nuisance_dim; ++q) {
          BtB[q * nuisance_dim + q] += eps_b;
        }

        double Binv[25] = {0};
        bool ok_b = true;
        for (int col = 0; col < nuisance_dim; ++col) {
          double Acol[25];
          double rhs_col[5] = {0};
          double x_col[5] = {0};
          rhs_col[col] = 1.0;
          for (int q = 0; q < nuisance_dim; ++q) {
            for (int qq = 0; qq < nuisance_dim; ++qq) {
              Acol[q * nuisance_dim + qq] = BtB[q * nuisance_dim + qq];
            }
          }
          if (!solve_linear_small(Acol, rhs_col, nuisance_dim, x_col)) {
            ok_b = false;
            break;
          }
          for (int q = 0; q < nuisance_dim; ++q) {
            Binv[q * nuisance_dim + col] = x_col[q];
          }
        }

        if (ok_b) {
          double Y[30] = {0};  // (J'WB) (B'WB)^-1
          for (int a = 0; a < 6; ++a) {
            for (int q = 0; q < nuisance_dim; ++q) {
              double s = 0.0;
              for (int qq = 0; qq < nuisance_dim; ++qq) {
                s += C[a * nuisance_dim + qq] * Binv[qq * nuisance_dim + q];
              }
              Y[a * nuisance_dim + q] = s;
            }
          }
          for (int a = 0; a < 6; ++a) {
            for (int b = 0; b < 6; ++b) {
              double sub = 0.0;
              for (int q = 0; q < nuisance_dim; ++q) {
                sub += Y[a * nuisance_dim + q] * C[b * nuisance_dim + q];
              }
              Hsolve[a * 6 + b] -= sub;
            }
            double subg = 0.0;
            for (int q = 0; q < nuisance_dim; ++q) {
              subg += Y[a * nuisance_dim + q] * bt[q];
            }
            gsolve[a] -= subg;
          }
        }
      }

      double A[36];
      double rhs[6];
      double delta_v[6] = {0, 0, 0, 0, 0, 0};
      for (int a = 0; a < 6; ++a) {
        for (int b = 0; b < 6; ++b) {
          A[a * 6 + b] = Hsolve[a * 6 + b];
        }
        if (use_whiten_modes) {
          A[a * 6 + a] += lambda;
        } else {
          A[a * 6 + a] += lambda * std::max(Hsolve[a * 6 + a], 1e-6);
        }
        rhs[a] = -gsolve[a];
      }

      if (!solve_linear6(A, rhs, delta_v)) {
        lambda *= 4.0;
        continue;
      }

      double dxi[6];
      if (use_whiten_modes && Linv_u2xi != nullptr) {
        for (int a = 0; a < 6; ++a) {
          double s = 0.0;
          for (int b = 0; b < 6; ++b) {
            s += Linv_u2xi[a * 6 + b] * delta_v[b];
          }
          dxi[a] = s;
        }
      } else {
        for (int a = 0; a < 6; ++a) {
          dxi[a] = delta_v[a];
        }
      }
      Pose6 delta = make_pose(dxi[0], dxi[1], dxi[2], dxi[3], dxi[4], dxi[5]);
      double dnorm = std::sqrt(delta.tx * delta.tx + delta.ty * delta.ty + delta.tz * delta.tz +
                               delta.rx * delta.rx + delta.ry * delta.ry + delta.rz * delta.rz);

      bool accepted = false;
      double alpha = 1.0;
      double c_new = c_prev;

      for (int ls = 0; ls < 8; ++ls) {
        Pose6 cand = add_pose(cur, scale_pose(delta, alpha));
        c_new = eval_pose(cand).cost;
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

    CostResult evalf = eval_pose(cur);
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
                             bool adaptive_huber,
                             int n_restarts,
                             bool use_overlap_taper,
                             bool capture_boost,
                             double capture_trans_mm,
                             double capture_rot_deg,
                             int capture_topk) {
  std::vector<Pose6> seeds;
  int n_keep = std::max(n_restarts, capture_boost ? capture_topk : n_restarts);
  if (n_keep > 0 || capture_boost) {
    std::vector<Pose6> cand =
        build_restart_candidates(init_pose, capture_boost, capture_trans_mm, capture_rot_deg);
    seeds = select_best_restart_seeds(
        init_pose,
        cand,
        n_keep,
        [&](const Pose6& p) {
          return evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, p, huber_k, use_overlap_taper).cost;
        });
  } else {
    seeds.push_back(init_pose);
  }

  Pose6 best_pose = init_pose;
  CostResult best_eval =
      evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, init_pose, huber_k, use_overlap_taper);
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
    double huber_eff = huber_k;
    if (adaptive_huber) {
      double sigma = estimate_residual_sigma_dense(
          data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, cur);
      huber_eff = std::max(1e-6, huber_k * sigma);
    }
    double lambda = lambda0;
    CostResult eval0 =
        evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, cur, huber_eff, use_overlap_taper);
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

        double cp =
            evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, p_plus, huber_eff, use_overlap_taper)
                .cost;
        double cm =
            evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, p_minus, huber_eff, use_overlap_taper)
                .cost;

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
        c_new = evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, cand, huber_eff, use_overlap_taper)
                    .cost;
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

    CostResult evalf =
        evaluate_cost(data, tmpl, nx, ny, nz, nxyz, t, pixdim, samples, cur, huber_eff, use_overlap_taper);

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
                                             double trim_frac,
                                             bool valid_only,
                                             int nthreads = 1) {
  NumericVector out(nxyz);
  int trim_n = static_cast<int>(std::floor(trim_frac * nt));
  trim_n = std::max(0, std::min(trim_n, (nt - 1) / 2));

  parallel_for_range(static_cast<size_t>(nxyz), nthreads, [&](size_t begin, size_t end) {
    std::vector<double> vals(static_cast<size_t>(nt));
    double src[3];
    for (size_t vv = begin; vv < end; ++vv) {
      int v = static_cast<int>(vv);
      int k = v / (nx * ny);
      int rem = v - k * nx * ny;
      int j = rem / nx;
      int i = rem - j * nx;

      int n_valid = 0;
      for (int t = 0; t < nt; ++t) {
        bool in_bounds = false;
        fixed_to_moving_vox(i, j, k, nx, ny, nz, pixdim, poses[static_cast<size_t>(t)], src, &in_bounds);
        if (in_bounds) {
          double mv = sample_trilinear_frame(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2]);
          if (R_finite(mv)) {
            vals[static_cast<size_t>(n_valid)] = mv;
            ++n_valid;
          }
        } else if (!valid_only) {
          vals[static_cast<size_t>(n_valid)] = 0.0;
          ++n_valid;
        }
      }

      if (n_valid <= 0) {
        out[v] = 0.0;
        continue;
      }

      std::sort(vals.begin(), vals.begin() + n_valid);
      int trim_here = std::min(trim_n, (n_valid - 1) / 2);
      int lo = trim_here;
      int hi = n_valid - trim_here;
      double s = 0.0;
      for (int tt = lo; tt < hi; ++tt) {
        s += vals[static_cast<size_t>(tt)];
      }
      out[v] = s / static_cast<double>(hi - lo);
    }
  });

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

std::vector<int> get_int_vector_ctrl(const List& control,
                                     const std::string& name,
                                     const std::vector<int>& default_values) {
  if (!control.containsElementNamed(name.c_str())) {
    return default_values;
  }
  SEXP v = control[name];
  std::vector<int> out;
  if (Rf_isInteger(v)) {
    IntegerVector iv(v);
    out.reserve(iv.size());
    for (int i = 0; i < iv.size(); ++i) {
      out.push_back(iv[i]);
    }
  } else if (Rf_isReal(v)) {
    NumericVector nv(v);
    out.reserve(nv.size());
    for (int i = 0; i < nv.size(); ++i) {
      out.push_back(static_cast<int>(std::round(nv[i])));
    }
  }
  if (out.empty()) {
    return default_values;
  }
  return out;
}

int level_value(const std::vector<int>& values, size_t level, int fallback) {
  if (values.empty()) {
    return fallback;
  }
  size_t idx = std::min(level, values.size() - 1);
  return values[idx];
}

std::vector<int> sanitize_downsample_levels(std::vector<int> levels) {
  std::vector<int> out;
  out.reserve(levels.size() + 1);
  for (size_t i = 0; i < levels.size(); ++i) {
    int ds = std::max(1, levels[i]);
    if (out.empty() || out.back() != ds) {
      out.push_back(ds);
    }
  }
  if (out.empty()) {
    out.push_back(1);
  }
  if (out.back() != 1) {
    out.push_back(1);
  }
  return out;
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
  std::string ref_mode = get_str_ctrl(control, "ref_mode", "middle");
  if (control.containsElementNamed("robust_template")) {
    robust_template = as<bool>(control["robust_template"]);
  }
  if (control.containsElementNamed("template_trim")) {
    trim_frac = as<double>(control["template_trim"]);
  }

  std::vector<int> downsample_levels =
      sanitize_downsample_levels(get_int_vector_ctrl(control, "downsample", std::vector<int>{1}));
  bool use_pyramid = false;
  if (control.containsElementNamed("use_pyramid")) {
    use_pyramid = as<bool>(control["use_pyramid"]);
  }
  if (!use_pyramid) {
    downsample_levels = std::vector<int>{1};
  }
  std::vector<int> max_iter_levels =
      get_int_vector_ctrl(control, "max_iter", std::vector<int>{20});
  std::vector<int> n_samples_levels =
      get_int_vector_ctrl(control, "n_samples", std::vector<int>{20000});
  int max_iter = std::max(0, level_value(max_iter_levels, max_iter_levels.size() - 1, 20));
  int n_samples = std::max(0, level_value(n_samples_levels, n_samples_levels.size() - 1, 20000));
  if (!use_pyramid) {
    max_iter_levels = std::vector<int>{max_iter};
    n_samples_levels = std::vector<int>{n_samples};
  }
  int n_restarts = get_int_ctrl(control, "n_restarts", 0);
  double huber_k = get_dbl_ctrl(control, "huber_k", 1.5);
  bool adaptive_huber = false;
  if (control.containsElementNamed("adaptive_huber")) {
    adaptive_huber = as<bool>(control["adaptive_huber"]);
  }
  bool adaptive_huber_guard = false;
  if (control.containsElementNamed("adaptive_huber_guard")) {
    adaptive_huber_guard = as<bool>(control["adaptive_huber_guard"]);
  }
  double lambda0 = get_dbl_ctrl(control, "lambda0", 1e-2);
  double tol_step = get_dbl_ctrl(control, "tol_step", 1e-5);
  double tol_cost = get_dbl_ctrl(control, "tol_cost", 1e-6);
  bool use_overlap_taper = true;
  if (control.containsElementNamed("use_overlap_taper")) {
    use_overlap_taper = as<bool>(control["use_overlap_taper"]);
  }
  bool use_ic_stencil = false;
  if (control.containsElementNamed("use_ic_stencil")) {
    use_ic_stencil = as<bool>(control["use_ic_stencil"]);
  }
  bool whiten_modes = false;
  if (control.containsElementNamed("whiten_modes")) {
    whiten_modes = as<bool>(control["whiten_modes"]);
  }
  double whiten_reg = get_dbl_ctrl(control, "whiten_reg", 1e-6);
  bool use_nuisance_marginalization = false;
  if (control.containsElementNamed("use_nuisance_marginalization")) {
    use_nuisance_marginalization = as<bool>(control["use_nuisance_marginalization"]);
  }
  bool nuisance_include_template = false;
  if (control.containsElementNamed("nuisance_include_template")) {
    nuisance_include_template = as<bool>(control["nuisance_include_template"]);
  }
  bool nuisance_orthogonalize = true;
  if (control.containsElementNamed("nuisance_orthogonalize")) {
    nuisance_orthogonalize = as<bool>(control["nuisance_orthogonalize"]);
  }
  double nuisance_ridge = get_dbl_ctrl(control, "nuisance_ridge", 1e-4);
  double nuisance_orth_reg = get_dbl_ctrl(control, "nuisance_orth_reg", 1e-6);
  bool capture_boost = true;
  if (control.containsElementNamed("capture_boost")) {
    capture_boost = as<bool>(control["capture_boost"]);
  }
  double capture_trans_mm = get_dbl_ctrl(control, "capture_trans_mm", 8.0);
  double capture_rot_deg = get_dbl_ctrl(control, "capture_rot_deg", 8.0);
  int capture_topk = get_int_ctrl(control, "capture_topk", 6);
  int stencil_n = get_int_ctrl(control, "stencil_n", n_samples);
  int stencil_bx = get_int_ctrl(control, "stencil_bx", 4);
  int stencil_by = get_int_ctrl(control, "stencil_by", 4);
  int stencil_bz = get_int_ctrl(control, "stencil_bz", 2);
  double stencil_gamma = get_dbl_ctrl(control, "stencil_gamma", 0.5);
  bool use_dense_polish = false;
  if (control.containsElementNamed("use_dense_polish")) {
    use_dense_polish = as<bool>(control["use_dense_polish"]);
  }
  bool use_ic_polish = false;
  if (control.containsElementNamed("use_ic_polish")) {
    use_ic_polish = as<bool>(control["use_ic_polish"]);
  }
  int polish_iter = get_int_ctrl(control, "polish_iter", 5);
  int polish_n_samples = get_int_ctrl(control, "polish_n_samples", 0);
  bool parallel_frames = false;
  if (control.containsElementNamed("parallel_frames")) {
    parallel_frames = as<bool>(control["parallel_frames"]);
  }
  double edge_exclude_frac = get_dbl_ctrl(control, "edge_exclude_frac", 0.0);
  int edge_exclude_vox = get_int_ctrl(control, "edge_exclude_vox", 0);
  int nthreads = std::max(1, get_int_ctrl(control, "nthreads", 1));

  NumericMatrix par(nt, 6);
  NumericMatrix mats = identity_mats(nt);

  bool robust_mean_ref = (ref_mode == "robust_mean") && nt > 1;
  NumericVector template_vol = robust_mean_ref
                                   ? trimmed_mean_template(data, nx, ny, nz, nt, trim_frac, nthreads)
                                   : extract_ref_volume(data, nx, ny, nz, ref_index);

  LogicalVector fit_mask = build_mask(mask, nxyz);
  if (mask.isNull()) {
    fit_mask = auto_mask_from_template(template_vol, nxyz);
  }
  fit_mask = exclude_mask_edges(fit_mask, nx, ny, nz, edge_exclude_frac, edge_exclude_vox);
  std::vector<int> samples = make_sample_points(fit_mask, nxyz, n_samples);
  if (samples.empty()) {
    stop("mask produced zero fitting samples");
  }
  std::vector<int> polish_samples;
  if (use_dense_polish && polish_iter > 0) {
    int polish_target = polish_n_samples > 0 ? polish_n_samples : n_samples;
    polish_samples = make_sample_points(fit_mask, nxyz, polish_target);
    if (polish_samples.empty()) {
      stop("mask produced zero polishing samples");
    }
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

  int ref_t = robust_mean_ref
                  ? select_reference_frame_from_template(data, template_vol, nt, nxyz, fit_mask)
                  : (ref_index - 1);
  Vec3 com_ref = center_of_mass_3d(template_vol, nx, ny, nz, pixdim, fit_mask);
  Vec3 com_init_ref = center_of_mass_frame(data, nx, ny, nz, nxyz, ref_t, pixdim, fit_mask);
  Pose6 ref_init = make_pose(com_ref.x - com_init_ref.x,
                             com_ref.y - com_init_ref.y,
                             com_ref.z - com_init_ref.z,
                             0.0,
                             0.0,
                             0.0);
  int warm_start_passes = 0;
  int momentum_backward = 0;
  int momentum_forward = 0;
  int com_seed_count = 0;
  int blended_seed_count = 0;
  int parallel_frame_fits = 0;
  auto solve_sequence = [&](const Pose6& ref_seed, bool blend_current_pose_seed) {
    warm_start_passes += 1;
    bool parallel_frame_pass = parallel_frames && blend_current_pose_seed && nthreads > 1 && nt > 2;
    int fit_nthreads = parallel_frame_pass ? 1 : nthreads;
    struct LevelContext {
      int downsample;
      int max_iter;
      std::vector<int> samples;
      std::vector<int> polish_samples;
      std::vector<StencilPoint> stencil;
      std::vector<StencilPoint> polish_stencil;
      bool use_ic;
      bool use_whiten;
      bool polish_use_whiten;
      int nuisance_dim;
      int polish_nuisance_dim;
      WhitenMetric wm;
      WhitenMetric polish_wm;
    };

    std::vector<LevelContext> levels;
    levels.reserve(downsample_levels.size());
    for (size_t li = 0; li < downsample_levels.size(); ++li) {
      int ds = downsample_levels[li];
      LogicalVector level_mask = make_downsampled_mask(fit_mask, nx, ny, nz, ds);
      int level_n_samples = std::max(0, level_value(n_samples_levels, li, n_samples));
      int level_max_iter = std::max(0, level_value(max_iter_levels, li, max_iter));

      LevelContext level;
      level.downsample = ds;
      level.max_iter = level_max_iter;
      level.samples = make_sample_points(level_mask, nxyz, level_n_samples);
      if (level.samples.empty()) {
        level.samples = samples;
      }
      level.use_ic = false;
      level.use_whiten = false;
      level.polish_use_whiten = false;
      level.nuisance_dim = 0;
      level.polish_nuisance_dim = 0;
      level.wm.ok = false;
      level.polish_wm.ok = false;
      set_identity6(level.wm.L);
      set_identity6(level.wm.Linv);
      set_identity6(level.wm.LinvT);
      set_identity6(level.polish_wm.L);
      set_identity6(level.polish_wm.Linv);
      set_identity6(level.polish_wm.LinvT);

      if (use_dense_polish && polish_iter > 0 && li + 1 == downsample_levels.size()) {
        int polish_target = polish_n_samples > 0 ? polish_n_samples : level_n_samples;
        level.polish_samples = make_sample_points(level_mask, nxyz, polish_target);
        if (level.polish_samples.empty()) {
          level.polish_samples = polish_samples;
        }
      }

      if (use_ic_stencil) {
        level.stencil = build_ic_stencil(template_vol,
                                         data,
                                         nx,
                                         ny,
                                         nz,
                                         nt,
                                         nxyz,
                                         pixdim,
                                         level_mask,
                                         stencil_n,
                                         stencil_bx,
                                         stencil_by,
                                         stencil_bz,
                                         stencil_gamma);
        if (!level.stencil.empty()) {
          level.use_ic = true;
          if (whiten_modes) {
            level.wm = build_whiten_metric(level.stencil, whiten_reg);
            if (level.wm.ok) {
              apply_whiten_metric_to_stencil(&level.stencil, level.wm);
              level.use_whiten = true;
            } else {
              copy_raw_jacobian_to_whitened(&level.stencil);
            }
          } else {
            copy_raw_jacobian_to_whitened(&level.stencil);
          }
          if (use_nuisance_marginalization) {
            level.nuisance_dim = prepare_nuisance_basis(&level.stencil,
                                                        nx,
                                                        ny,
                                                        nz,
                                                        nuisance_include_template,
                                                        nuisance_orthogonalize,
                                                        level.use_whiten,
                                                        nuisance_orth_reg);
          }
        }
      }
      if (use_ic_polish && polish_iter > 0 && li + 1 == downsample_levels.size()) {
        int polish_target = polish_n_samples > 0 ? polish_n_samples : std::max(stencil_n * 2, stencil_n);
        level.polish_stencil = build_ic_stencil(template_vol,
                                                data,
                                                nx,
                                                ny,
                                                nz,
                                                nt,
                                                nxyz,
                                                pixdim,
                                                level_mask,
                                                polish_target,
                                                stencil_bx,
                                                stencil_by,
                                                stencil_bz,
                                                stencil_gamma);
        if (!level.polish_stencil.empty()) {
          if (whiten_modes) {
            level.polish_wm = build_whiten_metric(level.polish_stencil, whiten_reg);
            if (level.polish_wm.ok) {
              apply_whiten_metric_to_stencil(&level.polish_stencil, level.polish_wm);
              level.polish_use_whiten = true;
            } else {
              copy_raw_jacobian_to_whitened(&level.polish_stencil);
            }
          } else {
            copy_raw_jacobian_to_whitened(&level.polish_stencil);
          }
          if (use_nuisance_marginalization) {
            level.polish_nuisance_dim = prepare_nuisance_basis(&level.polish_stencil,
                                                               nx,
                                                               ny,
                                                               nz,
                                                               nuisance_include_template,
                                                               nuisance_orthogonalize,
                                                               level.polish_use_whiten,
                                                               nuisance_orth_reg);
          }
        }
      }
      levels.push_back(level);
    }

    auto fit_level = [&](const LevelContext& level, int frame_t, const Pose6& init_pose, bool first_level, bool final_level) {
      if (level.use_ic) {
        auto run_ic = [&](bool use_adaptive) {
          return optimize_frame_pose_ic(data,
                                        nx,
                                        ny,
                                        nz,
                                        nxyz,
                                        frame_t,
                                        pixdim,
                                        level.stencil,
                                        init_pose,
                                        level.max_iter,
                                        lambda0,
                                        tol_step,
                                        tol_cost,
                                        huber_k,
                                        use_adaptive,
                                        first_level ? n_restarts : 0,
                                        use_overlap_taper,
                                        first_level && capture_boost,
                                        capture_trans_mm,
                                        capture_rot_deg,
                                        first_level ? capture_topk : 0,
                                        level.use_whiten,
                                        level.use_whiten ? level.wm.LinvT : nullptr,
                                        use_nuisance_marginalization,
                                        level.nuisance_dim,
                                        nuisance_ridge,
                                        fit_nthreads);
        };
        auto fixed_objective = [&](const Pose6& pose) {
          if (use_nuisance_marginalization && level.nuisance_dim > 0) {
            return evaluate_cost_stencil_nuisance(data,
                                                  nx,
                                                  ny,
                                                  nz,
                                                  nxyz,
                                                  frame_t,
                                                  pixdim,
                                                  level.stencil,
                                                  pose,
                                                  huber_k,
                                                  use_overlap_taper,
                                                  level.nuisance_dim,
                                                  nuisance_ridge,
                                                  fit_nthreads)
                .cost;
          }
          return evaluate_cost_stencil(data,
                                       nx,
                                       ny,
                                       nz,
                                       nxyz,
                                       frame_t,
                                       pixdim,
                                       level.stencil,
                                       pose,
                                       huber_k,
                                       use_overlap_taper,
                                       fit_nthreads)
              .cost;
        };

        FrameFit fit_ic = run_ic(adaptive_huber);
        if (adaptive_huber && adaptive_huber_guard) {
          FrameFit fit_fixed = run_ic(false);
          double c_adaptive = fixed_objective(fit_ic.pose);
          double c_fixed = fixed_objective(fit_fixed.pose);
          if (c_fixed <= c_adaptive) {
            fit_ic = fit_fixed;
          } else {
            fit_ic.cost_final = c_adaptive;
          }
        }
        if (final_level && use_ic_polish && polish_iter > 0 && !level.polish_stencil.empty()) {
          auto run_ic_polish = [&](bool use_adaptive) {
            return optimize_frame_pose_ic(data,
                                          nx,
                                          ny,
                                          nz,
                                          nxyz,
                                          frame_t,
                                          pixdim,
                                          level.polish_stencil,
                                          fit_ic.pose,
                                          polish_iter,
                                          lambda0,
                                          tol_step,
                                          tol_cost,
                                          huber_k,
                                          use_adaptive,
                                          0,
                                          use_overlap_taper,
                                          false,
                                          capture_trans_mm,
                                          capture_rot_deg,
                                          0,
                                          level.polish_use_whiten,
                                          level.polish_use_whiten ? level.polish_wm.LinvT : nullptr,
                                          use_nuisance_marginalization,
                                          level.polish_nuisance_dim,
                                          nuisance_ridge,
                                          fit_nthreads);
          };
          auto polish_fixed_objective = [&](const Pose6& pose) {
            if (use_nuisance_marginalization && level.polish_nuisance_dim > 0) {
              return evaluate_cost_stencil_nuisance(data,
                                                    nx,
                                                    ny,
                                                    nz,
                                                    nxyz,
                                                    frame_t,
                                                    pixdim,
                                                    level.polish_stencil,
                                                    pose,
                                                    huber_k,
                                                    use_overlap_taper,
                                                    level.polish_nuisance_dim,
                                                    nuisance_ridge,
                                                    fit_nthreads)
                  .cost;
            }
            return evaluate_cost_stencil(data,
                                         nx,
                                         ny,
                                         nz,
                                         nxyz,
                                         frame_t,
                                         pixdim,
                                         level.polish_stencil,
                                         pose,
                                         huber_k,
                                         use_overlap_taper,
                                         fit_nthreads)
                .cost;
          };
          FrameFit fit_polish = run_ic_polish(adaptive_huber);
          if (adaptive_huber && adaptive_huber_guard) {
            FrameFit fit_fixed = run_ic_polish(false);
            double c_adaptive = polish_fixed_objective(fit_polish.pose);
            double c_fixed = polish_fixed_objective(fit_fixed.pose);
            if (c_fixed <= c_adaptive) {
              fit_polish = fit_fixed;
            } else {
              fit_polish.cost_final = c_adaptive;
            }
          }
          double c_polish_init = polish_fixed_objective(fit_ic.pose);
          double c_polish_final = polish_fixed_objective(fit_polish.pose);
          double c_primary_init = fixed_objective(fit_ic.pose);
          double c_primary_final = fixed_objective(fit_polish.pose);
          if (c_polish_final <= c_polish_init && c_primary_final <= c_primary_init) {
            fit_polish.restart_used = fit_ic.restart_used;
            fit_polish.iters += fit_ic.iters;
            fit_polish.flags |= fit_ic.flags;
            return fit_polish;
          }
        }
        if (final_level && use_dense_polish && polish_iter > 0 && !level.polish_samples.empty()) {
          FrameFit fit_polish = optimize_frame_pose(data,
                                                    template_vol,
                                                    nx,
                                                    ny,
                                                    nz,
                                                    nxyz,
                                                    frame_t,
                                                    pixdim,
                                                    level.polish_samples,
                                                    fit_ic.pose,
                                                    polish_iter,
                                                    lambda0,
                                                    tol_step,
                                                    tol_cost,
                                                    huber_k,
                                                    adaptive_huber,
                                                    0,
                                                    use_overlap_taper,
                                                    false,
                                                    capture_trans_mm,
                                                    capture_rot_deg,
                                                    0);
          fit_polish.restart_used = fit_ic.restart_used;
          fit_polish.iters += fit_ic.iters;
          fit_polish.flags |= fit_ic.flags;
          return fit_polish;
        }
        return fit_ic;
      }
      return optimize_frame_pose(data,
                                 template_vol,
                                 nx,
                                 ny,
                                 nz,
                                 nxyz,
                                 frame_t,
                                 pixdim,
                                 level.samples,
                                 init_pose,
                                 level.max_iter,
                                 lambda0,
                                 tol_step,
                                 tol_cost,
                                 huber_k,
                                 adaptive_huber,
                                 first_level ? n_restarts : 0,
                                 use_overlap_taper,
                                 first_level && capture_boost,
                                 capture_trans_mm,
                                 capture_rot_deg,
                                 first_level ? capture_topk : 0);
    };

    auto fit_one = [&](int frame_t, const Pose6& init_pose) {
      Pose6 cur = init_pose;
      FrameFit out = FrameFit{init_pose, NA_REAL, NA_REAL, NA_REAL, 0, false, 0};
      bool have_fit = false;
      for (size_t li = 0; li < levels.size(); ++li) {
        FrameFit fit = fit_level(levels[li], frame_t, cur, li == 0, li + 1 == levels.size());
        if (!have_fit) {
          out = fit;
          have_fit = true;
        } else {
          out.pose = fit.pose;
          out.cost_final = fit.cost_final;
          out.overlap = fit.overlap;
          out.iters += fit.iters;
          out.restart_used = out.restart_used || fit.restart_used;
          out.flags |= fit.flags;
        }
        cur = fit.pose;
      }
      return out;
    };

    std::vector<Pose6> seed_poses = poses;
    store_fit(ref_t, fit_one(ref_t, ref_seed));

    auto backward_seed = [&](int t) {
      Pose6 init = seed_poses[t + 1];
      if (t + 2 <= ref_t) {
        init = add_pose(seed_poses[t + 1], sub_pose(seed_poses[t + 1], seed_poses[t + 2]));
        momentum_backward += 1;
      }
      Vec3 com_t = center_of_mass_frame(data, nx, ny, nz, nxyz, t, pixdim, fit_mask);
      Pose6 com_seed = make_pose(com_ref.x - com_t.x,
                                 com_ref.y - com_t.y,
                                 com_ref.z - com_t.z,
                                 0.0,
                                 0.0,
                                 0.0);
      com_seed_count += 1;
      init = scale_pose(add_pose(init, com_seed), 0.5);
      if (blend_current_pose_seed) {
        init = scale_pose(add_pose(init, seed_poses[t]), 0.5);
        blended_seed_count += 1;
      }
      return init;
    };

    auto forward_seed = [&](int t) {
      Pose6 init = seed_poses[t - 1];
      if (t - 2 >= ref_t) {
        init = add_pose(seed_poses[t - 1], sub_pose(seed_poses[t - 1], seed_poses[t - 2]));
        momentum_forward += 1;
      }
      Vec3 com_t = center_of_mass_frame(data, nx, ny, nz, nxyz, t, pixdim, fit_mask);
      Pose6 com_seed = make_pose(com_ref.x - com_t.x,
                                 com_ref.y - com_t.y,
                                 com_ref.z - com_t.z,
                                 0.0,
                                 0.0,
                                 0.0);
      com_seed_count += 1;
      init = scale_pose(add_pose(init, com_seed), 0.5);
      if (blend_current_pose_seed) {
        init = scale_pose(add_pose(init, seed_poses[t]), 0.5);
        blended_seed_count += 1;
      }
      return init;
    };

    if (parallel_frame_pass) {
      std::vector<Pose6> seeds(nt, make_pose());
      std::vector<unsigned char> do_fit(nt, 0);
      for (int t = ref_t - 1; t >= 0; --t) {
        seeds[t] = backward_seed(t);
        do_fit[t] = 1;
      }
      for (int t = ref_t + 1; t < nt; ++t) {
        seeds[t] = forward_seed(t);
        do_fit[t] = 1;
      }

      std::vector<FrameFit> fits(nt, FrameFit{make_pose(), NA_REAL, NA_REAL, NA_REAL, 0, false, 0});
      parallel_for_range(static_cast<size_t>(nt), nthreads, [&](size_t begin, size_t end) {
        for (size_t tt = begin; tt < end; ++tt) {
          int t = static_cast<int>(tt);
          if (do_fit[t]) {
            fits[t] = fit_one(t, seeds[t]);
          }
        }
      });
      for (int t = 0; t < nt; ++t) {
        if (do_fit[t]) {
          store_fit(t, fits[t]);
          parallel_frame_fits += 1;
        }
      }
    } else {
      // Backward from reference
      for (int t = ref_t - 1; t >= 0; --t) {
        store_fit(t, fit_one(t, backward_seed(t)));
      }

      // Forward from reference
      for (int t = ref_t + 1; t < nt; ++t) {
        store_fit(t, fit_one(t, forward_seed(t)));
      }
    }
  };

  solve_sequence(ref_init, false);

  bool template_refresh = robust_template && nt > 1;
  if (control.containsElementNamed("template_refresh")) {
    template_refresh = as<bool>(control["template_refresh"]);
  }
  bool template_refresh_valid_only = false;
  if (control.containsElementNamed("template_refresh_valid_only")) {
    template_refresh_valid_only = as<bool>(control["template_refresh_valid_only"]);
  }

  if (template_refresh) {
        template_vol = build_aligned_trimmed_template(
            data, poses, nx, ny, nz, nt, nxyz, pixdim, trim_frac, template_refresh_valid_only, nthreads);
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
                      _["flags"] = flags,
                      _["warm_start"] = List::create(
                          _["passes"] = warm_start_passes,
                          _["momentum_backward"] = momentum_backward,
                          _["momentum_forward"] = momentum_forward,
                          _["momentum_seeds"] = momentum_backward + momentum_forward,
                          _["center_of_mass_seeds"] = com_seed_count,
                          _["blended_current_pose_seeds"] = blended_seed_count,
                          _["parallel_frame_fits"] = parallel_frame_fits),
                      _["pyramid"] = List::create(
                          _["downsample"] = wrap(downsample_levels),
                          _["max_iter"] = wrap(max_iter_levels),
                          _["n_samples"] = wrap(n_samples_levels)));
}
