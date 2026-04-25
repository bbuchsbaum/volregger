#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <thread>
using namespace Rcpp;

namespace {

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

inline double sample_at_clamp(const double* frame,
                              int nx,
                              int ny,
                              int nz,
                              int i,
                              int j,
                              int k) {
  i = clamp_int(i, 0, nx - 1);
  j = clamp_int(j, 0, ny - 1);
  k = clamp_int(k, 0, nz - 1);
  size_t idx = static_cast<size_t>(i + nx * (j + ny * k));
  return frame[idx];
}

inline double sample_at_raw(const double* frame,
                            int nx,
                            int ny,
                            int i,
                            int j,
                            int k) {
  size_t idx = static_cast<size_t>(i + nx * (j + ny * k));
  return frame[idx];
}

inline double sample_at_zero(const double* frame,
                             int nx,
                             int ny,
                             int nz,
                             int i,
                             int j,
                             int k) {
  if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) {
    return 0.0;
  }
  size_t idx = static_cast<size_t>(i + nx * (j + ny * k));
  return frame[idx];
}

template <bool ZeroPad>
inline double sample_at_padded(const double* frame,
                               int nx,
                               int ny,
                               int nz,
                               int i,
                               int j,
                               int k) {
  if (ZeroPad) {
    return sample_at_zero(frame, nx, ny, nz, i, j, k);
  }
  return sample_at_clamp(frame, nx, ny, nz, i, j, k);
}

template <bool ZeroPad>
inline double sample_trilinear(const double* frame,
                               int nx,
                               int ny,
                               int nz,
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

  bool inside = x0 >= 0 && x1 < nx && y0 >= 0 && y1 < ny && z0 >= 0 && z1 < nz;

  double c000, c100, c010, c110, c001, c101, c011, c111;
  if (inside) {
    c000 = sample_at_raw(frame, nx, ny, x0, y0, z0);
    c100 = sample_at_raw(frame, nx, ny, x1, y0, z0);
    c010 = sample_at_raw(frame, nx, ny, x0, y1, z0);
    c110 = sample_at_raw(frame, nx, ny, x1, y1, z0);
    c001 = sample_at_raw(frame, nx, ny, x0, y0, z1);
    c101 = sample_at_raw(frame, nx, ny, x1, y0, z1);
    c011 = sample_at_raw(frame, nx, ny, x0, y1, z1);
    c111 = sample_at_raw(frame, nx, ny, x1, y1, z1);
  } else {
    c000 = sample_at_padded<ZeroPad>(frame, nx, ny, nz, x0, y0, z0);
    c100 = sample_at_padded<ZeroPad>(frame, nx, ny, nz, x1, y0, z0);
    c010 = sample_at_padded<ZeroPad>(frame, nx, ny, nz, x0, y1, z0);
    c110 = sample_at_padded<ZeroPad>(frame, nx, ny, nz, x1, y1, z0);
    c001 = sample_at_padded<ZeroPad>(frame, nx, ny, nz, x0, y0, z1);
    c101 = sample_at_padded<ZeroPad>(frame, nx, ny, nz, x1, y0, z1);
    c011 = sample_at_padded<ZeroPad>(frame, nx, ny, nz, x0, y1, z1);
    c111 = sample_at_padded<ZeroPad>(frame, nx, ny, nz, x1, y1, z1);
  }

  double c00 = c000 * (1.0 - ax) + c100 * ax;
  double c10 = c010 * (1.0 - ax) + c110 * ax;
  double c01 = c001 * (1.0 - ax) + c101 * ax;
  double c11 = c011 * (1.0 - ax) + c111 * ax;

  double c0 = c00 * (1.0 - ay) + c10 * ay;
  double c1 = c01 * (1.0 - ay) + c11 * ay;

  return c0 * (1.0 - az) + c1 * az;
}

inline double cubic_weight(double x) {
  // Keys cubic kernel with a = -0.5 (Catmull-Rom style)
  const double a = -0.5;
  x = std::abs(x);
  if (x < 1.0) {
    return (a + 2.0) * x * x * x - (a + 3.0) * x * x + 1.0;
  }
  if (x < 2.0) {
    return a * x * x * x - 5.0 * a * x * x + 8.0 * a * x - 4.0 * a;
  }
  return 0.0;
}

template <bool ZeroPad>
inline double sample_cubic(const double* frame,
                           int nx,
                           int ny,
                           int nz,
                           double x,
                           double y,
                           double z) {
  int x0 = static_cast<int>(std::floor(x)) - 1;
  int y0 = static_cast<int>(std::floor(y)) - 1;
  int z0 = static_cast<int>(std::floor(z)) - 1;

  int xs[4], ys[4], zs[4];
  double wxs[4], wys[4], wzs[4];
  for (int ii = 0; ii < 4; ++ii) {
    xs[ii] = x0 + ii;
    wxs[ii] = cubic_weight(x - xs[ii]);
  }
  for (int jj = 0; jj < 4; ++jj) {
    ys[jj] = y0 + jj;
    wys[jj] = cubic_weight(y - ys[jj]);
  }
  for (int kk = 0; kk < 4; ++kk) {
    zs[kk] = z0 + kk;
    wzs[kk] = cubic_weight(z - zs[kk]);
  }

  bool inside = xs[0] >= 0 && xs[3] < nx &&
                ys[0] >= 0 && ys[3] < ny &&
                zs[0] >= 0 && zs[3] < nz;
  double out = 0.0;
  double wsum = 0.0;
  if (inside) {
    for (int kk = 0; kk < 4; ++kk) {
      int zk = zs[kk];
      double wz = wzs[kk];
      for (int jj = 0; jj < 4; ++jj) {
        int yj = ys[jj];
        double wy = wys[jj];
        for (int ii = 0; ii < 4; ++ii) {
          int xi = xs[ii];
          double wx = wxs[ii];
          double w = wx * wy * wz;
          out += w * sample_at_raw(frame, nx, ny, xi, yj, zk);
          wsum += w;
        }
      }
    }
  } else {
    for (int kk = 0; kk < 4; ++kk) {
      int zk = zs[kk];
      double wz = wzs[kk];
      for (int jj = 0; jj < 4; ++jj) {
        int yj = ys[jj];
        double wy = wys[jj];
        for (int ii = 0; ii < 4; ++ii) {
          int xi = xs[ii];
          double wx = wxs[ii];
          double w = wx * wy * wz;
          out += w * sample_at_padded<ZeroPad>(frame, nx, ny, nz, xi, yj, zk);
          wsum += w;
        }
      }
    }
  }
  if (wsum == 0.0) {
    return 0.0;
  }
  return out / wsum;
}

inline double sinc_pi(double x) {
  if (std::abs(x) < 1e-8) {
    return 1.0;
  }
  double px = M_PI * x;
  return std::sin(px) / px;
}

inline double lanczos_weight(double x, int a) {
  x = std::abs(x);
  if (x >= a) {
    return 0.0;
  }
  return sinc_pi(x) * sinc_pi(x / static_cast<double>(a));
}

template <bool ZeroPad>
inline double sample_lanczos(const double* frame,
                             int nx,
                             int ny,
                             int nz,
                             double x,
                             double y,
                             double z,
                             int a = 3) {
  int x0 = static_cast<int>(std::floor(x)) - a + 1;
  int y0 = static_cast<int>(std::floor(y)) - a + 1;
  int z0 = static_cast<int>(std::floor(z)) - a + 1;
  int width = 2 * a;

  int xs[14], ys[14], zs[14];
  double wxs[14], wys[14], wzs[14];
  for (int ii = 0; ii < width; ++ii) {
    xs[ii] = x0 + ii;
    wxs[ii] = lanczos_weight(x - xs[ii], a);
  }
  for (int jj = 0; jj < width; ++jj) {
    ys[jj] = y0 + jj;
    wys[jj] = lanczos_weight(y - ys[jj], a);
  }
  for (int kk = 0; kk < width; ++kk) {
    zs[kk] = z0 + kk;
    wzs[kk] = lanczos_weight(z - zs[kk], a);
  }

  bool inside = xs[0] >= 0 && xs[width - 1] < nx &&
                ys[0] >= 0 && ys[width - 1] < ny &&
                zs[0] >= 0 && zs[width - 1] < nz;
  double out = 0.0;
  double wsum = 0.0;
  if (inside) {
    for (int kk = 0; kk < width; ++kk) {
      int zk = zs[kk];
      double wz = wzs[kk];
      if (wz == 0.0) continue;
      for (int jj = 0; jj < width; ++jj) {
        int yj = ys[jj];
        double wy = wys[jj];
        if (wy == 0.0) continue;
        for (int ii = 0; ii < width; ++ii) {
          int xi = xs[ii];
          double wx = wxs[ii];
          if (wx == 0.0) continue;
          double w = wx * wy * wz;
          out += w * sample_at_raw(frame, nx, ny, xi, yj, zk);
          wsum += w;
        }
      }
    }
  } else {
    for (int kk = 0; kk < width; ++kk) {
      int zk = zs[kk];
      double wz = wzs[kk];
      if (wz == 0.0) continue;
      for (int jj = 0; jj < width; ++jj) {
        int yj = ys[jj];
        double wy = wys[jj];
        if (wy == 0.0) continue;
        for (int ii = 0; ii < width; ++ii) {
          int xi = xs[ii];
          double wx = wxs[ii];
          if (wx == 0.0) continue;
          double w = wx * wy * wz;
          out += w * sample_at_padded<ZeroPad>(frame, nx, ny, nz, xi, yj, zk);
          wsum += w;
        }
      }
    }
  }
  if (wsum == 0.0) {
    return 0.0;
  }
  return out / wsum;
}

inline void fixed_to_moving_vox(int i,
                                int j,
                                int k,
                                int nx,
                                int ny,
                                int nz,
                                const NumericVector& pixdim,
                                double tx,
                                double ty,
                                double tz,
                                double rx,
                                double ry,
                                double rz,
                                double* src_xyz) {
  double cx = 0.5 * (nx - 1.0);
  double cy = 0.5 * (ny - 1.0);
  double cz = 0.5 * (nz - 1.0);

  double xw = (i - cx) * pixdim[0];
  double yw = (j - cy) * pixdim[1];
  double zw = (k - cz) * pixdim[2];

  double xf = xw - tx;
  double yf = yw - ty;
  double zf = zw - tz;

  double R[9];
  rotation_matrix(rx, ry, rz, R);

  // Inverse map by R^T
  double xm = R[0] * xf + R[3] * yf + R[6] * zf;
  double ym = R[1] * xf + R[4] * yf + R[7] * zf;
  double zm = R[2] * xf + R[5] * yf + R[8] * zf;

  src_xyz[0] = xm / pixdim[0] + cx;
  src_xyz[1] = ym / pixdim[1] + cy;
  src_xyz[2] = zm / pixdim[2] + cz;
}

inline void fixed_to_moving_vox_zpad(int i,
                                     int j,
                                     int k,
                                     int nx,
                                     int ny,
                                     int nz,
                                     int zpad,
                                     const NumericVector& pixdim,
                                     double tx,
                                     double ty,
                                     double tz,
                                     double rx,
                                     double ry,
                                     double rz,
                                     double* src_xyz) {
  fixed_to_moving_vox(i + zpad,
                      j + zpad,
                      k + zpad,
                      nx + 2 * zpad,
                      ny + 2 * zpad,
                      nz + 2 * zpad,
                      pixdim,
                      tx,
                      ty,
                      tz,
                      rx,
                      ry,
                      rz,
                      src_xyz);
  src_xyz[0] -= zpad;
  src_xyz[1] -= zpad;
  src_xyz[2] -= zpad;
}

}  // namespace

// [[Rcpp::export]]
SEXP cpp_volreg_apply(
    NumericVector data,
    IntegerVector dim,
    NumericVector pixdim,
    NumericMatrix par,
    Nullable<NumericVector> slice_times = R_NilValue,
    Nullable<IntegerVector> mb_groups = R_NilValue,
    std::string interp = "bspline",
    int pad_mode = 0,
    int zpad = 0,
    int nthreads = 1) {
  if (dim.size() != 4) {
    stop("dim must be length 4");
  }
  if (par.nrow() != dim[3] || par.ncol() != 6) {
    stop("par must be nt x 6");
  }

  int nx = dim[0];
  int ny = dim[1];
  int nz = dim[2];
  int nt = dim[3];
  int nxyz = nx * ny * nz;
  zpad = std::max(0, zpad);
  if (zpad > 0) {
    // AFNI-style zpad means resampling a symmetrically zero-padded volume
    // and cropping back to the original grid.
    pad_mode = 1;
  }

  enum InterpMode { LINEAR = 0, BSPLINE = 1, LANCZOS = 2 };
  InterpMode mode = BSPLINE;
  if (interp == "linear") {
    mode = LINEAR;
  } else if (interp == "bspline") {
    mode = BSPLINE;
  } else if (interp == "lanczos") {
    mode = LANCZOS;
  } else {
    stop("interp must be one of: linear, bspline, lanczos");
  }

  NumericVector out(data.size());
  out.attr("dim") = dim;
  const double* data_ptr = data.begin();
  double* out_ptr = out.begin();
  bool zero_pad = pad_mode == 1;

  std::vector<double> slice_offsets(nz, 0.0);
  bool use_packet_pose = false;

  if (slice_times.isNotNull()) {
    NumericVector st(slice_times);
    if (st.size() > 0) {
      use_packet_pose = true;
      if (st.size() == nz) {
        for (int k = 0; k < nz; ++k) {
          slice_offsets[k] = st[k];
        }
      } else if (mb_groups.isNotNull()) {
        IntegerVector mg(mb_groups);
        if (mg.size() == nz) {
          for (int k = 0; k < nz; ++k) {
            int g = mg[k];
            int idx = g;
            if (idx < 0) idx = 0;
            if (idx >= st.size()) idx = st.size() - 1;
            slice_offsets[k] = st[idx];
          }
        }
      } else {
        for (int k = 0; k < nz; ++k) {
          slice_offsets[k] = st[k % st.size()];
        }
      }
    }
  } else if (mb_groups.isNotNull()) {
    IntegerVector mg(mb_groups);
    if (mg.size() == nz) {
      use_packet_pose = true;
      std::set<int> uniq;
      for (int k = 0; k < nz; ++k) {
        uniq.insert(mg[k]);
      }
      std::vector<int> gids(uniq.begin(), uniq.end());
      std::map<int, double> rank;
      for (size_t i = 0; i < gids.size(); ++i) {
        double denom = std::max(1.0, static_cast<double>(gids.size() - 1));
        rank[gids[i]] = static_cast<double>(i) / denom;
      }
      for (int k = 0; k < nz; ++k) {
        slice_offsets[k] = rank[mg[k]];
      }
    }
  }

  if (use_packet_pose) {
    double mn = *std::min_element(slice_offsets.begin(), slice_offsets.end());
    double mx = *std::max_element(slice_offsets.begin(), slice_offsets.end());
    double den = mx - mn;
    if (den > 0.0) {
      for (int k = 0; k < nz; ++k) {
        slice_offsets[k] = (slice_offsets[k] - mn) / den;
      }
    } else {
      std::fill(slice_offsets.begin(), slice_offsets.end(), 0.0);
    }
  }

  size_t tasks = static_cast<size_t>(nt) * static_cast<size_t>(nz);
  parallel_for_range(tasks, nthreads, [&](size_t begin, size_t end) {
    double src[3];
    for (size_t task = begin; task < end; ++task) {
      int t = static_cast<int>(task / static_cast<size_t>(nz));
      int k = static_cast<int>(task % static_cast<size_t>(nz));
      size_t off = static_cast<size_t>(nxyz) * static_cast<size_t>(t);
      const double* frame = data_ptr + off;

      double tx = par(t, 0);
      double ty = par(t, 1);
      double tz = par(t, 2);
      double rx = par(t, 3);
      double ry = par(t, 4);
      double rz = par(t, 5);

      if (use_packet_pose && nt > 1) {
        double u = static_cast<double>(t) + slice_offsets[k];
        int t0 = static_cast<int>(std::floor(u));
        t0 = clamp_int(t0, 0, nt - 1);
        int t1 = clamp_int(t0 + 1, 0, nt - 1);
        double a = std::min(1.0, std::max(0.0, u - t0));

        tx = (1.0 - a) * par(t0, 0) + a * par(t1, 0);
        ty = (1.0 - a) * par(t0, 1) + a * par(t1, 1);
        tz = (1.0 - a) * par(t0, 2) + a * par(t1, 2);
        rx = (1.0 - a) * par(t0, 3) + a * par(t1, 3);
        ry = (1.0 - a) * par(t0, 4) + a * par(t1, 4);
        rz = (1.0 - a) * par(t0, 5) + a * par(t1, 5);
      }

      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          if (zpad > 0) {
            fixed_to_moving_vox_zpad(i, j, k, nx, ny, nz, zpad, pixdim, tx, ty, tz, rx, ry, rz, src);
          } else {
            fixed_to_moving_vox(i, j, k, nx, ny, nz, pixdim, tx, ty, tz, rx, ry, rz, src);
          }
          double v = 0.0;
          if (mode == LINEAR) {
            if (zero_pad) {
              v = sample_trilinear<true>(frame, nx, ny, nz, src[0], src[1], src[2]);
            } else {
              v = sample_trilinear<false>(frame, nx, ny, nz, src[0], src[1], src[2]);
            }
          } else if (mode == BSPLINE) {
            if (zero_pad) {
              v = sample_cubic<true>(frame, nx, ny, nz, src[0], src[1], src[2]);
            } else {
              v = sample_cubic<false>(frame, nx, ny, nz, src[0], src[1], src[2]);
            }
          } else {
            if (zero_pad) {
              v = sample_lanczos<true>(frame, nx, ny, nz, src[0], src[1], src[2], 3);
            } else {
              v = sample_lanczos<false>(frame, nx, ny, nz, src[0], src[1], src[2], 3);
            }
          }
          out_ptr[off + static_cast<size_t>(i + nx * (j + ny * k))] = v;
        }
      }
    }
  });

  out.attr("dim") = dim;
  return out;
}
