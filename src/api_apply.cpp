#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
using namespace Rcpp;

namespace {

inline int clamp_int(int x, int lo, int hi) {
  return std::max(lo, std::min(hi, x));
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

inline double sample_trilinear(const NumericVector& data,
                               int nx,
                               int ny,
                               int nz,
                               int nxyz,
                               int t,
                               double x,
                               double y,
                               double z,
                               int pad_mode) {
  auto sample_at = [&](int i, int j, int k) {
    if (pad_mode == 1) {
      // zero padding
      if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) {
        return 0.0;
      }
    }
    i = clamp_int(i, 0, nx - 1);
    j = clamp_int(j, 0, ny - 1);
    k = clamp_int(k, 0, nz - 1);
    size_t off = static_cast<size_t>(nxyz) * t;
    size_t idx = static_cast<size_t>(i + nx * (j + ny * k));
    return data[off + idx];
  };

  int x0 = static_cast<int>(std::floor(x));
  int y0 = static_cast<int>(std::floor(y));
  int z0 = static_cast<int>(std::floor(z));
  int x1 = x0 + 1;
  int y1 = y0 + 1;
  int z1 = z0 + 1;

  double ax = x - x0;
  double ay = y - y0;
  double az = z - z0;

  double c000 = sample_at(x0, y0, z0);
  double c100 = sample_at(x1, y0, z0);
  double c010 = sample_at(x0, y1, z0);
  double c110 = sample_at(x1, y1, z0);
  double c001 = sample_at(x0, y0, z1);
  double c101 = sample_at(x1, y0, z1);
  double c011 = sample_at(x0, y1, z1);
  double c111 = sample_at(x1, y1, z1);

  double c00 = c000 * (1.0 - ax) + c100 * ax;
  double c10 = c010 * (1.0 - ax) + c110 * ax;
  double c01 = c001 * (1.0 - ax) + c101 * ax;
  double c11 = c011 * (1.0 - ax) + c111 * ax;

  double c0 = c00 * (1.0 - ay) + c10 * ay;
  double c1 = c01 * (1.0 - ay) + c11 * ay;

  return c0 * (1.0 - az) + c1 * az;
}

inline double sample_at(const NumericVector& data,
                        int nx,
                        int ny,
                        int nz,
                        int nxyz,
                        int t,
                        int i,
                        int j,
                        int k,
                        int pad_mode) {
  if (pad_mode == 1) {
    if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) {
      return 0.0;
    }
  }
  i = clamp_int(i, 0, nx - 1);
  j = clamp_int(j, 0, ny - 1);
  k = clamp_int(k, 0, nz - 1);
  size_t off = static_cast<size_t>(nxyz) * t;
  size_t idx = static_cast<size_t>(i + nx * (j + ny * k));
  return data[off + idx];
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

inline double sample_cubic(const NumericVector& data,
                           int nx,
                           int ny,
                           int nz,
                           int nxyz,
                           int t,
                           double x,
                           double y,
                           double z,
                           int pad_mode) {
  int x0 = static_cast<int>(std::floor(x)) - 1;
  int y0 = static_cast<int>(std::floor(y)) - 1;
  int z0 = static_cast<int>(std::floor(z)) - 1;

  double out = 0.0;
  double wsum = 0.0;
  for (int kk = 0; kk < 4; ++kk) {
    int zk = z0 + kk;
    double wz = cubic_weight(z - zk);
    for (int jj = 0; jj < 4; ++jj) {
      int yj = y0 + jj;
      double wy = cubic_weight(y - yj);
      for (int ii = 0; ii < 4; ++ii) {
        int xi = x0 + ii;
        double wx = cubic_weight(x - xi);
        double w = wx * wy * wz;
        out += w * sample_at(data, nx, ny, nz, nxyz, t, xi, yj, zk, pad_mode);
        wsum += w;
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

inline double sample_lanczos(const NumericVector& data,
                             int nx,
                             int ny,
                             int nz,
                             int nxyz,
                             int t,
                             double x,
                             double y,
                             double z,
                             int pad_mode,
                             int a = 3) {
  int x0 = static_cast<int>(std::floor(x)) - a + 1;
  int y0 = static_cast<int>(std::floor(y)) - a + 1;
  int z0 = static_cast<int>(std::floor(z)) - a + 1;

  double out = 0.0;
  double wsum = 0.0;
  for (int kk = 0; kk < 2 * a; ++kk) {
    int zk = z0 + kk;
    double wz = lanczos_weight(z - zk, a);
    if (wz == 0.0) continue;
    for (int jj = 0; jj < 2 * a; ++jj) {
      int yj = y0 + jj;
      double wy = lanczos_weight(y - yj, a);
      if (wy == 0.0) continue;
      for (int ii = 0; ii < 2 * a; ++ii) {
        int xi = x0 + ii;
        double wx = lanczos_weight(x - xi, a);
        if (wx == 0.0) continue;
        double w = wx * wy * wz;
        out += w * sample_at(data, nx, ny, nz, nxyz, t, xi, yj, zk, pad_mode);
        wsum += w;
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

  for (int t = 0; t < nt; ++t) {
    size_t off = static_cast<size_t>(nxyz) * t;
    for (int k = 0; k < nz; ++k) {
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
          double src[3];
          fixed_to_moving_vox(i, j, k, nx, ny, nz, pixdim, tx, ty, tz, rx, ry, rz, src);
          double v = 0.0;
          if (mode == LINEAR) {
            v = sample_trilinear(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2], pad_mode);
          } else if (mode == BSPLINE) {
            v = sample_cubic(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2], pad_mode);
          } else {
            v = sample_lanczos(data, nx, ny, nz, nxyz, t, src[0], src[1], src[2], pad_mode, 3);
          }
          out[off + static_cast<size_t>(i + nx * (j + ny * k))] = v;
        }
      }
    }
  }

  // Accepted in API but not yet used for threading.
  (void) nthreads;

  out.attr("dim") = dim;
  return out;
}
