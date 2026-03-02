#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace {

inline size_t idx4(int i, int j, int k, int t, int nx, int ny, int nz) {
  return static_cast<size_t>(i + nx * (j + ny * (k + nz * t)));
}

NumericVector compute_fd(const NumericMatrix& par, double radius = 50.0) {
  int nt = par.nrow();
  NumericVector fd(nt, 0.0);
  for (int t = 1; t < nt; ++t) {
    double v = 0.0;
    for (int c = 0; c < 3; ++c) {
      v += std::abs(par(t, c) - par(t - 1, c));
    }
    for (int c = 3; c < 6; ++c) {
      v += radius * std::abs(par(t, c) - par(t - 1, c));
    }
    fd[t] = v;
  }
  return fd;
}

NumericVector compute_dvars(const NumericVector& vol4d,
                            int nx,
                            int ny,
                            int nz,
                            int nt,
                            const std::vector<uint8_t>& maskv) {
  NumericVector out(nt, NA_REAL);
  if (nt < 2) {
    return out;
  }

  for (int t = 1; t < nt; ++t) {
    double ss = 0.0;
    double n = 0.0;
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          size_t vidx = static_cast<size_t>(i + nx * (j + ny * k));
          if (!maskv[vidx]) {
            continue;
          }
          double a = vol4d[idx4(i, j, k, t - 1, nx, ny, nz)];
          double b = vol4d[idx4(i, j, k, t, nx, ny, nz)];
          double d = b - a;
          ss += d * d;
          n += 1.0;
        }
      }
    }
    if (n > 0.0) {
      out[t] = std::sqrt(ss / n);
    }
  }

  return out;
}

NumericVector robust_clip_3x_median(const NumericVector& x) {
  std::vector<double> vals;
  vals.reserve(x.size());
  for (int i = 0; i < x.size(); ++i) {
    if (R_finite(x[i])) {
      vals.push_back(x[i]);
    }
  }
  if (vals.empty()) {
    return clone(x);
  }

  std::nth_element(vals.begin(), vals.begin() + vals.size() / 2, vals.end());
  double med = vals[vals.size() / 2];
  double cap = 3.0 * med;

  NumericVector out = clone(x);
  for (int i = 0; i < out.size(); ++i) {
    if (R_finite(out[i]) && out[i] > cap) {
      out[i] = cap;
    }
  }
  return out;
}

std::vector<uint8_t> as_mask(const Nullable<LogicalVector>& mask, int nxyz) {
  std::vector<uint8_t> m(nxyz, 1);
  if (mask.isNotNull()) {
    LogicalVector in(mask);
    if (in.size() != nxyz) {
      stop("mask size mismatch");
    }
    for (int i = 0; i < nxyz; ++i) {
      m[i] = in[i] == TRUE ? 1 : 0;
    }
  }
  return m;
}

}  // namespace

// [[Rcpp::export]]
List cpp_volreg_metrics(
    NumericVector original,
    Nullable<NumericVector> corrected,
    IntegerVector dim,
    Nullable<LogicalVector> mask,
    NumericMatrix par,
    NumericVector cost_init,
    NumericVector cost_final) {
  if (dim.size() != 4) {
    stop("dim must be length 4");
  }
  int nx = dim[0];
  int ny = dim[1];
  int nz = dim[2];
  int nt = dim[3];
  int nxyz = nx * ny * nz;

  if (par.nrow() != nt || par.ncol() != 6) {
    stop("par must be nt x 6");
  }

  NumericVector vol4d = original;
  if (corrected.isNotNull()) {
    vol4d = NumericVector(corrected);
  }
  if (vol4d.size() != static_cast<R_xlen_t>(nxyz) * nt) {
    stop("volume size mismatch");
  }

  if (cost_init.size() != nt || cost_final.size() != nt) {
    stop("cost vectors must have length nt");
  }

  std::vector<uint8_t> maskv = as_mask(mask, nxyz);

  NumericVector fd = compute_fd(par, 50.0);
  NumericVector dvars = compute_dvars(vol4d, nx, ny, nz, nt, maskv);
  NumericVector robust_dvars = robust_clip_3x_median(dvars);

  LogicalVector motion_spike(nt, false);
  LogicalVector fit_failure(nt, false);
  LogicalVector censor_suggest(nt, false);
  NumericVector cost_drop(nt, 0.0);

  for (int t = 0; t < nt; ++t) {
    motion_spike[t] = (fd[t] > 0.5);
    fit_failure[t] = (cost_final[t] > cost_init[t] + 1e-8);
    censor_suggest[t] = (motion_spike[t] || fit_failure[t]);
    cost_drop[t] = cost_init[t] - cost_final[t];
  }

  return List::create(
      _["fd"] = fd,
      _["dvars"] = dvars,
      _["robust_dvars"] = robust_dvars,
      _["motion_spike"] = motion_spike,
      _["fit_failure"] = fit_failure,
      _["censor_suggest"] = censor_suggest,
      _["cost_drop"] = cost_drop);
}
