#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List cpp_fit_pose_spline(
    NumericMatrix par,
    double tr,
    Nullable<NumericVector> slice_times = R_NilValue,
    Nullable<IntegerVector> mb_groups = R_NilValue,
    List control = List::create()) {
  int nt = par.nrow();
  if (par.ncol() != 6) {
    stop("par must be nt x 6");
  }

  int smooth_iter = 2;
  if (control.containsElementNamed("spline_smooth_iter")) {
    smooth_iter = as<int>(control["spline_smooth_iter"]);
    smooth_iter = std::max(0, smooth_iter);
  }

  NumericMatrix sm = clone(par);
  NumericMatrix tmp(nt, 6);

  for (int iter = 0; iter < smooth_iter; ++iter) {
    for (int c = 0; c < 6; ++c) {
      tmp(0, c) = sm(0, c);
      for (int t = 1; t < nt - 1; ++t) {
        // cubic B-spline style 1D smoothing stencil
        tmp(t, c) = (sm(t - 1, c) + 4.0 * sm(t, c) + sm(t + 1, c)) / 6.0;
      }
      if (nt > 1) {
        tmp(nt - 1, c) = sm(nt - 1, c);
      }
    }
    sm = clone(tmp);
  }

  NumericVector times(nt);
  for (int t = 0; t < nt; ++t) {
    times[t] = t * tr;
  }

  NumericVector packet_offsets;
  if (slice_times.isNotNull()) {
    packet_offsets = as<NumericVector>(slice_times);
  } else if (mb_groups.isNotNull()) {
    IntegerVector mg = as<IntegerVector>(mb_groups);
    packet_offsets = NumericVector(mg.size());
    std::vector<int> uniq = as<std::vector<int>>(sort_unique(mg));
    for (int k = 0; k < mg.size(); ++k) {
      auto it = std::find(uniq.begin(), uniq.end(), mg[k]);
      double rank = static_cast<double>(std::distance(uniq.begin(), it));
      double denom = std::max(1.0, static_cast<double>(uniq.size() - 1));
      packet_offsets[k] = rank / denom;
    }
  } else {
    packet_offsets = NumericVector::create();
  }

  return List::create(
      _["par_smooth"] = sm,
      _["times"] = times,
      _["tr"] = tr,
      _["packet_offsets"] = packet_offsets,
      _["method"] = "cubic_bspline_smoother");
}
