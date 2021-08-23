#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector IntersectSize(NumericVector x_start, NumericVector x_end, NumericVector y_start, NumericVector y_end) {
  int n = x_start.size();
  NumericVector out(n);

  for (int i = 0; i < n; i++) {
    if (y_start[i] > x_end[i] || x_start[i] > y_end[i]) {
      out[i] = 0;
    } else {
      out[i] = std::min(x_end[i], y_end[i]) - std::max(x_start[i], y_start[i]) + 1;
    }
  }

  return out;
}
