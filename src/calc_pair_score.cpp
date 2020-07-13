#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int pairScoreVector(NumericVector x, NumericVector y, int x_max, int y_max) {
  // x[0] and y[0] represent length value
  // x[1] and y[1] represent segment copy number value
  if (x[0] == y[0] && x[1] == y[1]) {
    return x_max * y_max;
  } else {
    int flag = 1;
    if (((x[1] - 2) > 0) ^ ((y[1] - 2) > 0)) {
      flag = -1;
    }
    return flag * (x_max - abs(x[0] - y[0])) * (y_max - abs(x[1] - y[1]));
  }
}

// [[Rcpp::export]]
NumericMatrix pairScoreMatrix(NumericMatrix x, NumericMatrix y, int x_max, int y_max) {
  int nrow = x.nrow(), ncol = y.nrow();
  NumericMatrix out(nrow, ncol);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      out(i, j) = pairScoreVector(x(i, _), y(j, _), x_max, y_max);
    }
  }

  return out;
}

// [[Rcpp::export]]
IntegerMatrix getScoreMatrix(IntegerMatrix indexMat, IntegerMatrix subMat) {
  // indexMat: each row represent the index in subMat (0-based)
  // subMat: a matrix stores match score
  int n = indexMat.nrow(), size = indexMat.ncol();
  IntegerMatrix out(n);
  int score = 0;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      score = 0;
      Rcpp::Rcout << "i:" << i << " j:" << j << std::endl;
      for (int k = 0; k < size; k++) {
        Rcpp::Rcout << "  score index to plus:" << indexMat(i, k) << "," << indexMat(j, k) << std::endl;
        Rcpp::Rcout << "  score to plus:" << subMat(indexMat(i, k), indexMat(j, k)) << std::endl;
        score += subMat(indexMat(i, k), indexMat(j, k));
      }
      out(i, j) = out(j, i) = score;
    }
  }

  return out;
}