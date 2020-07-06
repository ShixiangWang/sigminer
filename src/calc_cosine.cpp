#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cosineVector(NumericVector x, NumericVector y) {
  return sum(x * y) / (sqrt(sum(pow(x, 2))) * sqrt(sum(pow(y, 2))));
}

// [[Rcpp::export]]
NumericMatrix cosineMatrix(NumericMatrix x, NumericMatrix y) {
  int nrow = x.ncol(), ncol = y.ncol();
  NumericMatrix out(nrow, ncol); // output a cosine value matrix

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      out(i, j) = cosineVector(x(_, i), y(_, j));
    }
  }

  return out;
}
