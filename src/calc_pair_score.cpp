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
      // In this case, switch copy number 0 and 2 to
      // avoid result of copy number 0 is greater than
      // copy number 2
      //Rcpp::Rcout << x[1] << "," << y[1] << ";";
      if (x[1] == 0) {
        x[1] = 2;
      } else if (x[1] == 2) {
        x[1] = 0;
      }
      if (y[1] == 0) {
        y[1] = 2;
      } else if (y[1] == 2) {
        y[1] = 0;
      }
      //Rcpp::Rcout << x[1] << "," << y[1] << std::endl;
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
IntegerMatrix getScoreMatrix(IntegerMatrix indexMat, IntegerMatrix subMat, int bSize, bool verbose) {
  // indexMat: each row represent the index in subMat (0-based)
  // subMat: a matrix stores match score
  // bSize: block size to aggregrate
  int n = indexMat.nrow(), size = indexMat.ncol();
  int score = 0;

  if (bSize == 1) {
    IntegerMatrix out(n);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= i; j++) {
        score = 0;
        if (verbose) {
          Rcpp::Rcout << "Handling index pair (" << i+1 << "," << j+1 << ")." << std::endl;
        }
        for (int k = 0; k < size; k++) {
          // if (verbose) {
          //   Rcpp::Rcout << "  substitution matrix index to plus:" << indexMat(i, k)+1 << "," << indexMat(j, k)+1 << std::endl;
          //   Rcpp::Rcout << "  score to plus:" << subMat(indexMat(i, k), indexMat(j, k)) << std::endl;
          // }
          score += subMat(indexMat(i, k), indexMat(j, k));
        }
        out(i, j) = score;
        if (i != j) {
          out(j, i) = score;
        }
      }
    }

    return out;
  } else {
    int chunkSize = (n / bSize) + 1;
    IntegerMatrix out(chunkSize);
    int blockScore = 0;
    int eCounter = 0; // element counter to aggregrate

    if (verbose) {
      Rcpp::Rcout << "Running with block size: " << bSize << std::endl;
    }
    // Use ii & jj to represent block index
    // Use i & j to represent matrix (indexMat) index
    for (int ii = 0; ii < chunkSize; ii++) {
      for (int jj = 0; jj <= ii; jj++) {
        if (verbose) {
          Rcpp::Rcout << "Handling block pair (" << ii+1 << "," << jj+1 << ")." << std::endl;
        }
        // Each block
        // NOTE the indices of last block
        blockScore = 0;
        eCounter = 0;
        for (int i = 0 + ii * bSize; i < std::min((ii + 1) * bSize, n); i++) {
          for (int j = 0 + jj * bSize; j < std::min((jj + 1) * bSize, n); j++) {
            // Each element
            score = 0;
            for (int k = 0; k < size; k++) {
              score += subMat(indexMat(i, k), indexMat(j, k));
            }
            blockScore += score;
            eCounter++;
          }
        }

        // Calculate the mean
        out(ii, jj) = std::round(blockScore / float(eCounter));
        if (ii != jj) {
          out(jj, ii) = std::round(blockScore / float(eCounter));
        }
      }
    }

    return out;
  }
}

// [[Rcpp::export]]
IntegerMatrix getScoreMatrixRect(IntegerMatrix indexMat1, IntegerMatrix indexMat2, IntegerMatrix subMat, bool verbose) {
  // For not equal matrices, indexMat1 should be a bigger matrix than indexMat2
  // indexMat: each row represent the index in subMat (0-based)
  // subMat: a matrix stores match score
  int m = indexMat1.nrow(), n = indexMat2.nrow();
  int size = indexMat1.ncol();
  int score = 0;
  IntegerMatrix out(m, n);

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      score = 0;
      if (verbose) {
        Rcpp::Rcout << "Handling index pair (" << i+1 << "," << j+1 << ")." << std::endl;
      }
      for (int k = 0; k < size; k++) {
        score += subMat(indexMat1(i, k), indexMat2(j, k));
      }
      out(i, j) = score;
    }
  }

  return out;
}
