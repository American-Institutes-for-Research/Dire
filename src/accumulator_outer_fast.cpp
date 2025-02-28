#include <Rcpp.h>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericMatrix accumulator_outer_fast(const Rcpp::NumericMatrix& mat) {
  int n = mat.nrow();
  int k = mat.ncol();
  std::vector<double> result(k * k, 0.0);
  
  for (int row = 0; row < n; ++row) {
    for (int i = 0; i < k; ++i) {
      for (int j = 0; j < k; ++j) {
        result[i * k + j] += mat(row, i) * mat(row, j);
      }
    }
  }
  
  // Convert the std::vector back to NumericMatrix
  Rcpp::NumericMatrix result_mat(k, k);
  std::copy(result.begin(), result.end(), result_mat.begin());
  
  return result_mat;
}
