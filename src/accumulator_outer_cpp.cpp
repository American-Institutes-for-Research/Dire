#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix accumulator_outer_cpp(NumericMatrix mat) {
  int n = mat.nrow();
  int k = mat.ncol();
  NumericMatrix result(k, k);
  
  for(int row = 0; row < n; row++) {
    for(int i = 0; i < k; i++) {
      for(int j = 0; j < k; j++) {
        result(i, j) += mat(row, i) * mat(row, j);
      }
    }
  }
  
  return result;
}
