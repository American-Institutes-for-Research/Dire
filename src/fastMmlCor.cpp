#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calcRrij(int i, int j, NumericMatrix rr1fi, NumericMatrix rr2fj, double detSigma, NumericVector mvnResid) {
  return (rr1fi(i,_) * rr2fj(j,_) * (2/3.141592653589793) * exp(-0.5 * (log(detSigma) + mvnResid))); 
}
