#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]] 
arma::mat matTimesVec(arma::mat mat, arma::vec v) {
  mat.each_col() %= v;
  return mat;
}

#include <thread>

// [[Rcpp::export]] 
arma::mat calcHess(int K, arma::mat rr2, arma::mat rr2_, arma::mat trr2mxb, arma::mat X_, arma::mat nodesminusXB, arma::vec w, double s2, double s_) {
  arma::mat H;
  H.zeros(K+1, K+1);
  arma::mat rr2T = rr2.t();
  arma::mat rr2T_s2 = rr2T / s2;
  arma::rowvec denom = sum(rr2, 0);
  arma::rowvec denom2 = pow(denom, 2);
  arma::rowvec colSums_rr2_ = arma::sum(rr2_, 0);
  arma::mat rr2_nodsminusXB_t = (rr2_ % nodesminusXB).t();
  arma::mat trr2mxb_s2 = trr2mxb / s2;

  arma::cube fj_cube;
  
  for (int j=0; j < K; j++) {
    arma::mat Xj = X_.col(j);
    arma::mat fj_temp = matTimesVec(trr2mxb_s2, Xj).t();
    if (fj_cube.is_empty()) {
      fj_cube.set_size(fj_temp.n_rows, fj_temp.n_cols, 1);
      fj_cube.slice(0) = fj_temp;
    } else {
      fj_cube.insert_slices(fj_cube.n_slices, fj_temp);
    }
    
  }

  for (int i=0; i < K; i++) {
    arma::mat Xi = X_.col(i);
    arma::mat fi = matTimesVec(trr2mxb_s2, Xi).t();
    arma::mat rr2T_s2_Xi = matTimesVec(rr2T_s2, Xi);
    arma::rowvec num = arma::sum(fi, 0);
    for (int j=i; j < K; j++) {
      arma::mat Xj = X_.col(j);
      arma::mat rijT = matTimesVec(rr2T_s2_Xi, Xj).t();
      arma::mat fj = fj_cube.slice(j);
      arma::mat fjXBXi = matTimesVec((fj % nodesminusXB).t() / s2, Xi).t();
      arma::rowvec numPrime = arma::sum((rijT - fjXBXi), 0);
      arma::rowvec denomPrime = arma::sum(fj, 0);
      double ij = 2 * sum((w.t() % ((numPrime % denom) + (num % denomPrime))) / denom2);
      H(i, j) = ij;
      H(j, i) = ij; 
    }
    double gr0 = -2 * sum(((w.t() % arma::sum(fi, 0)) / denom));
    double gr_ = -2 * sum((w.t() % arma::sum((matTimesVec(rr2_nodsminusXB_t, Xi) / pow(s_, 2)).t(), 0)) / colSums_rr2_);
    double kv = (gr_ - gr0) / 1e-6;
    H(K, i) = kv;
    H(i, K) = kv;
  }
  return (H);
}



// [[Rcpp::export]] 
arma::mat calcHessOld(int K, arma::mat rr2, arma::mat rr2_, arma::mat trr2mxb, arma::mat X_, arma::mat nodesminusXB, arma::vec w, double s2, double s_) {
  arma::mat H;
  H.zeros(K+1, K+1);
  arma::mat rr2T = rr2.t();
  arma::rowvec denom = sum(rr2, 0);
  arma::rowvec denom2 = pow(denom, 2);
  arma::rowvec colSums_rr2_ = arma::sum(rr2_, 0);
  
  for (int i=0; i < K; i++) {
    arma::mat Xi = X_.col(i);
    arma::mat trr2mxb_s2 = trr2mxb / s2;
    arma::mat fi = matTimesVec(trr2mxb_s2, Xi).t();
    arma::rowvec num = arma::sum(fi, 0);
    for (int j=i; j < K; j++) {
      arma::mat Xj = X_.col(j);
      arma::mat fj = matTimesVec(trr2mxb_s2, Xj).t();
      arma::mat rijT = matTimesVec(matTimesVec(rr2T / s2, Xi), Xj).t();
      arma::mat fjXBXi = matTimesVec((fj % nodesminusXB).t() / s2, Xi).t();
      arma::rowvec numPrime = arma::sum((rijT - fjXBXi), 0);
      arma::rowvec denomPrime = arma::sum(fj, 0);
      double ij = 2 * sum((w.t() % ((numPrime % denom) + (num % denomPrime))) / denom2);
      H(i, j) = ij;
      H(j, i) = ij; 
    }
    double gr0 = -2 * sum(((w.t() % arma::sum(fi, 0)) / denom));
    double gr_ = -2 * sum((w.t() % arma::sum((matTimesVec((rr2_ % nodesminusXB).t(), Xi) / pow(s_, 2)).t(), 0)) / colSums_rr2_);
    double kv = (gr_ - gr0) / 1e-6;
    H(K, i) = kv;
    H(i, K) = kv;
  }
  return (H);
}
