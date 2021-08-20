#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cumsum1(NumericVector x){
  // initialize an accumulator variable
  double acc = 0;
  
  // initialize the result vector
  NumericVector res(x.size());
  
  for(int i = 0; i < x.size(); i++){
    acc += x[i];
    res[i] = acc;
  }
  return res;
}

// [[Rcpp::export]]
double GPCMC(NumericVector d, double a, double theta, double score, double D) {
  // catch cases 
  int size = d.size(); 
  // score is NA 
  if(std::isnan(score)){
    return score; 
  }
  // score is higher than max 
  if (score > size) {
    stop("Score is higher than maximum");
  }
  // score is lower than min 
  if (score <= 0) {
    stop("Score is lower than minimum");
  }
  // define helper variables and vector 
  double Da;
  double denom; 
  double numerator; 
  score = score - 1; // correction of indexing 
  NumericVector dSub = d[Rcpp::Range(0, score)]; 
  // calculations
  Da = D*a;
  // calculate numerator 
  numerator = exp(sum(Da * (theta - dSub))); 
  // calculate denominator 
  denom = sum(exp(cumsum1(Da * (theta - d)))); 
  // final 
  return numerator / denom; 
}

// [[Rcpp::export]]
NumericVector polyLvls(List d, NumericVector a, double theta, NumericVector score, NumericVector D){
  int size = d.size(); 
  NumericVector lvls(size); 
  // for every polytimous item
  for(int i = 0; i < size; i++){
    // apply gpcm
    lvls[i] = GPCMC(d[i], a[i], theta, score[i], D[i]);
  }
  return lvls; 
}

// [[Rcpp::export]]
NumericVector ansItems(List d, NumericVector a, NumericVector theta, NumericVector score, NumericVector D){
  // declare objects
  int size = theta.size(); 
  NumericVector polyItems(size); 
  // for every item student answered 
  for(int i = 0; i < size; i++){
    // for every polytimous item
    polyItems[i] = sum(log(polyLvls(d, a, theta[i], score, D)));
  }
  return polyItems; 
}
