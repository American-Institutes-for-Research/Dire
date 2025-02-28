#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector ldbinomC(NumericVector x, NumericVector pr){
  // declare container 
  NumericVector binom; 
  // apply ldbinom 
  binom = x*log(pr) + (1-x)*log(1-pr); 
  //replace NAs with zero 
  binom[is_nan(binom)] = 0; 
  // return 
  return(binom); 
}

// [[Rcpp::export]]
NumericVector multItems(NumericVector x1, NumericVector guess, NumericVector D, 
                        NumericVector slope, NumericVector nodes, NumericVector difficulty){
  // define helper variables
  int Q = nodes.size();
  NumericVector nodesNew(Q); //return container, we don't want to overwrite original nodes 
  // loop through nodes
  for(int i = 0; i < Q; i++){
    // sum ldbinom and append to return container 
    nodesNew[i] = sum(ldbinomC(x1, guess + ((1 - guess) / (1 + exp(-D * slope * (nodes[i] - difficulty)))))); 
  }
  // return 
  return nodesNew; 
}
