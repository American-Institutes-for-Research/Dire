// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cumsum1
NumericVector cumsum1(NumericVector x);
RcppExport SEXP _Dire_cumsum1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cumsum1(x));
    return rcpp_result_gen;
END_RCPP
}
// GPCMC
double GPCMC(NumericVector d, double a, double theta, double score, double D);
RcppExport SEXP _Dire_GPCMC(SEXP dSEXP, SEXP aSEXP, SEXP thetaSEXP, SEXP scoreSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(GPCMC(d, a, theta, score, D));
    return rcpp_result_gen;
END_RCPP
}
// polyLvls
NumericVector polyLvls(List d, NumericVector a, double theta, NumericVector score, NumericVector D);
RcppExport SEXP _Dire_polyLvls(SEXP dSEXP, SEXP aSEXP, SEXP thetaSEXP, SEXP scoreSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(polyLvls(d, a, theta, score, D));
    return rcpp_result_gen;
END_RCPP
}
// ansItems
NumericVector ansItems(List d, NumericVector a, NumericVector theta, NumericVector score, NumericVector D);
RcppExport SEXP _Dire_ansItems(SEXP dSEXP, SEXP aSEXP, SEXP thetaSEXP, SEXP scoreSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(ansItems(d, a, theta, score, D));
    return rcpp_result_gen;
END_RCPP
}
// ldbinomC
NumericVector ldbinomC(NumericVector x, NumericVector pr);
RcppExport SEXP _Dire_ldbinomC(SEXP xSEXP, SEXP prSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pr(prSEXP);
    rcpp_result_gen = Rcpp::wrap(ldbinomC(x, pr));
    return rcpp_result_gen;
END_RCPP
}
// multItems
NumericVector multItems(NumericVector x1, NumericVector guess, NumericVector D, NumericVector slope, NumericVector nodes, NumericVector difficulty);
RcppExport SEXP _Dire_multItems(SEXP x1SEXP, SEXP guessSEXP, SEXP DSEXP, SEXP slopeSEXP, SEXP nodesSEXP, SEXP difficultySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type guess(guessSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type D(DSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type slope(slopeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type difficulty(difficultySEXP);
    rcpp_result_gen = Rcpp::wrap(multItems(x1, guess, D, slope, nodes, difficulty));
    return rcpp_result_gen;
END_RCPP
}
// calcRrij
NumericVector calcRrij(int i, int j, NumericMatrix rr1fi, NumericMatrix rr2fj, double detSigma, NumericVector mvnResid);
RcppExport SEXP _Dire_calcRrij(SEXP iSEXP, SEXP jSEXP, SEXP rr1fiSEXP, SEXP rr2fjSEXP, SEXP detSigmaSEXP, SEXP mvnResidSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rr1fi(rr1fiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rr2fj(rr2fjSEXP);
    Rcpp::traits::input_parameter< double >::type detSigma(detSigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mvnResid(mvnResidSEXP);
    rcpp_result_gen = Rcpp::wrap(calcRrij(i, j, rr1fi, rr2fj, detSigma, mvnResid));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Dire_cumsum1", (DL_FUNC) &_Dire_cumsum1, 1},
    {"_Dire_GPCMC", (DL_FUNC) &_Dire_GPCMC, 5},
    {"_Dire_polyLvls", (DL_FUNC) &_Dire_polyLvls, 5},
    {"_Dire_ansItems", (DL_FUNC) &_Dire_ansItems, 5},
    {"_Dire_ldbinomC", (DL_FUNC) &_Dire_ldbinomC, 2},
    {"_Dire_multItems", (DL_FUNC) &_Dire_multItems, 6},
    {"_Dire_calcRrij", (DL_FUNC) &_Dire_calcRrij, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_Dire(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}