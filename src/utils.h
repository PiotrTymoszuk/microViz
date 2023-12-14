#ifndef __utils__
#define __utils__

#include <Rcpp.h>

using namespace Rcpp;

double Median(NumericVector x);
double geoMean(NumericVector x);
double harmMean(NumericVector x);
double Var(NumericVector x);
double SD(NumericVector x);
double GiniCpp(NumericVector x, bool unbiased);
double freqRatioCpp(NumericVector x);
double percUniqueCpp(NumericVector x);

NumericVector Quantile(NumericVector x, NumericVector probs);
NumericVector perci(NumericVector theta, double conf_level);
NumericVector bca(NumericVector theta, double conf_level);

IntegerVector Table(NumericVector x);

#endif // __utils__
