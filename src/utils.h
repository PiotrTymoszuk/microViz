#ifndef __utils__
#define __utils__

#include <Rcpp.h>

using namespace Rcpp;

// distribution statistics of numeric vector

double Mean(NumericVector x);
double Median(NumericVector x);
double geoMean(NumericVector x);
double harmMean(NumericVector x);
double Var(NumericVector x);
double SD(NumericVector x);
double SEM(NumericVector x);
double GiniCpp(NumericVector x, bool unbiased);
double freqRatioCpp(NumericVector x);
double percUniqueCpp(NumericVector x);

// variable missingness

double missingCpp(NumericVector x);

// arranging values of a numeric matrix and vector

NumericMatrix sortByFirstColumn(NumericMatrix mat, bool decreasing);
IntegerVector sortIndexes(NumericVector x, bool decreasing) ;

// integration

double trapezoidal_integration(NumericVector x, NumericVector y);

// quantiles and confidence intervals

NumericVector Quantile(NumericVector x, NumericVector probs);
NumericVector perci(NumericVector theta, double conf_level);
NumericVector bca(NumericVector theta, double conf_level);

// contingency table

IntegerVector Table(NumericVector x);

// matrix row- and colnames

LogicalVector checkNames(const NumericMatrix &x);

#endif // __utils__
