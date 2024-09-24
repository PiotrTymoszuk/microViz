#ifndef __normalization__
#define __normalization__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector minMaxVec(NumericVector x);
NumericVector zScoreVec(NumericVector x, String center, String dispersion);

#endif // __normalization__
