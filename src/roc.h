#ifndef __roc__
#define __roc__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector rocThreshold(IntegerVector outcome,
                           NumericVector marker,
                           double threshold,
                           String direction,
                           bool skip_control);

NumericMatrix rocMarker(IntegerVector outcome,
                        NumericVector marker,
                        String direction,
                        bool skip_control);

NumericVector aucVec(IntegerVector outcome,
                     NumericVector marker,
                     String direction,
                     bool skip_control);

#endif // __roc__
