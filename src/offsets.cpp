/*** Functions for offsetting vectors and matrices*/

#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::export]]

NumericVector offsetVector(NumericVector x) {

  // offsets values of a numeric vector so that the minimum is set to 0

  NumericVector result(x.size());
  double min_x = min(x);

  x = na_omit(x);

  if(min_x == 0) {

    result = x;

  } else if(min_x < 0) {

    result = x + abs(min_x);

  } else {

    result = x - min_x;

  }

  return result;

}

// [[Rcpp::export]]

NumericMatrix offsetMatrix(NumericMatrix x) {

  int n_rows = x.nrow();
  int n_cols = x.ncol();
  NumericVector col_inp(n_rows);

  NumericMatrix result(n_rows, n_cols);

  for(int i = 0; i < n_cols; i++) {

    col_inp = x(_, i);

    result(_, i) = offsetVector(col_inp);

  }

  return result;

}

// END
