/*** Functions for calculation of pooled variances and SD */

#include <Rcpp.h>
#include <Rmath.h>
#include "utils.h"

using namespace Rcpp;

// pooled variance for a paired T test

// [[Rcpp::export]]

double poolVarPaired(NumericVector x, NumericVector y) {

  // as per definition of the test, the vectors x and y are of equal lengths

  x = na_omit(x);
  y = na_omit(y);

  NumericVector res(x.size());

  for(int i = 0; i < x.size(); ++i) {

    res[i] = x[i] - y[i];

  }

  return Var(res);

}

// pooled variance for a Welch's T test

// [[Rcpp::export]]

double poolVarWelch(NumericVector x, NumericVector y) {

  x = na_omit(x);
  y = na_omit(y);

  return (Var(x) + Var(y))/2;

}

// pooled variance for a standard T test

// [[Rcpp::export]]

double poolVarStandard(NumericVector x, NumericVector y) {

  x = na_omit(x);
  y = na_omit(y);

  double mean_x = mean(x);
  double mean_y = mean(y);

  double n_diff = x.size() + y.size() - 2;

  if(n_diff <= 0) n_diff = 1;

  NumericVector squares_x(x.size());
  NumericVector squares_y(y.size());

  for(int i = 0; i < x.size(); ++i) {

    squares_x[i] = x[i] - mean_x;

    squares_x[i] = std::pow(squares_x[i], 2.0);

  }

  for(int j = 0; j < y.size(); ++j) {

    squares_y[j] = y[j] - mean_y;

    squares_y[j] = std::pow(squares_y[j], 2.0);

  }

  return (sum(squares_x) + sum(squares_y))/n_diff;


}

// pooled variance for linear modeling

// [[Rcpp::export]]

NumericVector poolVarLM(List x) {

  // takes a list of numeric vectors, each corresponding to one
  // level of the explanatory factor

  NumericVector x_vars(x.size());
  NumericVector avg_vars(x.size());
  NumericVector x_part;

  for(int i = 0; i < x.size(); ++i) {

    x_part = x[i];

    x_vars[i] = Var(na_omit(x_part));

  }

  for(int j = 0; j < x.size(); ++j) {

    avg_vars[j] = (x_vars[j] + x_vars[0])/2;

  }

  return(avg_vars);

}
