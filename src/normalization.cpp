/*** Rcpp functions for min-max scaling and Z-score calculation */

#include <Rcpp.h>
#include <Rmath.h>
#include "utils.h"
#include "normalization.h"

using namespace Rcpp;

// min-max scaling of a vector

// [[Rcpp::export]]

NumericVector minMaxVec(NumericVector x) {

  int n = x.size();

  // checking for NA in the vector: if there are
  // NA only, the vector is returned as it is

  LogicalVector na_check = is_na(x);

  if(na_omit(x).size() == 0) return x;

  // scaling

  double x_max = max(na_omit(x));
  double x_min = min(na_omit(x));

  NumericVector res(n, NA_REAL);

  for(int i = 0; i < n; ++i) {

    if(na_check[i]) continue;

    res[i] = (x[i] - x_min)/(x_max - x_min);

  }

  return res;

}

// column-wise min-max scaling of a matrix

// [[Rcpp::export]]

NumericMatrix minMaxCols(NumericMatrix x) {

  int n_cols = x.ncol();
  int n_rows = x.nrow();

  NumericVector colVec(n_rows);
  NumericMatrix res(n_rows, n_cols);

  for(int i = 0; i < n_cols; ++i) {

    colVec = x(_, i);

    res(_, i) = minMaxVec(colVec);

  }

  LogicalVector hasNames = checkNames(x);

  if(hasNames[1]) rownames(res) = rownames(x);
  if(hasNames[2]) colnames(res) = colnames(x);

  return res;

}

// row-wise min-max scaling of a matrix

// [[Rcpp::export]]

NumericMatrix minMaxRows(NumericMatrix x) {

  int n_cols = x.ncol();
  int n_rows = x.nrow();

  NumericVector rowVec(n_cols);
  NumericMatrix res(n_rows, n_cols);

  for(int i = 0; i < n_rows; ++i) {

    rowVec = x(i, _);

    res(i, _) = minMaxVec(rowVec);

  }

  LogicalVector hasNames = checkNames(x);

  if(hasNames[1]) rownames(res) = rownames(x);
  if(hasNames[2]) colnames(res) = colnames(x);

  return res;

}

// Z-scores of a vector

// [[Rcpp::export]]

NumericVector zScoreVec(NumericVector x,
                        String center = "mean",
                        String dispersion = "sd") {

  int n = x.size();

  // checking for NA in the vector: if there are
  // NA only, the vector is returned as it is

  LogicalVector na_check = is_na(x);

  if(na_omit(x).size() == 0) return x;

  // finding the functions for calculating the centrality
  // and dispersion stats

  double (*centerFun)(NumericVector);
  double (*dispFun)(NumericVector);

  if(center == "mean") {

    centerFun = Mean;

  } else if(center == "median") {

    centerFun = Median;

  } else if(center == "geo_mean") {

    centerFun = geoMean;

  } else {

    centerFun = harmMean;

  }

  if(dispersion == "sd") {

    dispFun = SD;

  } else {

    dispFun = SEM;

  }

  // calculation of Z-scores

  double centerVal = centerFun(na_omit(x));
  double dispVal = dispFun(na_omit(x));

  NumericVector res(n, NA_REAL);

  if(dispVal == 0) {

    warning("dispersion statistic is 0; returning NA");

    return res;

  }

  for(int i = 0; i < n; ++i) {

    if(na_check[i]) continue;

    res[i] = (x[i] - centerVal)/dispVal;

  }

  return res;

}

// Z-scores of matrix columns

// [[Rcpp::export]]

NumericMatrix zScoreCols(NumericMatrix x,
                         String center = "mean",
                         String dispersion = "sd") {

  int n_cols = x.ncol();
  int n_rows = x.nrow();

  NumericVector colVec(n_rows);
  NumericMatrix res(n_rows, n_cols);

  for(int i = 0; i < n_cols; ++i) {

    colVec = x(_, i);

    res(_, i) = zScoreVec(colVec, center, dispersion);

  }

  LogicalVector hasNames = checkNames(x);

  if(hasNames[1]) rownames(res) = rownames(x);
  if(hasNames[2]) colnames(res) = colnames(x);

  return res;

}

// Z-scores of matrix rows

// [[Rcpp::export]]

NumericMatrix zScoreRows(NumericMatrix x,
                         String center = "mean",
                         String dispersion = "sd") {

  int n_cols = x.ncol();
  int n_rows = x.nrow();

  NumericVector rowVec(n_cols);
  NumericMatrix res(n_rows, n_cols);

  for(int i = 0; i < n_rows; ++i) {

    rowVec = x(i, _);

    res(i, _) = zScoreVec(rowVec, center, dispersion);

  }

  LogicalVector hasNames = checkNames(x);

  if(hasNames[1]) rownames(res) = rownames(x);
  if(hasNames[2]) colnames(res) = colnames(x);

  return res;

}

// END
