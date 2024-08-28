/*** Row stats */

#include <Rcpp.h>
#include <Rmath.h>
#include "utils.h"

using namespace Rcpp;

// row medians

// [[Rcpp::export]]

NumericVector rowMed(NumericMatrix x, bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res[i] = Median(row_inp);

  }

  return res;

}

// geometric means of the rows

// [[Rcpp::export]]

NumericVector rowGeoMean(NumericMatrix x, bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res[i] = geoMean(row_inp);

  }

  return res;

}

// harmonic means of the rows

// [[Rcpp::export]]

NumericVector rowHarmMean(NumericMatrix x, bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res[i] = harmMean(row_inp);

  }

  return res;

}

// row variances

// [[Rcpp::export]]

NumericVector rowVariance(NumericMatrix x, bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res[i] = Var(row_inp);

  }

  return res;

}

// row standard deviations

// [[Rcpp::export]]

NumericVector rowSD(NumericMatrix x, bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res[i] = SD(row_inp);

  }

  return res;

}

// Gini coefficients of the rows

// [[Rcpp::export]]

NumericVector rowGi(NumericMatrix x, bool unbiased = true, bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res[i] = GiniCpp(row_inp, unbiased);

  }

  return res;

}

// row quantiles

// [[Rcpp::export]]

NumericMatrix rowQuant(NumericMatrix x,
                       NumericVector probs,
                       bool na_rm = true) {

  int x_nrows = x.nrow();
  int quant_size = probs.size();

  NumericVector row_inp;
  NumericMatrix res(x_nrows, quant_size);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res(i, _) = Quantile(row_inp, probs);

  }

  return res;

}

// row percentile confidence intervals

//[[Rcpp::export]]

NumericMatrix rowPerCi(NumericMatrix x,
                       double conf_level,
                       bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericMatrix res(x_nrows, 2);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res(i, _) = perci(row_inp, conf_level);

  }

  return res;

}

// BCA confidence intervals of the rows

//[[Rcpp::export]]

NumericMatrix rowBcaCi(NumericMatrix x,
                       double conf_level,
                       bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericMatrix res(x_nrows, 2);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res(i, _) = bca(row_inp, conf_level);

  }

  return res;

}

// row matrix counterparts of the percUnique() and freqRatio() functions

//[[Rcpp::export]]

NumericVector rowFreqRatio(NumericMatrix x) {

  int n_vars = x.nrow();
  double row_res;
  NumericVector res(n_vars);

  for(int i = 0; i < n_vars; ++i) {

    row_res = freqRatioCpp(x(i, _));

    res[i] = row_res;

  }

  return res;

}

//[[Rcpp::export]]

NumericVector rowPercUnique(NumericMatrix x) {

  int n_vars = x.nrow();
  double row_res;
  NumericVector res(n_vars);

  for(int i = 0; i < n_vars; ++i) {

    row_res = percUniqueCpp(x(i, _));

    res[i] = row_res;

  }

  return res;

}

// row minima and maxima

//[[Rcpp::export]]

NumericVector rowMi(NumericMatrix x, bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res[i] = Rcpp::min(row_inp);

  }

  return res;


}

//[[Rcpp::export]]

NumericVector rowMa(NumericMatrix x, bool na_rm = true) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    if(na_rm) {

      row_inp = Rcpp::na_omit(row_inp);

    }

    res[i] = Rcpp::max(row_inp);

  }

  return res;


}

// missing observations

//[[Rcpp::export]]

NumericVector rowMis(NumericMatrix x) {

  int x_nrows = x.nrow();

  NumericVector row_inp;
  NumericVector res(x_nrows);

  for(int i = 0; i < x_nrows; ++i) {

    row_inp = x(i, _);

    res[i] = missingCpp(row_inp);

  }

  return res;

}

// END
