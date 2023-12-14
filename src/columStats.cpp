/*** Columns stats */

#include <Rcpp.h>
#include <Rmath.h>
#include "utils.h"

using namespace Rcpp;

// column medians

// [[Rcpp::export]]

NumericVector colMed(NumericMatrix x, bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericVector res(x_ncols);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res[i] = Median(col_inp);

  }

  return res;

}

// column geometric means

// [[Rcpp::export]]

NumericVector colGeoMean(NumericMatrix x, bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericVector res(x_ncols);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res[i] = geoMean(col_inp);

  }

  return res;

}

// column harmonic means

// [[Rcpp::export]]

NumericVector colHarmMean(NumericMatrix x, bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericVector res(x_ncols);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res[i] = harmMean(col_inp);

  }

  return res;

}

// column variances

// [[Rcpp::export]]

NumericVector colVariance(NumericMatrix x, bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericVector res(x_ncols);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res[i] = Var(col_inp);

  }

  return res;

}

// column standard deviations

// [[Rcpp::export]]

NumericVector colSD(NumericMatrix x, bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericVector res(x_ncols);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res[i] = SD(col_inp);

  }

  return res;

}

// column Gini coefficients

// [[Rcpp::export]]

NumericVector colGi(NumericMatrix x, bool unbiased = true, bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericVector res(x_ncols);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res[i] = GiniCpp(col_inp, unbiased);

  }

  return res;

}

// column quantiles

// [[Rcpp::export]]

NumericMatrix colQuant(NumericMatrix x,
                       NumericVector probs,
                       bool na_rm = true) {

  int x_ncols = x.ncol();
  int quant_size = probs.size();

  NumericVector col_inp;
  NumericMatrix res(x_ncols, quant_size);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res(i, _) = Quantile(col_inp, probs);

  }

  return res;

}

// column percentile confidence intervals

//[[Rcpp::export]]

NumericMatrix colPerCi(NumericMatrix x,
                       double conf_level,
                       bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericMatrix res(x_ncols, 2);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res(i, _) = perci(col_inp, conf_level);

  }

  return res;

}

// column BCA confidence intervals

//[[Rcpp::export]]

NumericMatrix colBcaCi(NumericMatrix x,
                       double conf_level,
                       bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericMatrix res(x_ncols, 2);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res(i, _) = bca(col_inp, conf_level);

  }

  return res;

}

// matrix counterparts of the percUnique() and freqRatio() functions

//[[Rcpp::export]]

NumericVector colFreqRatio(NumericMatrix x) {

  int n_vars = x.ncol();
  double col_res;
  NumericVector res(n_vars);

  for(int i = 0; i < n_vars; ++i) {

    col_res = freqRatioCpp(x(_, i));

    res[i] = col_res;

  }

  return res;

}

//[[Rcpp::export]]

NumericVector colPercUnique(NumericMatrix x) {

  int n_vars = x.ncol();
  double col_res;
  NumericVector res(n_vars);

  for(int i = 0; i < n_vars; ++i) {

    col_res = percUniqueCpp(x(_, i));

    res[i] = col_res;

  }

  return res;

}

//deviations from data center

//[[Rcpp::export]]

NumericMatrix Delta(NumericMatrix x, NumericVector mu) {

  int n_vars = x.ncol();
  int n_obs = x.nrow();

  NumericVector col_res;

  NumericMatrix res(n_obs, n_vars);

  for(int i = 0; i < n_vars; ++i) {

    col_res = x(_, i) - mu[i];

    res(_, i) = col_res;

  }

  return res;

}

// column minima and maxima

//[[Rcpp::export]]

NumericVector colMi(NumericMatrix x, bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericVector res(x_ncols);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res[i] = Rcpp::min(col_inp);

  }

  return res;


}

//[[Rcpp::export]]

NumericVector colMa(NumericMatrix x, bool na_rm = true) {

  int x_ncols = x.ncol();

  NumericVector col_inp;
  NumericVector res(x_ncols);

  for(int i = 0; i < x_ncols; ++i) {

    col_inp = x(_, i);

    if(na_rm) {

      col_inp = Rcpp::na_omit(col_inp);

    }

    res[i] = Rcpp::max(col_inp);

  }

  return res;

}

// END


