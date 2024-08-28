/*** Calculation helpers to be use by other Cpp functions */

#include <Rcpp.h>
#include <Rmath.h>
#include "utils.h"

using namespace Rcpp;

// calculation of median

// [[Rcpp::export]]

double Median(NumericVector x) {

  int x_size = x.size();
  double med;

  x.sort();

  if(x_size % 2 != 0) {

    med = x[x_size/2];

  } else {

    med = (x[x_size/2 - 1] + x[x_size/2])/2;

  }

  return med;

}

// calculation of geometric mean

// [[Rcpp::export]]

double geoMean(NumericVector x) {

  int n = x.size();
  double product = 1.0;

  for (int i = 0; i < n; ++i) {

    product *= x[i];

  }

  return std::pow(product, 1.0 / n);

}

// calculation of harmonic mean

// [[Rcpp::export]]

double harmMean(NumericVector x) {

  int n = x.size();
  double sum = 0;

  for (int i = 0; i < n; ++i) {

    sum += 1 / x[i];

  }

  return n / sum;

}

// calculation of variance and SD

// [[Rcpp::export]]

double Var(NumericVector x) {

    int n = x.size();
    double x_mean = 0;
    double sum = 0;

    x_mean = mean(x);

    for (int i = 0; i < n; ++i) {

      sum += std::pow(x[i] - x_mean, 2);

    }

    return sum / (n - 1);

}

// [[Rcpp::export]]

double SD(NumericVector x) {

  return std::sqrt(Var(x));

}

// Gini coefficients

// [[Rcpp::export]]

double GiniCpp(NumericVector x, bool unbiased = true) {

  int n = x.size();
  double sum_x = sum(x);

  NumericVector x_sorted = clone(x).sort();

  NumericVector cum_sums = Rcpp::cumsum(x_sorted);

  double sum_i_x = sum(cum_sums);

  double gini = (n + 1 - 2 * sum_i_x / sum_x) / n;

  if(unbiased) {

    gini = gini * n/(n - 1);

  }

  return gini;

}

// computation of quantiles and confidence intervals

// [[Rcpp::export]]

NumericVector Quantile(NumericVector x, NumericVector probs) {

  // calculation of quantiles
  // many thanks to https://github.com/RcppCore/Rcpp/issues/967

  const size_t n = x.size(), np = probs.size();

  if (n == 0) return x;
  if (np == 0) return probs;

  NumericVector index = (n - 1.) * probs, y = x.sort(), x_hi(np), qs(np);
  NumericVector lo = floor(index), hi = ceiling(index);

  for (size_t i = 0; i < np; ++i) {

    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];

    if ((index[i] > lo[i]) && (x_hi[i] != qs[i])) {

      double h;
      h = index[i] - lo[i];
      qs[i] = (1.- h) * qs[i] + h * x_hi[i];

    }

  }

  return qs;

}

// [[Rcpp::export]]

NumericVector perci(NumericVector theta, double conf_level = 0.95) {

  // computes percentile confidence intervals

  NumericVector ci_probs{(1 - conf_level)/2, (1 + conf_level)/2};

  return Quantile(theta, ci_probs);

}

// [[Rcpp::export]]

NumericVector bca(NumericVector theta, double conf_level = 0.95) {

  // computes BCA confidence intervals based on the R code
  // provided by the coxed package
  // https://rdrr.io/cran/coxed/src/R/bca.R

  double low;
  double high;

  low = (1 - conf_level)/2;
  high = 1 - low;

  int sims = theta.size();

  NumericVector low_theta;

  low_theta = ifelse(theta < mean(theta), 1.0, 0.0);

  double z_inv;

  z_inv = sum(low_theta)/sims;

  double z;

  z = R::qnorm(z_inv, 0, 1, 1, 0);

  NumericVector U;

  U = (sims - 1) * (mean(theta) - theta);

  double top;
  double under;
  double a;

  top = sum(pow(U, 3));

  under = sum(pow(U, 2));

  under = 6 * std::pow(under, 1.5);

  a = top/under;

  double lower_inv;
  double upper_inv;

  double q_low;
  double q_high;

  q_low = R::qnorm(low, 0, 1, 1, 0);

  lower_inv = z + (z + q_low)/(1 - a * (z + q_low));

  lower_inv = R::pnorm(lower_inv, 0, 1, 1, 0);

  q_high = R::qnorm(high, 0, 1, 1, 0);

  upper_inv = z + (z + q_high)/(1 - a * (z + q_high));

  upper_inv = R::pnorm(upper_inv, 0, 1, 1, 0);

  NumericVector probs = {lower_inv, upper_inv};

  return Quantile(theta, probs);

}

// calculation of element frequency in a numeric vector
// NAs are skipped!

// [[Rcpp::export]]

IntegerVector Table(NumericVector x) {

  IntegerVector counts;

  x = Rcpp::na_omit(x);

  counts = Rcpp::table(x);

  counts.sort(true);

  return counts;

}

// frequency ratio in a numeric vector:
// ratio of the frequency of the most common element to the
// second most common one

// [[Rcpp::export]]

double freqRatioCpp(NumericVector x) {

  NumericVector counts;

  counts = Table(x);

  if(counts.size() == 1) return 0.0;

  return counts[0]/counts[1];

}

// percent of unique observations; NAs are skippd!

// [[Rcpp::export]]

double percUniqueCpp(NumericVector x) {

  double x_size;
  double unique_size;

  x = Rcpp::na_omit(x);

  if(x.size() == 0 || x.size() == 1) return 0.0;

  NumericVector unique_vals = Rcpp::unique(x);

  x_size = x.size();
  unique_size = unique_vals.size();

  return unique_size/x_size * 100.0;

}

// fill a matrix given a vector of numeric values and a two column
// matrix of x and y indexes

// [[Rcpp::export]]

NumericMatrix fillMat(NumericVector x,
                      IntegerMatrix ind,
                      int dim) {

  if(x.size() != ind.nrow()) stop("Improper x length.");

  NumericMatrix res(dim, dim);

  for(int i = 0; i < x.size(); ++i) {

    res(ind(i, 0), ind(i, 1)) = x[i];
    res(ind(i, 1), ind(i, 0)) = x[i];

  }

  return res;

}

// number of missing observations

//[[Rcpp::export]]

double missingCpp(NumericVector x) {

  double x_size = x.size();
  double x_complete_size = na_omit(x).size();

  return x_size - x_complete_size;

}

// END
