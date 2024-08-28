/*** Calculation of distances between elements of a list of vectors */

#include <Rcpp.h>
#include <Rmath.h>
#include "similUtils.h"

using namespace Rcpp;

// similarity of two character vectors

// [[Rcpp::export]]

double vecSimilarity(CharacterVector x,
                     CharacterVector y,
                     String method = "jaccard") {

  // finding the similarity-calculating function

  double (*fun)(CharacterVector, CharacterVector);

  x = unique(x);
  y = unique(y);

  if(method == "jaccard") {

    fun = jaccardCpp;

  } else if(method == "overlap") {

    fun = overlapCpp;

  } else if(method == "dice") {

    fun = diceCpp;

  }

  // computation of similarity

  double res = fun(x, y);

  return res;

}

// [[Rcpp::export]]

double vecTversky(CharacterVector x,
                  CharacterVector y,
                  double a = 1,
                  double b = 1) {

  x = unique(x);
  y = unique(y);

  double res = tverskyCpp(x, y, a, b);

  return res;

}

// similarity of list elements
// a separate function for calculation of the Tversky coefficient

// [[Rcpp::export]]

NumericMatrix setSimilarity(List x,
                            String method = "jaccard") {

  // finding the similarity-calculating function

  double (*fun)(CharacterVector, CharacterVector);

  if(method == "jaccard") {

    fun = jaccardCpp;

  } else if(method == "overlap") {

    fun = overlapCpp;

  } else if(method == "dice") {

    fun = diceCpp;

  }

  // computation of similarity

  int dim = x.length();
  NumericMatrix sim_mtx(dim);

  for(int i = 0; i < dim; ++i) {

    CharacterVector el = x[i];

    x[i] = unique(el);

  }

  for(int x_index = 0; x_index < dim; ++x_index) {

    for(int y_index = 0; y_index < dim; ++y_index) {

      sim_mtx(x_index, y_index) = fun(x[x_index], x[y_index]);

    }

  }

  return sim_mtx;

}

// [[Rcpp::export]]

NumericMatrix setTversky(List x,
                         double a,
                         double b) {

  int dim = x.length();
  NumericMatrix sim_mtx(dim);

  for(int i = 0; i < dim; ++i) {

    CharacterVector el = x[i];

    x[i] = unique(el);

  }

  for(int x_index = 0; x_index < dim; ++x_index) {

    for(int y_index = 0; y_index < dim; ++y_index) {

      sim_mtx(x_index, y_index) = tverskyCpp(x[x_index], x[y_index], a, b);

    }

  }

  return sim_mtx;

}

// END
