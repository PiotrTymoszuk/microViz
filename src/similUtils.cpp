/*** Calculation of distances between vectors */

#include <Rcpp.h>
#include <Rmath.h>
#include "similUtils.h"

using namespace Rcpp;

// calculates Jaccard similarity between two character vectors

double jaccardCpp(CharacterVector x,
                  CharacterVector y) {

  CharacterVector un;
  CharacterVector inter;

  un = union_(x, y);
  inter = intersect(x, y);

  double denom = un.length();

  if(denom == 0) return NA_REAL;

  return inter.length()/denom;

}

// calculates overlap or Szymkiewicz–Simpson similarity coefficient of
// two character vectors

double overlapCpp(CharacterVector x, CharacterVector y) {

  CharacterVector inter;

  double denom = x.length();

  if(y.length() < x.length()) denom = y.length();

  if(denom == 0) return NA_REAL;

  inter = intersect(x, y);

  return inter.length()/denom;

}

// calculates Dice or Dice-Sørensen similarity coefficient of
// two character vectors

double diceCpp(CharacterVector x, CharacterVector y) {

  CharacterVector inter;

  double denom = x.length() + y.length();

  if(denom == 0) return NA_REAL;

  inter = intersect(x, y);

  return 2 * inter.length()/denom;

}

// calculates Tversky coefficient of similarity between two character
// vectors

double tverskyCpp(CharacterVector x,
                  CharacterVector y,
                  double a = 1,
                  double b = 1) {

  CharacterVector inter;
  CharacterVector diff1;
  CharacterVector diff2;

  inter = intersect(x, y);
  diff1 = setdiff(x, y);
  diff2 = setdiff(y, x);

  double denom;
  double inter_len = inter.length();

  denom = inter_len + a * diff1.length() + b * diff1.length();

  if(denom == 0) return NA_REAL;

  return inter_len/denom;

}

// END
