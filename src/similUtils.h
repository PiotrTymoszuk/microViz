#ifndef __similUtils__
#define __similUtils__

#include <Rcpp.h>

using namespace Rcpp;

double jaccardCpp(CharacterVector x, CharacterVector y);
double overlapCpp(CharacterVector x, CharacterVector y);
double diceCpp(CharacterVector x, CharacterVector y);
double tverskyCpp(CharacterVector x, CharacterVector y, double a, double b);

#endif // __similUtils__
