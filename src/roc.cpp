/*** Rcpp functions for ROC analysis */

#include <Rcpp.h>
#include <Rmath.h>
#include "utils.h"
#include "roc.h"

using namespace Rcpp;

// sensitivity, specificity, Youden's J, TP, TN, FP and FN count calculation
// for a given threshold

// [[Rcpp::export]]

NumericVector rocThreshold(IntegerVector outcome,
                           NumericVector marker,
                           double threshold,
                           String direction = "<",
                           bool skip_control = false) {

  // entry control

  int n = outcome.size();

  if(!skip_control) {

    if(marker.size() != n) {

      stop("Incompatible sizes of outcome and marker vectors.");

    }

    LogicalVector non_binary = ((outcome > 1) | (outcome < 0));

    if(sum(non_binary) > 0) {

      stop("outcome is not binary!");

    }

    LogicalVector outcome_idx = (outcome > 0);
    LogicalVector control_idx = (outcome == 0);

    if(sum(outcome_idx) == 0) {

      stop("There are no outcome cases!");

    }

    if(sum(control_idx) == 0) {

      stop("There are no control cases!");

    }

  }

  LogicalVector outcome_na = is_na(outcome);
  LogicalVector marker_na = is_na(marker);

  // storage vectors

  double true_positive = 0;
  double false_positive = 0;
  double true_negative = 0;
  double false_negative = 0;

  double sensitivity = NA_REAL;
  double specificity = NA_REAL;

  // counting true and false hits

  if(direction == "<") {

    for (int i = 0; i < n; ++i) {

      if(outcome_na[i] | marker_na[i]) continue;

      if (marker[i] >= threshold) {

        if (outcome[i] == 1) {

          true_positive++;

        } else {

          false_positive++;

        }

      } else {

        if (outcome[i] == 1) {

          false_negative++;

        } else {

          true_negative++;

        }

      }

    }

  } else {

    for (int i = 0; i < n; ++i) {

      if(outcome_na[i] | marker_na[i]) continue;

      if (marker[i] <= threshold) {

        if (outcome[i] == 1) {

          true_positive++;

        } else {

          false_positive++;

        }

      } else {

        if (outcome[i] == 1) {

          false_negative++;

        } else {

          true_negative++;

        }

      }

    }

  }

  // computation of sensitivity and specificity

  if(true_positive + false_negative > 0.0) {

    sensitivity = true_positive/(true_positive + false_negative);

  }

  if(true_negative + false_positive) {

    specificity = true_negative/(true_negative + false_positive);

  }

  NumericVector res {threshold,
                     true_positive,
                     false_positive,
                     true_negative,
                     false_negative,
                     sensitivity,
                     specificity,
                     sensitivity + specificity};

  res.names() =
    CharacterVector({"threshold", "TP", "FP", "TN", "FN", "Se", "Sp", "J"});

  return res;

}

// sensitivity, specificity, Youden's J, TP, TN, FP and FN count calculation for
// all values of a vector

//[[Rcpp::export]]

NumericMatrix rocMarker(IntegerVector outcome,
                        NumericVector marker,
                        String direction = "auto",
                        bool skip_control = false) {

  // entry control: done only once for the entire vector

  int n = outcome.size();

  LogicalVector outcome_idx = (outcome > 0);
  LogicalVector control_idx = (outcome == 0);

  if(!skip_control) {

    if(marker.size() != n) {

      stop("Incompatible sizes of outcome and marker vectors.");

    }

    LogicalVector non_binary = ((outcome > 1) | (outcome < 0));

    if(sum(non_binary) > 0) {

      stop("outcome is not binary!");

    }

    if(sum(outcome_idx) == 0) {

      stop("There are no outcome cases!");

    }

    if(sum(control_idx) == 0) {

      stop("There are no control cases!");

    }

  }

  // determination of the comparison direction
  // if direction = "auto", the direction is determined by difference in
  // medians of the marker in the control and outcome groups

  if(direction == "auto") {

    NumericVector markerOutcome = marker[outcome_idx];
    NumericVector markerControl = marker[control_idx];

    double medianOutcome = Median(markerOutcome);
    double medianControl = Median(markerControl);

    if(medianOutcome >= medianControl) {

      direction = "<";

    } else {

      direction = ">";

    }

  }

  // ROC analysis

  NumericVector uniqueThrx = unique(marker);
  int nThrx = uniqueThrx.size();

  NumericVector rocThrx(8);

  NumericMatrix rocResults(nThrx, 8);

  colnames(rocResults) =
    CharacterVector({"marker", "TP", "FP", "TN", "FN", "Se", "Sp", "J"});

  for(int i = 0; i < nThrx; ++i) {

    rocThrx = rocThreshold(outcome, marker, uniqueThrx[i], direction, true);

    rocResults(i, _) = rocThrx;

  }

  return rocResults;

}

// Calculation of AUC by trapezoidal integration of Sensitivity
// over 1 - Specificity: version for one marker vector

//[[Rcpp::export]]

NumericVector aucVec(IntegerVector outcome,
                     NumericVector marker,
                     String direction = "auto",
                     bool skip_control = false) {

  // the input control is accomplished by the downstream functions

  // ROC analysis for all unique values of the marker vector
  // extraction of 1 - Specificity and Sensitivity

  NumericMatrix rocThrx =
    rocMarker(outcome, marker, direction, skip_control);

  // calculation of ROC AUC by trapezoidal integration

  NumericMatrix aucInputRaw(rocThrx.nrow(), 2);

  aucInputRaw(_, 0) = 1 - rocThrx(_, 6);
  aucInputRaw(_, 1) = rocThrx(_, 5);

  NumericMatrix aucInputSorted = sortByFirstColumn(aucInputRaw, false);

  double aucValue = trapezoidal_integration(aucInputSorted(_, 0),
                                            aucInputSorted(_, 1));

  // finding of the optimal cutoff of the marker by maximizing
  // the Youden's J
  // starting with a numeric matrix of indexes and Youden's J and finding
  // the index that corresponds to J's maximum

  IntegerVector j_indexes = sortIndexes(rocThrx(_, 7), true);

  int cutoffIdx = j_indexes[0];

  // the output numeric vector: AUC, values of the best cutoff, sensitivity,
  // specificity and Youden's J at the best cutoff

  NumericVector output {aucValue,
                        rocThrx(cutoffIdx, 0),
                        rocThrx(cutoffIdx, 5),
                        rocThrx(cutoffIdx, 6),
                        rocThrx(cutoffIdx, 7)};

  output.names() =
    CharacterVector {"auc", "cutoff", "Se", "Sp", "J"};

  return output;

}

// calculation of AUC for a matrix of markers

//[[Rcpp::export]]

NumericMatrix aucMtx(IntegerVector outcome,
                     NumericMatrix markers,
                     String direction = "auto") {

  // input control

  int n = outcome.size();

  LogicalVector outcome_idx = (outcome > 0);
  LogicalVector control_idx = (outcome == 0);

  if(markers.nrow() != n) {

    stop("Incompatible sizes of the outcome vector and marker matrix.");

  }

  LogicalVector non_binary = ((outcome > 1) | (outcome < 0));

  if(sum(non_binary) > 0) {

    stop("outcome is not binary!");

  }

  if(sum(outcome_idx) == 0) {

    stop("There are no outcome cases!");

  }

  if(sum(control_idx) == 0) {

    stop("There are no control cases!");

  }

  // calculation of AUC and ROC stats for the best cutoffs

  int n_markers = markers.ncol();

  NumericMatrix res(n_markers, 5);

  LogicalVector name_check = checkNames(markers);

  if(name_check[1]) rownames(res) = colnames(markers);

  colnames(res) = CharacterVector {"auc", "cutoff", "Se", "Sp", "J"};

  for(int i = 0; i < n_markers; ++i) {

    res(i, _) = aucVec(outcome, markers(_, i), direction, true);

  }

  return res;

}

// END
