# functions for calculation of numeric stats

#' @include imports.R

  NULL

# Numeric vectors ----------

#' Distribution statistics for numeric vectors and matrices.
#'
#' @description
#' A bunch of functions that compute extra numeric statistics not provided by
#' base R such as geometric and harmonic means or BCA confidence intervals.
#'
#' @details
#' `Gmean()`, `Hmean()` and `Gini()` compute geometric mean, harmonic mean,
#' and Gini coefficient, respectively.
#' `perCI()` and `bcaCI()` return confidence intervals calculated with the
#' percentile and BCA (bias corrected and accelerated) method.
#' `freqRatio()` computes the ratio of frequencies of the first most frequent
#' element of a vector to the second most frequent one.
#' `percUnique()` returns the percentage of unique observations.
#' The `freqRatio()` and `percUnique()` functions skip silently any `NA` in the
#' vector.
#'
#' @references
#' DiCiccio TJ, Efron B. Bootstrap confidence intervals.
#' https://doi.org/101214/ss/1032280214 (1996) 11:189–228.
#' doi:10.1214/SS/1032280214
#'
#' @param x a numeric vector or a numeric matrix.
#' @param na.rm logical, should `NAs` ne removed prior to the statistic
#' computation?
#' @param conf_level a single numeric with the width of the confidence interval.
#' @param unbiased logical, should estimated of the population Gini coefficient
#' be returned (multiplied by `n/(n - 1)`)?
#'
#' @return a single numeric or a numeric vector, see Details.
#'
#' @export

  Gmean <- function(x, na.rm = TRUE) {

    stopifnot(is.numeric(x))
    stopifnot(is.logical(na.rm))

    if(na.rm) x <- x[!is.na(x)]

    geoMean(x)

  }

#' @rdname Gmean
#' @export

  Hmean <- function(x, na.rm = TRUE) {

    stopifnot(is.numeric(x))
    stopifnot(is.logical(na.rm))

    if(na.rm) x <- x[!is.na(x)]

    harmMean(x)

  }

#' @rdname Gmean
#' @export

  Gini <- function(x, unbiased = TRUE, na.rm = TRUE) {

    stopifnot(is.numeric(x))
    stopifnot(is.logical(na.rm))
    stopifnot(is.logical(unbiased))

    if(na.rm) x <- x[!is.na(x)]

    GiniCpp(x, unbiased)

  }


#' @rdname Gmean
#' @export

  perCI <- function(x, conf_level = 0.95, na.rm = TRUE) {

    stopifnot(is.numeric(x))
    stopifnot(is.numeric(conf_level))
    stopifnot(is.logical(na.rm))

    conf_level <- conf_level[1]

    stopifnot(conf_level <= 1)
    stopifnot(conf_level >= 0)

    if(na.rm) x <- x[!is.na(x)]

    perci(x, conf_level)

  }

#' @rdname Gmean
#' @export

  bcaCI <- function(x, conf_level = 0.95, na.rm = TRUE) {

    stopifnot(is.numeric(x))
    stopifnot(is.numeric(conf_level))
    stopifnot(is.logical(na.rm))

    conf_level <- conf_level[1]

    stopifnot(conf_level <= 1)
    stopifnot(conf_level >= 0)

    if(na.rm) x <- x[!is.na(x)]

    bca(x, conf_level)

  }

#' @rdname Gmean
#' @export

  freqRatio <- function(x) {

    stopifnot(is.numeric(x))

    freqRatioCpp(x);

  }

#' @rdname Gmean
#' @export

  percUnique <- function(x) {

    stopifnot(is.numeric(x))

    percUniqueCpp(x);

  }

# Data frames and matrices: column distribution and variable selection ------

#' Compute statistics for numeric columns of a data frame or a matrix.
#'
#' @description
#' The functions compute column medians, variances, standard deviations,
#' Gini coefficients, quantiles and confidence intervals as well as metrics of
#' unique values.
#' `distr_stats` computes a bunch of distribution stats that
#' are helpful at selection of variant features for the further analysis, e.g.
#' differential gene expression.
#'
#' @references
#' DiCiccio TJ, Efron B. Bootstrap confidence intervals.
#' https://doi.org/101214/ss/1032280214 (1996) 11:189–228.
#' doi:10.1214/SS/1032280214
#'
#' @param x a numeric data frame or matrix.
#' @param na.rm logical: should NAs be removed?
#' @param probs a numeric vector of quantiles.
#' @param conf_level a single numeric with the width of the confidence interval.
#' @param method method of confidence interval calculation: percentile (default)
#' or BCA.
#' @param unbiased logical, should estimated of the population Gini coefficient
#' be returned (multiplied by `n/(n - 1)`)?
#' @param freqCut the cutoff of ratio of the most frequent element to the
#' second most frequent element. Used for identification of near-zero
#' variance variables.
#' @param uniqueCut the cutoff of percentage of unique observations. Used for
#' identification of near-zero variance variables.
#' @param ... extra arguments passed to methods.
#'
#' @return `colMedians()`, `colGmeans()`, `colHmeans()`, `colVars()`,
#' `colSDs()`, `colGini()` return numeric vectors with, respectively,
#' column medians, geometric means, harmonic means, variances, standard
#' deviations and Gini coefficients.
#' `colQuantiles()`, and `colCI()` return numeric matrices with variables
#' listed in rows.
#' `colFreqRatios()` returns ratios of frequencies of the first most common
#' element to the second one for each column; the function silently skips any
#' `NAs`.
#' `colPercUniques()` returns percentages of unique observations in each column;
#' `NAs` are silently removed.
#' `colMissing` returns numbers of missing observations.
#'
#' `distr_stats()` returns a data frame with the following columns:
#' * `variable`: variable names
#' * `median`: variable medians
#' * `mean`: variable means
#' * `var`: variances
#' * `gini_coef`: Gini coefficients
#' * `var_mean_ratio`: variance to mean ratios
#' * `freqRatio`: ratios of frequencies of the first most common element to the
#' second most common element.
#' * `percentUnique`: percentage of unique observations
#' * `zeroVar`: logical, indicates if the variable has only one unique value.
#' * `nzv`: logical, indicates if the variable passes the variability criteria
#' defined by the `freqCut` and `uniqueCut` arguments.
#' * `n`: total observation number.
#' * `n_missing`: number of missing observations
#' * `perc_missing`: percentage of missing observations
#' For the `distr_stats()` output, any missing observations are removed
#' silently prior to calculation of the numeric statistics.
#'
#' @export

  colMedians <- function(x, ...) UseMethod('colMedians')

#' @rdname colMedians
#' @export

  colMedians.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- colMed(x, na.rm)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colMedians.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    colMedians.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colMins <- function(x, ...) UseMethod('colMins')

#' @rdname colMedians
#' @export

  colMins.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- colMi(x, na.rm)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colMins.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    colMins.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colMax <- function(x, ...) UseMethod('colMax')

#' @rdname colMedians
#' @export

  colMax.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- colMa(x, na.rm)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colMax.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    colMax.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colGmeans <- function(x, ...) UseMethod('colGmeans')

#' @rdname colMedians
#' @export

  colGmeans.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- colGeoMean(x, na.rm)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colGmeans.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    colGmeans.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colHmeans <- function(x, ...) UseMethod('colHmeans')

#' @rdname colMedians
#' @export

  colHmeans.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- colHarmMean(x, na.rm)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colHmeans.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    colHmeans.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colVars <- function(x, ...) UseMethod('colVars')

#' @rdname colMedians
#' @export

  colVars.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- colVariance(x, na.rm)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colVars.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    colVars.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colSDs <- function(x, ...) UseMethod('colSDs')

#' @rdname colMedians
#' @export

  colSDs.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- colSD(x, na.rm)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colSDs.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    colSDs.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colGini <- function(x, unbiased, ...) UseMethod('colGini')

#' @rdname colMedians
#' @export

  colGini.matrix <- function(x, unbiased = TRUE, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- colGi(x, unbiased, na.rm)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colGini.data.frame <- function(x, unbiased = TRUE, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    colGini.matrix(as.matrix(x), na.rm = na.rm)


  }

#' @rdname colMedians
#' @export

  colFreqRatios <- function(x, ...) UseMethod('colFreqRatios')

#' @rdname colMedians
#' @export

  colFreqRatios.matrix <- function(x, ...) {

    check_mtx(x)

    res <- colFreqRatio(x)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colFreqRatios.data.frame <- function(x, ...) {

    check_df(x)

    colFreqRatios.matrix(as.matrix(x))

  }

#' @rdname colMedians
#' @export

  colPercUniques <- function(x, ...) UseMethod('colPercUniques')

#' @rdname colMedians
#' @export

  colPercUniques.matrix <- function(x, ...) {

    check_mtx(x)

    res <- colPercUnique(x)

    if(!is.null(colnames(x))) {

      res <- set_names(res, colnames(x))

    }

    res

  }

#' @rdname colMedians
#' @export

  colPercUniques.data.frame <- function(x, ...) {

    check_df(x)

    colPercUniques.matrix(as.matrix(x))

  }

#' @rdname colMedians
#' @export

  colQuantiles <- function(x, probs, ...) UseMethod('colQuantiles')

#' @rdname colMedians
#' @export

  colQuantiles.matrix <- function(x, probs, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    stopifnot(is.numeric(probs))

    if(any(probs > 1) | any(probs < 0)) {

      stop("'probs' must be numeric values between 0 and 1.", call. = FALSE)

    }

    res <- colQuant(x, probs, na.rm)

    if(!is.null(colnames(x))) {

      rownames(res) <- colnames(x)

    }

    colnames(res) <- probs

    res

  }

#' @rdname colMedians
#' @export

  colQuantiles.data.frame <- function(x, probs, na.rm = TRUE, ...) {

    check_df(x)

    colQuantiles.matrix(as.matrix(x), probs = probs, na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colCI <- function(x, conf_level, ...) UseMethod('colCI')

#' @rdname colMedians
#' @export

  colCI.matrix <- function(x,
                           conf_level = 0.95,
                           method = c('percentile', 'bca'),
                           na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    stopifnot(is.numeric(conf_level))

    conf_level <- conf_level[1]

    if(conf_level > 1 | conf_level < 0) {

      stop("'conf_level' must be a numeric value between 0 and 1.",
           call. = FALSE)

    }

    method <- match.arg(method[1], c('percentile', 'bca'))

    if(method == 'percentile') {

      res <- colPerCi(x, conf_level, na.rm)

    } else {

      res <- colBcaCi(x, conf_level, na.rm)

    }

    if(!is.null(colnames(x))) {

      rownames(res) <- colnames(x)

    }

    colnames(res) <- c('lower_ci', 'upper_ci')

    res

  }

#' @rdname colMedians
#' @export

  colCI.data.frame <- function(x,
                               conf_level = 0.95,
                               method = c('percentile', 'bca'),
                               na.rm = TRUE, ...) {

    check_df(x)

    colCI.matrix(as.matrix(x), conf_level = conf_level, na.rm = na.rm)

  }

#' @rdname colMedians
#' @export

  colMissing <- function(x, ...) UseMethod('colMissing')

#' @rdname colMedians
#' @export

  colMissing.matrix <- function(x, ...) {

    check_mtx(x)

    res <- colMis(x)

    if(!is.null(colnames(x))) res <- set_names(res, colnames(x))

    res

  }

#' @rdname colMedians
#' @export

  colMissing.data.frame <- function(x, ...) {

    check_df(x)

    colMissing.matrix(as.matrix(x))

  }

#' @rdname colMedians
#' @export

  distr_stats <- function(x,
                          unbiased = TRUE,
                          freqCut = 95/5,
                          uniqueCut = 10, ...) {

    ## the input check is done by downstream functions

    if(!is.data.frame(x) & !is.matrix(x)) {

      stop('Unsupported input.', call. = FALSE)

    } else if(is.data.frame(x)) {

      check_df(x)

    } else {

      check_mtx(x)

    }

    var <- NULL
    variable <- NULL
    gini_coef <- NULL
    percentUnique <- NULL
    zeroVar <- NULL
    n <- NULL
    n_missing <- NULL
    perc_missing <- NULL

    if(!is.null(colnames(x))) {

      var_names <- colnames(x)

    } else {

      var_names <- 1:ncol(x)

    }

    mean_variance <-
      tibble(variable = var_names,
             median = colMedians(x, na.rm = TRUE),
             mean = colMeans(x, na.rm = TRUE),
             var = colVars(x, na.rm = TRUE),
             gini_coef = colGini(x, unbiased = unbiased, na.rm = TRUE))

    mean_variance <- mutate(mean_variance,
                            var_mean_ratio = var/mean,
                            freqRatio = colFreqRatios(x, na.rm = TRUE),
                            percentUnique = colPercUniques(x, na.rm = TRUE),
                            zeroVar = percentUnique == 0,
                            nzv = ifelse(zeroVar,
                                         TRUE,
                                         freqRatio > freqCut | percentUnique <= uniqueCut),
                            n = nrow(x),
                            n_missing = colMissing(x),
                            perc_missing = n_missing/n * 100)

    mean_variance

  }

# Data frames and matrices: row distribution and variable selection -------

#' Compute distribution metrics for rows of a data frame or a matrix.
#'
#' @description
#' The functions compute row medians, variances, standard deviations, Gini
#' coefficients, quantiles, confidence intervals and metrics of uniqueness.
#' `row_stats` computes a bunch of distribution stats that
#' are helpful at selection of variant features for the further analysis, e.g.
#' differential gene expression.
#'
#' @references
#' DiCiccio TJ, Efron B. Bootstrap confidence intervals.
#' https://doi.org/101214/ss/1032280214 (1996) 11:189–228.
#' doi:10.1214/SS/1032280214
#'
#' @inheritParams colMedians
#'
#' @return `rowMedians()`, `rowMins()`, `rowMax()`, `rowGmeans()`,
#' `rowHmeans()`, `rowVars()`, `rowSDs()`, `rowGini()` return numeric vectors
#' with, respectively, row medians, row minima, row maxima, geometric means,
#' harmonic means, variances, standard deviations and Gini coefficients.
#' `rowQuantiles()`, and `rowCI()` return numeric matrices of quantiles and
#' confidence intervals of rows.
#' `rowFreqRatios()` returns ratios of frequencies of the first most common
#' element to the second one for each row; the function silently skips any
#' `NAs`.
#' `rowPercUniques()` returns percentages of unique observations in each row;
#' `NAs` are silently removed.
#'
#' `row_stats()` returns a data frame with the following columns:
#' * `variable`: if the input object has defined row names, they will be
#' displayed here, otherwise observation numbers.
#' * `median`: row medians
#' * `mean`: row means
#' * `var`: row variances
#' * `gini_coef`: Gini coefficients of the rows
#' * `var_mean_ratio`: variance to mean ratios
#' * `freqRatio`: ratios of frequencies of the first most common element to the
#' second most common element.
#' * `percentUnique`: percentage of unique observations
#' * `zeroVar`: logical, indicates if the row has only one unique value.
#' * `nzv`: logical, indicates if the row passes the variability criteria
#' defined by the `freqCut` and `uniqueCut` arguments.
#' * `n`: numbers of row entries
#' * `n_missing`: numbrers of missing row entries
#' * `perc_missing`: percentage of missing row entries
#' Note, that `row_stats()` silently removes missing values prior to
#' computation of the distribution statistics.
#'
#' @export

  rowMedians <- function(x, ...) UseMethod('rowMedians')

#' @rdname rowMedians
#' @export

  rowMedians.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- rowMed(x, na.rm)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowMedians.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    rowMedians.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowMins <- function(x, ...) UseMethod('rowMins')

#' @rdname rowMedians
#' @export

  rowMins.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- rowMi(x, na.rm)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowMins.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    rowMins.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowMax <- function(x, ...) UseMethod('rowMax')

#' @rdname rowMedians
#' @export

  rowMax.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- rowMa(x, na.rm)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowMax.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    rowMax.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowGmeans <- function(x, ...) UseMethod('rowGmeans')

#' @rdname rowMedians
#' @export

  rowGmeans.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- rowGeoMean(x, na.rm)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowGmeans.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    rowGmeans.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowHmeans <- function(x, ...) UseMethod('rowHmeans')

#' @rdname rowMedians
#' @export

  rowHmeans.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- rowHarmMean(x, na.rm)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowHmeans.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    rowHmeans.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowVars <- function(x, ...) UseMethod('rowVars')

#' @rdname rowMedians
#' @export

  rowVars.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- rowVariance(x, na.rm)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowVars.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    rowVars.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowSDs <- function(x, ...) UseMethod('rowSDs')

#' @rdname rowMedians
#' @export

  rowSDs.matrix <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- rowSD(x, na.rm)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowSDs.data.frame <- function(x, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    rowSDs.matrix(as.matrix(x), na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowGini <- function(x, unbiased, ...) UseMethod('rowGini')

#' @rdname rowMedians
#' @export

  rowGini.matrix <- function(x, unbiased = TRUE, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    res <- rowGi(x, unbiased, na.rm)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowGini.data.frame <- function(x, unbiased = TRUE, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_df(x)

    rowGini.matrix(as.matrix(x), na.rm = na.rm)


  }

#' @rdname rowMedians
#' @export

  rowFreqRatios <- function(x, ...) UseMethod('rowFreqRatios')

#' @rdname rowMedians
#' @export

  rowFreqRatios.matrix <- function(x, ...) {

    check_mtx(x)

    res <- rowFreqRatio(x)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowFreqRatios.data.frame <- function(x, ...) {

    check_df(x)

    rowFreqRatios.matrix(as.matrix(x))

  }

#' @rdname rowMedians
#' @export

  rowPercUniques <- function(x, ...) UseMethod('rowPercUniques')

#' @rdname rowMedians
#' @export

  rowPercUniques.matrix <- function(x, ...) {

    check_mtx(x)

    res <- rowPercUnique(x)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname rowMedians
#' @export

  rowPercUniques.data.frame <- function(x, ...) {

    check_df(x)

    rowPercUniques.matrix(as.matrix(x))

  }

#' @rdname rowMedians
#' @export

  rowQuantiles <- function(x, probs, ...) UseMethod('rowQuantiles')

#' @rdname rowMedians
#' @export

  rowQuantiles.matrix <- function(x, probs, na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    stopifnot(is.numeric(probs))

    if(any(probs > 1) | any(probs < 0)) {

      stop("'probs' must be numeric values between 0 and 1.", call. = FALSE)

    }

    res <- rowQuant(x, probs, na.rm)

    if(!is.null(rownames(x))) {

      rownames(res) <- rownames(x)

    }

    colnames(res) <- probs

    res

  }

#' @rdname rowMedians
#' @export

  rowQuantiles.data.frame <- function(x, probs, na.rm = TRUE, ...) {

    check_df(x)

    rowQuantiles.matrix(as.matrix(x), probs = probs, na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowCI <- function(x, conf_level, ...) UseMethod('rowCI')

#' @rdname rowMedians
#' @export

  rowCI.matrix <- function(x,
                           conf_level = 0.95,
                           method = c('percentile', 'bca'),
                           na.rm = TRUE, ...) {

    stopifnot(is.logical(na.rm))

    check_mtx(x)

    stopifnot(is.numeric(conf_level))

    conf_level <- conf_level[1]

    if(conf_level > 1 | conf_level < 0) {

      stop("'conf_level' must be a numeric value between 0 and 1.",
           call. = FALSE)

    }

    method <- match.arg(method[1], c('percentile', 'bca'))

    if(method == 'percentile') {

      res <- rowPerCi(x, conf_level, na.rm)

    } else {

      res <- rowBcaCi(x, conf_level, na.rm)

    }

    if(!is.null(rownames(x))) {

      rownames(res) <- rownames(x)

    }

    colnames(res) <- c('lower_ci', 'upper_ci')

    res

  }

#' @rdname rowMedians
#' @export

  rowCI.data.frame <- function(x,
                               conf_level = 0.95,
                               method = c('percentile', 'bca'),
                               na.rm = TRUE, ...) {

    check_df(x)

    rowCI.matrix(as.matrix(x), conf_level = conf_level, na.rm = na.rm)

  }

#' @rdname rowMedians
#' @export

  rowMissing <- function(x, ...) UseMethod('rowMissing')

#' @rdname rowMedians
#' @export

  rowMissing.matrix <- function(x, ...) {

    check_mtx(x)

    res <- rowMis(x)

    if(!is.null(rownames(x))) res <- set_names(res, rownames(x))

    res

  }

#' @rdname rowMedians
#' @export

  rowMissing.data.frame <- function(x, ...) {

    check_df(x)

    rowMissing.matrix(as.matrix(x))

  }

#' @rdname rowMedians
#' @export

  row_stats <- function(x,
                          unbiased = TRUE,
                          freqCut = 95/5,
                          uniqueCut = 10, ...) {

    ## the input check is done by downstream functions

    if(!is.data.frame(x) & !is.matrix(x)) {

      stop('Unsupported input.', call. = FALSE)

    } else if(is.data.frame(x)) {

      check_df(x)

    } else {

      check_mtx(x)

    }

    var <- NULL
    variable <- NULL
    gini_coef <- NULL
    percentUnique <- NULL
    zeroVar <- NULL
    n <- NULL
    n_missing <- NULL
    perc_missing <- NULL

    if(!is.null(rownames(x))) {

      var_names <- rownames(x)

    } else {

      var_names <- 1:nrow(x)

    }

    mean_variance <-
      tibble(variable = var_names,
             median = rowMedians(x, na.rm = TRUE),
             mean = rowMeans(x, na.rm = TRUE),
             var = rowVars(x, na.rm = TRUE),
             gini_coef = rowGini(x, unbiased = unbiased, na.rm = TRUE))

    mean_variance <- mutate(mean_variance,
                            var_mean_ratio = var/mean,
                            freqRatio = rowFreqRatios(x, na.rm = TRUE),
                            percentUnique = rowPercUniques(x, na.rm = TRUE),
                            zeroVar = percentUnique == 0,
                            nzv = ifelse(zeroVar,
                                         TRUE,
                                         freqRatio > freqCut | percentUnique <= uniqueCut),
                            n = ncol(x),
                            n_missing = rowMissing(x),
                            perc_missing = n_missing/n * 100)

    mean_variance

  }

# Minimum shift of numeric data ----------

#' Offset numeric vectors be their minimum.
#'
#' @description
#' Offsets a numeric vectors or columns of a numeric matrix or a data frame by
#' the respective minimal values.
#' This results in non-negative values with minimum corresponding to 0.
#' Any `NA`s are removed silently for minimum calculation but kept in the output.
#'
#' @return a numeric vector, matrix or data frame, depending on the input.
#'
#' @param x an object: numeric vector, matrix or data frame.
#' @param ... exra arguments passed to the methods.
#'
#' @export

  minimum_shift <- function(x, ...) UseMethod('minimum_shift')

#' @rdname minimum_shift
#' @export

  minimum_shift.default <- function(x, ...) {

    stopifnot(is.numeric(x))

    offsetVector(x)

  }

#' @rdname minimum_shift
#' @export

  minimum_shift.matrix <- function(x, ...) {

    check_mtx(x)

    res <- offsetMatrix(x)

    if(!is.null(rownames(x))) {

      res <- set_names(res, rownames(x))

    }

    res

  }

#' @rdname minimum_shift
#' @export

  minimum_shift.data.frame <- function(x, ...) {

    check_df(x)

    as.data.frame(minimum_shift(as.matrix(x)))

  }

# Numbers of complete cases in columns and rows ---------

#' Numbers of complete cases in rows and columns of a matrix or a data frame.
#'
#' @description
#' The functions calculate numbers of complete, i.e. non-`NA` cases in columns
#' or rows of a data frame or matrix.
#'
#' @return a numeric vector.
#'
#' @param x a data frame or a matrix.
#' @param ... extra arguments passed to methods; currently none.
#'
#' @export

  colComplete <- function(x, ...) UseMethod('colComplete')

#' @rdname colComplete
#' @export

  rowComplete <- function(x, ...) UseMethod('rowComplete')

#' @rdname colComplete
#' @export

  colComplete.data.frame <- function(x, ...) {

    stopifnot(is.data.frame(x))

    map_dbl(x, ~length(.x[!is.na(.x)]))

  }

#' @rdname colComplete
#' @export

  colComplete.matrix <- function(x, ...) {

    stopifnot(is.matrix(x))

    colComplete(as.data.frame(x))

  }

#' @rdname colComplete
#' @export

  rowComplete.data.frame <- function(x, ...) {

    stopifnot(is.data.frame(x))

    n_completes <- map(1:nrow(x), function(index) unlist(x[index, ]))

    map_dbl(n_completes, ~length(.x[!is.na(.x)]))

  }

#' @rdname colComplete
#' @export

  rowComplete.matrix <- function(x, ...) {

    stopifnot(is.matrix(x))

    n_completes <- map(1:nrow(x), function(index) x[index, ])

    map_dbl(n_completes, ~length(.x[!is.na(.x)]))

  }

# Fast computation of descriptive statistics of numeric variables --------

#' Fast computation of descriptive distribution statistics for numeric features.
#'
#' @description
#' The function takes a data frame with numeric variables and a splitting
#' variable in factor format, and returns a bunch of statistics of
#' distribution of the numeric variables such as median with interquartile
#' range and range, or mean, standard deviation and range.
#' Numbers of complete observations are returned as well.
#' The output is publication styled in a way that is accepted with most
#' scientific journals.
#'
#' @details
#' If any `NA` are detected, cells of the output data frame are appendend with
#' the information on numbers of complete observations.
#'
#' @return a data frame.
#'
#' @param data a data frame as describe above.
#' @param split_fct name of the factor variable used for splitting.
#' @param variables a vector of names of variables for which the statistics
#' will be computed. If `NULL`, all variables except of `split_factor` are
#' used.
#' @param format format of the output data frame cells. `median` (default):
#' median, interquartile range, and range; `mean`: mean, standard deviation,
#' and range; `full`: combination of the `mean` and `median` formats.
#' @param signif_digits significant digits used for the output.
#'
#' @export

  fast_num_stats <- function(data,
                             split_fct,
                             variables = NULL,
                             format = c('median', 'mean', 'full'),
                             signif_digits = 2) {

    ## entry control --------

    if(!is.data.frame(data)) {

      stop("'data' has to be a data frame.", call. = FALSE)

    }

    if(!split_fct %in% names(data)) {

      stop("'split_fct' is missing from the inpute data.", call. = FALSE)

    }

    if(is.null(variables)) {

      variables <- names(data)[names(data) != split_fct]

    }

    if(any(!variables %in% names(data))) {

      stop("Some of variables are missing from 'data'.", call. = FALSE)

    }

    data <- filter(data,
                   !is.na(.data[[split_fct]]))

    split_vec <- data[[split_fct]]

    data <- data[variables]

    err_txt <- "'split_fct' has to be a factor with two or more levels."

    if(!is.factor(split_vec)) stop(err_txt, call. = FALSE)

    if(length(levels(split_vec)) < 2) stop(err_txt, call. = FALSE)

    classes <- map_lgl(data, is.numeric)

    if(any(!classes)) {

      stop('At least one of the variables is not numeric.', call. = FALSE)

    }

    format <- match.arg(format[1], c('median', 'mean', 'full'))

    stopifnot(is.numeric(signif_digits))

    signif_digits <- as.integer(signif_digits)

    ## checking if there are any NAs in the data frame
    ## this will have an effect on presence of numbers of complete
    ## onservations in each cell of the output data frame

    na_index <- any(is.na(as.matrix(data)))

    ## stat functions ------

    tbl_cell <- NULL
    variable <- NULL

    medians <- NULL
    means <- NULL
    lower_quart <- NULL
    upper_quart <- NULL

    sds <- NULL

    median_stat <- function(x, sig) {

      stats <- tibble(variable = names(x),
                      medians = colMedians(x),
                      lower_quart = colQuantiles(x, 0.25),
                      upper_quart = colQuantiles(x, 0.75),
                      min = colMins(x),
                      max = colMax(x))

      stats <-
        mutate(stats,
               tbl_cell = paste0(signif(medians, sig),
                                 ' [IQR: ', signif(lower_quart, sig),
                                 ' - ', signif(upper_quart, sig),']\nrange: ',
                                 signif(min, sig), ' - ', signif(max, sig)))

      stats

    }

    mean_stat <- function(x, sig) {

      stats <- tibble(variable = names(x),
                      means = colMeans(x),
                      sds = colSDs(x),
                      min = colMins(x),
                      max = colMax(x))

      stats <-
        mutate(stats,
               tbl_cell = paste0(signif(means, sig),
                                 ' [SD: ', signif(sds, sig),
                                 ']n\nrange: ',
                                 signif(min, sig), ' - ', signif(max, sig)))

      stats

    }

    full_stat <- function(x, sig) {

      stats <- tibble(variable = names(x),
                      means = colMeans(x),
                      sds = colSDs(x),
                      medians = colMedians(x),
                      lower_quart = colQuantiles(x, 0.25),
                      upper_quart = colQuantiles(x, 0.75),
                      min = colMins(x),
                      max = colMax(x))

      stats <-
        mutate(stats,
               tbl_cell = paste0('mean: ', signif(means, sig),
                                 ' [SD: ', signif(sds, sig),
                                 ']\nmedian: ', signif(medians, sig),
                                 ' [IQR: ', signif(lower_quart, sig),
                                 ' - ', signif(upper_quart, sig),
                                 ']n\nrange: ',
                                 signif(min, sig), ' - ', signif(max, sig)))

      stats

    }

    stat_fun <- switch(format,
                       median = function(x) median_stat(x, signif_digits),
                       mean = function(x) mean_stat(x, signif_digits),
                       full = function(x) full_stat(x, signif_digits))

    ## computation -------

    data_splits <- split(data, f = split_vec, drop = TRUE)

    levs <- names(data_splits)

    stat_tbl <- map(data_splits, stat_fun)

    if(na_index) {

      n_completes <- map(data_splits, colComplete)

      stat_tbl <-
        map2(stat_tbl,
             n_completes,
             ~mutate(.x,
                     n = .y,
                     tbl_cell = paste(tbl_cell, n,
                                      sep = '\ncomplete: n = ')))

    }

    stat_tbl <- map(stat_tbl, ~.x[c('variable', 'tbl_cell')])

    stat_tbl <- reduce(stat_tbl, left_join, by = 'variable')

    stat_tbl <- set_names(stat_tbl, c('variable', levs))

    ## appending with the n numbers of samples ------

    n_numbers <- count(mutate(data,
                              !!split_fct := split_vec),
                       .data[[split_fct]])

    n_numbers <-
      as_tibble(t(column_to_rownames(n_numbers, split_fct)))

    n_numbers <- mutate(n_numbers,
                        variable = 'Samples, N')

    n_numbers <- relocate(n_numbers, variable)

    rbind(n_numbers, stat_tbl)

  }

# Differences of centrality statistics --------

#' Difference of averages.
#'
#' @description
#' Calculates differences between avarages of numeric variables split by a
#' factor  variable.
#'
#' @param data a data frame.
#' @param split_fct name of the factor variable used for splitting.
#' @param variables a vector of names of variables for which the statistics
#' will be computed. If `NULL`, all variables except of `split_factor` are
#' used.
#' @param average_fun a function for averaging numeric variables.
#' @param extreme logical, should difference between the last and first level
#' of the splitting factor be returned. If `TRUE` (default), a single numeric
#' vector of differences of averages between the last and fist levels of the
#' splitting factor is returned. If `FALSE`, a list of numeric vector is
#' returned; each of the vector contains differences between the particular
#' level and the first level of the splitting factor.
#'
#' @return a numeric vector or a list of numeric vectors.
#'
#' @export

  delta <- function(data,
                    split_fct,
                    variables = NULL,
                    average_fun = colMeans,
                    extreme = TRUE) {

    ## entry control -----

    stopifnot(is.data.frame(data))

    if(!split_fct %in% names(data)) {

      stop("The splitting factor is absent from data.", call. = FALSE)

    }

    err_txt <-
      "'split_fct' has to define a factor variable with at least two levels,"

    if(!is.factor(data[[split_fct]])) stop(err_txt, call. = FALSE)

    data <- mutate(data, !!split_fct := droplevels(.data[[split_fct]]))

    split_vec <- data[[split_fct]]

    n_levels <- length(levels(split_vec))

    if(n_levels < 2) stop(err_txt, call. = FALSE)

    levs <- levels(split_vec)

    if(extreme) {

      levs <- levs[c(1, length(levs))]

      data <- filter(data, .data[[split_fct]] %in% levs)

      data <- mutate(data, !!split_fct := droplevels(.data[[split_fct]]))

      split_vec <- data[[split_fct]]

    }

    if(any(!variables %in% names(data))) {

      stop('Some of the variables are missing from the data.', call. = FALSE)

    }

    data <- data[variables]

    classes <- map_lgl(data, is.numeric)

    if(any(!classes)) {

      stop('Some of the variables are not numeric.', call. = FALSE)

    }

    stopifnot(is.function(average_fun))
    stopifnot(is.logical(extreme))

    ## calculation --------

    data_splits <- split(data, split_vec, drop = TRUE)

    data_splits <- map(data_splits, average_fun)

    if(extreme | n_levels == 2) return(data_splits[[2]] - data_splits[[1]])

    diff_value <-
      map(2:n_levels,
          ~data_splits[[.x]] - data_splits[[1]])

    set_names(diff_value, levs[-1])

  }

# END -----
