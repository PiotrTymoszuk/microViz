# Functions for normalization of a vector, matrix and data frame

# Min-Max scaling -------

#' Minimum - maximum scaling.
#'
#' @description
#' Computes minimum - maximum scaled values of a numeric vectors, columns or
#' rows of a numeric matrix, or columns of a data frame.
#'
#' @details
#' The scaled values of a vector `x` are computed with the following formula:
#'
#' \deqn{Sc(x_i) = \frac{x_i - min(x)}{max(i) - min(i)}}
#'
#' @return a numeric vector, numeric matrix or a data frame. Row and column
#' names are preserved.
#'
#' @param x a numeric vector, numeric matrix or a data frame.
#' @param feature_type specifies if columns (default) or rows are subjected to
#' scaling.
#' @param variables a vector of names of numeric variables to be normalized.
#' If `NULL`, all variables of a data frame are scaled. Please note, that
#' variable selection makes the computation slower.
#' @param ... extra arguments passed to methods.
#'
#' @export

  minMax <- function(x, ...) UseMethod('minMax')

#' @rdname minMax
#' @export

  minMax.default <- function(x, ...) {

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric vector.", call. = FALSE)

    }

    if(length(x) <= 1) return(x)

    res <- minMaxVec(x)

    if(!is.null(names(x))) return(set_names(res, names(x)))

    res

  }

#' @rdname minMax
#' @export

  minMax.matrix <- function(x,
                            feature_type = c('columns', 'rows'), ...) {

    feature_type <- match.arg(feature_type[1], c('columns', 'rows'))

    check_mtx(x)

    if(feature_type == 'columns') {

      return(minMaxCols(x))

    } else {

      return(minMaxRows(x))

    }

  }

#' @rdname minMax
#' @export

  minMax.data.frame <- function(x,
                                variables = NULL, ...) {

    stopifnot(is.data.frame(x))

    ## no variable selection, the fastest option

    if(is.null(variables)) {

      check_df(x)

      if(!is_tibble(x)) {

        return(as.data.frame(minMaxCols(as.matrix(x))))

      } else {

        return(as_tibble(as.data.frame(minMaxCols(as.matrix(x)))))

      }

    }

    if(any(!variables %in% names(x))) {

      stop("Some variables are missing from the data frame 'x'.",
           call. = FALSE)

    }

    num_part <- x[, variables, drop = FALSE]

    check_df(num_part)

    scaled_part <- minMaxCols(as.matrix(num_part))

    x[variables] <- scaled_part

    x

  }

# Z-scores -----------

#' Normalization, standardization and Z-scores.
#'
#' @description
#' Calculates classical Z-scores or their variants employing median, geometric
#' mean, harmonic mean as centrality statistics, and SD or SEM as dispersion
#' statistics.
#'
#' @details
#' The Z-scores are calculated for a numeric vector `x` with a general formula:
#'
#' \deqn{Z(x_i) = \frac{x_i - centFunc(x)}{dispFunc(x)}}
#'
#' where \eqn{centFunc} stands for a function that returns a centrality
#' statistic such as mean or median, and \eqn{dispFunc(x)} represents a
#' function that returns a dispersion statistic such as standard deviation (SD)
#' or standard error of the mean (SEM).
#'
#' @return a numeric vector, numeric matrix or a data frame.
#'
#' @inheritParams minMax
#' @param center a string that specifies the centrality statistic: mean, median,
#' geometric mean or harmonic mean.
#' @param dispersion a string that specified the dispersion statistic:
#' SD or SEM.
#' @param ... extra arguments passed to methods.
#'
#' @export

  zScores <- function(x, ...) UseMethod('zScores')

#' @rdname zScores
#' @export

  zScores.default <- function(x,
                              center = c('mean',
                                         'median',
                                         'geo_mean',
                                         'harm_mean'),
                              dispersion = c('sd', 'sem'), ...) {

    if(!is.numeric(x)) {

      stop("'x' has to be a numeric vector.", call. = FALSE)

    }

    center <- match.arg(center[1],
                        c('mean',
                          'median',
                          'geo_mean',
                          'harm_mean'))

    dispersion <- match.arg(dispersion[1], c('sd', 'sem'))

    if(length(x) <= 1) return(x)

    res <- zScoreVec(x, center, dispersion)

    if(!is.null(names(x))) return(set_names(res, names(x)))

    res

  }

#' @rdname zScores
#' @export

  zScores.matrix <- function(x,
                            feature_type = c('columns', 'rows'),
                            center = c('mean',
                                       'median',
                                       'geo_mean',
                                       'harm_mean'),
                            dispersion = c('sd', 'sem'), ...) {

    feature_type <- match.arg(feature_type, c('columns', 'rows'))

    center <- match.arg(center[1],
                        c('mean',
                          'median',
                          'geo_mean',
                          'harm_mean'))

    dispersion <- match.arg(dispersion[1], c('sd', 'sem'))

    check_mtx(x)

    if(feature_type == 'columns') {

      return(zScoreCols(x, center, dispersion))

    } else {

      return(zScoreRows(x, center, dispersion))

    }

  }

#' @rdname zScores
#' @export

  zScores.data.frame <- function(x,
                                 variables = NULL,
                                 center = c('mean',
                                            'median',
                                            'geo_mean',
                                            'harm_mean'),
                                 dispersion = c('sd', 'sem'), ...) {

    stopifnot(is.data.frame(x))

    center <- match.arg(center[1],
                        c('mean',
                          'median',
                          'geo_mean',
                          'harm_mean'))

    dispersion <- match.arg(dispersion[1], c('sd', 'sem'))

    ## no variable selection, the fastest option

    if(is.null(variables)) {

      check_df(x)

      if(!is_tibble(x)) {

        return(as.data.frame(zScoreCols(as.matrix(x),
                                        center,
                                        dispersion)))

      } else {

        return(as_tibble(as.data.frame(zScoreCols(as.matrix(x),
                                                  center,
                                                  dispersion))))

      }

    }

    if(any(!variables %in% names(x))) {

      stop("Some variables are missing from the data frame 'x'.",
           call. = FALSE)

    }

    num_part <- x[, variables, drop = FALSE]

    check_df(num_part)

    scaled_part <- zScoreCols(as.matrix(num_part), center, dispersion)

    x[variables] <- scaled_part

    x

  }


# END -------
