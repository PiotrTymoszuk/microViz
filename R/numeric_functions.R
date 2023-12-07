# functions for calculation of numeric stats

# Variable distribution and variable selection ------

#' Compute statistics for numeric columns of a data frame.
#'
#' @description
#' The functions compute column medians, variances, standard deviations and
#' Gini coefficients. `distr_stats` computes a bunch of distribution stats that
#' are helpful at selection of variant features for the further analysis, e.g.
#' differential gene expression.
#'
#' @param x a numeric data frame.
#' @param na.rm logical: should NAs be removed?
#' @param .parallel logical, should the computation be done in parallel? Makes
#' sense for large datasets.
#' @param ... extra arguments passed to \code{\link[caret]{nearZeroVar}}.
#'
#' @return `colMedians()`, `colGmeans()`, `colHmeans()`, `colVars()`,
#' `colSDs()`, `colGini()` return numeric vectors with, respectively,
#' column medians, geometric means, harmonic means, variances, standard
#' deviations and Gini coefficients.
#' `distr_stats()` returns a data frame with variable names, their means,
#' variances, Gini coefficients, variance to mean ratios and the output of the
#' `caret` function \code{\link[caret]{nearZeroVar}}.
#'
#' @export

  colMedians <- function(x, na.rm = TRUE, .parallel = FALSE) {

    stopifnot(is.logical(na.rm))
    stopifnot(is.logical(.parallel))

    err_txt <- 'A numeric data frame is required'

    if(!is.data.frame(x)) stop(err_txt, call. = FALSE)

    col_check <- purrr::map_lgl(x, is.numeric)

    if(any(!col_check)) stop(err_txt, call. = FALSE)

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    furrr::future_map_dbl(x, stats::median, na.rm = na.rm,
                          .options = furrr::furrr_options(seed = TRUE))

  }

#' @rdname colMedians
#' @export

  colGmeans <- function(x, na.rm = TRUE, .parallel = FALSE) {

    stopifnot(is.logical(na.rm))
    stopifnot(is.logical(.parallel))

    err_txt <- 'A numeric data frame is required'

    if(!is.data.frame(x)) stop(err_txt, call. = FALSE)

    col_check <- purrr::map_lgl(x, is.numeric)

    if(any(!col_check)) stop(err_txt, call. = FALSE)

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    if(na.rm) {

      x_list <- purrr::map(x, ~.x[!is.na(.x)])

    } else {

      x_list <- as.list(x)

    }

    furrr::future_map_dbl(x_list,
                          ~prod(.x)^(1/length(.x)),
                          .options = furrr::furrr_options(seed = TRUE))

  }

#' @rdname colMedians
#' @export

  colHmeans <- function(x, na.rm = TRUE, .parallel = FALSE) {

    stopifnot(is.logical(na.rm))
    stopifnot(is.logical(.parallel))

    err_txt <- 'A numeric data frame is required'

    if(!is.data.frame(x)) stop(err_txt, call. = FALSE)

    col_check <- purrr::map_lgl(x, is.numeric)

    if(any(!col_check)) stop(err_txt, call. = FALSE)

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    furrr::future_map_dbl(x,
                          ~(1/mean(1/.x, na.rm = na.rm)),
                          .options = furrr::furrr_options(seed = TRUE))

  }

#' @rdname colMedians
#' @export

  colVars <- function(x, na.rm = TRUE, .parallel = FALSE) {

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    stopifnot(is.logical(na.rm))
    stopifnot(is.logical(.parallel))

    err_txt <- 'A numeric data frame is required'

    if(!is.data.frame(x)) stop(err_txt, call. = FALSE)

    col_check <- purrr::map_lgl(x, is.numeric)

    if(any(!col_check)) stop(err_txt, call. = FALSE)

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    furrr::future_map_dbl(x, stats::var, na.rm = na.rm,
                          .options = furrr::furrr_options(seed = TRUE))

  }

#' @rdname colMedians
#' @export

  colSDs <- function(x, na.rm = TRUE, .parallel = FALSE) {

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    stopifnot(is.logical(na.rm))
    stopifnot(is.logical(.parallel))

    err_txt <- 'A numeric data frame is required'

    if(!is.data.frame(x)) stop(err_txt, call. = FALSE)

    col_check <- purrr::map_lgl(x, is.numeric)

    if(any(!col_check)) stop(err_txt, call. = FALSE)

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    furrr::future_map_dbl(x, stats::sd, na.rm = na.rm,
                          .options = furrr::furrr_options(seed = TRUE))

  }

#' @rdname colMedians
#' @export

  colGini <- function(x, na.rm = TRUE, .parallel = FALSE) {

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    stopifnot(is.logical(na.rm))
    stopifnot(is.logical(.parallel))

    err_txt <- 'A numeric data frame is required'

    if(!is.data.frame(x)) stop(err_txt, call. = FALSE)

    col_check <- purrr::map_lgl(x, is.numeric)

    if(any(!col_check)) stop(err_txt, call. = FALSE)

    neg_check <- purrr::map_lgl(x, ~all(.x < 0))

    if(any(neg_check)) {

      warning(paste('Your data frame contains negative values, Gini',
                    'coefficient may not be computed correctly.'),
              call = FALSE)

    }

    if(na.rm) {

      x <- purrr::map(x, ~.x[!is.na(.x)])

    }

    if(.parallel) future::plan('multisession')

    on.exit(future::plan('sequential'))

    ginis <- furrr::future_map_dbl(x, DescTools::Gini, conf.level = NA,
                                   .options = furrr::furrr_options(seed = TRUE))

    purrr::map_dbl(ginis, ~.x[1])

  }

#' @rdname colMedians
#' @export

  distr_stats <- function(x, na.rm = TRUE, ..., .parallel = FALSE) {

    ## the input check is done by downstream functions

    var <- NULL

    mean_variance <-
      tibble::tibble(variable = names(x),
                     mean = colMeans(x),
                     var = colVars(x,
                                   na.rm = na.rm,
                                   .parallel = .parallel),
                     gini_coef = colGini(x,
                                         na.rm = na.rm,
                                         .parallel = .parallel))

    mean_variance <- dplyr::mutate(mean_variance,
                                   var_mean_ratio = var/mean)

    near_zero <- caret::nearZeroVar(x, saveMetrics = TRUE, ...)

    near_zero <- tibble::rownames_to_column(near_zero, 'variable')

    dplyr::left_join(mean_variance, near_zero, by = 'variable')

  }

# END -----
