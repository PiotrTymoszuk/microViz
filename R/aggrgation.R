# Functions used to aggregate rows and columns by a statistic

#' Aggregate rows or columns by a statistic.
#'
#' @description
#' The functions aggregate rows or columns of a data frame or matrix that
#' correspond to the same feature (e.g. microarray probes or transcript
#' variants that map to the same gene) by a statistic such as arithmetic mean,
#' geometric mean, harmonic mean, median, minimum, maximum or sum.
#'
#' @return a data frame or a matrix, dependent on the input. Rows are named
#' after unique elements of `f`
#'
#' @param x a data frame or a matrix.
#' @param f a vector that defines elements to be aggregated, e.g. a vector
#' of gene symbols.
#' @param fun a function that takes a matrix and returns a vector of column
#' (for `aggRows()`) or row (for `aggCols()`) statistics.
#' @param verbose logical, should messages and execution time be printed?
#' @param .parallel should the computation be run in parallel? experimental!
#' @param .paropts options of the parallel workers, as specified by
#' \code{\link[furrr]{furrr_options}}.
#' @param ... additional arguments passed to `fun`.
#'
#' @export

  aggRows <- function(x,
                      f,
                      fun = colMeans,
                      verbose = TRUE,
                      .parallel = FALSE,
                      .paropts = furrr::furrr_options(seed = TRUE),
                      ...) {

    ## entry control ----------

    if(!is.matrix(x) & !is.data.frame(x)) {

      stop('Unsupported input data.', call. = FALSE)

    } else if(is.matrix(x)) {

      check_mtx(x)

      init_mtx <- TRUE

    } else {

      check_df(x)

      init_mtx <- FALSE

    }

    x <- as.data.frame(x)

    if(!is.atomic(f)) {

      stop("'f' has to be an atomic vector.", call. = FALSE)

    }

    if(length(f) != nrow(x)) {

      stop("The length of 'f' must match the number of rows.",
           call. = FALSE)

    }

    stopifnot(is.function(fun))
    stopifnot(is.logical(.parallel))

    ## parallel backend ---------

    start_time <- Sys.time()

    if(.parallel) plan('multisession')

    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))
    on.exit(plan('sequential'), add = TRUE)

    ## splitting --------

    duplicated_f <- unique(f[duplicated(f)])

    unique_f <- f[!f %in% duplicated_f]

    if(length(duplicated_f) == 0) {

      if(init_mtx) return(as.matrix(x)) else return(x)

    }

    x_splits <- split(x, f, drop = TRUE)

    ## keeping separate the elements that do not need to be aggregated
    ## aggregation --------

    multi_splits <- future_map(x_splits[duplicated_f],
                               fun, ...,
                               .options = .paropts)

    multi_splits <- do.call('rbind', multi_splits[duplicated_f])

    rownames(multi_splits) <- duplicated_f

    if(length(unique_f) > 0) {

      single_splits <- do.call('rbind', x_splits[unique_f])

      rownames(single_splits) <- unique_f

      multi_splits <- rbind(single_splits, multi_splits)

    }

    if(init_mtx) {

      return(as.matrix(multi_splits))

    } else {

      return(multi_splits)

    }

  }

#' @rdname aggRows
#' @export

  aggCols <- function(x,
                      f,
                      fun = rowMeans,
                      verbose = TRUE,
                      .parallel = FALSE,
                      .paropts = furrr::furrr_options(seed = TRUE),
                      ...) {

    ## entry control ----------

    if(!is.matrix(x) & !is.data.frame(x)) {

      stop('Unsupported input data.', call. = FALSE)

    } else if(is.matrix(x)) {

      check_mtx(x)

      init_mtx <- TRUE

    } else {

      check_df(x)

      init_mtx <- FALSE

    }

    x <- as.data.frame(x)

    if(!is.atomic(f)) {

      stop("'f' has to be an atomic vector.", call. = FALSE)

    }

    if(length(f) != ncol(x)) {

      stop("The length of 'f' must match the number of columns.",
           call. = FALSE)

    }

    stopifnot(is.function(fun))
    stopifnot(is.logical(.parallel))

    ## parallel backend ---------

    start_time <- Sys.time()

    if(.parallel) plan('multisession')

    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))
    on.exit(plan('sequential'), add = TRUE)

    ## splitting -------

    duplicated_f <- unique(f[duplicated(f)])

    unique_f <- f[!f %in% duplicated_f]

    if(length(duplicated_f) == 0) {

      if(init_mtx) return(as.matrix(x)) else return(x)

    }

    split_lst <- split(1:length(f), f, drop = TRUE)

    split_x <- map(split_lst, ~x[, .x, drop = FALSE])

    multi_splits <- future_map(split_x[duplicated_f],
                               fun, ...,
                               .options = .paropts)

    multi_splits <- do.call('cbind', multi_splits[duplicated_f])

    colnames(multi_splits) <- duplicated_f

    if(length(unique_f) > 0) {

      single_splits <- do.call('cbind', split_x[unique_f])

      colnames(single_splits) <- unique_f

      multi_splits <- cbind(single_splits, multi_splits)

    }

    if(!is.null(rownames(x))) {

      rownames(multi_splits) <- rownames(x)

    }

    if(init_mtx) {

      return(as.matrix(multi_splits))

    } else {

      return(as.data.frame(multi_splits))

    }


  }

# END -----
