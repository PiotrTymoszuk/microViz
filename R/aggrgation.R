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

    multi_splits <- compact(x_splits[duplicated_f])

    duplicated_f <- names(multi_splits)

    multi_splits <- future_map(multi_splits[duplicated_f],
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

    multi_splits <- compact(split_x[duplicated_f])

    duplicated_f <- names(multi_splits)

    multi_splits <- future_map(multi_splits[duplicated_f],
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

# Collapsing expression of duplicated features -------

#' Collapse values of duplicated features by a function.
#'
#' @description
#' The function collapses duplicated features such as genes or proteins by a
#' user-provided function that computes returns a vector of numeric values
#' for a matrix or data frame columns (such as \code{\link[base]{colMeans}} or
#' \code{\link{colMedians}}).
#'
#' @details
#' The function takes a matrix or a data frames with defined row and column
#' names.
#' Features such as genes of proteins are placed in rows, samples or
#' observations are provided in columns. Note: row names should contain unique
#' feature identifiers such as microarray probe IDs.
#' The assignment of feature IDs to gene symbols or other names is defined by
#' a data frame with two columns:
#'
#' * a columns with the feature identifiers corresponding to the row names of
#' the matrix.
#'
#' * a column with gene, protein or other names.
#'
#' This data frame is provided as the `annotation` argument.
#'
#' @return
#' If `transpose = FALSE`, a numeric matrix is returned with row names
#' corresponding to the second column of the `annotation` data frame, and
#' columns corresponding to the samples.
#' If `transpose = TRUE`, a data frame is returned with samples in the rows and
#' features in columns and named after the content of the second column of
#' `annotation`, sample identifiers are stored in the `sample_id` column of
#' the output.
#'
#' @param data a numeric matrix with features in rows and samples in columns.
#' @param annotation a data frames with assignment of the features in data
#' (first column) to gene, protein or other symbols/identifiers (second column).
#' @param fun a function to handle column of duplicated features,
#' `colMean` by default.
#' @param transpose logical, should the output be returned as a data frame with
#' observation in rows and features in columns?
#' @param trans_fun a function to transform the output numeric variables,
#' defaults to `identity`.
#' @param ... extra arguments passed to \code{\link{aggRows}} and `fun`.
#'
#' @export

  integrate_expression <- function(data,
                                   annotation,
                                   fun = colMeans,
                                   transpose = TRUE,
                                   trans_fun = identity, ...) {

    ## entry control ------

    err_txt <- "'data' has to be a numeric matrix or a data frame."

    if(!is.matrix(data) & !is.data.frame(data)) stop(err_txt, call. = FALSE)

    if(is.matrix(data)) {

      if(!is.numeric(data)) stop(err_txt, call. = FALSE)

    } else {

      col_check <- map_lgl(data, is.numeric)

      if(any(!col_check)) stop(err_txt, call. = FALSE)

      data <- as.matrix(data)

    }

    stopifnot(is.function(fun))
    stopifnot(is.function(trans_fun))

    stopifnot(is.logical(transpose))

    err_txt <- "'annotation' has to be a data frame with at least two columns."

    if(!is.data.frame(annotation)) stop(err_txt, call. = FALSE)

    if(ncol(annotation) < 2) stop(err_txt, call. = FALSE)

    if(all(c('probe_id', 'gene_symbol') %in% names(annotation))) {

      anotation <- annotation[, c('probe_id', 'gene_symbol')]

    } else {

      annotation <- set_names(annotation[, 1:2],
                              c('probe_id', 'gene_symbol'))

    }

    annotation <- filter(annotation, complete.cases(annotation))

    ## integration of values: only for the duplicated features ------

    gene_symbol <- NULL

    dupl_genes <-
      filter(annotation, duplicated(gene_symbol))

    dupl_genes <- unique(dupl_genes$gene_symbol)

    single_annot <-
      filter(annotation, !gene_symbol %in% dupl_genes)

    dupl_annot <-
      filter(annotation, gene_symbol %in% dupl_genes)

    if(nrow(dupl_annot) == 0) {

      data <- data[annotation$probe_id, ]

      rownames(data) <- annotation$gene_symbol

    } else {

      agg_data <- aggRows(x = data[dupl_annot$probe_id, ],
                          f = dupl_annot$gene_symbol,
                          fun = fun,
                          na.rm = TRUE, ...)

      if(nrow(single_annot) == 0) {

        data <- agg_data

      } else {

        single_data <- data[single_annot$probe_id, ]

        rownames(single_data) <- single_annot$gene_symbol

        data <- rbind(single_data, agg_data)

      }

    }

    data <- trans_fun(data)

    if(!transpose) return(data)

    data <- as.data.frame(t(data))

    data <- rownames_to_column(data, 'sample_id')

    as_tibble(data)

  }

# END -----
