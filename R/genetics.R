# Transformation of genetic data

# Implementation of genetic information in R ---------

#' Binarize somatic mutation data.
#'
#' @description
#' Converts a data frame with detailed mutation information
#' (e.g. nucleotide sequence, variant classification) into a simplified
#' binarized form.
#' In the output, presence of at least one mutation in a particular gene is
#' coded as 1 or 'mutated'. Absence of mutations in a gene is
#' coded as 0 or 'WT'.
#' Gene identifiers are stored in columns, sample identifiers in rows.
#'
#' @details
#' The function bottlenecks will run in parallel if a parallel backend is
#' declared with \code{\link[future]{plan}}.
#'
#'
#' @return a data frame. The variable `sample_id` that stores sample
#' identifiers. The variables named after gene identifiers retrieved from the
#' `gene_id` variable of data store binarized indexes of presence of at least
#' one somatic mutation.
#'
#' @param data a data frame with detailed mutation information obtained e.g.
#' from cBioportal.
#' @param sample_id name of the variable in `data` that stores sample
#' identifiers.
#' @param gene_id name of the variable in `data` that stores gene identifiers
#' such as official HUGO symbols or Entrez IDs.
#' @param as_numeric logical, should absence/presence of mutations be 0/1 coded?
#'
#' @export

  extract_mutations <- function(data,
                                sample_id = 'Tumor_Sample_Barcode',
                                gene_id = 'Hugo_Symbol',
                                as_numeric = TRUE) {

    ## entry control ---------

    stopifnot(is.character(sample_id))
    stopifnot(is.character(gene_id))
    stopifnot(is.logical(as_numeric))

    if(!is.data.frame(data)) {

      stop("'data' has to be a data frame.", call. = FALSE)

    }

    if(!sample_id %in% names(data)) {

      stop("'sample_id' absent from the data frame.", call. = FALSE)

    }

    if(!gene_id %in% names(data)) {

      stop("'gene_id' absent from the data frame.", call. = FALSE)

    }

    ## counts mutations per gene and sample -------

    data[[gene_id]] <- factor(data[[gene_id]])

    data_splits <- split(data[[gene_id]], data[[sample_id]])

    mut_counts <- map(data_splits, table)

    sample_ids <- names(mut_counts)

    mut_counts <- do.call('rbind', mut_counts)

    mut_counts <- ifelse(mut_counts > 0, 1, mut_counts)

    mut_counts <- as_tibble(mut_counts)

    if(!as_numeric) {

      mut_counts <-
        future_map_dfc(mut_counts,
                       ~car::recode(as.character(.x),
                                    "'0' = 'WT'; '1' = 'mutated'"),
                       .options = furrr_options(seed = TRUE))

      mut_counts <-
        future_map_dfc(mut_counts,
                       factor,
                       c('WT', 'mutated'),
                       .options = furrr_options(seed = TRUE))

    }

    sample_id <- NULL

    mut_counts <- mutate(mut_counts,
                         sample_id = sample_ids)

    relocate(mut_counts, sample_id)

  }

# Frequency of binary data ------

#' Frequency of binary data.
#'
#' @description
#' Computes frequencies and percentages for binary data split by a factor.
#'
#' @details
#' The binary data can be provided in a 0/1 form or as two-level factors.
#'
#' @return a data frame with the variable names (`variable`), numbers of
#' 1-coded events (`n`), total numbers of observations (`n_total`), and
#' percentage of events within all observations (`percent`).
#'
#' @param data a data frame with the splitting factor and binary variables of
#' interest.
#' @param variables a character vector with names of the binary variables.
#' If `NULL`, the argument defaults to all variables in the data frame
#' except of `split_fct`.
#' @param split_fct name of the splitting factor variable. If `NULL`, global
#' frequencies in the data set will be quantified.
#'
#' @export

  count_binary <- function(data,
                           variables = NULL,
                           split_fct = NULL) {

    ## entry control --------

    if(!is.data.frame(data)) {

      stop("'data' has to be a data frame.", call. = FALSE)

    }

    if(is.null(variables)) variables <- names(data)

    if(!is.null(split_fct)) {

      stopifnot(is.character(split_fct))

      if(!split_fct %in% names(data)) {

        stop("'split_fct' missing from the input data frame.", call. = FALSE)

      }

      variables <- variables[variables != split_fct]

    }

    if(any(!variables %in% names(data))) {

      stop('Some of variables are absent from the input data frame.',
           call. = FALSE)

    }

    ## conversion to a numeric and checking for level numbers

    conv_numeric <- function(x) {

      if(is.numeric(x)) return(x)

      if(is.factor(x)) return(as.numeric(x) - 1)

      return(as.numeric(factor(x)) - 1)

    }

    data[variables] <- map_dfc(data[variables], conv_numeric)

    check_binary_df(data[variables])

    ## descriptive stats calculation -------

    variable <- NULL
    n <- NULL
    n_total <- NULL
    percent <- NULL

    if(is.null(split_fct)) {

      stats <- colSums(data[variables])

      stats <- tibble(variable = variables,
                      n = unname(stats),
                      n_total = nrow(data),
                      percent = unname(stats)/nrow(data) * 100)

    } else {

      data <- split(data[variables], f = data[[split_fct]])

      stats <- map(data, colSums)

      totals <- map_dbl(data, nrow)

      stats <-
        map2(stats, totals,
             ~tibble(variable = variables,
                     n = unname(.x),
                     n_total = unname(.y)))

      stats <-
        map2_dfr(stats, names(stats),
                 ~mutate(.x,
                         percent = n/n_total * 100,
                         !!split_fct := .y))

      stats <- relocate(stats, .data[[split_fct]])

    }

    arrange(stats, -percent)

  }

# END ------
