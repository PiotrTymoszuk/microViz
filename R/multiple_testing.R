# Re-adjustment of testing results and identification of significant effects

#' @include imports.R

  NULL

# Multiple testing adjustment ---------

#' Adjust for multiple testing.
#'
#' @description
#' Adjusts for multiple comparison. The input is a result of a statistic test
#' in a data frame format.
#'
#' @details
#' For available methods, please consult \code{\link[stats]{p.adjust}}.
#'
#' @param data a data frame with the testing results.
#' @param p_variable name of the variable that stores raw p values to be
#' adjusted.
#' @param method adjustment method.
#' @param simplify_p logical, should adjusted p values below 0.001 be presented
#' in a p < 0.001 form?
#'
#' @return
#' The genuine data frame appended with two columns:
#' * `p_adjusted` with adjusted p values
#' * `significance` with text indicating significance and p value, e.g.
#' 'ns(p = 0.12).
#'
#' @export

  re_adjust <- function(data,
                        p_variable = 'p_value',
                        method = 'BH',
                        simplify_p = TRUE) {

    ## input control -------

    if(!is.data.frame(data)) {

      stop("'x' has to be a dta frame.", call. = FALSE)

    }

    if(!p_variable %in% names(data)) {

      stop("'p_variable' is missing from the data.", call. = FALSE)

    }

    stopifnot(is.logical(simplify_p))

    ## adjustment ------

    data[['p_adjusted']] <- p.adjust(data[[p_variable]], method = method)

    significance <- NULL
    p_adjusted <- NULL

    if(simplify_p) {

      data <-
        mutate(data,
               significance = ifelse(p_adjusted < 0.001,
                                     'p < 0.001',
                                     ifelse(p_adjusted >= 0.05,
                                            paste0('ns (p = ', signif(p_adjusted, 2), ')'),
                                            paste('p =', signif(p_adjusted, 2)))))

    } else {

      data <-
        mutate(data,
               significance = ifelse(p_adjusted >= 0.05,
                                     paste0('ns (p = ', signif(p_adjusted, 2), ')'),
                                     paste('p =', signif(p_adjusted, 2))))

    }

    data

  }

# Identification of significant effects ----------

#' Identify significant effects.
#'
#' @description
#' Identifies significant effects with a given significance and regulation
#' or effect size threshold.
#'
#' @return
#' If `return_data = FALSE`: a character vector or a list of character vectors
#' with names of significant variables derived from `label_variable`.
#' Otherwise, the function returns the input data frame with the following
#' addtional columns:
#' * `.significant`: indicates if the effect reaches statistical significance
#' * `.regulation`: indicates the sign of significant regulation
#'
#' @param data a data frame with the testing results.
#' @param label_variable name of the variable storing names of the tested
#' variables.
#' @param p_variable name of the variable storing p values.
#' @param regulation_variable name of the variable storing regulation or
#' effect size figures. If `NULL`, significance is determined solely by the
#' p value threshold.
#' @param signif_level p value threshold.
#' @param regulation_level regulation or effect size threshold.
#' @param return_data if `TRUE`, the input data frame is returned with additional
#' columns indicating significance and regulation sign.
#'
#' @export

  identify_significant <- function(data,
                                   label_variable,
                                   p_variable,
                                   regulation_variable = NULL,
                                   signif_level = 0.05,
                                   regulation_level = 0,
                                   return_data = FALSE) {

    ## input control --------

    if(!is.data.frame(data)) {

      stop('A data frame is required.', call. = FALSE)

    }

    vars_to_check <- c(label_variable, p_variable)

    if(!is.null(regulation_variable)) {

      vars_to_check <- c(vars_to_check, regulation_variable)

    }

    err_txt <-
      paste0('At least one of ',
             paste(vars_to_check, collapse = ' or '),
             ' is missing from the data.')

    if(!any(vars_to_check %in% names(data))) stop(err_txt, call. = FALSE)

    class_check <- map_lgl(data[vars_to_check[-1]], is.numeric)

    err_txt <-
      paste0('All of ',
             paste(vars_to_check[-1], collapse = ' and '),
             ' variables have to be numeric.')

    if(any(!class_check)) stop(err_txt, call. = FALSE)

    stopifnot(is.numeric(signif_level))
    stopifnot(is.numeric(regulation_level))

    ## identification of significant effects -------

    .significant <- NULL

    data <- mutate(data,
                   .significant = ifelse(.data[[p_variable]] < signif_level,
                                        'significant', 'ns'),
                   .significant = factor(.significant, c('significant', 'ns')))

    if(is.null(regulation_variable)) {

      if(return_data) return(data)

      data <- filter(data, .significant == 'significant')

      return(data[[label_variable]])

    }

    .regulation <- NULL

    data <-
      mutate(data,
             .regulation = ifelse(.significant == 'ns',
                                  'ns',
                                  ifelse(.data[[regulation_variable]] == 0,
                                         'ns',
                                         ifelse(.data[[regulation_variable]] >= regulation_level,
                                                'upregulated',
                                                ifelse(.data[[regulation_variable]] <= -regulation_level,
                                                       'downregulated', 'ns')))),
             .regulation = factor(.regulation,
                                  c('upregulated', 'downregulated', 'ns')))

    if(return_data) return(data)

    data <- filter(data,
                   .regulation %in% c('upregulated', 'downregulated'))

    data <- mutate(data, .regulation = droplevels(.regulation))

    signif_vars <- split(data[[label_variable]],
                         data[['.regulation']])

    signif_vars <-
      map(signif_vars, function(x) if(length(x) == 0) NULL else x)

    signif_vars <- compact(signif_vars)

    if(length(signif_vars) == 1) {

      return(signif_vars[[1]])

    } else {

      return(signif_vars)

    }

  }

# Shared elements ----------

#' Tally and identify shared elements in sets.
#'
#' @description
#' The functions count and select elements shared by at least n elements of
#' a list.
#' As such they may be useful at identifying e.g. significantly regulated genes
#' or pathways in several data sets.
#'
#' @return `shared_features()`: a vector with the features shared by the given
#' number of elements.
#' `count_features()`: a data frame with the following columns:
#'
#' * `element`: element names
#'
#' * `n`: number of sets, where the element was present
#'
#' * `n_total`: total number of sets
#'
#' * `percent`: percentage of sets, where the given element was present
#'
#' @param x a list of character or numeric vectors.
#' @param m number of sets which share the common features.
#'
#' @export

  count_features <- function(x) {

    stopifnot(is.list(x))

    all_lst <- reduce(compact(x), c)

    freq_tbl <- as.data.frame(table(all_lst))

    freq_tbl <- set_names(freq_tbl, c('element', 'n'))

    n <- NULL
    n_total <- NULL
    percent <- NULL

    freq_tbl <- mutate(freq_tbl,
                       n_total = length(x),
                       percent = n/n_total * 100)

    freq_tbl <- arrange(freq_tbl, -percent)

    as_tibble(freq_tbl)

  }

#' @rdname count_features
#' @export

  shared_features <- function(x, m = 2) {

    ft_tbl <- count_features(x)

    n <- NULL

    filter(ft_tbl, n >= m)$element

  }

# END ------
