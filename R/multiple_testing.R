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

# Identification and counting of significant effects ----------

#' Identify and count significant effects.
#'
#' @description
#' Functions for identification and tallying of significant effects with a given
#' significance and regulation or effect size threshold.
#'
#' @return
#' `identify_significant()`, if `return_data = FALSE`:
#' a character vector or a list of character vectors
#' with names of significant variables derived from `label_variable`.
#' Otherwise, the function returns the input data frame with the following
#' additional columns:
#' * `.significant`: indicates if the effect reaches statistical significance
#' * `.regulation`: indicates the sign of significant regulation
#' `count_significant()`: a data frame with numbers and percentages of up- and
#' downregulated features.
#'
#' @param data a data frame with the testing results.
#' @param label_variable name of the variable storing names of the tested
#' variables.
#' @param p_variable name of the variable storing p values.
#' @param regulation_variable name of the variable storing regulation or
#' effect size figures. If `NULL`, significance is determined solely by the
#' p value threshold.
#' @param status_variable name of the variable storing information on
#' significance/regulation status. If not `NULL`, `p_variable`
#' and `regulation_variable` are ignored.
#' @param signif_level p value threshold.
#' @param regulation_level regulation or effect size threshold.
#' @param return_data if `TRUE`, the input data frame is returned with additional
#' columns indicating significance and regulation sign.
#'
#' @export

  identify_significant <- function(data,
                                   p_variable = "p_value",
                                   regulation_variable = NULL,
                                   signif_level = 0.05,
                                   regulation_level = 0,
                                   return_data = FALSE,
                                   label_variable = NULL) {

    ## input control --------

    if(!is.data.frame(data)) {

      stop('A data frame is required.', call. = FALSE)

    }

    stopifnot(is.character(p_variable))
    if(!p_variable %in% names(data)) stop("`p_variable` is missing from `data`.", call. = FALSE)
    if(!is.numeric(data[[p_variable]])) stop("`p_variable must be numeric.`", call. = FALSE)

    if(!is.null(label_variable)) {

      stopifnot(is.character(label_variable))

      if(!label_variable %in% names(data)) {

        stop("`label_variable` is missing from `data`.", call. = FALSE)

      }

    }

    if(!is.null(regulation_variable)) {

      stopifnot(is.character(regulation_variable))

      if(!regulation_variable %in% names(data)) {

        stop("`regulation_variable` is missing from `data`.",
             call. = FALSE)

      }

      if(!is.numeric(data[[regulation_variable]])) {

        stop("`regulation_variable must be numeric.`", call. = FALSE)

      }

    }

    if(is.null(label_variable) & !return_data) {

      stop("To return variables vectors (`return_data = TRUE`), please specify `label_variable`.",
           call. = FALSE)

    }

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

#' @rdname identify_significant
#' @export

  count_significant <- function(data,
                                p_variable = "p_value",
                                regulation_variable = NULL,
                                status_variable = NULL,
                                signif_level = 0.05,
                                regulation_level = 0) {

    ## input control --------

    if(!is.null(status_variable)) {

      if(!is.data.frame(data)) stop("A data frame is required.", call. = FALSE)

      if(!status_variable %in% names(data)) {

        stop("`status_variable` is missing from `data`.", call. = FALSE)

      }

    } else {

      data <- identify_significant(data,
                                   p_variable = p_variable,
                                   regulation_variable = regulation_variable,
                                   signif_level = signif_level,
                                   regulation_level = regulation_level,
                                   return_data = TRUE)

      if(is.null(regulation_variable)) {

        status_variable <- ".significant"

      } else {

        tatus_variable <- ".regulation"

      }

    }

    count_data <- count(data, .data[[status_variable]])

    count_data[["n_total"]] <- sum(count_data[["n"]])
    count_data[["percent"]] <- count_data[["n"]]/count_data[["n_total"]] * 100

    count_data

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
