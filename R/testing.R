# Functions for testing gene expression differences between two groups.

# Helpers -----

#' T test.
#'
#' @description Test for the difference with a T test.
#' @param data a data frame.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variable dependent variable name.
#' @param ... additional arguments passed to \code{\link[stats]{t.test}}.
#' @return a data frame with the regulation statistic, confidence intervals,
#' t and p values

  t_tester <- function(data, split_fct, variable, ...) {

    ## the entry control is accomplished by the upstream function

    tst_vec <- split(data[[variable]], data[[split_fct]])

    res <- tryCatch(t.test(x = tst_vec[[1]],
                           y = tst_vec[[2]], ...),
                    error = function(e) NULL)

    if(!is.null(res)) {

      tibble::tibble(test = 'T test',
                     response = variable,
                     grouping = split_fct,
                     stat_name = 't',
                     stat = res$statistic[1],
                     p_value = res$p.value[1],
                     expect_x = res$estimate[1],
                     expect_y = res$estimate[2],
                     estimate = res$estimate[2] - res$estimate[1],
                     lower_ci = -res$conf.int[2],
                     upper_ci = -res$conf.int[1])

    } else {

      NULL

    }

  }

#' Wilcoxon test.
#'
#' @description Test for the difference with a T test.
#' @param data a data frame.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variable dependent variable name.
#' @param conf.int logical, should confidence intervals be calculated?
#' @param ... additional arguments passed to \code{\link[stats]{wilcox.test}}.
#' @return a data frame with the regulation statistic, confidence intervals,
#' t and p values

  wilcox_tester <- function(data,
                            split_fct,
                            variable,
                            conf.int = FALSE, ...) {

    ## the entry control is accomplished by the upstream function

    tst_vec <- split(data[[variable]], data[[split_fct]])

    res <- tryCatch(wilcox.test(x = tst_vec[[1]],
                                y = tst_vec[[2]],
                                conf.int = conf.int, ...),
                    error = function(e) NULL)

    if(!is.null(res)) {

      tibble::tibble(test = 'Wilcoxon test',
                     response = variable,
                     grouping = split_fct,
                     stat_name = 'w',
                     stat = res$statistic[1],
                     p_value = res$p.value[1],
                     expect_x = if(conf.int) res$estimate[1] else NA,
                     expect_y = if(conf.int) res$estimate[2] else NA,
                     estimate = if(conf.int) res$estimate[1] else NA,
                     lower_ci = if(conf.int) res$conf.int[1] else NA,
                     upper_ci = if(conf.int) res$conf.int[2] else NA)

    } else {

      NULL

    }

  }

#' One-way ANOVA.
#'
#' @description A simplified version of one-way ANOVA suitable for whole-genome
#' testing.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variable dependent variable name.
#' @param ... additional arguments passed to \code{\link[stats]{aov}}.
#' @return a data frame with the F statistic, degrees of freedom and p value.

  anova_tester <- function(data,
                           split_fct,
                           variable, ...) {

    ## the entry control is done by the upstream function

    form <- as.formula(paste0('`',
                              variable,
                              '` ~ `',
                              split_fct,
                              '`'))

    aov_model <- tryCatch(aov(formula = form, data = data, ...),
                          error = function(e) NULL)

    if(is.null(aov_model)) return(NULL)

    aov_summary <- as.data.frame(summary(aov_model)[[1]])

    aov_summary <- tibble::tibble(test = 'one-way ANOVA',
                                  response = variable,
                                  grouping = split_fct,
                                  stat_name = 'F',
                                  stat = aov_summary[['F value']][1],
                                  df1 = aov_summary$Df[1],
                                  df2 = aov_summary$Df[2],
                                  p_value = aov_summary[['Pr(>F)']][1])

    lm_coefs <- summary.lm(aov_model)$coefficients

    lm_ci <- confint(aov_model)

    lm_summary <- tibble::tibble(test = 'lm',
                                 response = variable,
                                 grouping = split_fct,
                                 level = rownames(lm_coefs),
                                 stat_name = 't',
                                 stat = lm_coefs[, 't value'],
                                 estimate = lm_coefs[, 'Estimate'],
                                 se = lm_coefs[, 'Std. Error'],
                                 lower_ci = lm_ci[, 1],
                                 upper_ci = lm_ci[, 2],
                                 p_value = lm_coefs[, 'Pr(>|t|)'])

    list(anova = aov_summary,
         lm = lm_summary)

  }


# Tester: two groups ------

#' Test gene expression differences between two groups.
#'
#' @description Compares gene expression differences between two groups with
#' the T test.
#' @param data a data frame with the expression values.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variables a vector with dependent variable names.
#' @param type test type, currently 't' for T test (default) or 'wilcox' for
#' the Wilcoxon/Mann-Whitney test.
#' @param adj_method p value adjustment method,
#' see: \code{\link[stats]{p.adjust}}.
#' @param .parallel should the testing be run in parallel?
#' @param ... additional arguments passed to \code{\link[stats]{t.test}} or
#' \code{\link[stats]{wilcox.test}}.
#' @return a data frame with the regulation statistic, confidence intervals,
#' t and p values (raw and corrected).
#' @export

  test_two_groups <- function(data,
                              split_fct,
                              variables,
                              type = c('t', 'wilcox'),
                              adj_method = 'BH',
                              .parallel = FALSE, ...) {

    ## entry control

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(split_fct, variables) %in% names(data))) {

      stop('Some variables are missing from the data.', call. = FALSE)

    }

    if(!is.factor(data[[split_fct]])) {

      stop('split_fct needs to be a factor.', call. = FALSE)

    }

    if(length(levels(data[[split_fct]])) != 2) {

      stop('split_factor variable needs to have exactly two levels.',
           call. = FALSE)

    }

    type <- match.arg(type, c('t', 'wilcox'))

    test_fun <- switch(type,
                       t = microViz:::t_tester,
                       wilcox = microViz:::wilcox_tester)

    stopifnot(is.logical(.parallel))

    ## benchmarking

    start_time <- Sys.time()
    message(paste('Testing for n =', length(variables), 'variables'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## data split and testing

    data_lst <- purrr::map(variables,
                           ~data[c(split_fct, .x)])

    if(.parallel) {

      future::plan('multisession')

      tst_res <- furrr::future_map2(variables,
                                    data_lst,
                                    function(x, y) test_fun(data = y,
                                                            variable = x,
                                                            split_fct = split_fct, ...),
                                    .options = furrr::furrr_options(seed = TRUE,
                                                                    packages = c('tibble')))

      future::plan('sequential')

    } else {

      tst_res <- purrr::map2(variables,
                             data_lst,
                             function(x, y) test_fun(data = y,
                                                     variable = x,
                                                     split_fct = split_fct, ...))

    }

    tst_res <- purrr::compact(tst_res)

    tst_res <- do.call('rbind', tst_res)

    dplyr::mutate(tst_res,
                  p_adjusted = p.adjust(p_value, method = adj_method),
                  significant = ifelse(p_adjusted < 0.05, 'yes', 'no'))


  }

# Tester: one-way ANOVA with LM -----

#' Test gene expression differences between multiple groups with one-way ANOVA.
#'
#' @description Compares expression differences between multiple groups with
#' one-way ANOVA and linear modeling.
#' @param data a data frame with the expression values.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variables a vector with dependent variable names.
#' @param adj_method p value adjustment method,
#' see: \code{\link[stats]{p.adjust}}.
#' @param .parallel should the testing be run in parallel?
#' @param ... additional arguments passed to \code{\link[stats]{aov}}.
#' @return a list of data frames containing the ANOVA and linear model
#' estimates, respectively.
#' @export

  test_anova <- function(data,
                         split_fct,
                         variables,
                         adj_method = 'BH',
                         .parallel = FALSE, ...) {

    ## entry control

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(split_fct, variables) %in% names(data))) {

      stop('Some variables are missing from the data.', call. = FALSE)

    }

    if(!is.factor(data[[split_fct]])) {

      stop('split_fct needs to be a factor.', call. = FALSE)

    }

    stopifnot(is.logical(.parallel))

    ## benchmarking

    start_time <- Sys.time()
    message(paste('Testing for n =', length(variables), 'variables'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## data split and testing

    data_lst <- purrr::map(variables,
                           ~data[c(split_fct, .x)])

    if(.parallel) {

      future::plan('multisession')

      tst_res <- furrr::future_map2(variables,
                                    data_lst,
                                    function(x, y) microViz:::anova_tester(data = y,
                                                                           split_fct = split_fct,
                                                                           variable = x, ...),
                                    .options = furrr::furrr_options(seed = TRUE,
                                                                    packages = c('tibble')))

      future::plan('sequential')

    } else {

      tst_res <- purrr::map2(variables,
                             data_lst,
                             function(x, y) microViz:::anova_tester(data = y,
                                                                    split_fct = split_fct,
                                                                    variable = x, ...))

    }

    tst_res <- purrr::transpose(purrr::compact(tst_res))

    tst_res <- purrr::map(tst_res, ~do.call('rbind', .x))

    purrr::map(tst_res,
               dplyr::mutate,
               p_adjusted = p.adjust(p_value, method = adj_method),
               significant = ifelse(p_adjusted < 0.05, 'yes', 'no'))

  }



