# Functions for testing gene expression differences between two groups.

#' @include imports.R

  NULL

# Helpers -----

#' T and Wilcoxon test.
#'
#' @description Test for the difference with a T test.
#'
#' @param data a data frame.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variable dependent variable name.
#' @param conf.int logical, should confidence intervals be calculated?
#' @param ... additional arguments passed to \code{\link[stats]{t.test}} or
#' \code{\link[stats]{wilcox.test}}.
#'
#' @return a data frame with the regulation statistic, confidence intervals,
#' t and p values as well as effect size (Cohen's d for T tests and r for
#' Wilcoxon tests).

  t_tester <- function(data, split_fct, variable, ...) {

    ## the entry control is accomplished by the upstream function

    tst_vec <- split(data[[variable]], data[[split_fct]])

    res <- tryCatch(stats::t.test(x = tst_vec[[1]],
                                  y = tst_vec[[2]], ...),
                    error = function(e) NULL)

    if(!is.null(res)) {

      ## difference between the means

      if(is.na(res$estimate[2])) {

        estimate <- res$estimate[1]

      } else {

        estimate <- res$estimate[2] - res$estimate[1]

      }

      ## pooled SD

      method <- res$method

      if(stri_detect(method, fixed = 'Paired')) {

        if(length(tst_vec[[1]]) != length(tst_vec[[2]])) return(NULL)

        sd <- sqrt(var(tst_vec[[1]] - tst_vec[[2]]))

      } else if(stri_detect(method, fixed = 'Welch')) {

        sd <- sqrt((var(tst_vec[[1]]) + var(tst_vec[[2]]))/2)

      } else {

        sd_sq1 <- (tst_vec[[1]] - mean(tst_vec[[1]], na.rm = TRUE))^2
        sd_sq2 <- (tst_vec[[2]] - mean(tst_vec[[2]], na.rm = TRUE))^2

        n_diff <- length(tst_vec[[1]] + length(tst_vec[[2]]) - 2)

        if(n_diff == 0) n_diff <- 1

        sd <- sqrt((sum(sd_sq1) + sum(sd_sq2))/n_diff)

      }

      ## output

      tibble(test = 'T test',
             method = method,
             response = variable,
             grouping = split_fct,
             stat_name = 't',
             stat = res$statistic[1],
             p_value = res$p.value[1],
             expect_x = res$estimate[1],
             expect_y = res$estimate[2],
             estimate = estimate,
             error_name = 'SEM',
             error = res$stderr,
             lower_ci = -res$conf.int[2],
             upper_ci = -res$conf.int[1],
             sd_name = 'pooled SD',
             sd = sd,
             effect_size_name = 'd',
             effect_size = estimate/sd)

    } else {

      NULL

    }

  }

#' @rdname t_tester

  wilcox_tester <- function(data,
                            split_fct,
                            variable,
                            conf.int = FALSE, ...) {

    ## the entry control is accomplished by the upstream function

    tst_vec <- split(data[[variable]], data[[split_fct]])

    res <- tryCatch(stats::wilcox.test(x = tst_vec[[1]],
                                       y = tst_vec[[2]],
                                       conf.int = conf.int, ...),
                    error = function(e) NULL)

    if(!is.null(res)) {

      ## checking if the test was paired

      method <- res$method

      extra_args <- rlang::list2(...)

      if(!is.null(extra_args$paired)) {

        paired <- extra_args$paired

      } else {

        paired <- FALSE

      }

      ## effect size calculation

      if(paired){

        method <- paste('Paired', method)

        wilcox_r <- rcompanion::wilcoxonPairedR(x = data[[variable]],
                                                g = data[[split_fct]])

      } else {

        wilcox_r <- rcompanion::wilcoxonR(x = data[[variable]],
                                          g = data[[split_fct]])

      }

      tibble(test = 'Wilcoxon test',
             method = method,
             response = variable,
             grouping = split_fct,
             stat_name = 'w',
             stat = res$statistic[1],
             p_value = res$p.value[1],
             expect_x = if(conf.int) res$estimate[1] else NA,
             expect_y = if(conf.int) res$estimate[2] else NA,
             estimate = if(conf.int) res$estimate[1] else NA,
             error_name = 'none',
             error = NA,
             lower_ci = if(conf.int) res$conf.int[1] else NA,
             upper_ci = if(conf.int) res$conf.int[2] else NA,
             sd_name = '1 - R-square',
             sd = 1 - wilcox_r^2,
             effect_size_name = 'r',
             ## reversing effect size, the first level assumed control
             effect_size = -unname(wilcox_r))

    } else {

      NULL

    }

  }

#' One-way ANOVA and negative binomial testing.
#'
#' @description
#' Simplified versions of one-way ANOVA and negative binomial (NB)
#' modeling suitable for whole-genome testing.
#'
#' @param data a data frame.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param confounder name of the confounding variable, NULL by default.
#' @param variable dependent variable name.
#' @param dev_test name of the test used for analysis of deviance. One of:
#' 'Rao', 'LRT', 'Chisq' (default), 'F', 'Cp'.
#' @param ... additional arguments passed to \code{\link[stats]{aov}} or
#' \code{\link[MASS]{glm.nb}}.
#'
#' @return a list of data frames with ANOVA results or results of analysis of
#' deviance (F/test statistic, degrees of freedom and p value) and
#' linear model coefficient estimates (estimates, test statistic, errors and
#' confidence intervals, pooled standard deviation and Cohen's d).

  anova_tester <- function(data,
                           split_fct,
                           variable,
                           confounder = NULL, ...) {

    ## the entry control is done by the upstream function

    ## one-way ANOVA ------

    if(is.null(confounder)) {

      form <- as.formula(paste0('`',
                                variable,
                                '` ~ `',
                                split_fct,
                                '`'))

    } else {

      form <- as.formula(paste0('`',
                                variable,
                                '` ~ `',
                                confounder,
                                '` + `',
                                split_fct,
                                '`'))

    }

    aov_model <- tryCatch(stats::aov(formula = form, data = data, ...),
                          error = function(e) NULL)

    if(is.null(aov_model)) return(NULL)

    aov_summary <- as.data.frame(summary(aov_model)[[1]])

    aov_summary$eta_sq <-
      aov_summary[['Sum Sq']]/sum(aov_summary[['Sum Sq']])

    split_rex <- paste0('(', split_fct, ')|(Residuals)')

    trunk_rows <- stri_detect(rownames(aov_summary), regex = split_rex)

    aov_summary_trunk <- aov_summary[trunk_rows, ]

    if(is.null(confounder)) {

      aov_method <- 'One-way ANOVA'
      lm_method <- 'Linear modeling'

    } else {

      aov_method <- 'One way ANOVA with a confounder'
      lm_method <- 'Linear modeling with a confounder'

    }

    aov_summary <-
      tibble(test = 'one-way ANOVA',
             method = aov_method,
             response = variable,
             grouping = split_fct,
             confounder = if(is.null(confounder)) NA else confounder,
             n = nobs(aov_model),
             stat_name = 'F',
             stat = aov_summary_trunk[['F value']][1],
             df1 = aov_summary_trunk$Df[1],
             df2 = aov_summary_trunk$Df[1],
             p_value = aov_summary_trunk[['Pr(>F)']][1],
             effect_size_name = '\u03B7\u00B2',
             effect_size = aov_summary_trunk[['eta_sq']][1])

    ## linear modeling ------

    lm_coefs <- stats::summary.lm(aov_model)$coefficients

    coef_rex <- paste0('(', split_fct, ')|(Intercept)')

    trunk_rows <- stri_detect(rownames(lm_coefs), regex = coef_rex)

    lm_coefs_trunk <- lm_coefs[trunk_rows, ]

    lm_ci <- confint(aov_model)

    trunk_rows <- stri_detect(rownames(lm_ci), regex = coef_rex)

    lm_ci_trunk <- lm_ci[trunk_rows, ]

    ## effect size computation for betas ------

    ## pooled SD: computed as a square-root of average variance
    ## (baseline and the current level)

    split_levs <- levels(data[[split_fct]])

    split_levs[1] <- '(Intercept)'

    resp_splits <- split(data[[variable]], data[[split_fct]])

    split_variances <- map(resp_splits, var)

    sd <- map_dbl(split_variances,
                  ~sqrt((.x + split_variances[[1]])/2))

    level <- NULL
    n <- NULL

    sd_table <- tibble(level = split_levs,
                       n = map_dbl(resp_splits, length),
                       sd = sd)

    ## lm output -------

    estimate <- NULL

    lm_summary <-
      tibble(test = 'lm',
             method = lm_method,
             response = variable,
             grouping = split_fct,
             confounder = if(is.null(confounder)) NA else confounder,
             level = stri_replace(rownames(lm_coefs_trunk),
                                  fixed = split_fct,
                                  replacement = ''),
             stat_name = 't',
             stat = lm_coefs_trunk[, 't value'],
             estimate = lm_coefs_trunk[, 'Estimate'],
             error_name = 'SEM',
             error = lm_coefs_trunk[, 'Std. Error'],
             lower_ci = lm_ci_trunk[, 1],
             upper_ci = lm_ci_trunk[, 2],
             p_value = lm_coefs_trunk[, 'Pr(>|t|)'],
             sd_name = 'pooled SD')

    lm_summary <- left_join(lm_summary, sd_table, by = 'level')

    effect_size_name <- NULL
    effect_size <- NULL

    lm_summary <- mutate(lm_summary,
                         effect_size_name = 'd',
                         effect_size = estimate/sd)

    list(anova = aov_summary,
         lm = lm_summary[c('test', 'method', 'response', 'grouping',
                           'confounder', 'level', 'n', 'stat_name',
                           'stat', 'estimate', 'error_name', 'error',
                           'lower_ci', 'upper_ci', 'p_value',
                           'sd_name', 'sd', 'effect_size_name', 'effect_size')])

  }

#' @rdname anova_tester

  nb_tester <- function(data,
                        split_fct,
                        variable,
                        confounder = NULL,
                        dev_test = 'Chisq', ...) {

    ## the entry control is done by the upstream function

    ## model formulas -------

    if(is.null(confounder)) {

      form <- as.formula(paste0('`',
                                variable,
                                '` ~ `',
                                split_fct,
                                '`'))

    } else {

      form <- as.formula(paste0('`',
                                variable,
                                '` ~ `',
                                confounder,
                                '` + `',
                                split_fct,
                                '`'))

    }

    ## construction of an NB model, extraction of coefs and CI -------

    nb_model <- tryCatch(MASS::glm.nb(formula = form, data = data, ...),
                         error = function(e) NULL)

    if(is.null(nb_model)) return(NULL)

    coef_rex <- paste0('(', split_fct, ')|(Intercept)')

    lm_coefs <- tryCatch(summary(nb_model)$coefficients,
                         error = function(e) NULL)

    if(is.null(lm_coefs)) return(NULL)

    lm_ci <- tryCatch(suppressMessages(confint(nb_model)),
                      error = function(e) NULL)

    if(is.null(lm_ci)) return(NULL)

    trunk_rows <- stri_detect(rownames(lm_coefs), regex =  coef_rex)

    lm_coefs_trunk <- lm_coefs[trunk_rows, , drop = FALSE]

    if(nrow(lm_coefs_trunk) == 1) {

      warning("An intercept only-model, redefine the split factor?",
              call. = FALSE)

      return(NULL)

    }

    trunk_rows <- stri_detect(rownames(lm_ci), regex =  coef_rex)

    lm_ci_trunk <- lm_ci[trunk_rows, , drop = FALSE]

    if(!is.null(confounder)) {

      lm_method <- 'NB modeling with a confounder'
      aod_method <- 'One-way AOD with a confounder'

    } else {

      lm_method <- 'NB modeling'
      aod_method <- 'One-way AOD'

    }

    ## pooled SD: equal to the beta error, the model works with the Z statistic

    lm_summary <-
      tibble(test = 'nb lm',
             method = lm_method,
             response = variable,
             grouping = split_fct,
             confounder = if(is.null(confounder)) NA else confounder,
             level = stri_replace(rownames(lm_coefs_trunk),
                                           fixed = split_fct,
                                           replacement = ''),
             stat_name = 'z',
             stat = lm_coefs_trunk[, 'z value'],
             estimate = lm_coefs_trunk[, 'Estimate'],
             error_name = 'SEM',
             error = lm_coefs_trunk[, 'Std. Error'],
             lower_ci = lm_ci_trunk[, 1],
             upper_ci = lm_ci_trunk[, 2],
             p_value = lm_coefs_trunk[, 'Pr(>|z|)'],
             sd_name = 'pooled SD',
             sd = lm_coefs_trunk[, 'Std. Error'],
             effect_size_name = 'd',
             effect_size = lm_coefs_trunk[, 'Estimate']/lm_coefs_trunk[, 'Std. Error'])

    ## extracting the ANOVA table ------

    ## warnings concerning re-estimation of theta are suppressed

    aov_summary <- suppressWarnings(anova(nb_model, test = dev_test))

    aov_summary <- as.data.frame(aov_summary)

    aov_summary$eta_sq <-
      aov_summary[['Deviance']]/aov_summary[['Resid. Dev']][1]

    split_rex <- paste0('(', split_fct, ')|(NULL)')

    trunk_rows <-
      stri_detect(rownames(aov_summary), regex = split_rex)

    aov_summary <- aov_summary[trunk_rows, ]

    aov_summary <-
      tibble(test = 'analysis of deviance',
             response = variable,
             grouping = split_fct,
             counfounder = if(is.null(confounder)) NA else confounder,
             deviance = aov_summary[['Deviance']][2],
             residual_deviance = aov_summary[['Resid. Dev']][2],
             df = aov_summary$Df[2],
             p_value = aov_summary[[ncol(aov_summary) - 1]][2],
             effect_size_name = '\u03B7\u00B2',
             effect_size = aov_summary[['eta_sq']][2])

    list(aod = aov_summary,
         lm = lm_summary)

  }


# Testers for two groups ------

#' Test gene expression differences between two groups.
#'
#' @description
#' Compares gene expression differences between two groups with
#' T tests or Wilcoxon tests.
#'
#' @details
#' The function returns Cohen's d as an effect size metric for T tests and
#' r as an effect size statistic for Wilcoxon tests. The later is computed
#' with \code{\link[rcompanion]{wilcoxonR}}.
#'
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
#'
#' @return a data frame with the regulation statistic, confidence intervals,
#' t and p values (raw and corrected) as well as effect size measures.
#'
#' @export

  test_two_groups <- function(data,
                              split_fct,
                              variables,
                              type = c('t', 'wilcox'),
                              adj_method = 'BH',
                              .parallel = FALSE, ...) {

    ## entry control -------

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
                       t = t_tester,
                       wilcox = wilcox_tester)

    stopifnot(is.logical(.parallel))

    p_value <- NULL
    p_adjusted <- NULL

    ## benchmarking ------

    start_time <- Sys.time()
    message(paste('Testing for n =', length(variables), 'variables'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## data split and testing --------

    data_lst <- map(variables,
                    ~data[c(split_fct, .x)])

    if(.parallel) {

      plan('multisession')

      tst_res <- future_map2(variables,
                             data_lst,
                             function(x, y) test_fun(data = y,
                                                     variable = x,
                                                     split_fct = split_fct, ...),
                             .options = furrr_options(seed = TRUE,
                                                      packages = c('tibble')))

      plan('sequential')

    } else {

      tst_res <- map2(variables,
                      data_lst,
                      function(x, y) test_fun(data = y,
                                              variable = x,
                                              split_fct = split_fct, ...))

    }

    tst_res <- compact(tst_res)

    tst_res <- do.call('rbind', tst_res)

    mutate(tst_res,
           p_adjusted = p.adjust(p_value, method = adj_method),
           significant = ifelse(p_adjusted < 0.05, 'yes', 'no'))


  }

# One-way ANOVA with LM -----

#' Test gene expression differences between multiple groups with one-way ANOVA.
#'
#' @description
#' Compares expression differences between multiple groups with
#' one-way ANOVA and linear modeling.
#'
#' @details If a confounding variable is specified it will be included in the
#' model as a leading fixed term. Statistics for the confounding term are
#' not shown in the output. Note: effect size statistic, Cohen's d, returned for
#' the betas of linear modeling is computed as the ratio of beta to square root
#' of variance means (baseline and the i-th variable level of interest):
#'
#' \deqn{d_i = \frac{\beta_{i}}{\sqrt{0.5 \times (Var_{i} + Var_{baseline})}}}
#'
#' As such, the Cohen's d values are only very coarse estimates of effect size,
#' when a confounder is included in the model. The effect size of ANOVA is
#' eta-squared.
#'
#' @param data a data frame with the expression values.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variables a vector with dependent variable names.
#' @param confounder name of the confounding variable. Defaults to NULL, which
#' means that no confounder is included in the model.
#' @param adj_method p value adjustment method,
#' see: \code{\link[stats]{p.adjust}}.
#' @param .parallel should the testing be run in parallel?
#' @param ... additional arguments passed to \code{\link[stats]{aov}}.
#'
#' @return a list of data frames containing the ANOVA and linear model
#' estimates.
#'
#' @export

  test_anova <- function(data,
                         split_fct,
                         variables,
                         confounder = NULL,
                         adj_method = 'BH',
                         .parallel = FALSE, ...) {

    ## entry control ------

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(split_fct, variables) %in% names(data))) {

      stop('Some variables are missing from the data.', call. = FALSE)

    }

    if(!is.null(confounder)) {

      if(!confounder %in% names(data)) {

        stop('The confounder variable absent from data.', call. = FALSE)

      }

    }

    if(!is.factor(data[[split_fct]])) {

      stop('split_fct needs to be a factor.', call. = FALSE)

    }

    stopifnot(is.logical(.parallel))

    p_value <- NULL
    p_adjusted <- NULL

    ## benchmarking -------

    start_time <- Sys.time()
    message(paste('Testing for n =', length(variables), 'variables'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## data split and testing --------

    if(is.null(confounder)) {

      data_lst <- map(variables, ~data[c(split_fct, .x)])

    } else {

      data_lst <- map(variables, ~data[c(split_fct, confounder, .x)])

    }

    if(.parallel) {

      plan('multisession')

      tst_res <-
        future_map2(variables,
                    data_lst,
                    function(x, y) anova_tester(data = y,
                                                split_fct = split_fct,
                                                confounder = confounder,
                                                variable = x, ...),
                    .options = furrr_options(seed = TRUE,
                                             packages = c('tibble')))

      plan('sequential')

    } else {

      tst_res <-
        map2(variables,
             data_lst,
             function(x, y) anova_tester(data = y,
                                         split_fct = split_fct,
                                         confounder = confounder,
                                         variable = x, ...))

    }

    tst_res <- transpose(compact(tst_res))

    tst_res <- map(tst_res, ~do.call('rbind', .x))

    map(tst_res,
        mutate,
        p_adjusted = p.adjust(p_value, method = adj_method),
        significant = ifelse(p_adjusted < 0.05, 'yes', 'no'))

  }

# Negative binomial modeling ------

#' Test gene expression differences with negative binomial modeling.
#'
#' @description
#' Compares expression differences between multiple groups with
#' negative binomial modeling. The method is particularly useful in finding
# differentially expressed genes with untransformed RNAseq counts.
#'
#' @details
#' If a confounding variable is specified it will be included in the
#' model as a leading fixed term. Statistics for the confounding term are
#' not shown in the output. Effect size for the analysis of deviance (AOD),
#' eta-squared statistic, is computed as the ratio of deviance to total deviance.
#' Effect sizes for NB model coefficients, Cohen's d metrics, are calculated as
#' ratios of betas to their errors:
#'
#' \deqn{d_{i} = \frac{\beta{1}}{SEM_{i}}}
#'
#' As such, these effects size statistics are likely too optimistic, and in
#' case of models with a confounder, potentially biased. Please treat them with
#' caution as raw estimates of the effect size.
#'
#' @param data a data frame with the expression values.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variables a vector with dependent variable names.
#' @param confounder name of the confounding variable. Defaults to NULL, which
#' means that no confounder is include in the model.
#' @param dev_test name of the test used for analysis of deviance. One of:
#' 'Rao', 'LRT', 'Chisq' (default), 'F', 'Cp'.
#' @param adj_method p value adjustment method,
#' see: \code{\link[stats]{p.adjust}}.
#' @param .parallel should the testing be run in parallel?
#' @param ... additional arguments passed to \code{\link[MASS]{glm.nb}}.
#'
#' @return a list of data frames containing the analysis of deviance (AOD)
#' and linear model estimates.
#'
#' @export

  test_nb <- function(data,
                      split_fct,
                      variables,
                      confounder = NULL,
                      adj_method = 'BH',
                      dev_test = 'Chisq',
                      .parallel = FALSE, ...) {

    ## entry control --------

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(split_fct, variables) %in% names(data))) {

      stop('Some variables are missing from the data.', call. = FALSE)

    }

    if(!is.null(confounder)) {

      if(!confounder %in% names(data)) {

        stop('The confounder variable absent from data.', call. = FALSE)

      }

    }

    if(!is.factor(data[[split_fct]])) {

      stop('split_fct needs to be a factor.', call. = FALSE)

    }

    stopifnot(is.logical(.parallel))

    p_value <- NULL
    p_adjusted <- NULL

    ## benchmarking ---------

    start_time <- Sys.time()
    message(paste('Testing for n =', length(variables), 'variables'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## data split and testing ---------

    if(is.null(confounder)) {

      data_lst <- map(variables, ~data[c(split_fct, .x)])

    } else {

      data_lst <- map(variables, ~data[c(split_fct, confounder, .x)])

    }

    if(.parallel) {

      plan('multisession')

      tst_res <-
        future_map2(variables,
                    data_lst,
                    function(x, y) nb_tester(data = y,
                                             split_fct = split_fct,
                                             confounder = confounder,
                                             dev_test = dev_test,
                                             variable = x, ...),
                    .options = furrr_options(seed = TRUE,
                                             packages = c('tibble')))

      plan('sequential')

    } else {

      tst_res <-
        map2(variables,
             data_lst,
             function(x, y) nb_tester(data = y,
                                      split_fct = split_fct,
                                      confounder = confounder,
                                      dev_test = dev_test,
                                      variable = x, ...))

    }

    tst_res <- transpose(compact(tst_res))

    tst_res <- map(tst_res, ~do.call('rbind', .x))

    map(tst_res,
        mutate,
        p_adjusted = p.adjust(p_value, method = adj_method),
        significant = ifelse(p_adjusted < 0.05, 'yes', 'no'))

  }


# END -----
