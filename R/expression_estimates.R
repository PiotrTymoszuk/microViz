# Helper functions which kelp at obtaining expression regulation estimates
# for different contrast settings.

#' @include imports.R

  NULL

# Deviation from mean --------

#' Average deviation from the mean or median.
#'
#' @description
#' The function computes deviations from the grand measure of central tendency
#' and spits them by an user defined factor.
#'
#' @details
#' A bunch of metrics of the central tendency is available, both for
#' computation of the grand statistic (argument `grand_center`) and averaging
#' of the the deviations (argument `split_center`). They are: arithmetic mean
#' ('mean', default), median ('median'), geometric mean ('gmean') and
#' harmonic mean ('hmean'). Note that, not all of those statistic are
#' appropriate for all types of numeric data: e.g. geometric mean for a mixture
#' of positive and negative values.
#' The function computes also the Student's t statistic for the difference
#' of the averaged deviations and the grand mean as well as p-value of a single
#' sample t test. Prior to using those values, checking for the normality
#' assumption is highly recommended.
#' Finally, 95% confidence intervals are returned for three scenarios. For
#' `ci = 'distr'`, confidence intervals of the average deviation from the grand
#' mean are computed based on the critical values of the t statistic and
#' standard errors. By contrast, `ci = 'percentile'` and `ci = 'bca'` calculate
#' 95% CI for single deviations from the grand means with the percentile and
#' bias-corrected method, respectively.
#'
#' @return
#' A data frame with the grand averages (`grand_center`), average deviations
#' from the grand average (`deviation_center`), standard deviations and
#' standard errors (`deviation_sd` and `deviation-sem`), t statistic, degrees of
#' freedom (`df`), limits of the 95% confidence interval (`lower_ci` and
#' `upper_ci`), raw and multiple testing-adjusted one-sample T test p values
#' (`p_value` and `p_adjusted`) for variables split by the user-defined
#' splitting feature.
#'
#' @param data a data frame.
#' @param split_fct name of the splitting factor.
#' @param variables names of numeric variables, for which deviations will be
#' computed.
#' @param grand_center statistic of the grand central tendency, see Details.
#' @param split_center statistic of the central tendency in subsets defined
#' by `split_fct`.
#' @param ci type of 95% confidence intervals, as specified in Details.
#' @param adj_method method of multiple testing adjustment.
#' See: \code{\link[stats]{p.adjust}}.
#' @param .parallel logical: should the function be run in parallel?
#'
#' @export

  avg_deviation <- function(data,
                            split_fct,
                            variables,
                            grand_center = c('mean', 'median', 'gmean', 'hmean'),
                            split_center = c('mean', 'median', 'gmean', 'hmean'),
                            ci = c('distr', 'percentile', 'bca'),
                            adj_method = 'BH',
                            .parallel = FALSE) {

    ## input control --------

    stopifnot(is.data.frame(data))

    if(any(!variables %in% names(data))) {

      stop('Some variables are absent from the data.', call. = FALSE)

    }

    if(!split_fct %in% names(data)) {

      stop('The splitting factor is missing from the data.', call. = FALSE)

    }

    if(!is.factor(data[[split_fct]])) {

      stop('The splitting variable has to be a factor.', call. = FALSE)

    }

    data[[split_fct]] <- droplevels(data[[split_fct]])

    grand_center <-
      match.arg(grand_center[1],
                c('mean', 'median', 'gmean', 'hmean'))

    if(!is.null(split_center)) {

      split_center <- grand_center

    }

    split_center <-
      match.arg(split_center[1],
                c('mean', 'median', 'gmean', 'hmean'))

    ci <- match.arg(ci[1], c('distr', 'percentile', 'bca'))

    stopifnot(is.logical(.parallel))

    ## data, central tendency function and CI function-----

    grand_data <- data[variables]

    data_check <- purrr::map_lgl(grand_data, is.numeric)

    if(any(!data_check)) {

      stop('Numeric data is required.', call. = FALSE)

    }

    grand_fun <-
      switch(grand_center,
             mean = function(x) colMeans(x, na.rm = TRUE),
             median = function(x) colMedians(x, na.rm = TRUE, .parallel = .parallel),
             gmean = function(x) colHmeans(x, na.rm = TRUE, .parallel = .parallel),
             hmean = function(x) colHmeans(x, na.rm = TRUE, .parallel = .parallel))

    split_fun <-
      switch(split_center,
             mean = function(x) colMeans(x, na.rm = TRUE),
             median = function(x) colMedians(x, na.rm = TRUE, .parallel = .parallel),
             gmean = function(x) colHmeans(x, na.rm = TRUE, .parallel = .parallel),
             hmean = function(x) colHmeans(x, na.rm = TRUE, .parallel = .parallel))

    ci_fun <-
      switch(ci,
             percentile = function(x) stats::quantile(x,
                                                      c(0.025, 0.975),
                                                      na.rm = TRUE),
             bca = function(x) coxed::bca(x[!is.na(x)]))

    ## grand stats and single deviations from the grand stat -------

    grand_centers <- grand_fun(grand_data)

    diffs <- purrr::map2(grand_data, grand_centers, `-`)

    split_diffs <- purrr::map(diffs, split, f = data[[split_fct]])

    split_diffs <- purrr::transpose(split_diffs)

    split_diffs <- purrr::map(split_diffs, as.data.frame)

    ## average, SD and confidence intervals for the deviations ------

    split_centers <- purrr::map(split_diffs, split_fun)

    split_variance <- list()

    split_ci <- list()

    for(i in names(split_centers)) {

      split_variance[[i]] <-
        purrr::map2_dbl(split_diffs[[i]], split_centers[[i]],
                        ~sum((.x - .y)^2, na.rm = TRUE)/(length(.x) - 1))

      if(ci != 'distr') {

        split_ci[[i]] <- purrr::map(split_diffs[[i]], ci_fun)

        split_ci[[i]] <- do.call('rbind', split_ci[[i]])

        colnames(split_ci[[i]]) <- c('lower_ci', 'upper_ci')

      }

    }

    ## output ----------

    output <- list()

    split_counts <- table(data[[split_fct]])

    for(i in names(split_centers)) {

     output[[i]] <-
       tibble::tibble(!!split_fct := i,
                      n = split_counts[i],
                      variable = names(grand_centers),
                      grand_center = grand_centers,
                      deviation_center = split_centers[[i]],
                      deviation_sd = sqrt(split_variance[[i]]),
                      deviation_sem = sqrt(split_variance[[i]])/sqrt(split_counts[i]))

     ## t statistic, its degrees of freedom and one sample t test

     output[[i]][['t']] <-
       output[[i]][['deviation_center']]/output[[i]][['deviation_sem']]

     output[[i]][['df']] <- split_counts[i] - 1

     output[[i]][['p_value']] <-
       stats::pt(abs(output[[i]][['t']]),
                 df = output[[i]][['df']],
                 lower.tail = FALSE)

     if(ci != 'distr') {

       output[[i]] <- cbind(output[[i]],
                            tibble::as_tibble(split_ci[[i]]))

     } else {

      output[[i]][['lower_ci']] <- output[[i]][['deviation_center']] +
        stats::qt(0.025, df = output[[i]][['df']]) *
        output[[i]][['deviation_sem']]

      output[[i]][['upper_ci']] <- output[[i]][['deviation_center']] +
        stats::qt(0.975, df = output[[i]][['df']]) *
        output[[i]][['deviation_sem']]

     }

     output[[i]] <- tibble::as_tibble(output[[i]])

    }

    output <- do.call('rbind', output)

    output[[split_fct]] <- factor(output[[split_fct]],
                                  levels(data[[split_fct]]))

    output[['p_adjusted']] <-
      stats::p.adjust(output[['p_value']], method = adj_method)

    output

  }

# END ----
