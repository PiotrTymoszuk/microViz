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
#' @param conf_level width of the confidence interval.
#' @param adj_method method of multiple testing adjustment.
#' See: \code{\link[stats]{p.adjust}}.
#'
#' @export

  avg_deviation <- function(data,
                            split_fct,
                            variables,
                            grand_center = c('mean', 'median', 'gmean', 'hmean'),
                            split_center = c('mean', 'median', 'gmean', 'hmean'),
                            ci = c('distr', 'percentile', 'bca'),
                            conf_level = 0.95,
                            adj_method = 'BH') {

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

    ## data, central tendency function and CI function-----

    grand_data <- data[variables]

    check_df(grand_data)

    grand_fun <-
      switch(grand_center,
             mean = function(x) colMeans(x, na.rm = TRUE),
             median = function(x) colMedians(x, na.rm = TRUE),
             gmean = function(x) colGmeans(x, na.rm = TRUE),
             hmean = function(x) colHmeans(x, na.rm = TRUE))

    split_fun <-
      switch(split_center,
             mean = function(x) colMeans(x, na.rm = TRUE),
             median = function(x) colMedians(x, na.rm = TRUE),
             gmean = function(x) colGmeans(x, na.rm = TRUE),
             hmean = function(x) colHmeans(x, na.rm = TRUE))

    ## grand stats and single deviations from the grand stat -------

    grand_centers <- grand_fun(grand_data)

    diffs <- Delta(as.matrix(grand_data), grand_centers)

    split_diffs <- split(as.data.frame(diffs), f = data[[split_fct]])

    split_diffs <- map(split_diffs, as.data.frame)

    ## average, SD and confidence intervals for the deviations ------

    split_centers <- map(split_diffs, split_fun)

    split_variance <- map(split_diffs, colVars)

    if(ci != 'distr') {

      split_ci <- map(split_diffs,
                      colCI,
                      conf_level = conf_level,
                      method = ci,
                      na.rm = TRUE)

      split_ci <- do.call('rbind', split_ci)

    }

    ## output ----------

    output <- list()

    split_counts <- table(data[[split_fct]])

    for(i in names(split_centers)) {

      output[[i]] <-
        tibble(!!split_fct := i,
               n = split_counts[i],
               variable = names(grand_centers),
               grand_center = grand_centers,
               deviation_center = split_centers[[i]],
               deviation_sd = sqrt(split_variance[[i]]),
               deviation_sem = sqrt(split_variance[[i]])/sqrt(split_counts[i]))

    }

    output <- do.call('rbind', output)

    ## t statistic, its degrees of freedom and p values ------

   output[['t']] <-
     output[['deviation_center']]/output[['deviation_sem']]

   output[['df']] <- split_counts[i] - 1

   output[['p_value']] <-
     pt(abs(output[['t']]),
        df = output[['df']],
        lower.tail = FALSE)

   if(ci != 'distr') {

     output <- cbind(output, as_tibble(split_ci))

   } else {

     output[['lower_ci']] <- output[['deviation_center']] +
       qt(0.025, df = output[['df']]) *
       output[['deviation_sem']]

     output[['upper_ci']] <- output[['deviation_center']] +
       qt(0.975, df = output[['df']]) *
       output[['deviation_sem']]

   }

   output[[split_fct]] <- factor(output[[split_fct]],
                                 levels(data[[split_fct]]))

   output[['p_adjusted']] <-
     p.adjust(output[['p_value']], method = adj_method)

   as_tibble(output[c(split_fct,
                      'n', 'variable',
                      'grand_center',
                      'deviation_center',
                      'deviation_sd',
                      'deviation_sem',
                      't', 'df',
                      'lower_ci', 'upper_ci',
                      'p_value', 'p_adjusted')])

  }

# END ----
