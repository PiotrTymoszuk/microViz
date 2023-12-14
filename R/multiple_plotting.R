# visualization tools for multi-data set analysis results

#' @include imports.R

  NULL

# Scatter plots of estimates --------

#' Scatter/Forest plot of estimates.
#'
#' @description
#' Draws a scatter plot of estimates specified by the user. Point color codes
#' for significance and regulation status. As an option, a central statistic
#' over all estimates such as mean or median with an error measure
#' (SD, SEM, etc.) may be displayed.
#'
#' @return a `ggplot` object.
#'
#' @param data a data frame with analysis results.
#' @param label_variable name of the variable storing feature names.
#' @param regulation_variable name of the variable storing the estimate to
#' be plotted.
#' @param p_variable optional, name of the variable storing statistical test
#' results.
#' @param signif_level significance cutoff used for identification of
#' significant estimates, ignored if `p_variable = NULL`.
#' @param regulation_level estimate value cutoff used for identification of
#' significant estimates.
#' @param fill_scale colors coding for the estimate status.
#' @param central_stat measure of central tendency to be plotted: none, mean,
#' median, geometric mean ('gmean') or harmonic mean ('hmean').
#' @param error_stat measure of distribution to be plotted as a horizontal bar:
#' none, percentile confidence interval, BCA (bias-corrected and accelerated)
#' confidence interval , or interquartile range ('iqr').
#' @param conf_level width of the confidence interval.
#' @param point_size size of the data points.
#' @param point_alpha alpha of the data points.
#' @param point_wjitter data point jittering width.
#' @param point_hjitter data point jittering height.
#' @param center_symbol_size size of the symbol representing the central
#' tendency measure.
#' @param stat_color color of the symbols of the central tendency measure and
#' the error bars.
#' @param error_bar_size size of the error bar.
#' @param plot_title plot title.
#' @param plot_subtitle plot_subtitle.
#' @param x_lab X axis title.
#' @param y_lab Y axis title.
#' @param fill_lab fill scale title.
#' @param plot_tag plot tag.
#' @param cust_theme custom `ggplot` theme.#'
#'
#' @export

  draw_estimate_scatter <- function(data,
                                    label_variable,
                                    regulation_variable,
                                    p_variable = NULL,
                                    signif_level = 0.05,
                                    regulation_level = 0,
                                    fill_scale = c(upregulated = 'firebrick',
                                                    downregulated = 'steelblue',
                                                    ns = 'gray60'),
                                    central_stat = c('none',
                                                     'mean',
                                                     'median',
                                                     'gmean',
                                                     'hmean'),
                                    error_stat = c('none',
                                                   'percentile_ci',
                                                   'bca_ci',
                                                   'iqr'),
                                    conf_level = 0.95,
                                    point_size = 2,
                                    point_alpha = 0.75,
                                    point_wjitter = 0,
                                    point_hjitter = 0.1,
                                    center_symbol_size = 3,
                                    stat_color = 'orangered3',
                                    error_bar_size = 0.75,
                                    plot_title = NULL,
                                    plot_subtitle = NULL,
                                    x_lab = 'Regulation',
                                    y_lab = NULL,
                                    fill_lab = regulation_variable,
                                    plot_tag = NULL,
                                    cust_theme = microViz::theme_micro()) {

    ## input control -------

    if(!is.data.frame(data)) {

      stop('A data frame is required.', call. = FALSE)

    }

    vars_to_check <- c(label_variable, regulation_variable)

    if(!is.null(p_variable)) {

      vars_to_check <- c(vars_to_check, p_variable)

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
    stopifnot(is.numeric(point_alpha))

    central_stat <- match.arg(central_stat[1],
                              c('none',
                                'mean',
                                'median',
                                'gmean',
                                'hmean'))

    error_stat <- match.arg(error_stat[1],
                            c('none',
                              'percentile_ci',
                              'bca_ci',
                              'sd',
                              'sem',
                              '2sem',
                              'iqr'))

    stopifnot(inherits(cust_theme, 'theme'))

    ## base plotting data -------

    .regulation <- NULL
    .significant <- NULL

    if(is.null(p_variable)) {

      plot_tbl <-
        mutate(data,
               .regulation = ifelse(.data[[regulation_variable]] == 0,
                                    'ns',
                                    ifelse(.data[[regulation_variable]] > regulation_level,
                                           'upregulated',
                                           ifelse(.data[[regulation_variable]] < -regulation_level,
                                                  'downregulated', 'ns'))))

    } else {

      plot_tbl <-
        mutate(data,
               .significant = ifelse(.data[[p_variable]] < signif_level,
                                     'yes', 'no'),
               .regulation = ifelse(.data[[regulation_variable]] == 0 | .significant == 'no',
                                    'ns',
                                    ifelse(.data[[regulation_variable]] > regulation_level,
                                           'upregulated',
                                           ifelse(.data[[regulation_variable]] < -regulation_level,
                                                  'downregulated', 'ns'))))

    }

    plot_tbl <-
      mutate(plot_tbl,
             .regulation = factor(.regulation,
                                  c('upregulated', 'downregulated', 'ns')))

    ## base plot --------

    base_plot <- ggplot(plot_tbl,
                        aes(x = .data[[regulation_variable]],
                            y = reorder(.data[[label_variable]],
                                        .data[[regulation_variable]]),
                            fill = .regulation)) +
      geom_point(shape = 21,
                 size = point_size,
                 alpha = point_alpha,
                 position = position_jitter(width = point_wjitter,
                                            height = point_hjitter)) +
      scale_fill_manual(values = fill_scale) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab,
           fill = fill_lab,
           tag = plot_tag)

    if(central_stat == 'none' & error_stat == 'none') return(base_plot)

    ## appending with the central tendency stat and errors -----

    if(central_stat != 'none') {

      central_fun <-
        switch(central_stat,
               mean = function(x) mean(x, na.rm = TRUE),
               median = function(x) median(x, na.rm = TRUE),
               gmean = function(x) Gmean(x, na.rm = TRUE),
               hmean = function(x) Hmean(x, na.rm = TRUE))

      central_data <- group_by(plot_tbl, .data[[label_variable]])

      central_data <-
        summarise(central_data,
                  !!regulation_variable := central_fun(.data[[regulation_variable]]))

      central_data <- mutate(central_data, .regulation = 'ns')

      base_plot <- base_plot +
        geom_point(data = central_data,
                   shape = 23,
                   size = center_symbol_size,
                   color = stat_color,
                   fill = stat_color)

    }

    if(error_stat != 'none') {

      error_fun <-
        switch(error_stat,
               percentile_ci = function(x) perCI(x, conf_level),
               bca_ci = function(x) bcaCI(x, conf_level),
               iqr = function(x) quantile(x, c(0.25, 0.75)))

      error_data <-
        split(plot_tbl[[regulation_variable]], plot_tbl[[label_variable]])

      error_data <- map(error_data, error_fun)

      error_data <- map(error_data, reduce, cbind)

      error_data <- map(error_data, as.data.frame)

      lower_bound <- NULL
      upper_bound <- NULL

      error_data <- map(error_data,
                        set_names,
                        c('lower_bound', 'upper_bound'))

      error_data <- map2_dfr(error_data, names(error_data),
                             ~mutate(.x,
                                     !!label_variable := .y,
                                     !!regulation_variable := 1))

      base_plot <- base_plot +
        geom_errorbarh(data = error_data,
                       aes(xmin = lower_bound,
                           xmax = upper_bound),
                       color = stat_color,
                       height = 0,
                       linewidth = error_bar_size)

    }

    return(base_plot)

  }

# Ellipsoid plot ----------
