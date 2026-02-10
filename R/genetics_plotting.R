# Plots for genetics data

# Oncoplots for binary data -------

#' Oncoplots for binary data.
#'
#' @description
#' Generates oncoplots also known as waterfall plots for binary data.
#'
#' @details
#' The function takes a data frame that stores binary, i.e. 0/1 numeric data or
#' two-level factors, as the first argument and generates a two-color heat map
#' with observations in the X axis and features in the Y axis. The variables are
#' ordered by frequency of the event. The observations are arranged by
#' the frequency of event in the first variable with the highest frequency of
#' the  event, followed by the second one, till the n-th variable with the
#' lowest event frequency.
#'
#' The user may facet the plot by providing a name of the splitting variable
#' (`split_fct`, faceting of the columns) and a two-column data frame
#' as `variable_classification` (faceting of the rows). In the
#' `variable_classification` data frame, the first column stores variable names
#' and the second one codes for the variable subset assignment. Alternatively,
#' the user may provide a \code{\link[clustTools]{clust_analysis}} object with
#' results of a clustering analysis as the `variable_classification` argument.
#'
#' By default, i.e. when `one_plot = TRUE`, the function returns a single panel
#' of three plots: the oncoplot, an Y axis rug bar plot with overall
#' frequencies of the event per variable, and an X axis rug bar plot with
#' overall frequencies of the event per observation. Those frequencies may be
#' represented as absolute numbers or percentages as defined by `x_rug_scale`
#' and `y_rug_scale` arguments. If `one_plot = FALSE`, the function returns
#' a list of thoe three plots.
#'
#' @return a `ggplot` graphic object or, if `one_plot = FALSE`, a list of three
#' `ggplot` objects with the main plot and X and Y axis rug bar plots.
#'
#' @param data a data frame with binary 0/1 numeric data or two-level factors,
#' where 0 codes for absence and 1 codes for presence of the event of interest.
#' @param variables names of variables in `data` that store the binary event
#' indexes. If `NULL`, all variables will be used.
#' @param split_fct an optional name of a variable used for classification of
#' observation and column faceting of the plots.
#' @param variable_classification a two-variable data frame or a
#' \code{\link[clustTools]{clust_analysis}} object that defines categories of
#' the variables used for row faceting of the plots.
#' @param plot_title plot title
#' @param plot_subtitle plot subtitle.
#' @param x_lab X axis title.
#' @param y_lab Y axis title.
#' @param fill_lab title of the fill legend.
#' @param hide_x_axis_text logical, should the X axis text be hidden?
#' @param y_text_face font face of the Y axis text. An italic font face may be
#' useful e.g. for gene symbols.
#' @param cust_theme custom `ggplot` theme applied to the plots.
#' @param color_scale a two-element vector storing the colors for the
#' absence and presence of the event, respectively.
#' @param color_labels a two-element character vector with labels for the
#' event absence and presence, respectively.
#' @param x_rug_scale statistic of frequency of the event per observation,
#' presented in the X axis rug bar plot.
#' @param y_rug_scale statistic of frequency of the event per variable,
#' presented in the Y axis rug bar plot.
#' @param x_rug_fill fill color for bars of the X axis rug plot.
#' @param x_rug_color color for the bars of the X axis rug plot.
#' @param y_rug_fill fill color for bars of the Y axis rug plot.
#' @param y_rug_color color for the bars of the Y axis rug plot.
#' @param one_plot logical, should the oncoplot and the rug plots be stitched
#' into one panel?
#' @param rel_heights relative heights of the output plot panel. The defaults
#' are expected to work fine in most cases.
#' @param rel_widths relative widths of the output plot panel. The defaults
#' are expected to work fine in most cases.
#' @param legend_position placement of the oncoplot fill color legend in the
#' plot panel. By default it is placed in the right bottom corner (
#' `legend_position = 'inside'`).
#' @param plot_margins margins of the rug plots, which may be used for fine
#' tuning of their position in the output panel. The defaults are expected to
#' work fine in most cases.
#' @param title_panel_ratio ratio of the title/subtitle banner width to the
#' width of the output panel.
#' @param legend_panel_ratio ratio of the legend size to the size of the output
#' panel. Ignored if `legend_position = 'inside'`.
#'
#' @export

  plot_bionco <- function(data,
                          variables = NULL,
                          split_fct = NULL,
                          variable_classification = NULL,
                          plot_title = NULL,
                          plot_subtitle = NULL,
                          x_lab = 'observation',
                          y_lab = 'variable',
                          fill_lab = 'status',
                          hide_x_axis_text = FALSE,
                          y_text_face = 'plain',
                          cust_theme = theme_micro(),
                          color_scale = c('white', 'steelblue'),
                          color_labels = c('WT', 'altered'),
                          x_rug_scale = c('number', 'percentage'),
                          y_rug_scale = c('number', 'percentage'),
                          x_rug_fill = 'gray50',
                          x_rug_color = x_rug_fill,
                          y_rug_color = 'white',
                          y_rug_fill = color_scale[2],
                          one_plot = TRUE,
                          rel_heights = c(0.75, 0.25),
                          rel_widths = c(0.75, 0.25),
                          legend_position = c('inside', 'right',
                                              'bottom', 'top', 'left'),
                          plot_margins = margin(t = -35, r = -35,
                                                b = -10, l = 0,
                                                unit = 'pt'),
                          title_panel_ratio = c(0.1, 0.9),
                          legend_panel_ratio = c(0.1, 0.9)) {

    ## draws an oncoplot for binary data:
    ## features are presented in the X axis, samples are presented in the y axis

    ## entry control for the variables and the splitting factor -----

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

    ## entry control for the data, conversion, level numbers --------

    err_txt <- paste("'data' has to be a binary numeric data frame",
                     "or a data frame of two-level factors")

    if(!is.data.frame(data)) stop(err_txt, call. = FALSE)

    class_check <- map_lgl(data[variables], function(x) (is.numeric(x) | is.factor(x)))

    if(any(!class_check)) stop(err_txt, call. = FALSE)

    conv_numeric <- function(x) {

      if(is.numeric(x)) return(x)

      if(is.factor(x)) return(as.numeric(x) - 1)

      return(as.numeric(factor(x)) - 1)

    }

    data[variables] <- map_dfc(data[variables], conv_numeric)

    check_binary_df(data[variables])

    if(is.null(split_fct)) {

      data <- data[variables]

    } else {

      data <- data[c(split_fct, variables)]

    }

    data$sample_id <- paste0('obs_', 1:nrow(data))

    ## entry control for the variable classification -------

    if(!is.null(variable_classification)) {

      if(inherits(variable_classification, 'clust_analysis')) {

        variable_classification <-
          variable_classification$clust_assignment[, c('observation', 'clust_id')]

      } else {

        stopifnot(is.data.frame(variable_classification))

        if(ncol(variable_classification) == 1) {

          stop("'variable_classification' must have at least two columns.",
               call. = FALSE)

        }

      }

      variable_classification <-
        set_names(variable_classification[, 1:2],
                  c('variable', 'variable_subset'))

    }

    ## entry control for the remaining arguments ---------

    if(is.null(plot_subtitle)) {

      ## showing N numbers in the subtitle

      if(is.null(split_fct)) {

        plot_subtitle <- paste('n =', nrow(data))

      } else {

        split_counts <- table(data[[split_fct]])

        split_counts <- map2_chr(names(split_counts),
                                 split_counts,
                                 paste, sep = ': n = ')

        plot_subtitle <- paste(split_counts, collapse = ', ')

      }

    }

    stopifnot(is.logical(hide_x_axis_text))

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' has to be a valid 'ggplot' theme", call. = FALSE)

    }

    if(length(color_scale) < 2) {

      stop('At least two colors are required', call. = FALSE)

    }

    stopifnot(is.character(color_labels))

    if(length(color_labels) < 2) {

      stop('At least two  color labels are required', call. = FALSE)

    }

    x_rug_scale <- match.arg(x_rug_scale[1], c('number', 'percentage'))
    y_rug_scale <- match.arg(y_rug_scale[1], c('number', 'percentage'))

    stopifnot(is.logical(one_plot))

    if(one_plot) {

      stopifnot(is.numeric(rel_widths))
      stopifnot(is.numeric(rel_heights))

      if(length(rel_widths) < 2) {

        stop("'rel_widths' has to have at least two elements", call. = FALSE)

      }

      if(length(rel_heights) < 2) {

        stop("'rel_heights' has to have at least two elements", call. = FALSE)

      }

      rel_widths <- rel_widths[1:2]
      rel_heights <- rel_heights[1:2]

      legend_position <- match.arg(legend_position[1],
                                   c('inside', 'right',
                                     'bottom', 'top', 'left'))

    }

    ## plotting data, variable and sample order --------

    long_data <- pivot_longer(data,
                              cols = all_of(variables),
                              names_to = 'variable',
                              values_to = 'status')

    if(!is.null(variable_classification)) {

      long_data <- left_join(long_data,
                             variable_classification,
                             by = 'variable')

    }

    ## the variables will be ordered by the ascending frequency of events

    variable_order <- names(sort(colSums(data[variables])))

    long_data$variable <- factor(long_data$variable,  variable_order)

    ## ordering samples by an overall mutation frequency and frequency
    ## of mutations in particular genes

    sort_data <- data

    sort_data[variable_order] <-
      map_dfc(sort_data[variable_order], ~.x * -1)

    sample_order <- arrange(sort_data, pick(all_of(rev(variable_order))))

    sample_order <- sample_order$sample_id

    long_data$sample_id <- factor(long_data$sample_id, sample_order)

    ## base oncoplot -------

    sample_id <- NULL
    variable <- NULL
    status <- NULL

    onco_plot <- ggplot(long_data,
                        aes(x = sample_id,
                            y = variable,
                            fill = status)) +
      geom_tile() +
      scale_fill_gradient(low = color_scale[1],
                          high = color_scale[2],
                          limits = c(0, 1),
                          breaks = c(0, 1),
                          labels = color_labels[1:2],
                          name = fill_lab) +
      guides(fill = guide_legend()) +
      cust_theme +
      theme(axis.text.y = element_text(face = y_text_face)) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    if(hide_x_axis_text) {

      onco_plot <- onco_plot +
        theme(axis.text.x = element_blank())

    }

    ## rug plot data with overall frequencies of the event ---------

    variable <- NULL
    variable_subset <- NULL
    sample_id <- NULL

    ## X axis: frequency of the event per observation

    if(x_rug_scale == 'number') {

      x_rug_data <- rowSums(data[variables])

      x_rug_lab <- 'event, # variables'

    } else {

      x_rug_data <- rowMeans(data[variables]) * 100

      x_rug_lab <- 'event, % of variables'

    }

    ## Y axis: frequency of the event per variable

    if(y_rug_scale == 'number') {

      y_rug_data <- colSums(data[variables])

      y_rug_lab <- 'event, # observations'

    } else {

      y_rug_data <- colMeans(data[variables]) * 100

      y_rug_lab <- 'event, % of observations'

    }

    rug_data <- list()

    rug_data[['x']] <-
      tibble(sample_id = factor(paste0('obs_', 1:length(x_rug_data)),
                                levels(long_data[['sample_id']])),
             freq = unname(x_rug_data))

    rug_data[['y']] <-
      tibble(variable = factor(variables, levels(long_data[['variable']])),
             freq = unname(y_rug_data))

    if(!is.null(split_fct)) {

      rug_data[['x']] <-
        left_join(rug_data[['x']],
                  filter(long_data[c('sample_id', split_fct)],
                         !duplicated(sample_id)),
                  by = 'sample_id')

    }

    if(!is.null(variable_classification)) {

      rug_data[['y']] <-
        left_join(rug_data[['y']],
                  filter(long_data[c('variable', 'variable_subset')],
                         !duplicated(variable)),
                  by = 'variable')

    }

    ## rug plots ----------

    rug_plots <-
      pmap(list(x = rug_data,
                y = c('sample_id', 'freq'),
                z = c('freq', 'variable'),
                u = c(x_rug_color, y_rug_color),
                w = c(x_rug_fill, y_rug_fill)),
           function(x, y, z, u, w) ggplot(x,
                                          aes(x = .data[[y]],
                                              y = .data[[z]])) +
             geom_bar(stat = 'identity',
                      color = u,
                      fill = w) +
             cust_theme)

    rug_plots[['x']] <- rug_plots[['x']] +
      labs(x = x_lab,
           y = x_rug_lab) +
      theme(panel.grid.major.x = element_blank(),
            axis.text.y = element_text(face = 'plain'))

    rug_plots[['y']] <- rug_plots[['y']] +
      labs(x = y_rug_lab,
           y = y_lab) +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_text(face = y_text_face))

    if(hide_x_axis_text) {

      rug_plots[['x']] <- rug_plots[['x']] +
        theme(axis.text.x = element_blank())

    }

    ## no faceting: output ---------

    if(is.null(split_fct) & is.null(variable_classification) & !one_plot) {

      return(list(main = onco_plot,
                  observation = rug_plots[['x']],
                  variable = rug_plots[['y']]))

    }

    ## optional faceting of the rug plots --------

    if(!is.null(split_fct)) {

      facet_formula <- as.formula(paste('. ~', split_fct))

      rug_plots[['x']] <- rug_plots[['x']] +
        facet_grid(facet_formula,
                   scales = 'free',
                   space = 'free')

    }

    if(!is.null(variable_classification)) {

      rug_plots[['y']] <- rug_plots[['y']] +
        facet_grid(variable_subset ~ .,
                   scales = 'free',
                   space = 'free')

    }

    ## optional faceting of the main plot ------

    if(!is.null(split_fct) & is.null(variable_classification) & !one_plot) {

      facet_formula <- paste('. ~', split_fct)

    } else if(is.null(split_fct) & !is.null(variable_classification)) {

      facet_formula <- 'variable_subset ~ .'

    } else if(!is.null(split_fct) & !is.null(variable_classification)) {

      facet_formula <- paste('variable_subset ~', split_fct)

    }

    if(!is.null(split_fct) | !is.null(variable_classification)) {

      facet_formula <- as.formula(facet_formula)

      onco_plot <- onco_plot +
        facet_grid(facet_formula,
                   scales = 'free',
                   space = 'free')

    }

    onco_lst <- list(main = onco_plot,
                     observation = rug_plots[['x']],
                     variable = rug_plots[['y']])

    if(!one_plot) return(onco_lst)

    ## optional plot panel -------

    onco_lst[['main']] <- onco_lst[['main']] +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            strip.background.y = element_blank(),
            strip.text.y = element_blank(),
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            plot.margin = margin(t = 0, r = 0, l = 0, b = 0))

    onco_lst[['variable']] <- onco_lst[['variable']] +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank())

    onco_lst[['observation']] <- onco_lst[['observation']] +
      theme(strip.text.x = element_blank(),
            strip.background.x = element_blank())

    onco_lst <-
      map(onco_lst[c('main', 'variable', 'observation')],
          ~.x +
            theme(legend.position = 'none'))

    if(legend_position == 'inside') {

      onco_lst <- c(onco_lst, list(get_legend(onco_plot)))

    }

    onco_lst[2:3] <- map(onco_lst[2:3],
                         ~.x + theme(plot.margin = plot_margins))

    onco_panel <-
      plot_grid(plotlist = onco_lst,
                ncol = 2,
                rel_heights = rel_heights,
                rel_widths = rel_widths,
                align = 'hv',
                axis = 'tblr')

    ## appending with the title, subtitle and, optionally, the legend

    x_pos <- NULL
    y_pos <- NULL

    if(is.null(plot_subtitle)) {

      title_banner <-
        plot_grid(get_title(onco_plot))

    } else {

      title_banner <-
        plot_grid(get_title(onco_plot),
                  get_subtitle(onco_plot),
                  nrow = 2)

    }

    onco_panel <- plot_grid(title_banner +
                              theme(plot.margin = cust_theme$plot.margin),
                            onco_panel,
                            nrow = 2,
                            rel_heights = title_panel_ratio)

    if(legend_position == 'bottom') {

      onco_panel <-
        plot_grid(onco_panel,
                  get_legend(onco_plot +
                               theme(legend.position = 'bottom')),
                  nrow = 2,
                  rel_heights = rev(legend_panel_ratio))

    } else if(legend_position == 'left') {

      onco_panel <-
        plot_grid(get_legend(onco_plot +
                               theme(legend.position = 'right')),
                  onco_panel,
                  ncol = 2,
                  rel_widths = legend_panel_ratio)

    } else if(legend_position == 'right') {

      onco_panel <-
        plot_grid(onco_panel,
                  get_legend(onco_plot +
                               theme(legend.position = 'right')),
                  ncol = 2,
                  rel_widths = rev(legend_panel_ratio))

    } else if(legend_position == 'top') {

      onco_panel <-
        plot_grid(get_legend(onco_plot +
                               theme(legend.position = 'bottom')),
                  onco_panel,
                  nrow = 2,
                  rel_widths)

    }

    onco_panel +
      theme(plot.margin = cust_theme$plot.margin)

  }

# Stack plots for binary data -------

#' Stack plots for binary data.
#'
#' @description
#' The function draws a classical stack plot with counts or percentages of the
#' event in binary data.
#'
#' @details
#' The function takes a data frame that stores binary, i.e. 0/1 numeric data or
#' two-level factors, as the first argument and generates a two-color stack
#' plot. The variables are ordered by frequency of the event.
#'
#' The user may facet the plot by providing a name of the splitting variable
#' (`split_fct`, faceting of the columns) and a two-column data frame
#' as `variable_classification` (faceting of the rows). In the
#' `variable_classification` data frame, the first column stores variable names
#' and the second one codes for the variable subset assignment. Alternatively,
#' the user may provide a \code{\link[clustTools]{clust_analysis}} object with
#' results of a clustering analysis as the `variable_classification` argument.
#'
#' @inheritParams plot_bionco
#' @param scale statistic of frequency to be presented in the plot: absolute
#' number (default) or frequency.
#' @param reverse_status logical, should the levels of the status be reversed?
#' By default, i.e. when `reverse_status = FALSE`, the event coded as 1 is
#' plotted first.
#' @param show_frequencies logical, should labels with the frequencies be
#' displayed in the plot?
#' @param bar_rim_color color of rims of the bars.
#' @param signif_digits number of significant digits used for rounding of the
#' text labels. Ignored if `show_frequencies = FALSE` or `scale = 'number'`.
#' @param text_size size of the frequency labels.
#' @param text_color color of the text in the frequency labels.
#' @param facet_scales a string that indicates which scales can vary between
#' the plot facets. By default, `facet_scales = 'free_y'`, which means that
#' only the Y axis scale will differ. See \code{\link[ggplot2]{facet_grid}} for
#' details.
#' @param facet_space a string that indicates which spaces can vary between
#' the plot facets. By default, `facet_space = 'free_y'`, which means that
#' only the Y axis scale spaces will differ. See
#' \code{\link[ggplot2]{facet_grid}} for details.
#' @param ... extra arguments that control appearance of the frequency labels
#' passed to \code{\link[ggplot2]{geom_label}}.
#'
#' @export

  plot_bistack <- function(data,
                           variables = NULL,
                           split_fct = NULL,
                           variable_classification = NULL,
                           scale = c('number', 'percentage'),
                           reverse_status = FALSE,
                           show_frequencies = FALSE,
                           plot_title = NULL,
                           plot_subtitle = NULL,
                           x_lab = 'events, # observations',
                           y_lab = 'variable',
                           fill_lab = 'status',
                           cust_theme = theme_micro(),
                           color_scale = c('white', 'steelblue'),
                           color_labels = c('WT', 'altered'),
                           bar_rim_color = 'gray60',
                           signif_digits = 2,
                           text_size = 2.75,
                           text_color = 'black',
                           facet_scales = c('free_y', 'free_x', 'free'),
                           facet_space = c('free_y', 'free_x', 'free'),
                           ...) {

    ## entry control for the variables and the splitting factor -----

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

    ## entry control for the data, conversion, level numbers --------

    err_txt <- paste("'data' has to be a binary numeric data frame",
                     "or a data frame of two-level factors")

    if(!is.data.frame(data)) stop(err_txt, call. = FALSE)

    class_check <- map_lgl(data[variables], function(x) (is.numeric(x) | is.factor(x)))

    if(any(!class_check)) stop(err_txt, call. = FALSE)

    conv_numeric <- function(x) {

      if(is.numeric(x)) return(x)

      if(is.factor(x)) return(as.numeric(x) - 1)

      return(as.numeric(factor(x)) - 1)

    }

    data[variables] <- map_dfc(data[variables], conv_numeric)

    check_binary_df(data[variables])

    if(is.null(split_fct)) {

      data <- data[variables]

    } else {

      data <- data[c(split_fct, variables)]

    }

    data$sample_id <- paste0('obs_', 1:nrow(data))

    ## entry control for the variable classification -------

    if(!is.null(variable_classification)) {

      if(inherits(variable_classification, 'clust_analysis')) {

        variable_classification <-
          variable_classification$clust_assignment[, c('observation', 'clust_id')]

      } else {

        stopifnot(is.data.frame(variable_classification))

        if(ncol(variable_classification) == 1) {

          stop("'variable_classification' must have at least two columns.",
               call. = FALSE)

        }

      }

      variable_classification <-
        set_names(variable_classification[, 1:2],
                  c('variable', 'variable_subset'))

    }

    ## entry control for the remaining arguments --------

    scale <- match.arg(scale[1], c('number', 'percentage'))

    stopifnot(is.logical(reverse_status))
    stopifnot(is.logical(show_frequencies))

    if(is.null(plot_subtitle)) {

      ## showing N numbers in the subtitle

      if(is.null(split_fct)) {

        plot_subtitle <- paste('n =', nrow(data))

      } else {

        split_counts <- table(data[[split_fct]])

        split_counts <- map2_chr(names(split_counts),
                                 split_counts,
                                 paste, sep = ': n = ')

        plot_subtitle <- paste(split_counts, collapse = ', ')

      }

    }

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' has to be a valid 'ggplot' theme", call. = FALSE)

    }

    if(length(color_scale) < 2) {

      stop('At least two colors are required', call. = FALSE)

    }

    stopifnot(is.character(color_labels))

    if(length(color_labels) < 2) {

      stop('At least two  color labels are required', call. = FALSE)

    }

    stopifnot(is.numeric(text_size))
    stopifnot(is.numeric(signif_digits))

    signif_digits <- as.integer(signif_digits[1])

    facet_scales <- match.arg(facet_scales, c('free_y', 'free_x', 'free'))
    facet_space <- match.arg(facet_space, c('free_y', 'free_x', 'free'))

    ## plotting data ----------

    plot_tbl <- count_binary(data,
                             variables,
                             split_fct)

    if(scale == 'number') {

      plot_tbl[['present']] <- plot_tbl[['n']]

      plot_tbl[['absent']] <- plot_tbl[['n_total']] - plot_tbl[['n']]

      round_fun <- identity

    } else {

      plot_tbl[['present']] <- plot_tbl[['percent']]

      plot_tbl[['absent']] <- 100 - plot_tbl[['percent']]

      round_fun <- function(x) signif(x, signif_digits)

    }

    plot_tbl <- pivot_longer(plot_tbl,
                             cols = all_of(c('present', 'absent')),
                             values_to = 'freq',
                             names_to = 'status')

    if(!is.null(split_fct)) {

      plot_tbl <- plot_tbl[c('variable', split_fct, 'status', 'freq')]

    } else {

      plot_tbl <- plot_tbl[c('variable', 'status', 'freq')]

    }

    ## variable classification and ordering by frequency of events

    if(!is.null(variable_classification)) {

      plot_tbl <- left_join(plot_tbl,
                            variable_classification,
                            by = 'variable')

    }

    variable_order <- names(sort(colSums(data[variables])))

    plot_tbl[['variable']] <- factor(plot_tbl[['variable']], variable_order)

    if(is.factor(data[[split_fct]])) {

      plot_tbl[[split_fct]] <-
        factor(plot_tbl[[split_fct]], levels(data[[split_fct]]))

    }

    status_levs <- c('absent', 'present')

    if(reverse_status) {

      status_levs <- rev(status_levs)

    }

    plot_tbl[['status']] <- factor(plot_tbl[['status']], status_levs)

    ## position of numeric plot labels

    status <- NULL
    freq <- NULL
    plot_pos <- NULL
    variable <- NULL
    x <- NULL

    if(show_frequencies) {

      plot_tbl <- arrange(plot_tbl, desc(status))

      if(!is.null(split_fct)) {

        plot_tbl <- group_by(plot_tbl, variable, .data[[split_fct]])

      } else {

        plot_tbl <- group_by(plot_tbl, variable)

      }

      plot_tbl <- mutate(plot_tbl,
                         plot_pos = cumsum(freq) - 0.5 * freq)

      plot_tbl <- ungroup(plot_tbl)

    }

    ## basic stack plot -------

    stack_plot <- ggplot(plot_tbl,
                         aes(x = freq,
                             y = variable,
                             fill = status)) +
      geom_bar(stat = 'identity',
               color = bar_rim_color) +
      scale_fill_manual(values = color_scale,
                        labels = color_labels,
                        name = fill_lab) +
      guides(fill = guide_legend(reverse = reverse_status)) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    if(show_frequencies) {

      stack_plot <- stack_plot +
        geom_label(aes(label = round_fun(freq),
                       x = plot_pos),
                   size = text_size,
                   color = text_color,
                   show.legend = FALSE, ...)

    }

    if(is.null(variable_classification) & is.null(split_fct)) {

      return(stack_plot)

    }

    ## optional faceting of the plot -------

    if(!is.null(split_fct) & is.null(variable_classification)) {

      facet_formula <- paste('. ~', split_fct)

    } else if(is.null(split_fct) & !is.null(variable_classification)) {

      facet_formula <- 'variable_subset ~ .'

    } else if(!is.null(split_fct) & !is.null(variable_classification)) {

      facet_formula <- paste('variable_subset ~', split_fct)

    }

    stack_plot +
      facet_grid(facet_formula,
                 scales = facet_scales,
                 space = facet_space)

  }

# END -------
