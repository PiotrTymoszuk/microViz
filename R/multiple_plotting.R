# visualization tools for multi-data set analysis results

#' @include imports.R

  NULL

# Before-after plots for lists of data frames --------

#' Before-after or change plots for multiple data frames.
#'
#' @description
#' Draws a list of plots of average statistics (i.e. mean or median) of numeric
#' variables split by a factor splitting variable. Averages for the same data
#' set are connected by lines.
#'
#' @return a list of `ggplot` objects.
#'
#' @param data a list of data frames.
#' @param variables a vector with names of the variables of interest.
#' @param split_fct name of the splitting factor.
#' @param normalize logical, should the data frame variables be normalized prior
#' to plotting?
#' @param norm_center defines centering of the variable during scaling: mean
#' (default) or median. Ignored if `normalize = FALSE`.
#' @param average_fun a function used to calculate average values of the
#' variables within levels of the splitting factor.
#' @param show_whiskers logical, should whiskers representing a quantile range
#' be displayed in the plot? Defaults to false.
#' @param whisker_size size of the whiskers.
#' @param probs a numeric vector of two quantiles defining the range for the
#' whiskers. Ignored if `show_whiskers = FALSE`.
#' @param dodge a numeric specifying dodging of the points, whiskers, labels,
#' and lines.
#' @param point_size size of the data points.
#' @param point_alpha alpha of the points and, optionally, whiskers.
#' @param linetype type of the line connecting the points.
#' @param line_alpha alpha of the connecting lines.
#' @param line_color color of the connecting lines.
#' @param show_txt logical, should labels with names of the data sets be
#' displayed?
#' @param txt_size size of the label text.
#' @param txt_position placement of the text labels.
#' @param txt_color color of the text. If `NULL`, the color is specified by
#' levels of the splitting factor.
#' @param txt_face fornt face of the text labels.
#' @param labeller a function transforming the data set names to the text labels
#' in the plot. Defaults to `identity`, i.e. the names are displayed as the are.
#' @param cust_theme a custom `ggplot` theme.
#' @param plot_titles a character vector with plot titles.
#' @param plot_subtitles a character vector with lot subtitles.
#' @param x_lab X axis title.
#' @param y_lab Y axis title.
#' @param palette colors for levels of the splitting factor. If `NULL`, default
#' R color and fill scales will be used.
#'
#' @export

  plot_common_change <- function(data,
                                 variables,
                                 split_fct,
                                 normalize = TRUE,
                                 norm_center = c('mean', 'median'),
                                 average_fun = colMeans,
                                 show_whiskers = FALSE,
                                 whisker_size = 0.5,
                                 probs = c(0.025, 0.975),
                                 dodge = 0,
                                 point_size = 2,
                                 point_alpha = 1,
                                 linetype = 'solid',
                                 line_alpha = 0.75,
                                 line_color = 'gray70',
                                 show_txt = TRUE,
                                 txt_size = 2.75,
                                 txt_position = c('left', 'right'),
                                 txt_color = NULL,
                                 txt_face = 'plain',
                                 labeller = identity,
                                 cust_theme = microViz::theme_micro(),
                                 plot_titles = variables,
                                 plot_subtitles = NULL,
                                 x_lab = split_fct,
                                 y_lab = 'average Z-score',
                                 palette = NULL) {


    ## entry control for the data argument ---------

    stopifnot(is.list(data))

    classes <- map_lgl(data, is.data.frame)

    if(any(!classes)) {

      stop("At least element of 'data' is not a data frame.", call. = FALSE)

    }

    split_present <- map_lgl(data, ~split_fct %in% names(.x))

    if(any(!split_present)) {

      stop("'split_fct' is missing from at least one data set.", call. = FALSE)

    }

    split_is_factor <- map_lgl(data, ~is.factor(.x[[split_fct]]))

    if(any(!split_is_factor)) {

      stop('In at least one of data frames, splitting variable is not factor.',
           call. = FALSE)

    }

    variables_present <- map_lgl(data, ~all(variables %in% names(.x)))

    if(any(!variables_present)) {

      stop('Some of the variables are missing from the data.', call. = FALSE)

    }

    variables_numeric <-
      map(data,
          function(ds) map_lgl(ds[variables], is.numeric))

    variables_numeric <- map_lgl(variables_numeric, all)

    if(any(!variables_numeric)) {

      stop('Some of the variables are not numeric.', call. = FALSE)

    }

    ## entry control for remaining arguments --------

    stopifnot(is.logical(normalize))

    norm_center <- match.arg(norm_center[1], c('mean', 'median'))

    if(!is.function(average_fun)) {

      stop("'average_fun' has to be a function.", call. = FALSE)

    }

    stopifnot(is.logical(show_whiskers))

    stopifnot(is.numeric(probs))

    if(length(probs) < 2) {

      stop('At least two probs are required.', call. = FALSE)

    }

    probs <- probs[1:2]

    if(any(probs < 0) | any(probs > 1)) {

      stop('Probs have to be within the [0, 1] range.', call. = FALSE)

    }

    stopifnot(is.numeric(dodge))

    stopifnot(is.numeric(point_size))
    stopifnot(is.numeric(point_alpha))
    stopifnot(is.numeric(line_alpha))

    stopifnot(is.logical(show_txt))
    stopifnot(is.numeric(txt_size))

    txt_position <- match.arg(txt_position[1], c('left', 'right'))

    if(!is.function(labeller)) {

      stop("'labeller' has to be a function.", call. = FALSE)

    }

    stopifnot(inherits(cust_theme, 'theme'))

    ## plotting data -------

    if(is.null(plot_subtitles)) {

      plot_subtitles <- paste('data sets: n =', length(data))

    }

    data_names <- names(data)

    data <- map(data, ~.x[c(split_fct, variables)])

    if(normalize) {

      center_fun <- switch(norm_center,
                           mean = function(x) mean(x, na.rm = TRUE),
                           median = function(x) median(x, na.rm = TRUE))

      for(i in seq_along(data)) {

        data[[i]][variables] <-
          map_dfc(data[[i]][variables],
                  ~scale(.x, center = center_fun(.x))[, 1])

      }

    }

    split_vec <- map(data, ~.x[[split_fct]])

    levs <- map(split_vec, levels)

    cmm_levs <- levels(droplevels(reduce(split_vec, c)))

    data_splits <-
      map2(data, split_vec,
           ~split(.x[variables], .y, drop = TRUE))

    averages <- map(data_splits, map, average_fun)
    errors <- map(data_splits, map, colQuantiles, probs)

    average_val <- NULL
    lower_quant <- NULL
    upper_quant <- NULL
    data_set <- NULL
    variable <- NULL
    shift_x <- NULL

    averages <- map(averages,
                    function(sp) map2_dfr(names(sp), sp,
                                          ~tibble(!!split_fct := .x,
                                                  variable = names(.y),
                                                  average_val = .y)))

    errors <- map(errors, map, as.data.frame)

    errors <- map(errors, map, set_names, c('lower_quant', 'upper_quant'))

    errors <- map(errors, map, rownames_to_column, 'variable')

    errors <- map(errors,
                  function(sp) map2_dfr(sp, names(sp),
                                        ~mutate(.x,
                                                !!split_fct := .y)))

    plot_tbl <-
      map2(averages, errors,
           left_join, by = c('variable', split_fct))

    plot_tbl <- map2_dfr(plot_tbl, names(plot_tbl),
                         ~mutate(.x, data_set = .y))

    ## overlap prevention by random dodging

    if(dodge != 0) {

      norm_dodge <- rnorm(n = nrow(plot_tbl), mean = 0, sd = dodge)

    } else {

      norm_dodge <- 0

    }

    plot_tbl <- mutate(plot_tbl,
                       variable = factor(variable, variables),
                       !!split_fct := factor(.data[[split_fct]], cmm_levs),
                       data_set = factor(data_set, data_names),
                       x_numeric = as.numeric(.data[[split_fct]]),
                       shift_x = norm_dodge)

    ## X scale breaks and labels: we're plotting on a numeric scale for
    ## possible overlap prevention, buth displaying the discrete values.

    x_breaks <- sort(unique(plot_tbl$x_numeric))
    x_labels <- levels(plot_tbl[[split_fct]])

    plot_tbl <- split(plot_tbl, plot_tbl[['variable']])

    ## plotting --------

    ch_plots <- map(plot_tbl,
                    ~ggplot(.x,
                            aes(x = x_numeric + shift_x,
                                y = average_val,
                                color = .data[[split_fct]],
                                fill = .data[[split_fct]])) +
                      geom_line(aes(group = data_set),
                                linetype = linetype,
                                alpha = line_alpha,
                                color = line_color))

    if(show_whiskers) {

      ch_plots <-
        map(ch_plots,
            ~.x +
              geom_errorbar(aes(ymin = lower_quant,
                                ymax = upper_quant),
                            linewidth = whisker_size,
                            width = 0,
                            alpha = point_alpha))

    }

    ch_plots <- map(ch_plots,
                    ~.x +
                      geom_point(shape = 21,
                                 size = point_size,
                                 alpha = point_alpha,
                                 color = 'black') +
                      scale_x_continuous(breaks = x_breaks,
                                         labels = x_labels))

    ch_plots <-
      pmap(list(x = ch_plots,
                y = plot_titles,
                z = plot_subtitles),
           function(x, y, z) x +
             labs(title = y,
                  subtitle = z,
                  x = x_lab,
                  y = y_lab) +
             cust_theme)

    if(!is.null(palette)) {

      ch_plots <- map(ch_plots,
                      ~.x +
                        scale_fill_manual(values = palette,
                                          name = x_lab) +
                        scale_color_manual(values = palette,
                                           name = x_lab))

    }

    if(!show_txt) return(ch_plots)

    ## data set labels --------

    if(txt_position == 'left') {

      txt_data <- map(plot_tbl,
                      filter,
                      .data[[split_fct]] == cmm_levs[1])

    } else {

      txt_data <- map(plot_tbl,
                      filter,
                      .data[[split_fct]] == cmm_levs[length(cmm_levs)])

    }

    if(is.null(txt_color)) {

      ch_plots <-
        map2(ch_plots,
             txt_data,
             ~.x +
               geom_text_repel(data = .y,
                               aes(label = labeller(data_set)),
                               size = txt_size,
                               fontface = txt_face,
                               show.legend = FALSE))

    } else {

      ch_plots <-
        map2(ch_plots,
             txt_data,
             ~.x +
               geom_text_repel(data = .y,
                               aes(label = labeller(data_set)),
                               size = txt_size,
                               color = txt_color,
                               fontface = txt_face,
                               show.legend = FALSE))

    }

    ch_plots

  }


# Scatter plots of estimates --------

#' Scatter/Forest plot of estimates in a multi-data set analysis.
#'
#' @description
#' Draws a scatter plot of estimates specified by the user. Point color codes
#' for significance and regulation status. As an option, a central statistic
#' over all estimates such as mean or median with a quantile range may be
#' displayed.
#'
#' @return a `ggplot` object.
#'
#' @param data a list of data frames with analysis results.
#' @param label_variable name of the variable storing feature names.
#' @param regulation_variable name of the variable storing the estimate to
#' be plotted.
#' @param regulation_status name of the variable storing regulation status index
#' (e.g. activated/inhibited/ns or regulated/ns) in factor form. This variable
#' levels will be represented by color of the data points.
#' @param palette colors coding for the regulation status. Most common
#' options are provided as a default.
#' @param show_average logical, should averages over all data srts be displayed?
#' @param average_fun function used for calculation of the average
#' regulation estimates over the data sets.
#' @param show_whiskers logical, should whiskers representing a quantile range
#' be displayed in the plot? Defaults to false.
#' @param probs a numeric vector of two quantiles defining the range for the
#' whiskers. Ignored if `show_whiskers = FALSE`.
#' @param point_size size of the data points.
#' @param point_alpha alpha of the data points.
#' @param point_wjitter data point jittering width.
#' @param point_hjitter data point jittering height.
#' @param average_shape code of shape of the average symbol.
#' @param average_size size of the symbol representing the central
#' tendency measure.
#' @param stat_color color of the symbols of the central tendency measure and
#' the error bars.
#' @param stat_rim color of the rim of the average statistic symbol.
#' @param whisker_size size of the percentile whiskers.
#' @param plot_title plot title.
#' @param plot_subtitle plot_subtitle.
#' @param x_lab X axis title.
#' @param y_lab Y axis title.
#' @param status_lab fill and color scale title.
#' @param plot_tag plot tag.
#' @param cust_theme custom `ggplot` theme.
#'
#' @export

  plot_common_estimates <- function(data,
                                    label_variable,
                                    regulation_variable,
                                    regulation_status,
                                    palette = c(upregulated = 'firebrick',
                                                activated = 'firebrick',
                                                positive = 'firebrick',
                                                downregulated = 'steelblue',
                                                inhibited = 'steelblue',
                                                negative = 'steelblue',
                                                regulated = 'plum4',
                                                significant = 'plum4',
                                                ns = 'gray80'),
                                    show_average = TRUE,
                                    average_fun = colMeans,
                                    show_whiskers = TRUE,
                                    whisker_size = 0.75,
                                    probs = c(0.025, 0.975),
                                    point_size = 2,
                                    point_alpha = 0.75,
                                    point_wjitter = 0,
                                    point_hjitter = 0.1,
                                    average_shape = 23,
                                    average_size = 3,
                                    stat_color = 'orangered3',
                                    stat_rim = stat_color,
                                    plot_title = NULL,
                                    plot_subtitle = NULL,
                                    x_lab = regulation_variable,
                                    y_lab = NULL,
                                    status_lab = regulation_status,
                                    plot_tag = NULL,
                                    cust_theme = microViz::theme_micro()) {

    ## input control for the data argument -------

    stopifnot(is.list(data))

    classes <- map_lgl(data, is.data.frame)

    if(any(!classes)) {

      stop("At least one element of 'data' is not a data frame.",
           call. = FALSE)

    }

    vars_to_check <- c(regulation_status, label_variable, regulation_variable)

    var_present <- map_lgl(data, ~all(vars_to_check %in% names(.x)))

    if(any(!var_present)) {

      stop(paste("At least one of 'regulation_status', 'label_variable', or",
                 "'regulation_variable' is missing."),
           call. = FALSE)

    }

    regulation_format <- map_lgl(data, ~is.numeric(.x[[regulation_variable]]))

    if(any(!regulation_format)) {

      stop("'regulation_variable' has to be numeric.", call. = FALSE)

    }

    index_format <- map_lgl(data, ~is.factor(.x[[regulation_status]]))

    if(any(!index_format)) {

      stop("'regulation_status' has to be a factor.", call. = FALSE)

    }

    ## input control for the remaining arguments -------

    stopifnot(is.logical(show_average))

    if(!is.function(average_fun)) {

      stop("'average_fun' has to be a function.", call. = FALSE)

    }

    stopifnot(is.logical(show_whiskers))
    stopifnot(is.numeric(probs))

    if(length(probs) < 2) {

      stop("Two 'probs' are required.", call. = FALSE)

    }
    probs <- probs[1:2]

    if(any(probs < 0) | any(probs > 1)) {

      stop("'probs' must be within the [0, 1] range.", call. = FALSE)

    }

    stopifnot(is.numeric(point_size))
    stopifnot(is.numeric(point_alpha))
    stopifnot(is.numeric(point_wjitter))
    stopifnot(is.numeric(point_hjitter))

    stopifnot(is.numeric(average_size))
    stopifnot(is.numeric(whisker_size))

    stopifnot(inherits(cust_theme, 'theme'))

    ## plotting data -------

    data <-
      map(data, ~.x[c(label_variable, regulation_status, regulation_variable)])

    data_names <- names(data)

    data_set <- NULL

    plot_tbl <- map2_dfr(data, names(data),
                         ~mutate(.x, data_set = .y))

    if(is.null(plot_subtitle)) {

      plot_subtitle <- paste('data sets: n =', length(data))

    }

    ## statistics to be plotted ---------

    estimate_tbl <- map(data, ~.x[c(label_variable, regulation_variable)])

    estimate_tbl <- reduce(estimate_tbl, full_join, by = label_variable)

    estimate_tbl <- t(column_to_rownames(estimate_tbl, label_variable))

    averages <- average_fun(estimate_tbl)
    errors <- as.data.frame(colQuantiles(estimate_tbl, probs = probs))

    errors <- set_names(errors, c('lower_quant', 'upper_quant'))

    errors <- rownames_to_column(errors, label_variable)

    stat_tbl <- mutate(errors, !!regulation_variable := averages)

    ## plotting ---------

    est_plot <-
      ggplot(plot_tbl,
             aes(x = .data[[regulation_variable]],
                 y = reorder(.data[[label_variable]],
                             .data[[regulation_variable]]))) +
      geom_point(aes(fill = .data[[regulation_status]]),
                 shape = 21,
                 size = point_size,
                 alpha = point_alpha,
                 position = position_jitter(width = point_wjitter,
                                            height = point_hjitter)) +
      scale_fill_manual(values = palette,
                        name = status_lab) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    if(show_whiskers) {

      lower_quant <- NULL
      upper_quant <- NULL

      est_plot <- est_plot +
        geom_errorbarh(data = stat_tbl,
                       aes(xmin = lower_quant,
                           xmax = upper_quant),
                       height = 0,
                       linewidth = whisker_size,
                       color = stat_color)

    }

    if(show_average) {

      est_plot <- est_plot +
        geom_point(data = stat_tbl,
                   shape = average_shape,
                   size = average_size,
                   fill = stat_color,
                   color = stat_rim)

    }

    est_plot

  }

# Counts of significant effects and average effect sizes --------

#' Counts of significant effects and average effect sizes in multi-data set
#' analysis.
#'
#' @description
#' Draws a scatter plot of average effect sizes (or regulation; X axis) and
#' counts of data sets with significant regulation (Y axis).
#'
#' @return a `ggplot` object.
#'
#' @param data a list of data frames with analysis results.
#' @param label_variable name of the variable storing feature names.
#' @param regulation_variable name of the variable storing the estimate to
#' be plotted.
#' @param regulation_status name of the variable storing regulation status index
#' (e.g. activated/inhibited/ns or regulated/ns) in factor form. This variable
#' levels will be represented by color of the data points.
#' @param ns_level level of the `regulation_status` variable that represent
#' non-significant regulation. `ns` by default.
#' @param average_fun function used for calculation of the average
#' regulation estimates over the data sets.
#' @param cutoff cutoff value of the count of data sets with significant
#' effects. Useful for identification of common significant effects shared by
#' multiple data sets.
#' @param palette colors coding for the regulation status. Most common
#' options are provided as a default.
#' @param average_fun function used for calculation of the average
#' regulation estimates over the data sets.
#' @param show_whiskers logical, should whiskers representing a quantile range
#' be displayed in the plot? Defaults to false.
#' @param probs a numeric vector of two quantiles defining the range for the
#' whiskers. Ignored if `show_whiskers = FALSE`.
#' @param point_size size of the data points.
#' @param point_alpha alpha of the data points.
#' @param point_wjitter data point jittering width.
#' @param point_hjitter data point jittering height.
#' @param whisker_size size of the percentile whiskers.
#' @param show_txt logical, should labels with names of the data sets be
#' displayed?
#' @param txt_size size of the label text.
#' @param txt_color color of the text. If `NULL`, the color is specified by
#' levels of the splitting factor.
#' @param txt_face font face of the text labels.
#' @param labeller a function transforming the data set names to the text labels
#' in the plot. Defaults to `identity`, i.e. the names are displayed as the are.
#' @param plot_title plot title.
#' @param plot_subtitle plot_subtitle.
#' @param x_lab X axis title.
#' @param y_lab Y axis title.
#' @param status_lab fill and color scale title.
#' @param plot_tag plot tag.
#' @param cust_theme custom `ggplot` theme.
#'
#' @export

  plot_common_effect <- function(data,
                                 label_variable,
                                 regulation_variable,
                                 regulation_status,
                                 ns_level = 'ns',
                                 average_fun = colMeans,
                                 cutoff = 0,
                                 palette = c(upregulated = 'firebrick',
                                             activated = 'firebrick',
                                             positive = 'firebrick',
                                             downregulated = 'steelblue',
                                             inhibited = 'steelblue',
                                             negative = 'steelblue',
                                             regulated = 'plum4',
                                             significant = 'plum4',
                                             ns = 'gray70'),
                                 show_whiskers = TRUE,
                                 whisker_size = 0.75,
                                 probs = c(0.025, 0.975),
                                 point_size = 2,
                                 point_alpha = 0.75,
                                 point_wjitter = 0,
                                 point_hjitter = 0,
                                 show_txt = TRUE,
                                 txt_size = 2.75,
                                 txt_color = NULL,
                                 txt_face = 'plain',
                                 labeller = identity,
                                 plot_title = NULL,
                                 plot_subtitle = NULL,
                                 x_lab = regulation_variable,
                                 y_lab = NULL,
                                 status_lab = regulation_status,
                                 plot_tag = NULL,
                                 cust_theme = microViz::theme_micro()) {

    ## input control for the data argument -------

    stopifnot(is.list(data))

    classes <- map_lgl(data, is.data.frame)

    if(any(!classes)) {

      stop("At least one element of 'data' is not a data frame.",
           call. = FALSE)

    }

    vars_to_check <- c(regulation_status, label_variable, regulation_variable)

    var_present <- map_lgl(data, ~all(vars_to_check %in% names(.x)))

    if(any(!var_present)) {

      stop(paste("At least one of 'regulation_status', 'label_variable', or",
                 "'regulation_variable' is missing."),
           call. = FALSE)

    }

    regulation_format <- map_lgl(data, ~is.numeric(.x[[regulation_variable]]))

    if(any(!regulation_format)) {

      stop("'regulation_variable' has to be numeric.", call. = FALSE)

    }

    index_format <- map_lgl(data, ~is.factor(.x[[regulation_status]]))

    if(any(!index_format)) {

      stop("'regulation_status' has to be a factor.", call. = FALSE)

    }

    n_levels <- map_dbl(data, ~length(levels(.x[[regulation_status]])))

    if(any(n_levels < 2)) {

      stop("'regulation_status' hast to have at least two levels.",
           call. = FALSE)

    }

    ns_present <-
      map_lgl(data, ~ns_level %in% levels(.x[[regulation_status]]))

    if(any(!ns_present)) {

      stop("'ns_level' not detected in the regulation status variable.",
           call. = TRUE)

    }

    ## input control for the remaining arguments ---------

    stopifnot(is.numeric(cutoff))

    if(!is.function(average_fun)) {

      stop("'average_fun' has to be a function.", call. = FALSE)

    }

    stopifnot(is.logical(show_whiskers))
    stopifnot(is.numeric(whisker_size))

    stopifnot(is.logical(show_whiskers))
    stopifnot(is.numeric(probs))

    if(length(probs) < 2) {

      stop("Two 'probs' are required.", call. = FALSE)

    }

    probs <- probs[1:2]

    if(any(probs < 0) | any(probs > 1)) {

      stop("'probs' must be within the [0, 1] range.", call. = FALSE)

    }

    stopifnot(is.numeric(point_size))
    stopifnot(is.numeric(point_alpha))
    stopifnot(is.numeric(point_wjitter))
    stopifnot(is.numeric(point_hjitter))

    stopifnot(inherits(cust_theme, 'theme'))

    ## plotting data: average estimates and their quantile ranges --------

    estimate_tbl <- map(data, ~.x[c(label_variable, regulation_variable)])

    estimate_tbl <- reduce(estimate_tbl, full_join, by = label_variable)

    estimate_tbl <- t(column_to_rownames(estimate_tbl, label_variable))

    averages <- average_fun(estimate_tbl)
    errors <- as.data.frame(colQuantiles(estimate_tbl, probs = probs))

    errors <- set_names(errors, c('lower_quant', 'upper_quant'))

    errors <- rownames_to_column(errors, label_variable)

    stat_tbl <- mutate(errors, !!regulation_variable := averages)

    ## numbers of significant effects --------

    ## if there are multiple modes of regulation, e.g.
    ## activation and inhibition of gene expression, the effects should be
    ## counted separately for each regulation mode

    all_variables <- map(data, ~.x[label_variable])

    all_variables <- unique(unname(unlist(all_variables)))

    count_tbl <- map_dfr(data,
                         filter,
                         !duplicated(.data[[label_variable]]))

    count_tbl <-
      split(count_tbl[[label_variable]], count_tbl[[regulation_status]])

    count_tbl <- map(count_tbl, table)

    n <- NULL

    count_tbl <- map(count_tbl,
                     ~tibble(!!label_variable := names(.x),
                             n = as.numeric(.x)))

    count_tbl <-
      map(count_tbl,
          function(x) if(nrow(x) == 0) {

            tibble(!!label_variable := all_variables,
                   n = 0)

            } else x)


    count_tbl <- count_tbl[names(count_tbl) != ns_level]

    count_tbl <- map2_dfr(count_tbl, names(count_tbl),
                          ~mutate(.x,
                                  !!regulation_status := factor(.y, .y)))

    .reg_index <- NULL

    count_tbl <-
      mutate(count_tbl,
             .reg_index = ifelse(n >= cutoff,
                                 as.character(.data[[regulation_status]]),
                                 ns_level),
             .reg_index = factor(.reg_index,
                                 c(levels(.data[[regulation_status]]),
                                   ns_level)))

    count_tbl <- count_tbl[c(label_variable, 'n', '.reg_index')]

    ## handling special cases: no significant effects for a variable

    count_tbl <- split(count_tbl, count_tbl[[label_variable]])

    ns_handler <- function(x) {

      if(all(x$.reg_index == ns_level)) {

        if(all(x$n == 0)) {

          return(x[1, ])

        } else {

          return(filter(x, n != 0))

        }

      } else {

        return(filter(x, .reg_index != 'ns'))

      }

    }

    count_tbl <-
      map_dfr(count_tbl, ns_handler)

    ## plotting data frame an management of possible overlaps -------

    plot_tbl <- left_join(count_tbl, stat_tbl, by = label_variable)

    if(point_wjitter != 0) {

      x_dodge <- rnorm(nrow(plot_tbl), mean = 0, sd = point_wjitter)

    } else {

      x_dodge <- 0

    }

    if(point_hjitter != 0) {

      y_dodge <- rnorm(nrow(plot_tbl), mean = 0, sd = point_hjitter)

    } else {

      y_dodge <- 0

    }

    x_shift <- NULL
    y_shift <- NULL
    txt_label <- NULL

    plot_tbl <- mutate(plot_tbl,
                       x_shift = x_dodge,
                       y_shift = y_dodge,
                       txt_label = ifelse(.reg_index == ns_level,
                                          NA, .data[[label_variable]]))

    ## plotting -----

    eff_plot <-
      ggplot(plot_tbl,
             aes(x = .data[[regulation_variable]] + x_shift,
                 y = n + y_shift,
                 fill = .reg_index))

    if(show_whiskers) {

      lower_quant <- NULL
      upper_quant <- NULL

      eff_plot <- eff_plot +
        geom_errorbarh(aes(xmin = lower_quant + x_shift,
                           xmax = upper_quant + x_shift,
                           color = .reg_index),
                       height = 0,
                       linewidth = whisker_size,
                       alpha = point_alpha)

    }

    if(cutoff != 0) {

      eff_plot <- eff_plot +
        geom_hline(yintercept = cutoff,
                   linetype = 'dashed')

    }

    eff_plot <- eff_plot +
      geom_point(shape = 21,
                 size = point_size,
                 alpha = point_alpha) +
      scale_fill_manual(values = palette,
                        name = status_lab) +
      scale_color_manual(values = palette,
                         name = status_lab) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    if(!show_txt) return(eff_plot)

    ## labels for selected effects ------

    txt_data <- filter(plot_tbl, !is.na(txt_label))

    if(!is.null(txt_color)) {

      eff_plot <- eff_plot +
        geom_text_repel(aes(label = labeller(txt_label)),
                        size = txt_size,
                        fontface = txt_face,
                        color = 'black')

    } else {

      eff_plot <- eff_plot +
        geom_text_repel(aes(label = labeller(txt_label),
                            color = .reg_index),
                        size = txt_size,
                        fontface = txt_face)

    }

    eff_plot

  }

# END ------
