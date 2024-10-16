# Plotting functions.

#' @include imports.R

  NULL

# Volcano plot -----

#' Draw a volcano plot.
#'
#' @description
#' Generates a Volcano plot with the effect size/regulation
#' statistic value on the X axis and -log10 p value on the y axis. Regulation
#' and significance is coded by the point fill color. Top n significant
#' genes/objects may be labeled with their names.
#'
#' @param data a data frame.
#' @param regulation_variable name of the variable storing the regulation/effect
#' size data.
#' @param p_variable name of the variable storing p values.
#' @param signif_level significance threshold, p = 0.05 be default.
#' @param regulation_level the regulation cutoff for identification of the
#' significant effects, 1.5 difference versus control by default.
#' @param fill_scale color palette for coding significance/regulation.
#' Three colors required: for upregulation, downregulation and non-significance.
#' @param x_lab x axis title.
#' @param y_lab y axis title.
#' @param top_significant top n significant genes/objects to be labeled in
#' the plot with their names.
#' @param top_regulated top n most regulated genes/objects to be labeled in
#' the plot with their names. Ignored if `top_significant` is provided.
#' @param label_variable variable storing the gene/object names.
#' @param label_type type of the gene/object label: 'label' (default)
#' or plain 'text'.
#' @param txt_size size of the gene/object name text.
#' @param txt_color gene/object name text color.
#' @param txt_face face of the gene/object name text.
#' @param fill_title title of the point fill legend.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag text. If NULL, the numbers of significantly
#' regulated genes is displayed.
#' @param cust_theme custom ggplot theme.
#' @param point_alpha plot point alpha.
#' @param point_hjitter height jittering of the points, defaults to 0.
#' @param point_wjitter width jittering of the points, defaults to 0.
#' @param show_hline logical, should the horizontal line representing the
#' significance cutoff be presented in the plot?
#' @param show_vlines specifies which vertical lines representing the regulation
#' cutoffs will be plotted.
#'
#' @return a `ggplot` object.
#'
#' @export

  plot_volcano <- function(data,
                           regulation_variable,
                           p_variable,
                           signif_level = 0.05,
                           regulation_level = 1.5,
                           fill_scale = c(upregulated = 'firebrick',
                                          downregulated = 'steelblue',
                                          ns = 'gray60'),
                           x_lab = regulation_variable,
                           y_lab = expression('-log'[10]*' p'),
                           top_significant = 0,
                           top_regulated = 0,
                           label_variable = NULL,
                           label_type = c('label', 'text'),
                           txt_size = 2.75,
                           txt_color = 'black',
                           txt_face = 1,
                           fill_title = '',
                           plot_title = NULL,
                           plot_subtitle = NULL,
                           plot_tag = NULL,
                           cust_theme = microViz::theme_micro(),
                           point_alpha = 0.8,
                           point_hjitter = 0,
                           point_wjitter = 0,
                           show_hline = TRUE,
                           show_vlines = c('both', 'right', 'left')) {

    ## entry control ------

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(regulation_variable, p_variable) %in% names(data))) {

      stop('Regulation_variable or p_variable absent from the data.',
           call. = FALSE)

    }

    top_significant <- as.integer(top_significant)
    top_regulated <- as.integer(top_regulated)

    if(!is.null(label_variable)) {

      if(!label_variable %in% names(data)) {

        stop('label_variable missing from data.', call. = FALSE)

      }

    }

    label_type <- match.arg(label_type[1], c('label', 'text'))

    if(!inherits(cust_theme, 'theme')) {

      stop('Please provide a valid ggplot theme object.', call. = FALSE)

    }

    significant <- NULL
    regulation <- NULL

    stopifnot(is.logical(show_hline))

    show_vlines <- match.arg(show_vlines[1],
                             c('both', 'right', 'left'))

    ## plotting data ------

    plot_tbl <-
      filter(data,
             complete.cases(data[c(regulation_variable,
                                   p_variable)]))

    plot_tbl <-
      mutate(plot_tbl,
             significant = ifelse(.data[[p_variable]] < signif_level,
                                  'yes', 'no'),
             regulation = ifelse(significant == 'no',
                                 'ns',
                                 ifelse(.data[[regulation_variable]] > regulation_level,
                                        'upregulated',
                                        ifelse(.data[[regulation_variable]] < -regulation_level,
                                               'downregulated', 'ns'))),
             regulation = factor(regulation, c('upregulated',
                                               'downregulated',
                                               'ns')))

    ## numbers of regulated genes -------

    if(is.null(plot_tag)) {

      n_genes <- count(plot_tbl, regulation, .drop = FALSE)

      plot_tag <- paste0('upregulated: n = ', n_genes$n[1],
                         ', downregulated: n = ', n_genes$n[2])

    }

    ## plotting ------

    volc <- ggplot(plot_tbl,
                   aes(x = .data[[regulation_variable]],
                       y = -log10(.data[[p_variable]]),
                       fill = regulation)) +
      geom_point(size = 2,
                 shape = 21,
                 alpha = point_alpha,
                 position = ggplot2::position_jitter(width = point_wjitter,
                                                     height = point_hjitter)) +
      scale_fill_manual(values = fill_scale,
                        name = fill_title) +
      cust_theme +
      labs(x = x_lab,
           y = y_lab,
           title = plot_title,
           subtitle = plot_subtitle,
           tag = plot_tag)

    if(show_hline) {

      volc <- volc +
        geom_hline(yintercept = -log10(signif_level), lty = 2)

    }

    if(regulation_level > 0) {

      if(show_vlines == 'both') {

        volc <- volc  +
          geom_vline(xintercept = -regulation_level, lty = 2) +
          geom_vline(xintercept = regulation_level, lty = 2)

      } else if(show_vlines == 'left') {

        volc <- volc  +
          geom_vline(xintercept = -regulation_level, lty = 2)

      } else {

        volc <- volc  +
          geom_vline(xintercept = regulation_level, lty = 2)

      }

    } else if(regulation_level == 0) {

      volc <- volc +
        geom_vline(xintercept = 0, lty = 2)

    }

    if(!is.null(label_variable)) {

      desc_tbl <- map(c('upregulated', 'downregulated'),
                      ~filter(plot_tbl, regulation == .x))

      if(top_significant > 0) {

        desc_tbl <-
          map_dfr(desc_tbl,
                  ~top_n(.x,
                         n = top_significant,
                         -.data[[p_variable]]))

      } else if(top_regulated > 0) {

        desc_tbl <-
          map_dfr(desc_tbl,
                  ~top_n(.x,
                         n = top_regulated,
                         abs(.data[[regulation_variable]])))

      } else {

        return(volc)

      }

      if(nrow(desc_tbl) == 0) return(volc)

      if(label_type == 'label') {

        volc <- volc +
          geom_label_repel(data = desc_tbl,
                           aes(label = .data[[label_variable]]),
                           size = txt_size,
                           color = txt_color,
                           label.padding = 0.1,
                           box.padding = 0.1,
                           show.legend = FALSE,
                           fontface = txt_face)

      } else {

        volc <- volc +
          geom_text_repel(data = desc_tbl,
                          aes(label = .data[[label_variable]]),
                          size = txt_size,
                          color = txt_color,
                          box.padding = 0.1,
                          show.legend = FALSE,
                          fontface = txt_face)

      }

    }

    return(volc)

  }

# Cascade plot -----

#' Regulation bar plot.
#'
#' @description Draws a bar plot with the regulation/effect size values.
#' Bar color codes for the regulation sign and significance.
#' @param data a data frame.
#' @param regulation_variable name of the variable storing the regulation/effect
#' size data.
#' @param p_variable name of the variable storing p values.
#' @param signif_level significance threshold, p = 0.05 be default.
#' @param regulation_level the regulation cutoff for identification of the
#' significant effects, 1.5 difference versus control by default.
#' @param fill_scale color palette for coding significance/regulation.
#' Three colors required: for upregulation, downregulation and non-significance.
#' @param fill_title fill legend title.
#' @param x_lab x axis title.
#' @param y_lab y axis title.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag text. If NULL, the numbers of significantly
#' regulated genes is displayed.
#' @param cust_theme custom ggplot theme.
#' @param bar_alpha plot bar alpha.
#' @param show_trend logical, should a LOESS trend for
#' regulation be displayed in the plot?
#' @param ... extra arguments passed to \code{\link[ggplot2]{geom_smooth}}.
#' @return a ggplot object.
#' @export

  plot_sign <- function(data,
                        regulation_variable,
                        p_variable,
                        signif_level = 0.05,
                        regulation_level = 1.5,
                        fill_scale = c(upregulated = 'firebrick',
                                       downregulated = 'steelblue',
                                       ns = 'gray60'),
                        fill_title = '',
                        bar_alpha = 0.75,
                        plot_title = NULL,
                        plot_subtitle = NULL,
                        plot_tag = NULL,
                        x_lab = 'Gene',
                        y_lab = 'Regulation',
                        show_trend = TRUE,
                        cust_theme = microViz::theme_micro(), ...) {

    ## entry control --------

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(regulation_variable, p_variable) %in% names(data))) {

      stop('Regulation_variable or p_variable absent from the data.',
           call. = FALSE)

    }

    if(!any(class(cust_theme) != 'theme')) {

      stop('Please provide a valid ggplot theme object.', call. = FALSE)

    }

    stopifnot(is.logical(show_trend))

    significant <- NULL
    regulation <- NULL
    plot_order <- NULL

    ## plotting data ---------

    plot_tbl <-
      filter(data,
             complete.cases(data[c(p_variable,
                                   regulation_variable)]))

    plot_tbl <-
      mutate(plot_tbl,
             significant = ifelse(.data[[p_variable]] < signif_level,
                                  'yes', 'no'),
             regulation = ifelse(significant == 'no',
                                 'ns',
                                 ifelse(.data[[regulation_variable]] > regulation_level,
                                        'upregulated',
                                        ifelse(.data[[regulation_variable]] < -regulation_level,
                                               'downregulated', 'ns'))),
             regulation = factor(regulation, c('upregulated',
                                               'downregulated',
                                               'ns')))

    plot_tbl <- arrange(plot_tbl, -.data[[regulation_variable]])

    plot_tbl <- mutate(plot_tbl, plot_order = 1:nrow(plot_tbl))

    ## numbers of regulated items

    if(is.null(plot_tag)) {

      n_genes <- count(plot_tbl, regulation)

      plot_tag <- paste0('upregulated: n = ', n_genes$n[1],
                         ', downregulated: n = ', n_genes$n[2])

    }

    ## plot --------

    bar <- ggplot(plot_tbl,
                  aes(x = plot_order,
                      y = .data[[regulation_variable]])) +
      geom_bar(stat = 'identity',
                        alpha = bar_alpha,
                        aes(fill = regulation)) +
      scale_fill_manual(values = fill_scale,
                        name = fill_title) +
      cust_theme +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           tag = plot_tag,
           x = x_lab,
           y = y_lab)

    if(show_trend) {

      bar <- bar +
        ggplot2::geom_smooth(show.legend = FALSE, ...)

    }

    return(bar)

  }

# General Forest plot ------

#' Draw regulation statistics for selected genes.
#'
#' @description
#' Draws a Forest plot with the regulation estimates and,
#' optionally confidence intervals for all items presented in a data frame
#' (`plot_forest`) or top up- and downregulated fratures (`plot_top`).
#'
#' @param data a data frame.
#' @param regulation_variable name of the variable storing the regulation/effect
#' size data.
#' @param label_variable variable storing the gene/object names.
#' @param p_variable variable storing the p values.
#' @param signif_level significance threshold, p = 0.05 be default.
#' @param regulation_level the regulation cutoff for identification of the
#' significant effects, 0 difference versus control by default.
#' @param lower_ci_variable variable storing the lower CI values.
#' @param upper_ci_variable variable storing the upper CI values.
#' @param top_regulated top n regulated genes/objects to be presented in the plot.
#' @param fill_scale regulation colors, a vector of three elements: for
#' upregulated, downregulated and non-significant items, respectively.
#' @param x_lab x axis title.
#' @param fill_title title of the point fill legend.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag text.
#' @param cust_theme custom ggplot theme.
#' @param signif_digits significant digits used for rounding of the
#' displayed numeric variables.
#' @param show_txt logical, label the points with the estimate value?
#' @param show_ci_txt logical, should the CI be included in the text label?
#' @param txt_size text size.
#' @param txt_hjust horizontal justification of the text label.
#' @param txt_vjust vertical justification of the text label.
#'
#' @return a ggplot object.
#'
#' @export

  plot_forest <- function(data,
                          regulation_variable,
                          label_variable,
                          p_variable,
                          signif_level = 0.05,
                          regulation_level = 0,
                          lower_ci_variable = NULL,
                          upper_ci_variable = NULL,
                          fill_scale = c(upregulated = 'firebrick',
                                         downregulated = 'steelblue',
                                         ns = 'gray60'),
                          fill_title = '',
                          plot_title = NULL,
                          plot_subtitle = NULL,
                          plot_tag = NULL,
                          x_lab = 'Regulation',
                          cust_theme = microViz::theme_micro(),
                          signif_digits = 2,
                          show_txt = FALSE,
                          show_ci_txt = FALSE,
                          txt_size = 2.75,
                          txt_hjust = 0.5,
                          txt_vjust = -0.8) {

    ## entry control -------

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(regulation_variable,
              label_variable,
              p_variable) %in% names(data))) {

      stop('Regulation_variable or label_variable absent from the data.',
           call. = FALSE)

    }

    if(!any(class(cust_theme) != 'theme')) {

      stop('Please provide a valid ggplot theme object.', call. = FALSE)

    }

    stopifnot(is.logical(show_txt))
    stopifnot(is.logical(show_ci_txt))
    stopifnot(is.numeric(signif_level))
    stopifnot(is.numeric(regulation_level))
    stopifnot(is.numeric(txt_size))
    stopifnot(is.numeric(txt_hjust))
    stopifnot(is.numeric(txt_vjust))

    if(!is.null(lower_ci_variable)) {

      if(!lower_ci_variable %in% names(data)) {

        stop('Lower and upper CI variables absent from the data.',
             call. = FALSE)

      }

      stopifnot(is.numeric(data[[lower_ci_variable]]))

    }

    if(!is.null(upper_ci_variable)) {

      if(!upper_ci_variable %in% names(data)) {

        stop('Lower and upper CI variables absent from the data.',
             call. = FALSE)

      }

      stopifnot(is.numeric(data[[upper_ci_variable]]))

    }

    stopifnot(is.numeric(data[[regulation_variable]]))
    stopifnot(is.numeric(data[[p_variable]]))

    significant <- NULL
    regulation <- NULL
    plot_lab <- NULL

    ## plotting table --------

    plot_tbl <-
      mutate(data,
             significant = ifelse(.data[[p_variable]] < signif_level,
                                  'yes', 'no'),
             regulation = ifelse(significant == 'no',
                                 'ns',
                                 ifelse(.data[[regulation_variable]] > regulation_level,
                                        'upregulated',
                                        ifelse(.data[[regulation_variable]] < -regulation_level,
                                               'downregulated', 'ns'))),
             regulation = factor(regulation, c('upregulated',
                                               'downregulated',
                                               'ns')))

    ## plotting -------

    forest <-
      ggplot(plot_tbl,
             aes(x = .data[[regulation_variable]],
                 y = reorder(.data[[label_variable]],
                             .data[[regulation_variable]]),
                 color = regulation,
                 fill = regulation)) +
      scale_fill_manual(values = fill_scale,
                        name = fill_title) +
      scale_color_manual(values = fill_scale,
                         name = fill_title) +
      cust_theme +
      theme(axis.title.y = element_blank()) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           tag = plot_tag,
           x = x_lab) +
      geom_vline(xintercept = 0,
                 linetype = 'dashed')

    if(!is.null(lower_ci_variable) & !is.null(upper_ci_variable)) {

      forest <- forest +
        geom_errorbarh(aes(xmin = .data[[lower_ci_variable]],
                           xmax = .data[[upper_ci_variable]]),
                       height = 0)

    }

    forest <- forest +
      geom_point(size = 2,
                          shape = 16)

    ## optional labeling

    if(!show_txt) {

      return(forest)

    }

    desc_tbl <- mutate(plot_tbl,
                       plot_lab = signif(.data[[regulation_variable]],
                                         signif_digits))

    if(show_ci_txt & all(!is.null(c(lower_ci_variable, upper_ci_variable)))) {

      desc_tbl <-
        mutate(desc_tbl,
               plot_lab = paste0(plot_lab,
                                 ' [', signif(.data[[lower_ci_variable]],
                                              signif_digits),
                                 ' - ', signif(.data[[upper_ci_variable]],
                                               signif_digits), ']'))

    }

    forest <- forest +
      geom_text(data = desc_tbl,
                aes(label = plot_lab),
                size = txt_size,
                hjust = txt_hjust,
                vjust = txt_vjust)

    return(forest)

  }

#' @rdname plot_forest
#' @export

  plot_top <- function(data,
                       regulation_variable,
                       label_variable,
                       p_variable,
                       signif_level = 0.05,
                       regulation_level = 0,
                       lower_ci_variable = NULL,
                       upper_ci_variable = NULL,
                       top_regulated = 10,
                       fill_scale = c(upregulated = 'firebrick',
                                      downregulated = 'steelblue',
                                      ns = 'gray60'),
                       fill_title = '',
                       plot_title = NULL,
                       plot_subtitle = NULL,
                       plot_tag = NULL,
                       x_lab = 'Regulation',
                       cust_theme = microViz::theme_micro(),
                       signif_digits = 2,
                       show_txt = FALSE,
                       show_ci_txt = FALSE,
                       txt_size = 2.75,
                       txt_hjust = 0.5,
                       txt_vjust = -0.8) {

    ## entry control -------

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(regulation_variable,
              label_variable,
              p_variable) %in% names(data))) {

      stop('Regulation, label or p variable absent from the data.',
           call. = FALSE)

    }

    if(!any(class(cust_theme) != 'theme')) {

      stop('Please provide a valid ggplot theme object.', call. = FALSE)

    }

    stopifnot(is.logical(show_txt))
    stopifnot(is.logical(show_ci_txt))
    stopifnot(is.numeric(signif_level))
    stopifnot(is.numeric(regulation_level))
    stopifnot(is.numeric(txt_size))
    stopifnot(is.numeric(txt_hjust))
    stopifnot(is.numeric(txt_vjust))

    if(!is.null(lower_ci_variable)) {

      if(!lower_ci_variable %in% names(data)) {

        stop('Lower and upper CI variables absent from the data.',
             call. = FALSE)

      }

      stopifnot(is.numeric(data[[lower_ci_variable]]))

    }

    if(!is.null(upper_ci_variable)) {

      if(!upper_ci_variable %in% names(data)) {

        stop('Lower and upper CI variables absent from the data.',
             call. = FALSE)

      }

      stopifnot(is.numeric(data[[upper_ci_variable]]))

    }

    stopifnot(is.numeric(data[[regulation_variable]]))
    stopifnot(is.numeric(data[[p_variable]]))

    reg_sign <- NULL

    ## plotting data ---------

    top_regulated <- as.integer(top_regulated)

    plot_tbl <-
      mutate(data,
             reg_sign = ifelse(.data[[regulation_variable]] > 0,
                               'up', 'down'),
             reg_sign = factor(reg_sign, c('up', 'down')))

    plot_tbl <- group_by(plot_tbl, reg_sign)

    plot_tbl <- top_n(plot_tbl,
                      n = top_regulated,
                      abs(.data[[regulation_variable]]))

    plot_tbl <- ungroup(plot_tbl)

    ## plotting ------

    plot_forest(data = plot_tbl,
                regulation_variable = regulation_variable,
                label_variable = label_variable,
                p_variable = p_variable,
                signif_level = signif_level,
                regulation_level = regulation_level,
                lower_ci_variable = lower_ci_variable,
                upper_ci_variable = upper_ci_variable,
                fill_scale = fill_scale,
                fill_title = fill_title,
                plot_title = plot_title,
                plot_subtitle = plot_subtitle,
                plot_tag = plot_tag,
                x_lab = x_lab,
                cust_theme = cust_theme,
                signif_digits = signif_digits,
                show_txt = show_txt,
                show_ci_txt = show_ci_txt,
                txt_size = txt_size,
                txt_hjust = txt_hjust,
                txt_vjust = txt_vjust)

  }

# P value and regulation plotting -------

#' Plot top p values or regulation estimates as bar plots.
#'
#' @description
#' Plots top n p values or regulation estimates as a bar plot.
#'
#' @param data a data frame.
#' @param p_variable variable storing the p values.
#' @param regulation_variable variable storing the regulation estimates.
#' @param label_variable variable storing the gene/object names.
#' @param signif_level significance threshold, p = 0.05 be default.
#' @param regulation_level regulation cutoff.
#' @param top_significant top n significant genes/objects to be presented in
#' the plot.
#' @param top_regulated top n strongest regulated genes/objects to be presented
#' in the plot.
#' @param fill_scale significance or regulation colors, a vector of two or
#' three elements.
#' @param x_lab x axis title.
#' @param fill_title title of the point fill legend.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag text.
#' @param cust_theme custom ggplot theme.
#' @param show_txt logical, should values of the regulation estimates be
#' presented in the plot?
#' @param signif_digits significant digits for regulation estimates presented
#' in the plot.
#' @param txt_color text color.
#' @param txt_size text size.
#' @param txt_offset horizontal offsetting of the text relative to the
#' regulation variable values.
#' @param txt_vjust vertical justification of the text.
#'
#' @return a `ggplot` object.
#'
#' @export

  plot_significant <- function(data,
                               p_variable,
                               label_variable,
                               signif_level = 0.05,
                               top_significant = 10,
                               fill_scale = c(significant = 'coral3',
                                              ns = 'gray60'),
                               x_lab = expression('-log'[10]*' p'),
                               fill_title = '',
                               plot_title = NULL,
                               plot_subtitle = NULL,
                               plot_tag = NULL,
                               cust_theme = microViz::theme_micro()) {

    ## entry control -------

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(label_variable,
              p_variable) %in% names(data))) {

      stop('p_variable or label_variable absent from the data.',
           call. = FALSE)

    }

    if(inherits(cust_theme, 'theme')) {

      stop('Please provide a valid ggplot theme object.', call. = FALSE)

    }

    top_significant <- as.integer(top_significant)

    significant <- NULL

    ## plotting table -------

    plot_tbl <- top_n(data, n = top_significant, -.data[[p_variable]])

    plot_tbl <- mutate(plot_tbl,
                       significant = ifelse(.data[[p_variable]] < signif_level,
                                            'significant',
                                            'ns'),
                       significant = factor(significant,
                                            c('significant', 'ns')))

    ## plotting ------

    ggplot(plot_tbl,
           aes(x = -log10(.data[[p_variable]]),
               y = reorder(.data[[label_variable]],
                           -.data[[p_variable]]),
               fill = significant)) +
      geom_bar(stat = 'identity',
               color = 'black') +
      geom_vline(xintercept = -log10(signif_level),
                 linetype = 'dashed') +
      scale_fill_manual(values = fill_scale,
                        name = fill_title) +
      cust_theme +
      theme(axis.title.y = element_blank()) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           tag = plot_tag,
           x = x_lab)

  }

#' @rdname plot_significant
#' @export

  plot_regulated <- function(data,
                             regulation_variable,
                             label_variable,
                             p_variable,
                             signif_level = 0.05,
                             regulation_level = 0,
                             top_regulated = 10,
                             fill_scale = c(upregulated = 'firebrick',
                                            downregulated = 'steelblue',
                                            ns = 'gray60'),
                             fill_title = '',
                             plot_title = NULL,
                             plot_subtitle = NULL,
                             plot_tag = NULL,
                             x_lab = 'Regulation',
                             cust_theme = microViz::theme_micro(),
                             show_txt = FALSE,
                             signif_digits = 2,
                             txt_color = 'black',
                             txt_size = 2.75,
                             txt_offset = 0.05,
                             txt_vjust = 0.5) {

    ## entry control -------

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(label_variable,
              p_variable,
              regulation_variable) %in% names(data))) {

      stop(paste("'p_variable', 'label_variable' or 'regulation_variable'",
                 "absent from the data.'"),
           call. = FALSE)

    }

    if(!inherits(cust_theme, 'theme')) {

      stop('Please provide a valid ggplot theme object.', call. = FALSE)

    }

    stopifnot(is.numeric(signif_level))
    stopifnot(is.numeric(regulation_level))
    stopifnot(is.numeric(top_regulated))

    top_regulated <- as.integer(top_regulated)

    stopifnot(is.logical(show_txt))
    stopifnot(is.numeric(signif_digits))

    signif_digits <- as.integer(signif_digits)

    stopifnot(is.character(txt_color))
    stopifnot(is.numeric(txt_size))
    stopifnot(is.numeric(txt_offset))
    stopifnot(is.numeric(txt_vjust))

    ## plotting data --------

    significant <- NULL
    regulation <- NULL
    regulation_raw <- NULL

    plot_data <-
      mutate(data,
             significant = ifelse(.data[[p_variable]] < signif_level,
                                  'yes', 'no'),
             regulation = ifelse(significant == 'no',
                                 'ns',
                                 ifelse(.data[[regulation_variable]] > regulation_level,
                                        'upregulated',
                                        ifelse(.data[[regulation_variable]] < -regulation_level,
                                               'downregulated', 'ns'))),
             regulation = factor(regulation,
                                 c('upregulated',
                                   'downregulated',
                                   'ns')),
             regulation_raw = ifelse(.data[[regulation_variable]] > 0,
                                     'upregulated', 'downregulated'),
             regulation_raw = factor(regulation_raw,
                                     c('upregulated', 'downregulated')))

    plot_data <- group_by(plot_data, regulation_raw)

    plot_data <-
      top_n(plot_data,
            n = top_regulated,
            wt = abs(.data[[regulation_variable]]))

    plot_data <- ungroup(plot_data)

    ## label position

    txt_pos <- NULL

    plot_data <-
      mutate(plot_data,
             txt_pos = (1 - txt_offset) * .data[[regulation_variable]])

    ## plotting --------

    bar_plot <- ggplot(plot_data,
                       aes(x = .data[[regulation_variable]],
                           y = stats::reorder(.data[[label_variable]],
                                              .data[[regulation_variable]]),
                           fill = regulation)) +
      geom_bar(color = 'black',
               stat = 'identity') +
      scale_fill_manual(values = fill_scale,
                        name = fill_title) +
      cust_theme +
      theme(axis.title.y = element_blank()) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           tag = plot_tag,
           x = x_lab)

    if(show_txt) {

      bar_plot <- bar_plot +
        geom_text(aes(label = signif(.data[[regulation_variable]],
                                     signif_digits),
                      x = txt_pos),
                  vjust = txt_vjust,
                  size = txt_size,
                  color = txt_color)

    }

    bar_plot

  }

# Heat map -------

#' Heat map of levels of variables of interest.
#'
#' @description
#' The functions draw a heat map of levels of the variables of interest in
#' data set subsets defined by a splitting factor.
#' While `heat_map()` operates with a single data frame, `common_data_frame()`
#' takes a list of data frames as the first argument.
#' In the former case, observations or samples are presented in the X axis and
#' variables are shown in the Y axis, color intensity codes for the variable
#' value.
#' In the later case, levels of variables of interest are averaged within
#' levels of the splitting factor (argument `spit_fct`); data sets are presented
#' in the X axis, variables are presented in the Y axis, and tile color codes
#' for the average.
#'
#' @details
#' Normalization of variables id done by \code{\link{zScores}}.
#' The `variable_classification` argument may be provided in several forms:
#'
#' * as a data frame with at least two columns. The first variable stores the
#' variable name and the second indicates the variable subset. This input is
#' accepted both by `heat_map()` and `common_heat_map()`.
#'
#' * as a list of data frames like the one described above. This form is
#' accepted only by `common_heat_map()`.
#'
#' * as a `clust_analysis` object generated with one of the tools of
#' `clustTools` package. The cluster assignment is extracted automatically.
#' This form works with both functions.
#'
#' * as a list of `clust_analysis` objects described above. This input is
#' accepted only by `common_heat_map()`. #'
#'
#' @return a `ggplot` graphics.
#'
#' @param data a data frame (`heat_map()`) or a list of data frames
#' (`common_heat_map()`).
#' @param variables a vector with names of the variables of interest.
#' @param split_fct name of the splitting factor.
#' @param normalize logical, should the data frame variables be normalized prior
#' to plotting?
#' @param norm_center defines centering of the variable during Z-score
#' calculation: mean (default), median, geometric mean or harmonic mean.
#' Ignored if `normalize = FALSE`. Please refer to \code{\link{zScores}}
#' for details.
#' @param norm_dispersion name of the dispersion statistic used for variable
#' normalization: standard deviation (SD) or standard error of the mean (SEM).
#' Please refer to \code{\link{zScores}} for details.
#' @param average_fun a function used to calculate average values of the
#' variables within levels of the splitting factor.
#' @param variable_classification optional, a data frame, `clust_analysis`
#' object or a list thereof as outlined in Details.
#' If not provided, the variables are classified with by their
#' specificity for the data set subsets defined by `split_fct` with
#' \code{\link{classify}}.
#' @param direction determines direction of the comparison between the data set
#' subsets, see: \code{\link{classify}}. Ignored if `variable_classification`
#' is provided.
#' @param facet logical, should the variable classification scheme be presented
#' by horizontal facets of the plot?
#' @param plot_title plot title.
#' @param plot_subtitle plot_subtitle. If not provided, numbers of observations
#' in the data set subsets will be displayed here.
#' @param x_lab X axis label.
#' @param y_lab Y axis label.
#' @param fill_lab title of the fill scale.
#' @param hide_x_axis_text logical, hide the X axis text?
#' @param cust_theme a custom `ggplot` theme.
#' @param color_scale a character vector of length 3, which defines the lower,
#' middle and upper point of the color scale.
#' @param midpoint optional, the middle point of the fill scale.
#' @param ... additional arguments passed to
#' \code{\link[ggplot2]{scale_fill_gradient2}} (`heat_map()`), or additional
#' arguments passed to `heat_map()` (`common_heat_map()`).
#'
#' @export

  heat_map <- function(data,
                       variables,
                       split_fct,
                       normalize = TRUE,
                       norm_center = c('mean',
                                       'median',
                                       'geo_mean',
                                       'harm_mean'),
                       norm_dispersion = c('sd', 'sem'),
                       variable_classification = NULL,
                       direction = '<',
                       facet = TRUE,
                       plot_title = NULL,
                       plot_subtitle = NULL,
                       x_lab = 'observation',
                       y_lab = 'variable',
                       fill_lab = 'Z-score',
                       hide_x_axis_text = FALSE,
                       cust_theme = theme_micro(),
                       color_scale = c('steelblue', 'black', 'firebrick'),
                       midpoint = NULL, ...) {

    ## input control ---------

    stopifnot(is.data.frame(data))

    if(any(!variables %in% names(data))) {

      stop('At least one variable is missing from the data.',
           call. = FALSE)

    }

    classes <- map_lgl(data[variables], is.numeric)

    if(any(!classes)) {

      stop("Some of the variables are not numeric.", call. = FALSE)

    }

    if(!split_fct %in% names(data)) {

      stop('The splitting factor is missing from the data.',
           call. = FALSE)

    }

    if(!is.factor(data[[split_fct]])) {

      stop('The splitting factor has to be a vector.',
           call. = FALSE)

    }

    stopifnot(is.logical(normalize))

    norm_center <- match.arg(norm_center[1],
                             c('mean', 'median', 'geo_mean', 'harm_mean'))

    norm_dispersion <- match.arg(norm_dispersion[1], c('sd', 'sem'))

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

    } else {

      variable_classification <-
        classify(data = data,
                 variables = variables,
                 split_fct = split_fct,
                 direction = direction)$classification

    }

    variable_classification <-
      set_names(variable_classification[, 1:2],
                c('variable', 'variable_subset'))

    stopifnot(is.logical(facet))
    stopifnot(inherits(cust_theme, 'theme'))

    if(length(color_scale) < 3) {

      stop('Not enough colors.', call. = FALSE)

    }

    ## plotting data ------

    observation <- NULL

    data <- mutate(data[c(split_fct, variables)],
                   observation = paste0('obs_', 1:nrow(data)))

    if(is.null(plot_subtitle)) {

      n_numbers <- count(data, .data[[split_fct]])

      plot_subtitle <-
        map2_chr(n_numbers[[1]], n_numbers[[2]],
                 paste, sep = ': n = ')

      plot_subtitle <- paste(plot_subtitle, collapse = ', ')

    }

    if(normalize) {

      data[variables] <- zScores(data[variables],
                                 center = norm_center,
                                 dispersion = norm_dispersion)

    }

    variable <- NULL
    value <- NULL

    data <-
      pivot_longer(data,
                   cols = all_of(variables),
                   names_to = 'variable',
                   values_to = 'value')

    data <- left_join(data,
                      variable_classification,
                      by = 'variable')

    plot_order <- unique(variable_classification$variable)

    data <-
      mutate(data,
             variable = factor(variable, plot_order))

    ## plotting -------

    if(is.null(midpoint)) {

      midpoint <- mean(range(data[['value']], na.rm = TRUE))

    }

    if(facet) {

      facet_formula <-
        as.formula(paste('variable_subset ~', split_fct))

    } else {

      facet_formula <-
        as.formula(paste('. ~', split_fct))

    }

    hm_plot <-
      ggplot(data,
             aes(x = observation,
                 y = variable,
                 fill = value)) +
      geom_tile() +
      facet_grid(facet_formula,
                 scales = 'free',
                 space = 'free') +
      scale_fill_gradient2(low = color_scale[[1]],
                           mid = color_scale[[2]],
                           high = color_scale[[3]],
                           midpoint = midpoint,
                           name = fill_lab, ...) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    if(hide_x_axis_text) {

      hm_plot <- hm_plot +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_blank())

    }

    hm_plot

  }

#' @rdname heat_map
#' @export

  common_heat_map <- function(data,
                              variables,
                              split_fct,
                              normalize = TRUE,
                              norm_center = c('mean',
                                              'median',
                                              'geo_mean',
                                              'harm_mean'),
                              norm_dispersion = c('sd', 'sem'),
                              variable_classification = NULL,
                              average_fun = colMeans,
                              plot_title = NULL,
                              plot_subtitle = NULL,
                              x_lab = 'data set',
                              y_lab = 'variable',
                              fill_lab = 'average Z-score',
                              hide_x_axis_text = FALSE,
                              cust_theme = theme_micro(),
                              color_scale = c("steelblue", "black", "firebrick"),
                              midpoint = 0, ...) {

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

    ## entry control for the remaining arguments -------

    stopifnot(is.logical(normalize))

    norm_center <- match.arg(norm_center[1],
                             c('mean',
                               'median',
                               'geo_mean',
                               'harm_mean'))

    norm_dispersion <- match.arg(norm_dispersion[1], c('sd', 'sem'))

    if(!is.function(average_fun)) {

      stop("'average_fun' has to be a function.", call. = FALSE)

    }

    stopifnot(inherits(cust_theme, 'theme'))

    if(length(color_scale) < 3) {

      stop('Not enough colors.', call. = FALSE)

    }

    ## plotting data -------

    if(is.null(plot_subtitle)) {

      plot_subtitle <- paste('data sets: n =', length(data))

    }

    data_names <- names(data)

    data <- map(data, ~.x[c(split_fct, variables)])

    if(normalize) {

      data <- map(data,
                  zScores,
                  variables = variables,
                  center = norm_center,
                  dispersion = norm_dispersion)

    }

    split_vec <- map(data, ~.x[[split_fct]])

    levs <- map(split_vec, levels)

    data_splits <-
      map2(data, split_vec,
           ~split(.x[variables], .y, drop = TRUE))

    data_splits <- map(data_splits, map, average_fun)

    data_splits <- map(data_splits, reduce, rbind)

    data_splits <- map(data_splits, as.data.frame)

    data_splits <-
      map2(data_splits, levs,
           ~mutate(.x, !!split_fct := factor(.y, .y)))

    data_set <- NULL

    plot_tbl <- map2_dfr(data_splits, names(data_splits),
                         ~mutate(.x, data_set = .y))

    ax_labs <-
      set_names(plot_tbl$data_set,
                paste0('obs_', 1:nrow(plot_tbl)))

    ## variable classification ----------

    if(inherits(variable_classification, 'clust_analysis')) {

      variable_classification <-
        variable_classification$clust_assignment[, c('observation', 'clust_id')]

    }

    if(is.data.frame(variable_classification)) {

      if(ncol(variable_classification) < 2) {

        stop("'variable_classification' has to heve at least two columns.",
             call. = FALSE)

      }

      variable_classification <- variable_classification[, c(1:2)]

    }

    if(inherits(variable_classification, 'list')) {

      ## extraction of the subset assignment from a list

      feature_clust <- NULL
      variable <- NULL
      score <- NULL

      clust_assign <-
        map(variable_classification,
            function(x) {

          if(inherits(x, 'clust_analysis')) {

            return(x$clust_assignment[, c('observation', 'clust_id')])

          } else {

            return(x[, 1:2])

          }

        })

      clust_assign <-
        map_dfr(clust_assign,
                set_names,
                c('variable', 'feature_clust'))

      if(!is.factor(clust_assign$feature_clust)) {

        clust_assign$feature_clust <- factor(clust_assign$feature_clust)

      }

      feature_levs <- levels(clust_assign$feature_clust)

      ## counting occurrences in the feature cluster
      ## and voting for the consensus assignment

      clust_assign <- split(clust_assign$feature_clust,
                            clust_assign$variable)

      clust_assign <- map(clust_assign,
                          ~sort(table(.x), decreasing = TRUE))

      clust_assign <- map(clust_assign,
                          ~data.frame(feature_clust = names(.x)[1],
                                      score = unname(.x[1])))

      variable_classification <-
        map2_dfr(clust_assign, names(clust_assign),
                 ~mutate(.x, variable = .y))

      variable_classification$feature_clust <-
        factor(variable_classification$feature_clust, feature_levs)

      variable_classification <- arrange(variable_classification, score)

      variable_classification <-
        variable_classification[, c('variable', 'feature_clust')]

    }

    ## heat map -------

    hm_plot <- heat_map(plot_tbl,
                        variables = variables,
                        split_fct = split_fct,
                        normalize = FALSE,
                        variable_classification = variable_classification,
                        plot_title = plot_title,
                        plot_subtitle = plot_subtitle,
                        x_lab = x_lab,
                        y_lab = y_lab,
                        cust_theme = cust_theme,
                        color_scale = color_scale,
                        midpoint = midpoint, ...)


    hm_plot <- hm_plot +
      scale_x_discrete(labels = ax_labs) +
      guides(x = guide_axis(angle = 45))

    if(hide_x_axis_text) {

      hm_plot <- hm_plot +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

    }

    hm_plot

  }

# Word cloud ------

#' Plot a word cloud.
#'
#' @description
#' Plots a word cloud e.g. of gene or GO term names. Font color and size may
#' code for additional parameters such as p value or fold regulation.
#'
#' @details
#' The function internally employs
#' \code{\link[ggwordcloud]{geom_text_wordcloud}} and
#' \code{\link[ggwordcloud]{geom_text_wordcloud_area}} from the `wordcloud`
#' package.
#'
#' @return a `ggplot` object.
#'
#' @param data a data frame.
#' @param label_variable name of the variable storing the text to be plotted.
#' @param split_fct optional, name of a variable that defines dataset subsets.
#' Each subset will be presented as a separate facet of the plot.
#' @param size_variable optional, name of a variable numeric whose value will be
#' size coded.
#' @param color_variable optional, name of a variable numeric whose value will be
#' color coded.
#' @param size_type type of size scaling: `size` makes font size correspond to
#' the `size_variable`, `area` makes font text related to the `size_variable`.
#' @param size_range a numeric vector of size 2, which defines the minimum and
#' maximum of the size scale.
#' @param color_scale a character vector of length 3, which defines the lower,
#' middle and upper point color of the color scale.
#' @param midpoint optional, the middle point of the color scale.
#' @param size_lab title of the size scale.
#' @param color_lab title of the color scale.
#' @param cust_theme a custom `ggplot` theme.
#' @param wrap logical, should each text be split into multiple lines? May
#' be useful in case of long labels, e.g. descriptions of biological processes.
#' @param len number of words per line. Relevant only if `wrap_text = TRUE`.
#' @param fraction_rotated fraction of words to be rotated with a 90 degree
#' angle. May improve visibility.
#' @param nrow number of rows in a faceted plot. Ignored if `split_fct = NULL`.
#' @param ... extra arguments passed to
#' \code{\link[ggwordcloud]{geom_text_wordcloud}} or
#' \code{\link[ggwordcloud]{geom_text_wordcloud_area}}
#'
#'@export

  plot_wordcloud <- function(data,
                             label_variable,
                             split_fct = NULL,
                             size_variable = NULL,
                             color_variable = NULL,
                             size_type = c('size', 'area'),
                             size_range = c(2, 6),
                             color_scale = c('steelblue', 'black', 'firebrick'),
                             midpoint = NULL,
                             size_lab = size_variable,
                             color_lab = color_variable,
                             cust_theme = ggplot2::theme_void(),
                             wrap = FALSE,
                             len = 3,
                             fraction_rotated = 0,
                             nrow = NULL, ...) {


    ## input control -----

    stopifnot(is.data.frame(data))

    if(!label_variable %in% names(data)) {

      stop("'label_variable' is missing from the data.",
           call. = FALSE)

    }

    if(!is.null(split_fct)) {

      if(!split_fct %in% names(data)) {

        stop("'split_fct' missing from the data.",
             call. = FALSE)

      }

      if(!is.factor(data[[split_fct]])) {

        stop("'split_fct' has to define a factor variable.",
             call. = FALSE)

      }

    }

    if(!is.null(size_variable)) {

      if(!size_variable %in% names(data)) {

        stop("'size_variable' is missing from the data.",
             call. = FALSE)

      }

    }

    if(!is.null(color_variable)) {

      if(!color_variable %in% names(data)) {

        stop("'size_variable' is missing from the data.",
             call. = FALSE)

      }

    }

    size_type <- match.arg(size_type[1],
                           c('size', 'area'))

    stopifnot(is.numeric(size_range))

    if(length(size_range) < 2) {

      stop("Improper length of 'size_range'",
           call. = FALSE)

    }

    size_range <- size_range[1:2]

    stopifnot(is.character(color_scale))

    if(length(color_scale) < 3) {

      stop("Improper length of 'color_scale'",
           call. = FALSE)

    }

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' has to be a ggplot theme object.",
           call. = FALSE)

    }

    stopifnot(is.logical(wrap))

    len <- as.integer(len)

    stopifnot(is.numeric(fraction_rotated))

    if(fraction_rotated < 0 | fraction_rotated > 1) {

      stop("Improper value of 'fraction_rotated'.", call. = FALSE)

    }

    ## plotting data ------

    if(wrap) {

      data[[label_variable]] <-
        map_chr(data[[label_variable]],
                wrap_text,
                len = len)

    }

    if(fraction_rotated > 0) {

      rot_angle <- NULL

      angles <- sample(c(0, 90),
                       size = nrow(data),
                       replace = TRUE,
                       prob = c(1 - fraction_rotated,
                                fraction_rotated))

      data <- mutate(data,
                     rot_angle = angles)

    } else {

      data <- mutate(data, rot_angle = 0)

    }

    ## geoms and metadata -------

    text_geom <-
      switch(size_type,
             size = ggwordcloud::geom_text_wordcloud(...),
             ares = ggwordcloud::geom_text_wordcloud_area(...))

    size_scale <-
      switch(size_type,
             size = ggplot2::scale_size_continuous(range = size_range),
             area = ggplot2::scale_size_area(max_size = size_range[2]))

    if(!is.null(color_variable)) {

      if(is.null(midpoint)) {

        midpoint <- mean(range(data[[color_variable]], na.rm = TRUE))

      }

    } else {

      midpoint <- 0

    }

    ## plotting -------

    if(is.null(size_variable) & is.null(color_variable)) {

      word_plot <-
        ggplot(data,
               aes(label = .data[[label_variable]],
                   angle = rot_angle))

    } else if(!is.null(size_variable) & is.null(color_variable)) {

      word_plot <-
        ggplot(data,
               aes(label = .data[[label_variable]],
                   size = .data[[size_variable]],
                   angle = rot_angle))

    } else if(is.null(size_variable) & !is.null(color_variable)) {

      word_plot <-
        ggplot(data,
               aes(label = .data[[label_variable]],
                   color = .data[[color_variable]],
                   angle = rot_angle))

    } else {

      word_plot <-
        ggplot(data,
               aes(label = .data[[label_variable]],
                   size = .data[[size_variable]],
                   color = .data[[color_variable]],
                   angle = rot_angle))

    }

    if(!is.null(split_fct)) {

      wrap_formula <- as.formula(paste('~', split_fct))

      word_plot <- word_plot +
        ggplot2::facet_wrap(wrap_formula, nrow = nrow)

    }

    word_plot +
      text_geom +
      size_scale +
      ggplot2::scale_color_gradient2(low = color_scale[[1]],
                                     mid = color_scale[[2]],
                                     high = color_scale[[3]],
                                     midpoint = midpoint) +
      cust_theme +
      theme(legend.position = 'right') +
      labs(color = color_lab,
           size = size_lab)

  }

# Bubble plot -----

#' Draw a bubble plot.
#'
#' @description
#' A general form of a bubble plot: variables are presented in the Y axis,
#' splitting factor levels are shown in the x axis, size and fill of the bubbles
#' code for additional variables and, as an option, can be labeled with a text.
#'
#' @return a `ggplot` object.
#'
#' @param data a data frame.
#' @param x_variable variable to be presented in the X axis of the plot.
#' @param y_variable variable to be presented in the Y axis of the plot.
#' @param size_variable name of the variable defining the bubble size.
#' @param color_variable name of the variable defining the bubble color.
#' @param label_variable variable to be presented as text next to the data
#' points.
#' @param point_color color of the points, relevant only if color variable is
#' not provided.
#' @param color_scale a character vector of length 3, which defines the lower,
#' middle and upper point color of the fill scale.
#' @param midpoint the middle point of the color scale, optional.
#' @param size_type type of size scaling: `size` makes bubble radius correspond
#' to the `size_variable`, `area` makes bubble area relates to the
#' `size_variable`.
#' @param size_range a numeric vector of size 2, which defines the minimum and
#' maximum of the size scale.
#' @param plot_title plot_title.
#' @param plot_subtitle plot_subtitle.
#' @param x_lab X axis title.
#' @param y_lab Y axis title.
#' @param size_lab title of the size scale.
#' @param color_lab title of the color scale.
#' @param cust_theme a custom `ggplot` theme.
#' @param txt_hjust horizontal justification of the text.
#' @param txt_vjust vertical justification of the text.
#' @param txt_size text size.
#' @param ... extra arguments passed to the size and color scales, e.g.
#' \code{\link[ggplot2]{scale_fill_gradient2}}.
#'
#' @export

  plot_bubble <- function(data,
                          x_variable,
                          y_variable,
                          size_variable,
                          color_variable = size_variable,
                          label_variable = NULL,
                          point_color = 'steelblue',
                          color_scale = c(low = 'steelblue',
                                          mid = 'white',
                                          high = 'firebrick'),
                          midpoint = NULL,
                          size_type = c('size', 'area'),
                          size_range = c(0.2, 4),
                          plot_title = NULL,
                          plot_subtitle = NULL,
                          x_lab = x_variable,
                          y_lab = y_variable,
                          size_lab = NULL,
                          color_lab = NULL,
                          cust_theme = microViz::theme_micro(),
                          txt_hjust = -1.4,
                          txt_vjust = 0.5,
                          txt_size = 2.75, ...) {

    ## input check -------

    stopifnot(is.data.frame(data))

    if(any(!c(x_variable, y_variable, size_variable) %in% names(data))) {

      stop(paste("Data has to contain 'x_variable', 'y_variable' and",
                 "'size_variable' variables."),
           call. = FALSE)

    }

    if(!is.null(color_variable)) {

      if(!color_variable %in% names(data)) {

        stop("'color_variable' missing from the data.", call. = FALSE)

      }

    } else {

      color_var <- NULL

      data <- mutate(data, color_var = point_color)

      color_variable <- 'color_var'

    }

    if(!is.null(label_variable)) {

      if(!label_variable %in% names(data)) {

        stop("'label_variable' is missing from the data.",
             call. = FALSE)

      }

    }

    if(length(color_scale) < 3) {

      stop("'color_scale' needs to have at least 3 values.",
           call. = FALSE)

    }

    if(!is.null(midpoint)) stopifnot(is.numeric(midpoint))

    size_type <- match.arg(size_type[1],
                           c('size', 'area'))

    stopifnot(is.numeric(size_range))

    if(length(size_range) < 2) {

      stop("'size_range' has to have at least 2 values.",
           call. = FALSE)

    }

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' has to be a valid ggplot theme.",
           call. = FALSE)

    }

    stopifnot(is.numeric(txt_size))
    stopifnot(is.numeric(txt_hjust))
    stopifnot(is.numeric(txt_vjust))

    ## geoms and scales -------

    size_scale <-
      switch(size_type,
             size = ggplot2::scale_radius(range = size_range, ...),
             area = ggplot2::scale_size_area(max_size = size_range[2], ...))

    if(!is.na(color_variable)) {

      if(is.null(midpoint)) {

        midpoint <- mean(range(data[[color_variable]], na.rm = TRUE))

      }

    }

    ## plotting -------

    if(is.null(color_variable)) {

      bubble_plot <-
        ggplot(data,
               aes(x = .data[[x_variable]],
                   y = .data[[y_variable]],
                   size = .data[[size_variable]]))

    } else {

      bubble_plot <-
        ggplot(data,
               aes(x = .data[[x_variable]],
                   y = .data[[y_variable]],
                   size = .data[[size_variable]],
                   fill = .data[[color_variable]]))

    }

    bubble_plot <- bubble_plot +
      geom_point(shape = 21) +
      size_scale +
      scale_fill_gradient2(low = color_scale[1],
                           mid = color_scale[2],
                           high = color_scale[3],
                           midpoint = midpoint, ...) +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab,
           fill = color_lab,
           size = size_lab)

    if(!is.null(label_variable)) {

      bubble_plot <- bubble_plot +
        geom_text(aes(label = .data[[label_variable]]),
                  size = txt_size,
                  hjust = txt_hjust,
                  vjust = txt_vjust)

    }

    bubble_plot

  }

# Box plot from distribution statistics -------

#' Plot a box plot given median, interquartile and 95% percentile range.
#'
#' @description
#' Generates a `ggplot` box plot from distribution statistics: median,
#' 25% and 75% percentiles, and 2.5% and 95% percentiles.
#' The function may be useful at plotting distribution statistic for very large
#' data sets (e.g. whole transcriptome). In this case, a box plot made with
#' pre-calculated distribution statistic is going to spare a lot of memory.
#'
#' @return a `ggplot` object.
#'
#' @param data a data frame with the following columns: `median`,
#' `perc25`, `perc75`, `perc025`, and `perc975` storing the stats.
#' @param x_variable name of the variable presented in the X axis.
#' @param fill_variable name of the variable defining color of the box plots.
#' @param fill_color fill color for the boxes, ignored if `fill_variable = NULL`.
#' @param width width of the boxes.
#' @param alpha alpha of the box plots.
#' @param plot_title plot title,
#' @param plot_subtitle plot subtitle.
#' @param x_lab x axis title.
#' @param y_lab y axis title.
#' @param cust_theme custom `ggplot` theme
#'
#' @export

  box_from_stats <- function(data,
                             x_variable,
                             fill_variable = NULL,
                             fill_color = 'steelblue',
                             width = 0.8,
                             alpha = 1,
                             plot_title = NULL,
                             plot_subtitle = NULL,
                             x_lab = NULL,
                             y_lab = NULL,
                             cust_theme = microViz::theme_micro()) {

    ## input control -------

    if(!is.data.frame(data)) {

      stop("'data' has to be a data frame.", call. = FALSE)

    }

    if(!x_variable %in% names(data)) {

      stop("'x_variable' missing from the input data.", call. = FALSE)

    }

    if(!is.factor(data[[x_variable]])) {

      data[[x_variable]] <- factor(data[[x_variable]])

    }

    if(!is.null(fill_variable)) {

      if(!fill_variable %in% names(data)) {

        stop("'fill_variable' is missing from the input data.", call. = FALSE)

      }

    }

    if(any(!c('median', 'perc25', 'perc75', 'perc025', 'perc975') %in% names(data))) {

      error_txt <-
        paste("'data' has to have the following variables:",
              "'median', 'perc25', 'perc75', 'perc025', and 'perc975'")

     stop(error_txt, call. = FALSE)

    }

    stopifnot(is.numeric(width))
    stopifnot(is.numeric(alpha))

    if(!inherits(cust_theme, 'theme')) {

      stop("'cust_theme' has to be a valid ggplot theme.",
           call. = FALSE)

    }

    ## plotting ---------

    ## plots a box plot given a data frame with the columns
    ## `median`, `perc25`, `perc75`, `perc025` and `perc975`
    ## storing the median expression, interquartile and 95 percentile ranges

    center_pos <- NULL
    median <- NULL
    perc25 <- NULL
    perc75 <- NULL
    perc025 <- NULL
    perc975 <- NULL

    data <- mutate(data, center_pos = as.numeric(.data[[x_variable]]))

    if(is.null(fill_variable)) {

      box_plot <-
        ggplot(data,
               aes(x = .data[[x_variable]],
                   y = median)) +
        geom_errorbar(aes(ymin = perc025,
                          ymax = perc975),
                      width = 0) +
        geom_rect(aes(xmin = center_pos - 0.5 * width,
                      xmax = center_pos + 0.5 * width,
                      ymin = perc25,
                      ymax = perc75),
                  color = 'black',
                  fill = fill_color)

    } else {

      box_plot <-
        ggplot(data,
               aes(x = .data[[x_variable]],
                   y = median,
                   fill = .data[[fill_variable]])) +
        geom_errorbar(aes(ymin = perc025,
                          ymax = perc975),
                      width = 0) +
        geom_rect(aes(xmin = center_pos - 0.5 * width,
                      xmax = center_pos + 0.5 * width,
                      ymin = perc25,
                      ymax = perc75),
                  color = 'black')

    }

    box_plot +
      geom_segment(aes(x = center_pos - 0.5 * width,
                       xend = center_pos + 0.5 * width,
                       y = median,
                       yend = median),
                   color = 'black') +
      cust_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

  }

# END -----
