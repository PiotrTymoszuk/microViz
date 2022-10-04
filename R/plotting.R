# Plotting functions.

# Volcano plot -----

#' Draw a volcano plot.
#'
#' @description Generates a Volcano plot with the effect size/regulation
#' statistic value on the X axis and -log10 p value on the y axis. Regulation
#' and significance is coded by the point fill color. Top n significant
#' genes/objects may be labeled with their names.
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
#' @return a ggplot object.
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
                           point_wjitter = 0) {

    ## entry control

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(regulation_variable, p_variable) %in% names(data))) {

      stop('Regulation_variable or p_variable absent from the data.',
           call. = FALSE)

    }

    top_significant <- as.integer(top_significant)

    if(!is.null(label_variable)) {

      if(!label_variable %in% names(data)) {

        stop('label_variable missing from data.', call. = FALSE)

      }

    }

    label_type <- match.arg(label_type[1], c('label', 'text'))

    if(!any(class(cust_theme) != 'theme')) {

      stop('Please provide a valid ggplot theme object.', call. = FALSE)

    }

    ## plotting data

    plot_tbl <- dplyr::filter(data, complete.cases(data))

    plot_tbl <- dplyr::mutate(plot_tbl,
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

    ## numbers of regulated genes

    if(is.null(plot_tag)) {

      n_genes <- dplyr::count(plot_tbl, regulation, .drop = FALSE)

      plot_tag <- paste0('upregulated: n = ', n_genes$n[1],
                         ', downregulated: n = ', n_genes$n[2])

    }

    ## plotting

    volc <- ggplot2::ggplot(plot_tbl,
                            ggplot2::aes(x = .data[[regulation_variable]],
                                         y = -log10(.data[[p_variable]]),
                                         fill = regulation)) +
      ggplot2::geom_point(size = 2,
                          shape = 21,
                          alpha = point_alpha,
                          position = ggplot2::position_jitter(width = point_wjitter,
                                                              height = point_hjitter)) +
      ggplot2::geom_hline(yintercept = -log10(signif_level), lty = 2) +
      ggplot2::scale_fill_manual(values = fill_scale,
                                 name = fill_title) +
      cust_theme +
      ggplot2::labs(x = x_lab,
                    y = y_lab,
                    title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag)

    if(regulation_level > 0) {

      volc <- volc  +
        ggplot2::geom_vline(xintercept = -regulation_level, lty = 2) +
        ggplot2::geom_vline(xintercept = regulation_level, lty = 2)

    }

    if(top_significant > 0 & !is.null(label_variable)) {

      desc_tbl <- purrr::map(c('upregulated', 'downregulated'),
                             ~dplyr::filter(plot_tbl, regulation == .x))

      desc_tbl <- purrr::map_dfr(desc_tbl,
                                 ~dplyr::top_n(.x,
                                               n = top_significant,
                                               -.data[[p_variable]]))

      if(nrow(desc_tbl) == 0) {

        return(volc)

      }

      if(label_type == 'label') {

        volc <- volc +
          ggrepel::geom_label_repel(data = desc_tbl,
                                    ggplot2::aes(label = .data[[label_variable]]),
                                    size = txt_size,
                                    color = txt_color,
                                    label.padding = 0.1,
                                    box.padding = 0.1,
                                    show.legend = FALSE,
                                    fontface = txt_face)

      } else {

        volc <- volc +
          ggrepel::geom_text_repel(data = desc_tbl,
                                   ggplot2::aes(label = .data[[label_variable]]),
                                   size = txt_size,
                                   color = txt_color,
                                   box.padding = 0.1,
                                   show.legend = FALSE,
                                   fontface = txt_face)

      }

    }

    return(volc)

  }

# GSEA plot -----

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

    ## entry control

    ## entry control

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

    ## plotting data

    plot_tbl <- dplyr::filter(data, complete.cases(data))

    plot_tbl <- dplyr::mutate(plot_tbl,
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

    plot_tbl <- dplyr::arrange(plot_tbl, -.data[[regulation_variable]])

    plot_tbl <- dplyr::mutate(plot_tbl, plot_order = 1:nrow(plot_tbl))

    ## numbers of regulated items

    if(is.null(plot_tag)) {

      n_genes <- dplyr::count(plot_tbl, regulation)

      plot_tag <- paste0('upregulated: n = ', n_genes$n[1],
                         ', downregulated: n = ', n_genes$n[2])

    }

    ## plot

    bar <- ggplot2::ggplot(plot_tbl,
                           ggplot2::aes(x = plot_order,
                                        y = .data[[regulation_variable]])) +
      ggplot2::geom_bar(stat = 'identity',
                        alpha = bar_alpha,
                        ggplot2::aes(fill = regulation)) +
      ggplot2::scale_fill_manual(values = fill_scale,
                                 name = fill_title) +
      cust_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::labs(title = plot_title,
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
#' @description Draws a Forest plot with the regulation estimates and,
#' optionally confidence intervals for all items presente in a data frame.
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
#' @return a ggplot object.
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

    ## entry control

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

    ## plotting table

    plot_tbl <-
      dplyr::mutate(data,
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

    ## plotting

    forest <-
      ggplot2::ggplot(plot_tbl,
                      ggplot2::aes(x = .data[[regulation_variable]],
                                   y = reorder(.data[[label_variable]],
                                               .data[[regulation_variable]]),
                                   color = regulation,
                                   fill = regulation)) +
      ggplot2::scale_fill_manual(values = fill_scale,
                                 name = fill_title) +
      ggplot2::scale_color_manual(values = fill_scale,
                                  name = fill_title) +
      cust_theme +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag,
                    x = x_lab) +
      ggplot2::geom_vline(xintercept = 0,
                          linetype = 'dashed')

    if(!is.null(lower_ci_variable) & !is.null(upper_ci_variable)) {

      forest <- forest +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[lower_ci_variable]],
                                             xmax = .data[[upper_ci_variable]]),
                                height = 0)

    }

    forest <- forest +
      ggplot2::geom_point(size = 2,
                          shape = 16)

    ## optional labeling

    if(!show_txt) {

      return(forest)

    }

    desc_tbl <- dplyr::mutate(plot_tbl,
                              plot_lab = signif(.data[[regulation_variable]],
                                                signif_digits))

    if(show_ci_txt & all(!is.null(c(lower_ci_variable, upper_ci_variable)))) {

      desc_tbl <-
        dplyr::mutate(desc_tbl,
                      plot_lab = paste0(plot_lab,
                                        ' [', signif(.data[[lower_ci_variable]],
                                                     signif_digits),
                                        ' - ', signif(.data[[upper_ci_variable]],
                                                      signif_digits), ']'))

    }

    forest <- forest +
      ggplot2::geom_text(data = desc_tbl,
                         ggplot2::aes(label = plot_lab),
                         size = txt_size,
                         hjust = txt_hjust,
                         vjust = txt_vjust)

    return(forest)

  }

# Top regulated Forest plot ------

#' Draw regulation statistics for top regulated genes.
#'
#' @description Draws a Forest plot with the regulation estimates and,
#' optionally confidence intervals for the n most strongly regulated items.
#' @inheritParams plot_forest
#' @param top_regulated top n regulated genes/objects to be presented in
#' the plot.
#' @return a ggplot object.
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


    ## plotting data

    top_regulated <- as.integer(top_regulated)

    plot_tbl <-
      dplyr::mutate(data,
                    reg_sign = ifelse(.data[[regulation_variable]] > 0,
                                      'up', 'down'),
                    reg_sign = factor(reg_sign, c('up', 'down')))

    plot_tbl <- dplyr::group_by(plot_tbl, reg_sign)

    plot_tbl <- dplyr::top_n(plot_tbl,
                             n = top_regulated,
                             abs(.data[[regulation_variable]]))

    plot_tbl <- dplyr::ungroup(plot_tbl)

    ## plotting

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

# P value plotting -------

#' Plot top p values as a bar plot.
#'
#' @description Plots top n p values as a bar plot.
#' @param data a data frame.
#' @param p_variable variable storing the p values.
#' @param label_variable variable storing the gene/object names.
#' @param signif_level significance threshold, p = 0.05 be default.
#' @param top_significant top n significant genes/objects to be presented in
#' the plot.
#' @param fill_scale regulation colors, a vector of two elements: for
#' significant and non-significant items.
#' @param x_lab x axis title.
#' @param fill_title title of the point fill legend.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag text.
#' @param cust_theme custom ggplot theme.
#' @return a ggplot object.
#' @export

  plot_signifcant <- function(data,
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

    ## entry control

    if(!is.data.frame(data)) {

      stop('Please provide a data frame as data.', call. = FALSE)

    }

    if(any(!c(label_variable,
              p_variable) %in% names(data))) {

      stop('p_variable or label_variable absent from the data.',
           call. = FALSE)

    }

    if(!any(class(cust_theme) != 'theme')) {

      stop('Please provide a valid ggplot theme object.', call. = FALSE)

    }

    top_significant <- as.integer(top_significant)

    ## plotting table

    plot_tbl <- dplyr::top_n(data, n = top_significant, -.data[[p_variable]])

    plot_tbl <- dplyr::mutate(plot_tbl,
                              significant = ifelse(.data[[p_variable]] < signif_level,
                                                   'significant',
                                                   'ns'),
                              significant = factor(significant,
                                                   c('significant', 'ns')))

    ## plotting

    ggplot2::ggplot(plot_tbl,
                    ggplot2::aes(x = -log10(.data[[p_variable]]),
                                 y = reorder(.data[[label_variable]],
                                             -.data[[p_variable]]),
                                 fill = significant)) +
      ggplot2::geom_bar(stat = 'identity',
                        color = 'black') +
      ggplot2::geom_vline(xintercept = -log10(signif_level),
                          linetype = 'dashed') +
      ggplot2::scale_fill_manual(values = fill_scale,
                                 name = fill_title) +
      cust_theme +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag,
                    x = x_lab)

  }

# END -----
