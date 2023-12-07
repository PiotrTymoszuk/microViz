# Exported utilities.


#' @include imports.R

  NULL

# Graphics theme -----

#' Deafult graphics theme for microViz tools.
#'
#' @description The default ggplot graphics theme for graphical output of the
#' package functions.
#' @description A ggplot theme object.
#' @export

  theme_micro <- function() {

    ggplot2::theme_classic() +
      theme(legend.text = element_text(size = 8,
                                       face = 'plain',
                                       color = 'black'),
            legend.title = element_text(size = 8,
                                        face = 'plain',
                                        color = 'black'),
            axis.text = element_text(size = 8,
                                     face = 'plain',
                                     color = 'black'),
            axis.title = element_text(size = 8,
                                      face = 'plain',
                                      color = 'black'),
            strip.text = element_text(size = 8,
                                      face = 'plain',
                                      color = 'black'),
            strip.background = ggplot2::element_rect(color = 'black',
                                                     fill = 'gray95'),
            plot.tag = element_text(size = 8,
                                    face = 'plain',
                                    color = 'black',
                                    margin = ggplot2::margin(8, 5, 5, 5),
                                    hjust = 0),
            plot.tag.position = 'bottom',
            plot.title = element_text(size = 8,
                                      face = 'plain',
                                      color = 'black'),
            plot.subtitle = element_text(size = 8,
                                         face = 'plain',
                                         color = 'black'),
            plot.margin = ggplot2::margin(t = 3,
                                          l = 3,
                                          r = 3,
                                          b = 3,
                                          unit = 'mm'),
            panel.grid.major = ggplot2::element_line(color = 'gray90'))

  }

# END ----
