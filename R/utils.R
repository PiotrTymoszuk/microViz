# Exported utilities.


# Vector utils ------

#' Split a vector into chunks of a defined size.
#'
#' @description Splits the vector into chuns of a defined size: the last one
#' can be shorter ('ceiling').
#' @param x a vector.
#' @param len chunk length.
#' @return a list of vector chunks.
#' @export

  split_vec <- function(x, len) {

    len <- as.integer(len)

    if(len >= length(x)) return(x)

    split(x, ceiling(seq_along(x)/len))

  }


# Text utils -----

#' Wrap text with a defined word count per line.
#'
#' @description Converts a character vector into a text with the defined line
#' length (word count). If a non-character vector is provided, it will be
#' coerced to a character one.
#' @param x a vector.
#' @param len line length in words.
#' @param word_sep word separator, comma by default.
#' @param line_sep line separator, comma + newline by default.
#' @return a text value.
#' @export

  wrap_vector <- function(x, len = 5, word_sep = ', ', line_sep = ',\n') {

    len <- as.integer(len)

    if(len >= length(x)) {

      return(paste(x, collapse = word_sep))

    }

    split_txt <-microViz::split_vec(x, len)

    split_txt <- purrr::map_chr(split_txt, paste, collapse = word_sep)

    paste(split_txt, collapse = line_sep)

  }

# Graphics theme -----

#' Deafult graphics theme for microViz tools.
#'
#' @description The default ggplot graphics theme for graphical output of the
#' package functions.
#' @description A ggplot theme object.
#' @export

  theme_micro <- function() {

    ggplot2::theme_classic() +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 8,
                                                         face = 'plain',
                                                         color = 'black'),
                     legend.title = ggplot2::element_text(size = 8,
                                                          face = 'plain',
                                                          color = 'black'),
                     axis.text = ggplot2::element_text(size = 8,
                                                       face = 'plain',
                                                       color = 'black'),
                     axis.title = ggplot2::element_text(size = 8,
                                                        face = 'plain',
                                                        color = 'black'),
                     strip.text = ggplot2::element_text(size = 8,
                                                        face = 'plain',
                                                        color = 'black'),
                     strip.background = ggplot2::element_rect(color = 'black',
                                                              fill = 'gray95'),
                     plot.tag = ggplot2::element_text(size = 8,
                                                      face = 'plain',
                                                      color = 'black',
                                                      margin = ggplot2::margin(8, 5, 5, 5),
                                                      hjust = 0),
                     plot.tag.position = 'bottom',
                     plot.title = ggplot2::element_text(size = 8,
                                                        face = 'plain',
                                                        color = 'black'),
                     plot.subtitle = ggplot2::element_text(size = 8,
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
