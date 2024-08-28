# Exported and non-exported utilities.


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
            panel.grid.major = element_line(color = 'gray90'))

  }

# Data frame check ------

#' Check for a numeric data frame.
#'
#' @description
#' The functions checks if all columns of a data frame or matrix are numeric or
#' binary.
#'
#' @param x a data frame or a matrix.
#'
#' @return a logical value.

  check_df <- function(x) {

    stopifnot(is.data.frame(x))

    col_check <- map_lgl(x, is.numeric)

    if(any(!col_check)) {

      stop('A numeric data frame is required', call. = FALSE)

    }

  }

#' @rdname check_df

  check_mtx <- function(x) {

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) {

      stop('A numeric matrix is required.', call. = FALSE)

    }

  }

#' @rdname check_df

  check_binary_df <- function(x) {

    check_df(x)

    lev_check <- map(x, ~.x[!is.na(.x)])

    lev_check <- map_lgl(lev_check, function(var) all(var == 0 | var == 1))

    if(any(!lev_check)) {

      stop('At least one of variables is not binary.', call. = FALSE)

    }

  }

# END ----
