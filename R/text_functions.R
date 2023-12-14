# Functions for operations on text and character vectors

#' @include imports.R

  NULL

# Text vectors ------

#' Split a vector into chunks of a defined size.
#'
#' @description
#' Splits the vector into chunks of a defined size: the last one
#' can be shorter ('ceiling').
#' @param x a vector.
#' @param len chunk length.
#' @return a list of vector chunks.
#'
#' @export

  split_vec <- function(x, len) {

    len <- as.integer(len)

    if(len >= length(x)) return(x)

    split(x, ceiling(seq_along(x)/len))

  }

#' Wrap text with a defined word count per line.
#'
#' @description
#' Converts a character vector or a text chunk into a text with the defined line
#' length (word count). If a non-character vector is provided, it will be
#' coerced to a character vector.
#'
#' @param x a vector or a text
#' @param len line length in words.
#' @param word_sep word separator, comma (`wrap_vec()`) or
#' space (`wrap_text()`) by default.
#' @param line_sep line separator, comma + newline (`wrap_vector()`) or
#' newline (`wrap_text()`) by default.
#' @return a text value.
#'
#' @export

  wrap_vector <- function(x,
                          len = 5,
                          word_sep = ', ',
                          line_sep = ',\n') {

    len <- as.integer(len)

    if(len >= length(x)) {

      return(paste(x, collapse = word_sep))

    }

    split_txt <- split_vec(x, len)

    split_txt <- map_chr(split_txt, paste, collapse = word_sep)

    paste(split_txt, collapse = line_sep)

  }

#' @rdname wrap_vector
#' @export

  wrap_text <- function(x,
                        len = 5,
                        word_sep = ' ',
                        line_sep = '\n') {

    len <- as.integer(len)

    txt_vector <- stri_split(x, regex = '\\s+')[[1]]

    wrap_vector(txt_vector,
                len = len,
                word_sep = word_sep,
                line_sep = line_sep)

  }

# END -----
