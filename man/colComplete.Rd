% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/numeric_functions.R
\name{colComplete}
\alias{colComplete}
\alias{rowComplete}
\alias{colComplete.data.frame}
\alias{colComplete.matrix}
\alias{rowComplete.data.frame}
\alias{rowComplete.matrix}
\title{Numbers of complete cases in rows and columns of a matrix or a data frame.}
\usage{
colComplete(x, ...)

rowComplete(x, ...)

\method{colComplete}{data.frame}(x, ...)

\method{colComplete}{matrix}(x, ...)

\method{rowComplete}{data.frame}(x, ...)

\method{rowComplete}{matrix}(x, ...)
}
\arguments{
\item{x}{a data frame or a matrix.}

\item{...}{extra arguments passed to methods; currently none.}
}
\value{
a numeric vector.
}
\description{
The functions calculate numbers of complete, i.e. non-\code{NA} cases in columns
or rows of a data frame or matrix.
}
